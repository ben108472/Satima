"""掃描帶取樣：輸出像素 -> 掃描帶的小數座標 -> 高斯距離加權，或在掃描帶的「索引空間」內插

微波輻射計的取樣很粗 (5~25 km)，而且兩個方向的間距不同：圓錐掃描的掃描線間距 (10~15 km) 常是
沿掃描方向像素間距 (4~6 km) 的 2~3.5 倍。這裡先把每個輸出像素換算成掃描帶中的小數座標
(第幾條掃描, 第幾個像素)，之後有兩類取值的方法：

  gauss (預設；NRL GeoIPS 的做法，即 pyresample 的 resample_gauss)
      像素值 = 周圍取樣點的加權平均，權重 exp(-d² / σ²)，d 為取樣點到像素的實際距離。各方向同樣平滑，
      所以兩個方向的間距相差很多時，與掃描線斜交的眼牆、雨帶也不會變成鋸齒；代價是解析度降到較粗的方向。
      σ 取「該處較粗方向取樣間距的倍數」(隨儀器、頻道、掃描位置調整) 或固定的公里數。
      由小數座標直接取附近的取樣點 (不必對每個像素查 KD-tree)，σ 的 GAUSS_CUT 倍以內的取樣點全部計入
      (pyresample 預設只取最近的 8 個，鄰點的組合改變時會有細微的不連續)；多組掃描 (A / B) 的取樣點一起加權。
  nearest / bilinear / cubic (索引空間)
      兩個方向各自以真正相鄰的取樣點內插：保留沿掃描方向的細節，但斜向的結構在取樣不對稱時呈鋸齒狀。

步驟
  1. Locator：掃描帶經緯度轉成單位球上的 3-D 座標 (不受換日線、極區影響)，只取與地圖範圍相交的掃描列，
     建 KD-tree。每個輸出像素先找最近的取樣點，以該點沿軌 / 沿掃描兩個方向的局部向量解出小數座標，
     再以所在格子的雙線性曲面做牛頓修正 (與之後內插所用的幾何一致)。
  2. 取值 ([swath] interpolation)：
       gauss     高斯距離加權 (見上)
       nearest   關閉內插：每個取樣點畫成一塊 (原始取樣的樣子)
       bilinear  雙線性
       cubic     三次 Hermite，斜率加上限制 (單調)：平滑，但不會在眼牆、海岸這類銳利邊界產生過衝的假冷 / 暖環
     周圍有無效值時自動退到較低階的方法 (cubic -> bilinear -> nearest)。
  3. A / B 掃描合併 (AMSR2、AMSR3 的 89 GHz)：兩組掃描各自有小數座標；沿掃描方向先在各條掃描線上內插，
     沿軌方向再取 A、B 掃描線的「聯集」中最近的上下各兩條內插。掃描帶中央 A、B 相間 (約 5 km)，
     邊緣兩者幾乎重疊，這個做法都能正確處理 (不需假設兩者等間隔交錯)。
  4. 掃描帶外、缺掃描 (相鄰掃描的間距超過正常值 GAP 倍) 為無資料 (NaN)。
"""
import numpy as np
from scipy.spatial import cKDTree

from ..core import parallel
from ..core.i18n import tr

try:                       # 內插的 numba 核心 (nogil，可多執行緒)；沒有 numba 時用 numpy
    from . import _swath_nb as _nb
except ImportError:
    _nb = None

EARTH_KM = 6371.0
PAD = 4                     # 地圖範圍外多取的掃描列數 (三次內插與斜率要用到鄰列)
MARGIN_KM = 150.            # 判斷掃描列是否與地圖相交時，地圖範圍外加的距離
GAP = 2.5                   # 相鄰掃描間距超過中位數的幾倍視為缺掃描
FAR = 1.2                   # 離最近取樣點超過「半格沿軌 + 半格沿掃描」的幾倍視為掃描帶外
CHUNK = 1 << 18             # 取值時每次處理的像素數 (控制暫存記憶體)
LATTICE = 8                 # 小數座標最多每隔幾個輸出像素實際反算一次 (其餘由格子點內插)
EPS = 1e-6
METHODS = ('nearest', 'bilinear', 'cubic')      # 索引空間的方法 (由低階到高階)
ALL_METHODS = METHODS + ('gauss',)
GAUSS_FACTOR = 0.7          # gauss：σ 預設為較粗方向取樣間距的幾倍
GAUSS_CUT = 2.5             # gauss：只計入距離在 σ 的幾倍以內的取樣點 (權重 exp(-6.25) = 0.002)
GAUSS_FILL = 1.05           # gauss：離最近的有效取樣點超過較粗間距的幾倍視為無資料 (只缺一點 / 一條掃描時由周圍補上)
GAUSS_MAXW = 16             # gauss：每個方向最多往外找幾個取樣點


def unit(lat, lon):
    """經緯度 (度) -> 單位球上的 3-D 座標 (..., 3)，float64；無效值為 NaN。"""
    p, l = np.deg2rad(np.asarray(lat, dtype=np.float64)), np.deg2rad(np.asarray(lon, dtype=np.float64))
    cp = np.cos(p)
    return np.stack([cp * np.cos(l), cp * np.sin(l), np.sin(p)], axis=-1)


def _dot(a, b):
    return np.einsum('ij,ij->i', a, b)


class Locator:
    """一組地理定位在某個地圖範圍內的定位器：輸出像素的經緯度 -> 掃描帶小數座標。"""

    def __init__(self, geo, center, radius):
        """center：地圖範圍中心的單位向量；radius：涵蓋地圖範圍的弦長半徑 (單位球，已含邊界)。"""
        xyz = unit(geo.lat, geo.lon)
        with np.errstate(invalid='ignore'):
            inside = np.linalg.norm(xyz - center, axis=-1) <= radius
        rows = np.flatnonzero(inside.any(axis=1))
        self.empty = rows.size == 0 or geo.shape[0] < 2 or geo.shape[1] < 2
        if self.empty:
            return
        ny = geo.shape[0]
        self.r0, self.r1 = max(0, int(rows[0]) - PAD), min(ny, int(rows[-1]) + 1 + PAD)
        if self.r1 - self.r0 < 2:
            self.r0, self.r1 = max(0, self.r1 - 2), min(ny, max(self.r1, 2))
        self.P = np.ascontiguousarray(xyz[self.r0:self.r1])
        self.shape = self.P.shape[:2]
        self.time = geo.time[self.r0:self.r1]
        valid = np.isfinite(self.P).all(axis=-1)
        self.idx = np.flatnonzero(valid.ravel())
        if self.idx.size == 0:
            self.empty = True
            return
        self.tree = cKDTree(self.P.reshape(-1, 3)[self.idx])
        # 相鄰掃描的間距 (各列取跨整條掃描的中位數) -> 缺掃描的位置
        with np.errstate(invalid='ignore', all='ignore'):
            step = np.linalg.norm(self.P[1:] - self.P[:-1], axis=-1)
            row_step = np.nanmedian(np.where(np.isfinite(step), step, np.nan), axis=1)
            med = np.nanmedian(row_step)
            self.ratio = row_step / med
            cols = np.linalg.norm(self.P[:, 1:] - self.P[:, :-1], axis=-1)
            col_step, col_max = np.nanmedian(cols), np.nanmax(cols)
        self.gap = ~(self.ratio <= GAP)             # NaN 也算缺
        # KD-tree 查詢的距離上限：掃描帶內的像素到最近取樣點不會超過「半格沿軌 + 半格沿掃描」。
        # 沒有上限時，離掃描帶很遠的像素 (全球圖的大部分) 幾乎要走訪整棵樹，慢上千倍
        self.reach = FAR * 0.5 * (GAP * (med if np.isfinite(med) else 0.) + (col_max if np.isfinite(col_max) else 0.))
        if not self.reach > 0.:
            self.reach = np.inf
        steps = [v for v in (med, col_step) if np.isfinite(v) and v > 0]
        self.step_km = EARTH_KM * min(steps) if steps else 10.      # 較細方向的取樣間距
        self._steps = None

    def steps(self):
        """各取樣點沿掃描、沿軌方向到相鄰取樣點的距離 (單位球上的弦長；(列, 行) float64)：前後兩個間距的平均。
        缺掃描處的間距不算；算不出來的位置用整體的中位數。"""
        if self._steps is None:
            with np.errstate(invalid='ignore', all='ignore'):
                dc = np.linalg.norm(self.P[:, 1:] - self.P[:, :-1], axis=-1)
                dr = np.linalg.norm(self.P[1:] - self.P[:-1], axis=-1)
                dr[self.gap] = np.nan
                self._steps = (_centered(dc, 1), _centered(dr, 0))
        return self._steps

    # ---- 小數座標 ----------------------------------------------------------
    def locate(self, lon, lat, workers=1):
        """輸出像素的經緯度 -> (r, c)：相對於 self.r0 的小數列號與行號 (float64，攤平；算不出來為 NaN)。
        離最近取樣點超過約半格 (FAR 倍) 的像素在掃描帶外，為 NaN；其餘是否在掃描帶內由 inside() 判斷
        (邊緣外半格以內的值由邊緣外插)。缺掃描處尚未處理 (snap())。"""
        p = unit(lat, lon).reshape(-1, 3)
        n = p.shape[0]
        r, c = np.full(n, np.nan), np.full(n, np.nan)
        good = np.flatnonzero(np.isfinite(p).all(axis=1))
        if self.empty or good.size == 0:
            return r, c
        dist, k = self.tree.query(p[good], workers=workers, distance_upper_bound=self.reach)
        near = k < self.idx.size                    # 超過距離上限的像素沒有結果 (掃描帶外)
        good, dist, k = good[near], dist[near], k[near]
        if good.size == 0:
            return r, c
        q = p[good]
        i, j = np.divmod(self.idx[k], self.shape[1])
        er, ec = self._diff(i, j, 0), self._diff(i, j, 1)
        dr, dc = _solve(q - self.P[i, j], er, ec)
        ra, ca = i + dr, j + dc
        rn, cn = self._newton(q, ra, ca)
        # 所在格子有無效的角點時牛頓法得不到結果：最近點半格以內仍以局部解為準 (畫成該取樣點)
        keep = ~(np.isfinite(rn) & np.isfinite(cn)) & (np.abs(dr) <= 0.5) & (np.abs(dc) <= 0.5)
        rn[keep], cn[keep] = ra[keep], ca[keep]
        # 掃描帶內的像素到最近取樣點的距離不會超過「半格沿軌 + 半格沿掃描」
        with np.errstate(invalid='ignore'):
            far = ~(dist <= FAR * 0.5 * (np.linalg.norm(er, axis=1) + np.linalg.norm(ec, axis=1)))
        rn[far], cn[far] = np.nan, np.nan
        r[good], c[good] = rn, cn
        return r, c

    def inside(self, r, c):
        """小數座標是否在掃描帶內 (邊緣各多半個取樣點)。"""
        ny, nx = self.shape
        with np.errstate(invalid='ignore'):
            return (r >= -0.5) & (r <= ny - 0.5) & (c >= -0.5) & (c <= nx - 0.5)

    def _diff(self, i, j, axis):
        """最近點 (i, j) 沿 axis 方向每一格的位移向量：前後差分的平均，缺一邊時用另一邊。"""
        n = self.shape[axis]
        here = self.P[i, j]
        if axis == 0:
            fwd = self.P[np.minimum(i + 1, n - 1), j] - here
            bwd = here - self.P[np.maximum(i - 1, 0), j]
            hasf, hasb = i + 1 < n, i - 1 >= 0
        else:
            fwd = self.P[i, np.minimum(j + 1, n - 1)] - here
            bwd = here - self.P[i, np.maximum(j - 1, 0)]
            hasf, hasb = j + 1 < n, j - 1 >= 0
        hasf = hasf & np.isfinite(fwd).all(axis=1)
        hasb = hasb & np.isfinite(bwd).all(axis=1)
        fwd = np.where(hasf[:, None], fwd, 0.)
        bwd = np.where(hasb[:, None], bwd, 0.)
        cnt = hasf.astype(np.float64) + hasb
        with np.errstate(invalid='ignore', divide='ignore'):
            return (fwd + bwd) / cnt[:, None]

    def _newton(self, q, r, c, iters=2):
        ny, nx = self.shape
        r, c = r.copy(), c.copy()
        for _ in range(iters):
            with np.errstate(invalid='ignore'):
                i0 = np.clip(np.floor(np.nan_to_num(r)), 0, ny - 2).astype(np.intp)
                j0 = np.clip(np.floor(np.nan_to_num(c)), 0, nx - 2).astype(np.intp)
            fr, fc = (r - i0)[:, None], (c - j0)[:, None]
            p00, p01 = self.P[i0, j0], self.P[i0, j0 + 1]
            p10, p11 = self.P[i0 + 1, j0], self.P[i0 + 1, j0 + 1]
            g = (1. - fr) * ((1. - fc) * p00 + fc * p01) + fr * ((1. - fc) * p10 + fc * p11)
            er = (1. - fc) * (p10 - p00) + fc * (p11 - p01)
            ec = (1. - fr) * (p01 - p00) + fr * (p11 - p10)
            dr, dc = _solve(q - g, er, ec)
            r, c = r + dr, c + dc
        return r, c

    def snap(self, r):
        """缺掃描處不跨越內插：離某一條掃描半個正常間距以內的畫成該掃描，其餘為無資料 (NaN)。"""
        if not self.gap.any():
            return r
        ny = self.shape[0]
        with np.errstate(invalid='ignore'):
            i0 = np.clip(np.floor(np.nan_to_num(r)), 0, ny - 2).astype(np.intp)
            f = r - i0
            hit = self.gap[i0] & (f > 0.) & (f < 1.)
            lim = 0.5 / np.where(np.isfinite(self.ratio[i0]), self.ratio[i0], np.inf)
            lower, upper = hit & (f <= lim), hit & (1. - f <= lim)
        out = np.where(hit, np.nan, r)
        out[lower] = i0[lower]
        out[upper] = i0[upper] + 1.
        return out


def _centered(d, axis):
    """相鄰取樣點的間距 d (沿 axis 少一個) -> 各取樣點的間距：前後兩個的平均 (缺一邊用另一邊，都缺用中位數)。"""
    d = np.moveaxis(d, axis, 0)
    pad = np.full((1,) + d.shape[1:], np.nan)
    a, b = np.concatenate([d, pad]), np.concatenate([pad, d])
    both = np.isfinite(a) & np.isfinite(b)
    out = np.where(both, 0.5 * (a + b), np.where(np.isfinite(a), a, b))
    med = np.nanmedian(d) if np.isfinite(d).any() else 1e-3
    out = np.where(np.isfinite(out) & (out > 0.), out, med)
    return np.ascontiguousarray(np.moveaxis(out, 0, axis))


def _solve(d, er, ec):
    """最小平方解 d = dr * er + dc * ec (3-D 向量，2 個未知數)。"""
    a11, a12, a22 = _dot(er, er), _dot(er, ec), _dot(ec, ec)
    b1, b2 = _dot(d, er), _dot(d, ec)
    with np.errstate(invalid='ignore', divide='ignore'):
        det = a11 * a22 - a12 * a12
        return (a22 * b1 - a12 * b2) / det, (a11 * b2 - a12 * b1) / det


# ---- 內插 ------------------------------------------------------------------
def _slope(sa, sb, ha, hb):
    """節點處的斜率：兩側斜率的加權平均 (中央差分)，限制在兩側較小斜率的 2 倍以內，
    兩側斜率反號 (極值) 時為 0 —— 內插結果不會超出相鄰取樣值的範圍。"""
    with np.errstate(invalid='ignore'):
        central = (hb * sa + ha * sb) / (ha + hb)
        lim = 2. * np.minimum(np.abs(sa), np.abs(sb))
        m = np.sign(central) * np.minimum(np.abs(central), lim)
        return np.where(sa * sb <= 0., 0., m)


def _hermite(y0, y1, y2, y3, x0, x1, x2, x3, x=0.):
    """節點 x0 < x1 <= x < x2 < x3 上的值 y0..y3，求 x 處的值 (限制斜率的三次 Hermite)。
    y0 / y3 無效時該側用單邊斜率；y1 / y2 無效時結果為 NaN。"""
    with np.errstate(invalid='ignore', divide='ignore'):
        h = np.maximum(x2 - x1, EPS)
        t = np.clip((x - x1) / h, 0., 1.)
        s1 = (y2 - y1) / h
        hl, hr = np.maximum(x1 - x0, EPS), np.maximum(x3 - x2, EPS)
        m1 = np.where(np.isfinite(y0), _slope((y1 - y0) / hl, s1, hl, h), s1)
        m2 = np.where(np.isfinite(y3), _slope(s1, (y3 - y2) / hr, h, hr), s1)
        t2 = t * t
        t3 = t2 * t
        return ((2. * t3 - 3. * t2 + 1.) * y1 + (t3 - 2. * t2 + t) * h * m1 +
                (3. * t2 - 2. * t3) * y2 + (t3 - t2) * h * m2)


def _lines(tb, r, c, ok, method):
    """一組掃描在各輸出像素附近的掃描線：回傳 (t, v)，形狀 (N, K)。
    t = 掃描線相對於像素的沿軌位置 (以掃描間距為單位，<= 0 在像素之前)；v = 該掃描線在像素橫向位置的值
    (沿掃描方向內插)。K：nearest 1、bilinear 2、cubic 4。"""
    ny, nx = tb.shape
    rc, cc = np.clip(np.nan_to_num(r), 0., ny - 1.), np.clip(np.nan_to_num(c), 0., nx - 1.)
    if method == 'nearest':
        i, j = np.rint(rc).astype(np.intp), np.rint(cc).astype(np.intp)
        return (i - rc)[:, None], np.where(ok, tb[i, j], np.nan)[:, None]
    i0 = np.clip(np.floor(rc), 0, ny - 2).astype(np.intp)
    j0 = np.clip(np.floor(cc), 0, nx - 2).astype(np.intp)
    g = (cc - j0)[:, None]
    offs = np.array((0, 1) if method == 'bilinear' else (-1, 0, 1, 2), dtype=np.intp)
    rows = i0[:, None] + offs
    rvalid = (rows >= 0) & (rows < ny) & ok[:, None]
    rows = np.clip(rows, 0, ny - 1)
    a, b = tb[rows, j0[:, None]], tb[rows, j0[:, None] + 1]
    if method == 'bilinear':
        v = a + (b - a) * g
    else:
        jl, jr = j0[:, None] - 1, j0[:, None] + 2
        yl = np.where(jl >= 0, tb[rows, np.maximum(jl, 0)], np.nan)
        yr = np.where(jr < nx, tb[rows, np.minimum(jr, nx - 1)], np.nan)
        v = _hermite(yl, a, b, yr, -1., 0., 1., 2., g)
    v = np.where(rvalid, v, np.nan)
    return rows - rc[:, None], v


def _pick(key, k):
    """每列 key 最小的 k 個的索引 (N, k)。"""
    if key.shape[1] <= k:
        return np.argsort(key, axis=1)
    return np.argsort(key, axis=1)[:, :k]


def _combine(t, v, method):
    """各掃描線 (t, v) (N, K) -> 像素值 (N,)：沿軌方向內插 (多組掃描時取聯集中最近的掃描線)。"""
    valid = np.isfinite(v)
    if method == 'nearest':
        k = np.argmin(np.where(valid, np.abs(t), np.inf), axis=1)[:, None]
        return np.take_along_axis(v, k, 1)[:, 0]
    need = 1 if method == 'bilinear' else 2
    kb = _pick(np.where(valid & (t <= 0.), -t, np.inf), need)      # 像素之前 (含剛好在線上)：由近到遠
    ka = _pick(np.where(valid & (t > 0.), t, np.inf), need)        # 像素之後

    def take(k, before):
        x, y = np.take_along_axis(t, k, 1)[:, 0], np.take_along_axis(v, k, 1)[:, 0]
        with np.errstate(invalid='ignore'):
            good = np.isfinite(y) & ((x <= 0.) if before else (x > 0.))
        return np.where(good, x, np.nan), np.where(good, y, np.nan)

    x1, y1 = take(kb[:, :1], True)
    x2, y2 = take(ka[:, :1], False)
    if method == 'bilinear':
        with np.errstate(invalid='ignore', divide='ignore'):
            w = np.clip(-x1 / np.maximum(x2 - x1, EPS), 0., 1.)
            out = y1 + (y2 - y1) * w
        # 剛好在掃描線上 (含掃描帶最後一條)：不需要另一側
        return np.where(np.isfinite(out), out, np.where(x1 == 0., y1, np.nan))
    x0, y0 = take(kb[:, 1:2], True) if kb.shape[1] > 1 else (np.full_like(x1, np.nan),) * 2
    x3, y3 = take(ka[:, 1:2], False) if ka.shape[1] > 1 else (np.full_like(x1, np.nan),) * 2
    x0 = np.where(np.isfinite(x0), x0, x1 - 1.)
    x3 = np.where(np.isfinite(x3), x3, x2 + 1.)
    out = _hermite(y0, y1, y2, y3, x0, x1, x2, x3)
    return np.where(np.isfinite(out), out, np.where(x1 == 0., y1, np.nan))


def interpolate(parts, method, engine=None):
    """parts = [(亮溫陣列, r, c, ok), ...] (一或多組掃描) -> 各像素的值 (N,) float32。
    周圍有無效值而算不出來的像素，改用較低階的方法。engine = 'numpy' 時不用 numba (測試用)。"""
    n = parts[0][1].size
    if _nb is not None and engine != 'numpy' and len(parts) <= 2:
        out = np.empty(n, dtype=np.float32)
        a, b = parts[0], parts[-1]
        _nb.interp(np.ascontiguousarray(a[0], np.float32), a[1], a[2], a[3],
                   np.ascontiguousarray(b[0], np.float32), b[1], b[2], b[3], len(parts), METHODS.index(method), out)
        return out
    out = np.full(n, np.nan, dtype=np.float32)
    order = METHODS[:METHODS.index(method) + 1][::-1]          # 例如 cubic -> bilinear -> nearest
    todo = np.flatnonzero(np.any([p[3] for p in parts], axis=0))
    for m in order:
        if todo.size == 0:
            break
        left = []
        for s in range(0, todo.size, CHUNK):
            sel = todo[s:s + CHUNK]
            tv = [_lines(tb, r[sel], c[sel], ok[sel], m) for tb, r, c, ok in parts]
            val = _combine(np.concatenate([x[0] for x in tv], axis=1), np.concatenate([x[1] for x in tv], axis=1), m)
            good = np.isfinite(val)
            out[sel[good]] = val[good]
            left.append(sel[~good])
        todo = np.concatenate(left) if left else todo[:0]
    return out


# ---- 高斯距離加權 ----------------------------------------------------------
def _gauss_np(P, sx, sy, tb, r, c, ok, sigma, factor, nsw, num, den, near):
    """_swath_nb.gauss 的 numpy 版 (沒有 numba 時使用；算法相同)。"""
    ny, nx = tb.shape
    todo = np.flatnonzero(ok)
    for s in range(0, todo.size, CHUNK // 4):
        sel = todo[s:s + CHUNK // 4]
        rr, cc = r[sel], c[sel]
        rc, ccl = np.clip(rr, 0., ny - 1.), np.clip(cc, 0., nx - 1.)
        i_n, j_n = np.rint(rc).astype(np.intp), np.rint(ccl).astype(np.intp)
        i0 = np.clip(np.floor(rc), 0, ny - 2).astype(np.intp)
        j0 = np.clip(np.floor(ccl), 0, nx - 2).astype(np.intp)
        fr, fc = (rr - i0)[:, None], (cc - j0)[:, None]
        q = ((1. - fr) * ((1. - fc) * P[i0, j0] + fc * P[i0, j0 + 1]) +
             fr * ((1. - fc) * P[i0 + 1, j0] + fc * P[i0 + 1, j0 + 1]))
        bad = ~np.isfinite(q).all(axis=1)
        q[bad] = P[i_n[bad], j_n[bad]]
        sc, sr = sx[i_n, j_n], sy[i_n, j_n]
        coarse = np.maximum(sc, sr)
        sig = np.full(sel.size, sigma) if sigma > 0. else factor * np.maximum(sc, sr / nsw)
        reach = np.maximum(GAUSS_CUT * sig, GAUSS_FILL * coarse)
        nr = np.minimum(GAUSS_MAXW, (reach / sr).astype(np.intp) + 1)
        nc = np.minimum(GAUSS_MAXW, (reach / sc).astype(np.intp) + 1)
        lim = (GAUSS_FILL * coarse) ** 2
        sw, sv, hit = np.zeros(sel.size), np.zeros(sel.size), np.zeros(sel.size, bool)
        for di in range(-int(nr.max()), int(nr.max()) + 1):
            for dj in range(-int(nc.max()), int(nc.max()) + 1):
                i, j = i_n + di, j_n + dj
                use = (abs(di) <= nr) & (abs(dj) <= nc) & (i >= 0) & (i < ny) & (j >= 0) & (j < nx)
                i, j = np.clip(i, 0, ny - 1), np.clip(j, 0, nx - 1)
                val = tb[i, j]
                with np.errstate(invalid='ignore'):
                    d2 = ((P[i, j] - q) ** 2).sum(axis=1)
                    use &= np.isfinite(val) & (d2 <= reach * reach)
                    w = np.where(use, np.exp(-np.minimum(d2 / (sig * sig), 700.)), 0.)
                    hit |= use & (d2 <= lim)
                sw += w
                sv += w * np.where(use, val, 0.)
        num[sel] += sv
        den[sel] += sw
        near[sel] |= hit


def gauss(parts, sigma, factor, engine=None):
    """parts = [(P, sx, sy, 亮溫, r, c, ok), ...] (一或多組掃描；各組的取樣點一起加權) -> 各像素的值 (N,) float32。
    sigma：固定的 σ (單位球上的弦長)；<= 0 時 σ = factor x 該處較粗方向的取樣間距。"""
    n = parts[0][4].size
    num, den, near = np.zeros(n), np.zeros(n), np.zeros(n, bool)
    for P, sx, sy, tb, r, c, ok in parts:
        if _nb is not None and engine != 'numpy':
            _nb.gauss(P, sx, sy, np.ascontiguousarray(tb, np.float32), r, c, ok, float(sigma), float(factor),
                      len(parts), GAUSS_CUT, GAUSS_FILL, GAUSS_MAXW, num, den, near)
        else:
            _gauss_np(P, sx, sy, tb, r, c, ok, float(sigma), float(factor), len(parts), num, den, near)
    with np.errstate(invalid='ignore', divide='ignore'):
        return np.where(near & (den > 0.), num / den, np.nan).astype(np.float32)


# ---- 輸出格點的取樣器 -------------------------------------------------------
class Resampler:
    """一張輸出格點 (grid.OutputGrid) 的取樣器；各組幾何的 Locator 只建一次。"""

    def __init__(self, grid, method='cubic', sigma_km=None, factor=GAUSS_FACTOR):
        """sigma_km、factor：gauss 的 σ —— 固定的公里數，或 (sigma_km = None) 較粗方向取樣間距的 factor 倍。"""
        if method not in ALL_METHODS:
            raise ValueError(tr('interpolation 須為 %s', 'interpolation must be %s') % ' / '.join(ALL_METHODS))
        self.grid = grid
        self.method = method
        self.sigma = float(sigma_km) / EARTH_KM if sigma_km else 0.
        self.factor = float(factor)
        lon, lat = grid.lonlat_sparse(8)
        p = unit(lat, lon).reshape(-1, 3)
        p = p[np.isfinite(p).all(axis=1)]
        mean = p.mean(axis=0)
        norm = np.linalg.norm(mean)
        if norm < 1e-6:                                 # 全球圖：中心沒有意義，取整個球
            self.center, self.radius = np.array([1., 0., 0.]), 2.1
        else:
            self.center = mean / norm
            self.radius = min(2.1, float(np.linalg.norm(p - self.center, axis=1).max()) + MARGIN_KM / EARTH_KM)
        self._loc = {}

    def locator(self, geo):
        with parallel.LOCK:
            if geo.key not in self._loc:
                self._loc[geo.key] = Locator(geo, self.center, self.radius)
            return self._loc[geo.key]

    def coords(self, geo, i0, i1, cache=None, workers=1):
        """第 i0~i1 列各像素在該組幾何中的 (Locator, r, c, ok)；cache (dict) 讓同一條內的多個頻道共用。

        小數座標是輸出位置的平滑函數，所以只在每隔 S 個像素的格子點上實際反算 (KD-tree + 局部求解)，
        其餘像素由格子點雙線性內插；S 取輸出像素不超過半個取樣間距的最大值 (1~LATTICE)，
        輸出像素比取樣粗時 S = 1 (逐像素反算)。格子點取全圖固定的位置，分條處理的結果與整張相同。
        四個角有算不出來的值、或座標變化異常大的格子 (掃描帶自己重疊處)，格內像素改為逐點反算。"""
        key = (geo.key, i0, i1)
        if cache is not None and key in cache:
            return cache[key]
        loc = self.locator(geo)
        ny, nx = self.grid.shape
        n = (i1 - i0) * nx
        if loc.empty:
            res = (loc, np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, bool))
        else:
            step = int(np.clip(0.5 * loc.step_km / self.grid.pixel_km, 1, LATTICE))
            if step == 1:
                lon, lat = self.grid.lonlat_rows(i0, i1)
                r, c = loc.locate(lon, lat, workers)
            else:
                r, c = self._lattice(loc, i0, i1, step, workers)
            r = loc.snap(r)
            res = (loc, r, c, loc.inside(r, c))
        if cache is not None:
            cache[key] = res
        return res

    def _lattice(self, loc, i0, i1, step, workers):
        ny, nx = self.grid.shape
        rows = np.unique(np.append(np.arange(0, ny, step), ny - 1))
        cols = np.unique(np.append(np.arange(0, nx, step), nx - 1))
        a = max(0, np.searchsorted(rows, i0, 'right') - 1)
        b = min(rows.size - 1, np.searchsorted(rows, i1 - 1, 'left'))
        rows = rows[a:max(b, a + 1) + 1] if rows.size > 1 else rows
        lon, lat = self.grid.lonlat_lattice(rows, cols)
        shape = (rows.size, cols.size)
        rl, cl = (v.reshape(shape) for v in loc.locate(lon, lat, workers))

        def weights(nodes, pos):
            """各像素位置 pos 所在的格子 (左節點索引) 與到右節點的權重。"""
            if nodes.size < 2:
                return np.zeros(pos.size, np.intp), np.zeros(pos.size)
            k = np.clip(np.searchsorted(nodes, pos, 'right') - 1, 0, nodes.size - 2)
            return k, (pos - nodes[k]) / (nodes[k + 1] - nodes[k])

        ky, wy = weights(rows, np.arange(i0, i1))
        kx, wx = weights(cols, np.arange(nx))
        ky1, kx1 = np.minimum(ky + 1, rows.size - 1), np.minimum(kx + 1, cols.size - 1)
        iy, ix = np.ix_(ky, kx)
        iy1, ix1 = np.ix_(ky1, kx1)
        wy, wx = wy[:, None], wx[None, :]

        def blend(v):
            return ((1. - wy) * ((1. - wx) * v[iy, ix] + wx * v[iy, ix1]) +
                    wy * ((1. - wx) * v[iy1, ix] + wx * v[iy1, ix1]))

        r, c = blend(rl), blend(cl)
        # 不規則的格子：角點算不出來，或四角的座標相差太多 (正常時一格內的變化不到 1)
        with np.errstate(invalid='ignore'):
            spread = np.zeros(shape, bool)
            for v in (rl, cl):
                quad = np.stack([v[:-1, :-1], v[:-1, 1:], v[1:, :-1], v[1:, 1:]]) if min(shape) > 1 else v[None]
                bad = ~np.isfinite(quad).all(axis=0) | ((quad.max(axis=0) - quad.min(axis=0)) > 2.)
                spread[:bad.shape[0], :bad.shape[1]] |= bad
        if spread.any():
            pr, pc = np.nonzero(spread[iy, ix])
            if pr.size:
                lon, lat = self.grid.lonlat_points(pr + i0, pc)
                er, ec = loc.locate(lon, lat, workers)
                r[pr, pc], c[pr, pc] = er, ec
        return r.ravel(), c.ravel()

    def sample(self, swath, i0, i1, cache=None, workers=1, method=None):
        """swath = [(Geo, 亮溫), ...] -> 第 i0~i1 列的值 (列, 行) float32，無資料為 NaN。"""
        method = method or self.method
        parts = []
        for geo, tb in swath:
            loc, r, c, ok = self.coords(geo, i0, i1, cache, workers)
            # 範圍內整組沒有資料的掃描不算一組 (AMSR-E 的 89 GHz A 掃描 2004/11 起故障，只剩 B 掃描：
            # 算進去的話 gauss 會把取樣間距當成 A、B 合併後的一半)
            if not loc.empty and ok.any() and self._alive(loc, tb):
                parts.append((loc, tb[loc.r0:loc.r1], r, c, ok))
        shape = (i1 - i0, self.grid.shape[1])
        if not parts:
            return np.full(shape, np.nan, dtype=np.float32)
        if method == 'gauss':
            return gauss([(loc.P,) + loc.steps() + (tb, r, c, ok) for loc, tb, r, c, ok in parts],
                         self.sigma, self.factor).reshape(shape)
        return interpolate([p[1:] for p in parts], method).reshape(shape)

    def _alive(self, loc, tb):
        """這組掃描在地圖範圍內是否有任何有效的值 (每個陣列只檢查一次)。"""
        key = (id(loc), id(tb))
        with parallel.LOCK:
            cache = self.__dict__.setdefault('_alive_cache', {})
            if key not in cache:
                if len(cache) > 64:
                    cache.clear()
                cache[key] = (tb, bool(np.isfinite(tb[loc.r0:loc.r1]).any()))      # 留住 tb，id 才不會被重複使用
            return cache[key][1]

    def coverage(self, geo, step=8):
        """這組幾何在地圖範圍內的涵蓋情形 (每隔 step 個像素取樣)：
        回傳 (涵蓋比例 0~1, 通過時間 (1970 起算的秒數；範圍內各像素掃描時間的中位數，不明為 None))。"""
        loc = self.locator(geo)
        if loc.empty:
            return 0., None
        lon, lat = self.grid.lonlat_sparse(step)
        r, c = loc.locate(lon, lat, workers=-1)
        r = loc.snap(r)
        ok = loc.inside(r, c)
        if not ok.any():
            return 0., None
        rows = np.clip(np.rint(r[ok]), 0, loc.shape[0] - 1).astype(np.intp)
        t = loc.time[rows]
        t = t[np.isfinite(t)]
        return float(ok.mean()), (float(np.median(t)) if t.size else None)
