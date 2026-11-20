"""輸出格點與取樣

做法：每個輸出像素的經緯度 -> 導航 (nav.to_rowcol) -> 原始影像的行列號 -> 直接取值。
  * 只讀取涵蓋地圖範圍的視窗 (由輸出格點反算，邊界只多留幾個像素)，不需要原始經緯度、也不需要 KD-tree。
  * 輸出比原始解析度粗時，先做 factor x factor 區塊平均 (減少雜訊與鋸齒)；區塊起點對齊 factor 的倍數，
    所以不同解析度的波段 (例如 0.5 km 紅光的 2x2 區塊與 1 km 藍光像素) 自然對齊。
  * 輸出像素對應的原始像素存成整數索引 (最近的像素中心)：floor((r + 0.5) / f) = floor(floor(r + 0.5) / f)，
    與浮點行列號的結果完全相同，但記憶體只要一半。經緯度與行列號都分段計算 (全圓盤等大圖不需整張的經緯度陣列)。
  * 雙線性內插 (interp='bilinear'，只用於 factor = 1)：另存相對最近像素中心的偏移 (-0.5~0.5，量化成 int8，
    精度 1/254 像素)，最近像素索引不變。輸出像素比原始像素細 (放大) 時，最近鄰會不規則地重複部分原始像素，
    畫面出現階梯與拼貼感；雙線性則平順過渡。輸出與原始像素一一對應 (全圓盤原始解析度) 時偏移為 0，結果與最近鄰相同。
"""
import numpy as np

from ..core import parallel

try:                        # 區塊平均與取值的 numba 核心 (nogil，可多執行緒)；沒有 numba 時用 numpy
    from . import _resample_nb as _nb
except ImportError:
    _nb = None

MISSING = -(1 << 30)        # 看不到 (太空、地球背面、地圖範圍外) 的索引
CHUNK = 1 << 21             # 分段計算的像素數
STRIP = 1 << 22             # 節省記憶體模式：每一條的輸出像素數 (約 400 萬)
MIN_STRIP = 1 << 18         # 平行處理時每一條的最少像素數


def strip_pixels(size, low_memory, workers):
    """每一條的像素數：單核心時同以往 (節省記憶體 STRIP，否則整張一條)；
    多核心時切成約每個執行緒 4 條 (不超過 STRIP)，讓各執行緒的工作量平均。"""
    if workers <= 1:
        return STRIP if low_memory else None
    return min(STRIP, max(MIN_STRIP, -(-size // (4 * workers))))
COARSE = 16                 # 節省記憶體模式：估計讀取範圍時，每隔幾個輸出像素取樣一次


class OutputGrid:
    """輸出地圖 (satima/maps.py) 的像素中心格點；第一列為北邊 (imshow 請用 origin='upper')。"""

    def __init__(self, gmap, shape, pixel_km):
        ny, nx = int(shape[0]), int(shape[1])
        self.map = gmap
        self.shape = (ny, nx)
        self.size = ny * nx
        self.pixel_km = float(pixel_km)
        x0, x1, y0, y1 = gmap.extent
        dx, dy = (x1 - x0) / nx, (y1 - y0) / ny
        self.x = x0 + (np.arange(nx) + 0.5) * dx
        self.y = y1 - (np.arange(ny) + 0.5) * dy
        if gmap.separable:
            # 麥卡托、等距圓柱：經度只與 x、緯度只與 y 有關，反算一列與一行即可
            self._lon1, _ = gmap.inverse(self.x, np.full(nx, 0.5 * (y0 + y1)))
            _, self._lat1 = gmap.inverse(np.full(ny, 0.5 * (x0 + x1)), self.y)
        self._lonlat = None
        self._rowcol = {}

    def _chunks(self):
        step = max(1, CHUNK // self.shape[1])
        for i0 in range(0, self.shape[0], step):
            yield i0, min(self.shape[0], i0 + step)

    def strips(self, max_pixels=None):
        """逐條產生 (i0, i1, 子格點)；max_pixels 為 None 或整張不超過時只有一條 (本身，沿用快取)。
        子格點的行列號快取只在該條存在期間保留 (請勿保存整個列表)。"""
        for i0, i1 in self.strip_ranges(max_pixels):
            yield i0, i1, self.strip(i0, i1)

    def strip_ranges(self, max_pixels=None):
        """各條的列範圍 [(i0, i1), ...]。"""
        if not max_pixels or self.size <= max_pixels:
            return [(0, self.shape[0])]
        step = max(1, int(max_pixels) // self.shape[1])
        return [(i0, min(self.shape[0], i0 + step)) for i0 in range(0, self.shape[0], step)]

    def strip(self, i0, i1):
        """第 i0~i1 列的子格點 (整張時回傳本身)。平行處理時在各工作內才建立，處理完即釋放。"""
        if i0 == 0 and i1 == self.shape[0]:
            return self
        return self.subgrid(np.arange(i0, i1))

    def subgrid(self, rows, cols=None):
        """取出部分列 (與行) 的格點 (像素中心不變)。"""
        g = object.__new__(OutputGrid)
        cols = np.arange(self.shape[1]) if cols is None else np.asarray(cols)
        rows = np.asarray(rows)
        g.map, g.pixel_km = self.map, self.pixel_km
        g.shape = (rows.size, cols.size)
        g.size = rows.size * cols.size
        g.x, g.y = self.x[cols], self.y[rows]
        if self.map.separable:
            g._lon1, g._lat1 = self._lon1[cols], self._lat1[rows]
        g._lonlat = None if self._lonlat is None else (self._lonlat[0][np.ix_(rows, cols)],
                                                       self._lonlat[1][np.ix_(rows, cols)])
        g._rowcol = {}
        return g

    def lonlat_rows(self, i0, i1):
        """第 i0 ~ i1-1 列的經緯度 (看不到的位置為 NaN)。"""
        if self._lonlat is not None:
            return self._lonlat[0][i0:i1], self._lonlat[1][i0:i1]
        if self.map.separable:
            return np.meshgrid(self._lon1, self._lat1[i0:i1])
        x, y = np.meshgrid(self.x, self.y[i0:i1])
        return self.map.inverse(x, y)

    def lonlat_at(self, index):
        """攤平索引 (列 * 寬 + 行) 處的經緯度。"""
        i, j = np.divmod(np.asarray(index), self.shape[1])
        if self.map.separable:
            return self._lon1[j], self._lat1[i]
        geos = getattr(self.map, 'geos', None)
        if geos is not None and _nb is not None:        # 全圓盤：numba 反算 (與 PROJ geos 相同公式，快且可多執行緒)
            lon0, h, a, b, sweep = geos
            ox, oy = self.map.origin
            lon = np.empty(i.size)
            lat = np.empty(i.size)
            _nb.geos_inverse((self.x[j] - ox) / h, (self.y[i] - oy) / h, float(lon0), 1. + h / a, a / b,
                             sweep == 'x', lon, lat)
            return lon, lat
        return self.map.inverse(self.x[j], self.y[i])

    def _full(self):
        if self._lonlat is None:
            lon = np.empty(self.shape)
            lat = np.empty(self.shape)
            for i0, i1 in self._chunks():
                lon[i0:i1], lat[i0:i1] = self.lonlat_rows(i0, i1)
            self._lonlat = lon, lat
        return self._lonlat

    @property
    def lon(self):
        return self._full()[0]

    @property
    def lat(self):
        return self._full()[1]

    def rowcol(self, nav, frac=False):
        """輸出像素在該導航影像中最近的像素 (列, 行) 索引 (int32，看不到為 MISSING)；同一導航只算一次。
        frac=True 時回傳 (列, 行, 列偏移, 行偏移)，偏移為 int8 (偏移像素 x 254，見模組說明)。"""
        key = nav.key()
        hit = self._rowcol.get(key)
        if hit is None or (frac and hit[2] is None):
            rows = np.empty(self.shape, dtype=np.int32)
            cols = np.empty(self.shape, dtype=np.int32)
            fr = np.zeros(self.shape, dtype=np.int8) if frac else None
            fc = np.zeros(self.shape, dtype=np.int8) if frac else None
            if _same_geos(self.map, nav):
                self._geos_rowcol(nav, rows, cols, fr, fc)
            else:
                for i0, i1 in self._chunks():
                    r, c = nav.to_rowcol(*self.lonlat_rows(i0, i1))
                    rows[i0:i1], cols[i0:i1] = _nearest(r), _nearest(c)
                    if frac:
                        fr[i0:i1], fc[i0:i1] = _frac8(r, rows[i0:i1]), _frac8(c, cols[i0:i1])
                bad = (rows == MISSING) | (cols == MISSING)
                rows[bad] = cols[bad] = MISSING
            hit = self._rowcol[key] = rows, cols, fr, fc
        return hit if frac else hit[:2]

    def _geos_rowcol(self, nav, rows, cols, fr=None, fc=None):
        """全圓盤與該波段是同一個靜止衛星投影 (星下點、高度、橢球、sweep 相同)：掃描角與行列號為線性關係，
        不需經過經緯度 (結果與經緯度路徑相同，只是快很多)；看不到地球的位置用與 proj 相同的判斷式遮掉。"""
        _, h, a, b, sweep = self.map.geos
        ox, oy = self.map.origin
        sx = (self.x - ox) / h                  # 掃描角 (弧度)
        sy = (self.y - oy) / h
        r = (sy - nav.y0) / nav.dy
        c = (sx - nav.x0) / nav.dx
        r1, c1 = _nearest(r), _nearest(c)
        if fr is not None:                      # 偏移只與列 (或行) 有關
            fr[...] = _frac8(r, r1)[:, None]
            fc[...] = _frac8(c, c1)[None, :]
        if _nb is not None:
            _nb.geos_rowcol(np.tan(sx), np.tan(sy), r1, c1, 1. + h / a, a / b, sweep == 'x', MISSING, rows, cols)
            return
        for i0, i1 in self._chunks():
            vis = _geos_visible(sx[None, :], sy[i0:i1, None], h, a, b, sweep)
            ok = vis & (r1[i0:i1, None] != MISSING) & (c1[None, :] != MISSING)
            rows[i0:i1] = np.where(ok, r1[i0:i1, None], MISSING)
            cols[i0:i1] = np.where(ok, c1[None, :], MISSING)


def _same_geos(gmap, nav):
    from .nav import GeosNav
    g = getattr(gmap, 'geos', None)
    if g is None or not isinstance(nav, GeosNav):
        return False
    lon0, h, a, b, sweep = g
    return (abs((nav.sub_lon - lon0 + 180.) % 360. - 180.) < 1e-9 and abs(nav.height - h) < 1e-3 and
            abs(nav.a - a) < 1e-3 and abs(nav.b - b) < 1e-3 and nav.sweep == sweep)


def _geos_visible(x, y, h, a, b, sweep):
    """掃描角 (x, y) 的視線是否碰到地球 (proj geos 反算的判別式 >= 0)。"""
    rg = 1. + h / a
    tx, ty = np.tan(x), np.tan(y)
    if sweep == 'x':
        vz = ty
        vy = tx * np.hypot(1., vz)
    else:
        vy = tx
        vz = ty * np.hypot(1., vy)
    q = vy * vy + (vz * (a / b)) ** 2 + 1.
    return rg * rg >= q * (rg * rg - 1.)


def _nearest(v):
    """浮點行列號 -> 最近的像素索引 (NaN 或超出範圍為 MISSING)。"""
    v = np.asarray(v, dtype=np.float64)
    out = np.full(v.shape, MISSING, dtype=np.int32)
    with np.errstate(invalid='ignore'):
        ok = np.isfinite(v) & (np.abs(v) < -MISSING // 2)
    out[ok] = np.floor(v[ok] + 0.5)
    return out


FRAC = 254                  # 偏移的量化：偏移 (-0.5~0.5 像素) x 254 -> int8 (-127~127)


def _frac8(v, idx):
    """浮點行列號相對最近像素索引 idx 的偏移，量化成 int8 (看不到的位置為 0)。"""
    v = np.asarray(v, dtype=np.float64)
    q = np.zeros(v.shape, dtype=np.int8)
    ok = idx != MISSING
    q[ok] = np.clip(np.rint((v[ok] - idx[ok]) * FRAC), -127, 127)
    return q


def find_window(rows, cols, shape, align=1, margin=2):
    """涵蓋所有 (rows, cols) 索引的視窗 (r0, r1, c0, c1)，向外擴 margin 個區塊並對齊 align 的倍數；看不到時回傳 None。"""
    pad = margin * align
    if _nb is not None:
        rmin, rmax, cmin, cmax, n = _nb.bounds(np.asarray(rows).reshape(-1), np.asarray(cols).reshape(-1), MISSING)
        if n == 0:
            return None
    else:
        ok = rows != MISSING
        if not ok.any():
            return None
        r, c = rows[ok], cols[ok]
        rmin, rmax, cmin, cmax = r.min(), r.max(), c.min(), c.max()
    r0 = int(rmin) - pad
    r1 = int(rmax) + 1 + pad
    c0 = int(cmin) - pad
    c1 = int(cmax) + 1 + pad
    r0, c0 = max(0, r0 // align * align), max(0, c0 // align * align)
    r1 = min(shape[0], -(-r1 // align) * align)
    c1 = min(shape[1], -(-c1 // align) * align)
    if r1 <= r0 or c1 <= c0:
        return None
    return r0, r1, c0, c1


def blur_margin(blur, margin=2):
    """高斯模糊的半徑 (區塊數；scipy 預設 truncate = 4) 不超過視窗留邊，分條處理才與整張相同。"""
    if not blur:
        return margin
    return max(margin, int(4. * float(np.max(blur)) + 0.5))


class Source:
    """一個波段的原始資料：整個讀取範圍只讀一次，可對整張或分條的格點取樣。

    一般模式 (coarse = False)：讀取範圍由整張格點的行列號決定 (與逐張處理完全相同)。
    節省記憶體模式 (coarse = True)：每隔 COARSE 個輸出像素估計讀取範圍 (外加相鄰取樣點的最大間距)，
    不需要整張的行列號；某一條超出估計範圍時，該條直接另外讀取，所以結果一定與一般模式相同。
    align：讀取範圍對齊的倍數 (取樣時的 factor 須能整除)。
    interp：預設的取值方式 'nearest' / 'bilinear' (雙線性只用於 factor = 1；鄰點都在留邊 margin >= 1 之內)。"""

    def __init__(self, reader, scene, channel, grid, align=1, margin=2, coarse=False, interp='nearest'):
        self.reader, self.scene, self.channel = reader, scene, channel
        self.nav = reader.nav(scene, channel)
        self.align, self.margin = int(align), int(margin)
        self.interp = interp
        self._exact = None                  # 讀取範圍正好由這個格點決定時，取樣不必重算
        if coarse and grid.size > STRIP:
            self.window = _coarse_window(grid, self.nav, self.align, self.margin)
        else:
            rows, cols = grid.rowcol(self.nav, frac=_bilinear(interp, align))[:2]
            self.window = find_window(rows, cols, self.nav.shape, self.align, self.margin)
            self._exact = grid
        self.data = reader.read_cached(scene, channel, self.window) if self.window else None

    def sample(self, grid, factor=1, blur=None, interp=None):
        """取樣到 grid (整張或其中一條)：float32，無資料為 NaN。interp 省略時用建立時的設定。"""
        interp = self.interp if interp is None else interp
        if grid is self._exact:
            w = self.window
        else:
            rows, cols = grid.rowcol(self.nav, frac=_bilinear(interp, factor))[:2]
            w = find_window(rows, cols, self.nav.shape, self.align, self.margin)
        if w is None:
            return np.full(grid.shape, np.nan, dtype=np.float32)
        W = self.window
        if W is not None and W[0] <= w[0] and w[1] <= W[1] and W[2] <= w[2] and w[3] <= W[3]:
            data = self.data[w[0] - W[0]:w[1] - W[0], w[2] - W[2]:w[3] - W[2]]
        else:
            with parallel.LOCK:             # 讀取器的檔案存取不保證執行緒安全
                data = self.reader.read_cached(self.scene, self.channel, w)
        return to_grid(grid, self.nav, w, data, factor, blur, interp)


def _bilinear(interp, factor):
    return interp == 'bilinear' and int(factor) == 1


def _coarse_window(grid, nav, align, margin):
    """以稀疏取樣的輸出像素估計讀取範圍。"""
    ny, nx = grid.shape
    rows = np.unique(np.r_[np.arange(0, ny, COARSE), ny - 1])
    cols = np.unique(np.r_[np.arange(0, nx, COARSE), nx - 1])
    r, c = grid.subgrid(rows, cols).rowcol(nav)
    ok = r != MISSING
    if not ok.any():
        return None
    # 相鄰取樣點之間最多跨幾個原始像素 (兩點都看得到時)
    step = 0
    for a in (r, c):
        for d, okd in ((np.diff(a, axis=0), ok[1:] & ok[:-1]), (np.diff(a, axis=1), ok[:, 1:] & ok[:, :-1])):
            if okd.any():
                step = max(step, int(np.abs(d[okd]).max()))
    return find_window(r, c, nav.shape, align, margin + -(-(step + 1) // align))


def read_window(reader, scene, channel, grid, align=1, margin=2):
    """讀取涵蓋輸出格點的視窗。回傳 (nav, window, data)；地圖範圍看不到時 window、data 為 None。"""
    src = Source(reader, scene, channel, grid, align, margin)
    return src.nav, src.window, src.data


def to_grid(grid, nav, window, data, factor=1, blur=None, interp='nearest'):
    """把視窗資料取樣到輸出格點 (先做 factor x factor 區塊平均、可選高斯模糊)。
    interp='bilinear' 且 factor = 1 時雙線性內插，否則取最近的像素 (或區塊)。"""
    out = np.full(grid.shape, np.nan, dtype=np.float32)
    if window is None:
        return out
    r0, _, c0, _ = window
    if factor > 1:
        data = block_mean(data, factor)
    if blur:
        data = gaussian_nan(data, blur)
    if _bilinear(interp, factor):
        rows, cols, fr, fc = grid.rowcol(nav, frac=True)
        args = (np.asarray(data, dtype=np.float32), rows.reshape(-1), cols.reshape(-1), fr.reshape(-1),
                fc.reshape(-1), int(r0), int(c0), out.reshape(-1))
        (_nb.bilinear if _nb is not None else _bilinear_numpy)(*args)
        return out
    rows, cols = grid.rowcol(nav)
    if _nb is not None:
        _nb.gather(np.asarray(data, dtype=np.float32), rows.reshape(-1), cols.reshape(-1),
                   int(r0), int(c0), int(factor), out.reshape(-1))
        return out
    # 看不到的索引 (MISSING) 減去 r0 後仍為負數，自然落在範圍外
    br = (rows - r0) // factor
    bc = (cols - c0) // factor
    inside = (br >= 0) & (bc >= 0) & (br < data.shape[0]) & (bc < data.shape[1])
    out[inside] = data[br[inside], bc[inside]]
    return out


def _bilinear_numpy(data, rows, cols, fr, fc, r0, c0, out):
    """_resample_nb.bilinear 的 numpy 版 (沒有 numba 時)，規則相同。"""
    ny, nx = data.shape
    br, bc = rows - r0, cols - c0
    inside = np.flatnonzero((rows != MISSING) & (br >= 0) & (bc >= 0) & (br < ny) & (bc < nx))
    br, bc = br[inside].astype(np.int64), bc[inside].astype(np.int64)
    v = data[br, bc].astype(np.float64)

    def pair(b, q, n):
        t = np.where(q > 0, q / 254., np.where(q < 0, 1. + q / 254., 0.))
        a = np.where(q < 0, b - 1, b)
        bad = (q == 0) | (a < 0) | (a + 1 >= n)
        return np.where(bad, b, a), np.where(bad, 0., t), np.where(bad, b, a + 1)

    ra, t, rb = pair(br, fr[inside].astype(np.float64), ny)
    ca, u, cb = pair(bc, fc[inside].astype(np.float64), nx)
    a, b, c, d = (data[i, j].astype(np.float64) for i, j in ((ra, ca), (ra, cb), (rb, ca), (rb, cb)))
    with np.errstate(invalid='ignore'):
        mix = (1. - t) * ((1. - u) * a + u * b) + t * ((1. - u) * c + u * d)
    use = ((t != 0.) | (u != 0.)) & np.isfinite(a) & np.isfinite(b) & np.isfinite(c) & np.isfinite(d)
    out[inside] = np.where(use, mix, v)


def sample(reader, scene, channel, grid, factor=1, blur=None, margin=2):
    """讀取並取樣單一波段到輸出格點 (float32，無資料為 NaN)。"""
    return Source(reader, scene, channel, grid, factor, blur_margin(blur, margin)).sample(grid, factor, blur)


def block_mean(arr, k):
    """k x k 區塊平均，忽略 NaN (全為 NaN 的區塊為 NaN)。有 numba 時以單次掃描的核心計算 (nogil，可多執行緒)。"""
    if k <= 1:
        return arr
    if _nb is not None:
        return _nb.block_mean(np.asarray(arr, dtype=np.float32), int(k))
    ny, nx = arr.shape[0] // k * k, arr.shape[1] // k * k
    blocks = arr[:ny, :nx].reshape(ny // k, k, nx // k, k)
    valid = np.isfinite(blocks)
    total = np.where(valid, blocks, np.float32(0)).sum(axis=(1, 3), dtype=np.float32)
    count = valid.sum(axis=(1, 3), dtype=np.int16)
    with np.errstate(invalid='ignore', divide='ignore'):
        return (total / count).astype(np.float32)


def gaussian_nan(arr, sigma):
    """高斯模糊，忽略 NaN (正規化摺積)。"""
    from scipy.ndimage import gaussian_filter
    valid = np.isfinite(arr)
    num = gaussian_filter(np.where(valid, arr, 0.).astype(np.float32), sigma)
    den = gaussian_filter(valid.astype(np.float32), sigma)
    with np.errstate(invalid='ignore', divide='ignore'):
        out = num / den
    return np.where(valid, out, np.nan).astype(np.float32)
