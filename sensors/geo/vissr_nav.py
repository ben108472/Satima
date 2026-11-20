"""JMA VISSR 自旋掃描導航 (GMS-1~5、GOES-9) — GMS-VISSR 與 GMS-AREA 共用
(檔名沒有 '-'，不會被當成讀取器載入)

依 MSC Technical Note No.23 (JMA/MSC 1991) 的 MGIVSR (舊版 gms_nav.pyx 是它的 Cython 版、satpy 的
gms5_vissr_navigation 也是同一套)，每一點都完整計算 (不內插)：
  影像行號 / 像素號 (VISSR line / pixel，由 1 起算)
    -> 觀測時刻 t = t0 + (floor((行號-1)/感測器數) + 像素號 x 取樣角 / 2π) / (自旋速率 x 1440)   [MJD]
    -> 該時刻的衛星位置與姿態 (MGI100/MGI110：軌道、姿態預報表線性內插；歲差章動矩陣取區間起點)
    -> 自旋座標系 (自旋軸 SP、以太陽方向與太陽-地球夾角 β 定出的 SX、SY)
    -> 掃描角 (行步進角 x (行號-中心行號)、像素取樣角 x (像素號-中心像素號)) 與失準矩陣 ELMIS 得到視線
    -> 與地球橢球的交點 -> 經緯度 (大地緯度)
反算 (經緯度 -> 行列號)：直接解上面的旋轉，與正算互為精確反函數 (MGIVSR 的反算用 VMIS 小角度近似，
GMS-1 的失準約 0.3°，近似誤差可達 1 個可見光像素)；觀測時刻依算出的行列號反覆修正到收斂。
地球橢球同 JMA (赤道半徑 6378136 m、扁率 1/298.257)。
計算方式由 Satima-GEO.toml 的 [navigation] engine 選擇，三者演算法相同、結果一致 (差 < 1e-10 像素)：
  numba   預設；不需 C 編譯器，第一次使用時編譯並快取在 sensors/__pycache__ (satpy 的 GMS-5 導航也用 numba)
  cython  vissr_nav_cy.pyx (結構沿用舊版 gms_nav.pyx)，需先執行 python sensors/geo/build_vissr_nav_cy.py 編譯
  numpy   不需任何編譯，較慢
找不到 Cython 模組時改用 numba，沒有 numba 時改用 numpy。
另含地球邊緣修正 (MSC Technical Note No.26)：earth_edge_deviation()。
"""
import math

import numpy as np

from satima.core.i18n import tr
from satima.geo.nav import Nav

try:
    import numba
except ImportError:                 # 沒有 numba：用 numpy 計算
    numba = None

EA = 6378136.0                  # 地球赤道半徑 (m)
EF = 1.0 / 298.257              # 扁率
EB = EA * (1.0 - EF)            # 極半徑
_F2 = (1.0 - EF) ** 2
_EE = 2.0 * EF - EF * EF        # 第一離心率平方
_CHUNK = 1 << 18                # numpy 計算時一次的點數 (控制記憶體)
_ITER = 5                       # 反算時觀測時刻修正的最多次數 (通常 2~3 次收斂)


def mjd_to_datetime(mjd):
    import datetime
    return datetime.datetime(1858, 11, 17) + datetime.timedelta(days=float(mjd))


def _clean(mjd, *arrays):
    """去掉無效 (MJD <= 0) 與重複的預報點，依時間排序。"""
    mjd = np.asarray(mjd, dtype=np.float64)
    ok = np.flatnonzero(mjd > 0)
    mjd, first = np.unique(mjd[ok], return_index=True)
    sel = ok[first]
    if mjd.size < 2:
        raise ValueError(tr('導航預報資料不足 (有效點 %d 個)',
                            'not enough navigation prediction data (%d valid points)') % mjd.size)
    return (mjd,) + tuple(np.asarray(a, dtype=np.float64)[sel] for a in arrays)


# ==================================================================================================
#  numba 版 (逐點)
# ==================================================================================================
if numba is not None:
    njit, prange = numba.njit, numba.prange

    @njit(cache=True)
    def interval(x, xp):
        """MGI100 的區間 xp[i] <= x < xp[i+1] (超出範圍時用第一或最後一個區間外插)。"""
        n = xp.shape[0]
        for i in range(n - 2):
            if x < xp[i + 1]:
                return i
        return n - 2

    @njit(cache=True)
    def unit3(x, y, z):
        r = math.sqrt(x * x + y * y + z * z)
        return x / r, y / r, z / r

    @njit(cache=True)
    def frame(t, att_t, att, orb_t, orb, npa):
        """MGI100 + MGI110 + 自旋座標系：時刻 t -> 衛星位置 SAT、SP、SX、SY。"""
        i = interval(t, att_t)
        w = (t - att_t[i]) / (att_t[i + 1] - att_t[i])
        alp = att[i, 0] + w * (att[i + 1, 0] - att[i, 0])
        dlt = att[i, 1] + w * (att[i + 1, 1] - att[i, 1])
        bet = att[i, 2] + w * (att[i + 1, 2] - att[i, 2])
        j = interval(t, orb_t)
        v = (t - orb_t[j]) / (orb_t[j + 1] - orb_t[j])
        sx0 = orb[j, 0] + v * (orb[j + 1, 0] - orb[j, 0])
        sy0 = orb[j, 1] + v * (orb[j + 1, 1] - orb[j, 1])
        sz0 = orb[j, 2] + v * (orb[j + 1, 2] - orb[j, 2])
        gst = orb[j, 3] + v * (orb[j + 1, 3] - orb[j, 3])
        az = orb[j, 4] + v * (orb[j + 1, 4] - orb[j, 4])
        el = orb[j, 5] + v * (orb[j + 1, 5] - orb[j, 5])
        a1x, a1y, a1z = math.sin(dlt), -math.cos(dlt) * math.sin(alp), math.cos(dlt) * math.cos(alp)
        a2x = npa[j, 0, 0] * a1x + npa[j, 0, 1] * a1y + npa[j, 0, 2] * a1z
        a2y = npa[j, 1, 0] * a1x + npa[j, 1, 1] * a1y + npa[j, 1, 2] * a1z
        a2z = npa[j, 2, 0] * a1x + npa[j, 2, 1] * a1y + npa[j, 2, 2] * a1z
        cg, sg = math.cos(gst), math.sin(gst)
        spx, spy, spz = unit3(cg * a2x + sg * a2y, -sg * a2x + cg * a2y, a2z)
        ssx, ssy, ssz = math.cos(el) * math.cos(az), math.cos(el) * math.sin(az), math.sin(el)
        w1x, w1y, w1z = unit3(spy * ssz - spz * ssy, spz * ssx - spx * ssz, spx * ssy - spy * ssx)
        w2x, w2y, w2z = unit3(w1y * spz - w1z * spy, w1z * spx - w1x * spz, w1x * spy - w1y * spx)
        sb, cb = math.sin(bet), math.cos(bet)
        xx, xy, xz = unit3(w1x * sb + w2x * cb, w1y * sb + w2y * cb, w1z * sb + w2z * cb)
        yx, yy, yz = unit3(spy * xz - spz * xy, spz * xx - spx * xz, spx * xy - spy * xx)
        return sx0, sy0, sz0, spx, spy, spz, xx, xy, xz, yx, yy, yz

    @njit(cache=True)
    def obs_time(line, pixel, p):
        return p[5] + (math.floor((line - 1.) / p[4]) + pixel * p[1] / (2. * math.pi)) / (p[6] * 1440.)

    @njit(parallel=True, cache=True)
    def forward(line, pixel, p, e, att_t, att, orb_t, orb, npa, lon, lat):
        f2 = (1. - p[7]) ** 2
        ea = p[8]
        for k in prange(line.shape[0]):
            L, P = line[k], pixel[k]
            if not (math.isfinite(L) and math.isfinite(P)):
                lon[k], lat[k] = np.nan, np.nan
                continue
            satx, saty, satz, spx, spy, spz, xx, xy, xz, yx, yy, yz = frame(obs_time(L, P, p), att_t, att, orb_t, orb, npa)
            ang_p = p[0] * (L - p[2])
            ang_q = p[1] * (P - p[3])
            pc, ps, qc, qs = math.cos(ang_p), math.sin(ang_p), math.cos(ang_q), math.sin(ang_q)
            w1x = e[0, 0] * pc + e[0, 2] * ps
            w1y = e[1, 0] * pc + e[1, 2] * ps
            w1z = e[2, 0] * pc + e[2, 2] * ps
            w2x = qc * w1x - qs * w1y
            w2y = qs * w1x + qc * w1y
            lx, ly, lz = unit3(xx * w2x + yx * w2y + spx * w1z, xy * w2x + yy * w2y + spy * w1z,
                               xz * w2x + yz * w2y + spz * w1z)
            a = f2 * (lx * lx + ly * ly) + lz * lz
            b = f2 * (satx * lx + saty * ly) + satz * lz
            c = f2 * (satx * satx + saty * saty - ea * ea) + satz * satz
            disc = b * b - a * c
            if disc < 0.:
                lon[k], lat[k] = np.nan, np.nan
                continue
            d = (-b - math.sqrt(disc)) / a
            px, py, pz = satx + d * lx, saty + d * ly, satz + d * lz
            lat[k] = math.degrees(math.atan2(pz, f2 * math.sqrt(px * px + py * py)))
            lon[k] = math.degrees(math.atan2(py, px))

    @njit(parallel=True, cache=True)
    def inverse(lon, lat, p, e, att_t, att, orb_t, orb, npa, iters, out_line, out_pixel):
        ea = p[8]
        ee = 2. * p[7] - p[7] * p[7]
        eb = ea * (1. - p[7])
        a_, b_ = e[2, 0], e[2, 2]
        rr, base = math.sqrt(a_ * a_ + b_ * b_), math.atan2(b_, a_)
        for k in prange(lon.shape[0]):
            out_line[k], out_pixel[k] = np.nan, np.nan
            if not (math.isfinite(lon[k]) and math.isfinite(lat[k])):
                continue
            phi, lam = math.radians(lat[k]), math.radians(lon[k])
            en = ea / math.sqrt(1. - ee * math.sin(phi) ** 2)
            stx = en * math.cos(phi) * math.cos(lam)
            sty = en * math.cos(phi) * math.sin(lam)
            stz = en * (1. - ee) * math.sin(phi)
            # 初值同 MGIVSR：行號由緯度粗估、像素號取中心
            L = p[2] - math.atan(math.sin(phi) / (6.610689 - math.cos(phi))) / p[0]
            P = p[3]
            satx = saty = satz = 0.
            for _ in range(iters):
                satx, saty, satz, spx, spy, spz, xx, xy, xz, yx, yy, yz = frame(obs_time(L, P, p), att_t, att, orb_t,
                                                                               orb, npa)
                vx, vy, vz = unit3(stx - satx, sty - saty, stz - satz)
                wx = xx * vx + xy * vy + xz * vz            # 自旋座標系的視線
                wy = yx * vx + yy * vy + yz * vz
                wz = spx * vx + spy * vy + spz * vz
                # 解 Rz(q) · ELMIS · (cos p, 0, sin p) = w：z 分量與 q 無關，先解 p，再由 x、y 分量解 q
                ang_p = base - math.acos(min(1., max(-1., wz / rr)))
                ux = e[0, 0] * math.cos(ang_p) + e[0, 2] * math.sin(ang_p)
                uy = e[1, 0] * math.cos(ang_p) + e[1, 2] * math.sin(ang_p)
                ang_q = math.atan2(wy, wx) - math.atan2(uy, ux)
                ang_q = (ang_q + math.pi) % (2. * math.pi) - math.pi
                L_new, P_new = p[2] + ang_p / p[0], p[3] + ang_q / p[1]
                done = abs(L_new - L) < 1e-6 and abs(P_new - P) < 1e-6
                L, P = L_new, P_new
                if done:
                    break
            # 地球背面看不到：衛星在該點的地平面之下
            if (satx - stx) * stx / (ea * ea) + (saty - sty) * sty / (ea * ea) + (satz - stz) * stz / (eb * eb) > 0.:
                out_line[k], out_pixel[k] = L, P


ENGINES = ('numba', 'cython', 'numpy')
_NOTICE = set()


def _notice(msg):
    if msg not in _NOTICE:
        _NOTICE.add(msg)
        print(tr('  [導航] ', '  [navigation] ') + msg)


def kernels(engine='numba'):
    """[navigation] engine -> (正算, 反算) 編譯函式；numpy 時回傳 None。
    cython 模組不存在時改用 numba，沒有 numba 時改用 numpy (結果都相同)。"""
    if engine == 'cython':
        try:
            import vissr_nav_cy
            return vissr_nav_cy.forward, vissr_nav_cy.inverse
        except ImportError as err:
            _notice(tr('無法載入 Cython 模組 vissr_nav_cy (%s)，改用 numba；'
                       '編譯方式：python sensors/geo/build_vissr_nav_cy.py',
                       'cannot load the Cython module vissr_nav_cy (%s); using numba instead; '
                       'to build: python sensors/geo/build_vissr_nav_cy.py') % err)
            engine = 'numba'
    if engine == 'numba':
        if numba is not None:
            return forward, inverse
        _notice(tr('沒有安裝 numba，改用 numpy 計算 (結果相同，較慢)',
                   'numba is not installed; using numpy instead (same result, slower)'))
    return None


# ==================================================================================================
#  導航類別
# ==================================================================================================
class VissrNav(Nav):
    """JMA VISSR 導航。

    shape     (列數, 行數)；陣列第 r 列 = 行號 first_line + r x line_step，第 c 行 = 像素號 first_pixel + c x pixel_step
    scan      dict(step=行步進角 RSTEP, samp=像素取樣角 RSAMP (rad), line0=中心行號 RFCL, pixel0=中心像素號 RFCP,
                   sensors=感測器數 SENSSU (一圈同時掃描的行數), elmis=失準矩陣 ELMIS (3x3，[列][行]),
                   t0=第 1 圈的掃描時刻 DTIMS (MJD), spin=自旋速率 DSPIN (rpm))
    attitude  dict(mjd, alpha=自旋軸赤經, delta=自旋軸赤緯, beta=太陽-地球夾角)  (rad，1950 平均座標)
    orbit     dict(mjd, sat=(n,3) 地固座標的衛星位置 (m), gst=格林威治恆星時, sun_az=太陽方位, sun_el=太陽仰角 (度，地固座標),
                   npa=(n,3,3) 歲差章動矩陣 (A1~A9 依 JMA 次序))
    """

    def __init__(self, shape, scan, attitude, orbit, first_line=1, first_pixel=1, line_step=1, pixel_step=1, key=None,
                 engine='numba'):
        self.shape = (int(shape[0]), int(shape[1]))
        self._k = kernels(engine)
        self.first_line, self.first_pixel = float(first_line), float(first_pixel)
        self.line_step, self.pixel_step = float(line_step), float(pixel_step)
        self.scan = dict(scan)
        self.elmis = np.ascontiguousarray(np.asarray(scan['elmis'], dtype=np.float64).reshape(3, 3))
        self._key = key
        a_mjd, alpha, delta, beta = _clean(attitude['mjd'], attitude['alpha'], attitude['delta'], attitude['beta'])
        # 角度先展開 (unwrap)，線性內插才不會跨過 ±π (同 MGI100 對 β、恆星時、太陽方位的 360 度處理)
        self._att_t = a_mjd
        self._att = np.ascontiguousarray(np.stack([np.unwrap(alpha), np.unwrap(delta), np.unwrap(beta)], axis=-1))
        o_mjd, sat, gst, az, el, npa = _clean(orbit['mjd'], orbit['sat'], orbit['gst'], orbit['sun_az'],
                                              orbit['sun_el'], orbit['npa'])
        self._orb_t = o_mjd
        self._orb = np.ascontiguousarray(np.column_stack([sat, np.unwrap(np.deg2rad(gst)), np.unwrap(np.deg2rad(az)),
                                                          np.unwrap(np.deg2rad(el))]))
        self._npa = np.ascontiguousarray(npa.reshape(-1, 3, 3))
        s = self.scan
        self._p = np.array([s['step'], s['samp'], s['line0'], s['pixel0'], s['sensors'], s['t0'], s['spin'], EF, EA],
                           dtype=np.float64)

    def key(self):
        return ('vissr', self._key) if self._key is not None else id(self)

    # ---- 行列號 <-> 行號、像素號 ------------------------------------------------
    def _line_pixel(self, row, col):
        return self.first_line + row * self.line_step, self.first_pixel + col * self.pixel_step

    def _row_col(self, line, pixel):
        return (line - self.first_line) / self.line_step, (pixel - self.first_pixel) / self.pixel_step

    def obs_time(self, line, pixel):
        """行號、像素號 (由 1 起算) -> 觀測時刻 (MJD)。"""
        s = self.scan
        spins = np.floor((np.asarray(line, dtype=np.float64) - 1.) / s['sensors'])
        spins = spins + np.asarray(pixel, dtype=np.float64) * s['samp'] / (2. * np.pi)
        return s['t0'] + spins / (s['spin'] * 1440.)

    # ---- 正算 / 反算 ------------------------------------------------------------
    def to_lonlat(self, row, col):
        row, col = np.broadcast_arrays(np.asarray(row, dtype=np.float64), np.asarray(col, dtype=np.float64))
        shape = row.shape
        line, pixel = self._line_pixel(row.ravel(), col.ravel())
        lon, lat = np.empty(line.size), np.empty(line.size)
        k = self._k
        if k is not None:
            k[0](np.ascontiguousarray(line), np.ascontiguousarray(pixel), self._p, self.elmis, self._att_t, self._att,
                 self._orb_t, self._orb, self._npa, lon, lat)
        else:
            for s in range(0, line.size, _CHUNK):
                sl = slice(s, s + _CHUNK)
                lon[sl], lat[sl] = self._forward_np(line[sl], pixel[sl])
        return lon.reshape(shape), lat.reshape(shape)

    def to_rowcol(self, lon, lat):
        lon, lat = np.broadcast_arrays(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        shape = lon.shape
        lon, lat = np.ascontiguousarray(lon.ravel()), np.ascontiguousarray(lat.ravel())
        line, pixel = np.empty(lon.size), np.empty(lon.size)
        k = self._k
        if k is not None:
            k[1](lon, lat, self._p, self.elmis, self._att_t, self._att, self._orb_t, self._orb, self._npa, _ITER,
                 line, pixel)
        else:
            for s in range(0, lon.size, _CHUNK):
                sl = slice(s, s + _CHUNK)
                line[sl], pixel[sl] = self._inverse_np(lon[sl], lat[sl])
        rows, cols = self._row_col(line, pixel)
        return rows.reshape(shape), cols.reshape(shape)

    # ---- numpy 版 (沒有 numba 時；與 numba 版相同的計算) ---------------------------
    def _frames_np(self, t):
        att_t, att, orb_t, orb, npa = self._att_t, self._att, self._orb_t, self._orb, self._npa
        i = np.clip(np.searchsorted(att_t, t, side='right') - 1, 0, len(att_t) - 2)
        w = ((t - att_t[i]) / (att_t[i + 1] - att_t[i]))[:, None]
        alp, dlt, bet = (att[i] + w * (att[i + 1] - att[i])).T
        j = np.clip(np.searchsorted(orb_t, t, side='right') - 1, 0, len(orb_t) - 2)
        v = ((t - orb_t[j]) / (orb_t[j + 1] - orb_t[j]))[:, None]
        o = orb[j] + v * (orb[j + 1] - orb[j])
        sat, gst, az, el = o[:, :3], o[:, 3], o[:, 4], o[:, 5]
        att1 = np.stack([np.sin(dlt), -np.cos(dlt) * np.sin(alp), np.cos(dlt) * np.cos(alp)], axis=-1)
        att2 = np.einsum('nij,nj->ni', npa[j], att1)
        cg, sg = np.cos(gst), np.sin(gst)
        unit = lambda a: a / np.linalg.norm(a, axis=-1, keepdims=True)
        sp = unit(np.stack([cg * att2[:, 0] + sg * att2[:, 1], -sg * att2[:, 0] + cg * att2[:, 1], att2[:, 2]], axis=-1))
        ss = np.stack([np.cos(el) * np.cos(az), np.cos(el) * np.sin(az), np.sin(el)], axis=-1)
        sw1 = unit(np.cross(sp, ss))
        sw2 = unit(np.cross(sw1, sp))
        sx = unit(sw1 * np.sin(bet)[:, None] + sw2 * np.cos(bet)[:, None])
        sy = unit(np.cross(sp, sx))
        return sat, sp, sx, sy

    def _forward_np(self, line, pixel):
        s, e = self.scan, self.elmis
        sat, sp, sx, sy = self._frames_np(self.obs_time(line, pixel))
        p = s['step'] * (line - s['line0'])
        q = s['samp'] * (pixel - s['pixel0'])
        pc, ps, qc, qs = np.cos(p), np.sin(p), np.cos(q), np.sin(q)
        w1 = e[:, 0][None, :] * pc[:, None] + e[:, 2][None, :] * ps[:, None]
        w2x = qc * w1[:, 0] - qs * w1[:, 1]
        w2y = qs * w1[:, 0] + qc * w1[:, 1]
        sl = sx * w2x[:, None] + sy * w2y[:, None] + sp * w1[:, 2:3]
        sl = sl / np.linalg.norm(sl, axis=-1, keepdims=True)
        a = _F2 * (sl[:, 0] ** 2 + sl[:, 1] ** 2) + sl[:, 2] ** 2
        b = _F2 * (sat[:, 0] * sl[:, 0] + sat[:, 1] * sl[:, 1]) + sat[:, 2] * sl[:, 2]
        c = _F2 * (sat[:, 0] ** 2 + sat[:, 1] ** 2 - EA * EA) + sat[:, 2] ** 2
        with np.errstate(invalid='ignore'):
            d = (-b - np.sqrt(b * b - a * c)) / a
        pt = sat + d[:, None] * sl
        lat = np.rad2deg(np.arctan2(pt[:, 2], _F2 * np.hypot(pt[:, 0], pt[:, 1])))
        lon = np.rad2deg(np.arctan2(pt[:, 1], pt[:, 0]))
        return lon, lat

    def _inverse_np(self, lon, lat):
        s, e = self.scan, self.elmis
        ok = np.isfinite(lon) & np.isfinite(lat)
        phi, lam = np.deg2rad(np.where(ok, lat, 0.)), np.deg2rad(np.where(ok, lon, 0.))
        en = EA / np.sqrt(1. - _EE * np.sin(phi) ** 2)
        stn = np.stack([en * np.cos(phi) * np.cos(lam), en * np.cos(phi) * np.sin(lam),
                        en * (1. - _EE) * np.sin(phi)], axis=-1)
        line = s['line0'] - np.arctan(np.sin(phi) / (6.610689 - np.cos(phi))) / s['step']
        pixel = np.full(line.shape, float(s['pixel0']))
        a_, b_ = e[2, 0], e[2, 2]
        rr, base = np.hypot(a_, b_), np.arctan2(b_, a_)
        for _ in range(_ITER):
            sat, sp, sx, sy = self._frames_np(self.obs_time(line, pixel))
            v = stn - sat
            v = v / np.linalg.norm(v, axis=-1, keepdims=True)
            w = np.stack([np.sum(sx * v, 1), np.sum(sy * v, 1), np.sum(sp * v, 1)], axis=-1)
            p = base - np.arccos(np.clip(w[:, 2] / rr, -1., 1.))
            u = e[:, 0][None, :] * np.cos(p)[:, None] + e[:, 2][None, :] * np.sin(p)[:, None]
            q = np.arctan2(w[:, 1], w[:, 0]) - np.arctan2(u[:, 1], u[:, 0])
            q = (q + np.pi) % (2. * np.pi) - np.pi
            new_line, new_pixel = s['line0'] + p / s['step'], s['pixel0'] + q / s['samp']
            done = np.all(np.abs(new_line - line) < 1e-6) and np.all(np.abs(new_pixel - pixel) < 1e-6)
            line, pixel = new_line, new_pixel
            if done:
                break
        normal = stn / np.array([EA * EA, EA * EA, EB * EB])
        visible = ok & (np.sum((sat - stn) * normal, 1) > 0)
        return np.where(visible, line, np.nan), np.where(visible, pixel, np.nan)


# ==================================================================================================
#  地球邊緣修正 (MSC Technical Note No.26)
# ==================================================================================================
def earth_edges(nav, rows, coarse=16, iters=14):
    """導航預測的各列地球西、東邊緣 (陣列行號，由 0 起算，次像素)；看不到完整地球寬度的列為 NaN。"""
    rows = np.asarray(rows, dtype=np.float64)
    nx = nav.shape[1]
    cols = np.unique(np.r_[np.arange(0, nx, coarse), nx - 1]).astype(np.float64)
    lon, _ = nav.to_lonlat(*np.meshgrid(rows, cols, indexing='ij'))
    fin = np.isfinite(lon)
    has = fin.any(1) & ~fin[:, 0] & ~fin[:, -1]
    iw = np.argmax(fin, 1)
    ie = len(cols) - 1 - np.argmax(fin[:, ::-1], 1)
    out = []
    for space, earth in ((cols[np.maximum(iw - 1, 0)], cols[iw]), (cols[np.minimum(ie + 1, len(cols) - 1)], cols[ie])):
        lo, hi = space.copy(), earth.copy()
        for _ in range(iters):                         # 二分法找太空 / 地球的交界
            mid = 0.5 * (lo + hi)
            on_earth = np.isfinite(nav.to_lonlat(rows, mid)[0])
            hi, lo = np.where(on_earth, mid, hi), np.where(on_earth, lo, mid)
        out.append(np.where(has, 0.5 * (lo + hi), np.nan))
    return out[0], out[1]


def earth_edge_deviation(nav, rows, west, east, lat_band=(15., 65.)):
    """地球邊緣修正 (MSC Technical Note No.26，Kigawa 1993「Fine Tuning of Stretched-VISSR Image Mapping」)：
    以各掃描線實際偵測的地球東西邊緣 (LCW 的 West / East side earth edge) 與導航預測的邊緣比較，求影像的偏移。
      東西：實際中心 - 預測中心 (同原文 (4)~(6) 式)
      南北：原文以 45 度帶的地球寬度變化 (扣掉經驗常數 B 的大氣寬度) 估計；這裡改以兩個半球 15~65 度的各線寬度
            同時用最小平方解出偏移與大氣寬度 A：實際寬度(L) = 預測寬度(L - 偏移) + A
            (南北半球寬度隨行號的變化正負相反，所以偏移與 A 可以分開，不需要經驗常數)
    rows：陣列列號；west、east：該列實際的邊緣 (像素號，由 1 起算)。
    回傳 dict(line=南北修正 (行，實際 - 導航，向南為正), pixel=東西修正 (像素，向東為正), A=大氣寬度 (像素),
              rms=寬度殘差, n=使用的線數)；資料不足時回傳 None。"""
    rows = np.asarray(rows, dtype=np.float64)
    we, ee = earth_edges(nav, rows)
    ok = np.isfinite(we) & np.isfinite(ee)
    if ok.sum() < 50:
        return None
    rows, west, east, we, ee = rows[ok], np.asarray(west, float)[ok], np.asarray(east, float)[ok], we[ok] + 1., ee[ok] + 1.
    order = np.argsort(rows)
    rows, west, east, we, ee = rows[order], west[order], east[order], we[order], ee[order]
    width = ee - we
    slope = np.gradient(width, rows)
    _, lat = nav.to_lonlat(rows, 0.5 * (we + ee) - 1.)
    zone = (np.abs(lat) >= lat_band[0]) & (np.abs(lat) <= lat_band[1])
    if zone.sum() < 50 or not ((lat[zone] > 0).any() and (lat[zone] < 0).any()):
        return None
    y = (east - west)[zone] - width[zone]
    a = np.column_stack([-slope[zone], np.ones(zone.sum())])
    keep = np.ones(zone.sum(), bool)
    for _ in range(3):                                 # 去掉離群值 (邊緣被雲或雜訊干擾的線) 後重解
        coef = np.linalg.lstsq(a[keep], y[keep], rcond=None)[0]
        resid = y - a @ coef
        mad = np.median(np.abs(resid[keep] - np.median(resid[keep]))) * 1.4826 + 1e-6
        keep = np.abs(resid) < 3. * mad
    centre = 0.5 * (west + east) - 0.5 * (we + ee)
    near = np.abs(lat) < 65.
    return dict(line=float(coef[0]), pixel=float(np.median(centre[near])), A=float(coef[1]),
                rms=float(np.sqrt(np.mean(resid[keep] ** 2))), n=int(keep.sum()))
