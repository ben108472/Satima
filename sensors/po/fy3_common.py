"""風雲三號成像儀 (MERSI、VIRR) 讀取器共用的工具

  find_ds / attr_text / attr_array    HDF5 資料集與屬性 (FY-3A/B 的資料集在最上層，FY-3C 起分 Data、Calibration… 群組)
  scaled                              資料集的 Slope / Intercept 與填補值
  nsmc_times                          NSMC 的時間欄位 (2000-01-01 12:00 UTC 起算) -> datetime64
  planck_bt                           輻亮度 mW/(m2 sr cm-1) -> 亮溫 (K)，等效中心波數 cm-1
  scan_interp / geo_interp            逐掃描帶的線性內插 (不跨掃描帶，避免 bow-tie 處出錯)：
                                      tie point -> 全解析度，或 1 km -> 250 m；經緯度於地心直角座標內插
  solar_angles                        由時間與經緯度計算太陽天頂角 / 方位角 (定位檔沒有角度時)
  view_angles                         由星下點推算衛星天頂角 / 方位角 (定位檔沒有角度時)
  TrueColor                           真實色：太陽天頂角校正 + CREFL (係數借用 MODIS 相近波段) + 增強
NSMC 的格式說明：https://satellite.nsmc.org.cn/DataPortal/cn/support/document.html (各衛星的「一級產品使用說明」)
"""
import numpy as np

from satima.po import truecolor as tc

EPOCH = np.datetime64('2000-01-01T12:00:00', 'us')     # NSMC 時間的起點 (J2000.0，UTC 計數，不含閏秒)
R_EARTH = 6371.0                                         # km
C1 = 1.191042e-5        # 2hc^2   mW/(m2 sr cm-4)
C2 = 1.4387752          # hc/k    K cm


# ---- HDF5 -------------------------------------------------------------------
def find_ds(f, name, groups=('Data', '')):
    """依序在各群組找資料集；找不到回傳 None。"""
    for g in groups:
        path = '%s/%s' % (g, name) if g else name
        if path in f:
            return f[path]
    return None


def _text(v):
    if isinstance(v, np.ndarray):
        v = v.ravel()
        if v.size == 0:
            return ''
        if v.dtype.kind in 'SO':
            v = v[0]
        else:
            return ' '.join(str(x) for x in v)
    return v.decode('latin1').strip() if isinstance(v, bytes) else str(v).strip()


def attr_text(obj, key, default=''):
    return _text(obj.attrs[key]) if key in obj.attrs else default


def attr_array(obj, key):
    """數值屬性 -> float64 一維陣列；沒有時 None。"""
    if key not in obj.attrs:
        return None
    try:
        return np.asarray(obj.attrs[key], dtype=np.float64).ravel()
    except (TypeError, ValueError):
        return None


def scaled(ds, index=None):
    """資料集 (第 index 個波段) 的 (slope, intercept)。屬性有逐波段的陣列時取第 index 個。"""
    out = []
    for key, default in (('Slope', 1.), ('Intercept', 0.)):
        v = attr_array(ds, key)
        if v is None or v.size == 0:
            out.append(default)
        elif index is not None and v.size > index and v.size > 1:
            out.append(float(v[index]))
        else:
            out.append(float(v[0]))
    return tuple(out)


def nsmc_times(values, unit='s'):
    """2000-01-01 12:00 UTC 起算的秒 (或毫秒) -> datetime64[us]；無效 (<= 0、NaN) 為 NaT。"""
    v = np.asarray(values, dtype=np.float64)
    scale = {'s': 1e6, 'ms': 1e3}[unit]
    bad = ~np.isfinite(v) | (v <= 0)
    us = np.round(np.where(bad, 0., v) * scale).astype(np.int64)
    out = EPOCH + us.astype('timedelta64[us]')
    out[bad] = np.datetime64('NaT', 'us')
    return out


def day_ms_times(days, ms):
    """Day_Count (2000-01-01 12:00 起的日數) + Millisecond_Count (當日 12:00 起的毫秒) -> datetime64[us]。"""
    days = np.asarray(days, dtype=np.float64).ravel()
    ms = np.asarray(ms, dtype=np.float64).ravel()
    return nsmc_times(days * 86400. + ms / 1000.)


# ---- 輻射 -------------------------------------------------------------------
def planck_bt(rad, wn):
    """輻亮度 mW/(m2 sr cm-1) -> 等效黑體亮溫 Te (K)；wn：等效中心波數 cm-1。輻亮度 <= 0 為 NaN。"""
    rad = np.asarray(rad, dtype=np.float64)
    with np.errstate(invalid='ignore', divide='ignore'):
        t = C2 * wn / np.log1p(C1 * wn ** 3 / rad)
    t[~(rad > 0)] = np.nan
    return t


# ---- 逐掃描帶內插 -----------------------------------------------------------
def _weights(pos_in, pos_out):
    pos_in = np.asarray(pos_in, dtype=np.float64)
    pos_out = np.asarray(pos_out, dtype=np.float64)
    if pos_in.size == 1:
        return np.zeros(pos_out.size, int), np.zeros(pos_out.size)
    k = np.clip(np.searchsorted(pos_in, pos_out, side='right') - 1, 0, pos_in.size - 2)
    w = (pos_out - pos_in[k]) / (pos_in[k + 1] - pos_in[k])
    return k, w


def scan_interp(src, rows_in, rows_out, cols_in, cols_out):
    """逐掃描帶的雙線性內插 (範圍外線性外插)，不跨掃描帶。
    src      (掃描數 x len(rows_in), len(cols_in))：每個掃描帶 len(rows_in) 列
    rows_in  每個掃描帶內各輸入列在輸出網格上的位置 (輸出列號，0 起算)，例如 tie point [0, 5]、1 km -> 250 m 為 4k+1.5
    rows_out 每個掃描帶的輸出列數
    cols_in  各輸入行在輸出網格上的位置 (輸出行號，絕對位置)
    cols_out 要輸出的行號 (絕對位置)
    回傳 (掃描數 x rows_out, len(cols_out))，float64。"""
    ni = len(rows_in)
    src = np.asarray(src, dtype=np.float64)
    nscan = src.shape[0] // ni
    kr, wr = _weights(rows_in, np.arange(rows_out))
    kc, wc = _weights(cols_in, cols_out)
    s = src[:nscan * ni].reshape(nscan, ni, src.shape[1])
    a = s[:, :, kc] * (1. - wc) + s[:, :, np.minimum(kc + 1, s.shape[2] - 1)] * wc if len(cols_in) > 1 else s[:, :, kc]
    if ni > 1:
        b = a[:, kr, :] * (1. - wr)[None, :, None] + a[:, kr + 1, :] * wr[None, :, None]
    else:
        b = np.repeat(a, rows_out, axis=1)
    return b.reshape(nscan * rows_out, len(cols_out))


def geo_interp(lat, lon, rows_in, rows_out, cols_in, cols_out):
    """經緯度的逐掃描帶內插：轉成地心直角座標各自內插再轉回 (換日線、極區不會出錯)。回傳 float32 (lat, lon)。"""
    la, lo = np.deg2rad(np.asarray(lat, dtype=np.float64)), np.deg2rad(np.asarray(lon, dtype=np.float64))
    xyz = (np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la))
    x, y, z = (scan_interp(v, rows_in, rows_out, cols_in, cols_out) for v in xyz)
    out_lat = np.rad2deg(np.arctan2(z, np.hypot(x, y))).astype(np.float32)
    out_lon = np.rad2deg(np.arctan2(y, x)).astype(np.float32)
    return out_lat, out_lon


def angle_interp(az, rows_in, rows_out, cols_in, cols_out):
    """方位角 (度) 的內插：以 sin / cos 內插，避免 0 / 360 度的跳動。"""
    r = np.deg2rad(np.asarray(az, dtype=np.float64))
    s = scan_interp(np.sin(r), rows_in, rows_out, cols_in, cols_out)
    c = scan_interp(np.cos(r), rows_in, rows_out, cols_in, cols_out)
    return np.rad2deg(np.arctan2(s, c)) % 360.


def clean_latlon(lat, lon):
    lat = np.array(lat, dtype=np.float32)
    lon = np.array(lon, dtype=np.float32)
    with np.errstate(invalid='ignore'):
        bad = ~(np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90.) & (np.abs(lon) <= 360.))
    lat[bad] = np.nan
    lon[bad] = np.nan
    return lat, lon


# ---- 角度 (定位檔沒有角度時的備援) --------------------------------------------
def solar_angles(times, lats, lons):
    """太陽天頂角、方位角 (度；方位角由北順時針)。times：每一列的 datetime64 (長度 = 列數)。
    Astronomical Almanac 的低精度太陽位置 (誤差約 0.01 度)。"""
    t = ((np.asarray(times, dtype='datetime64[us]') - EPOCH) / np.timedelta64(1, 'D')).astype(np.float64)[:, None]
    g = np.deg2rad(357.529 + 0.98560028 * t)
    q = 280.459 + 0.98564736 * t
    ecl = np.deg2rad(q + 1.915 * np.sin(g) + 0.020 * np.sin(2. * g))
    eps = np.deg2rad(23.439 - 3.6e-7 * t)
    ra = np.arctan2(np.cos(eps) * np.sin(ecl), np.cos(ecl))
    dec = np.arcsin(np.sin(eps) * np.sin(ecl))
    gmst = np.deg2rad((18.697374558 + 24.06570982441908 * t) * 15.)
    h = gmst + np.deg2rad(np.asarray(lons, dtype=np.float64)) - ra
    la = np.deg2rad(np.asarray(lats, dtype=np.float64))
    cos_z = np.sin(la) * np.sin(dec) + np.cos(la) * np.cos(dec) * np.cos(h)
    sza = np.rad2deg(np.arccos(np.clip(cos_z, -1., 1.)))
    saz = np.rad2deg(np.arctan2(-np.sin(h), np.tan(dec) * np.cos(la) - np.sin(la) * np.cos(h))) % 360.
    return sza, saz


def nadir_track(tie_lat, tie_lon):
    """由整條掃描寬的 tie point 經緯度找每一列 tie 的星下點：沿掃描方向相鄰點距離最小處 (掃描角 0 度)。
    回傳 (lat, lon)，長度 = tie 列數。"""
    la, lo = np.deg2rad(tie_lat.astype(np.float64)), np.deg2rad(tie_lon.astype(np.float64))
    with np.errstate(invalid='ignore'):
        h = (np.sin((la[:, 1:] - la[:, :-1]) / 2.) ** 2
             + np.cos(la[:, 1:]) * np.cos(la[:, :-1]) * np.sin((lo[:, 1:] - lo[:, :-1]) / 2.) ** 2)
        d = np.where(np.isfinite(h), h, np.inf)
    j = np.argmin(np.nanmedian(np.where(np.isfinite(d), d, np.nan), axis=0))        # 整段同一個位置 (中位數)
    return 0.5 * (tie_lat[:, j] + tie_lat[:, j + 1]), 0.5 * (tie_lon[:, j] + tie_lon[:, j + 1])


def view_angles(lats, lons, nadir_lat, nadir_lon, sat_height_km):
    """衛星天頂角、方位角 (度)：由各列的星下點與衛星高度推算 (球面地球)。方位角為像元往衛星 (星下點) 的方向。"""
    la1, lo1 = np.deg2rad(np.asarray(lats, dtype=np.float64)), np.deg2rad(np.asarray(lons, dtype=np.float64))
    la2, lo2 = np.deg2rad(np.asarray(nadir_lat, dtype=np.float64))[:, None], np.deg2rad(
        np.asarray(nadir_lon, dtype=np.float64))[:, None]
    dlon = lo2 - lo1
    with np.errstate(invalid='ignore'):
        cg = np.clip(np.sin(la1) * np.sin(la2) + np.cos(la1) * np.cos(la2) * np.cos(dlon), -1., 1.)
        g = np.arccos(cg)
        rs = R_EARTH + sat_height_km
        vza = np.rad2deg(np.arctan2(rs * np.sin(g), rs * np.cos(g) - R_EARTH))
        vaz = np.rad2deg(np.arctan2(np.sin(dlon) * np.cos(la2),
                                    np.cos(la1) * np.sin(la2) - np.sin(la1) * np.cos(la2) * np.cos(dlon))) % 360.
    return vza, vaz


# ---- 真實色 -----------------------------------------------------------------
TC_STYLES = {1: 'General', 2: 'CIRA', 3: 'Worldview'}
TC_STYLE_DESC = {1: 'General (sunz only, no Rayleigh)', 2: 'CIRA (sunz + CREFL + CIRA stretch)',
                 3: 'NASA Worldview (sunz + CREFL + Worldview stretch)'}
TC_NIGHT_SZA = 88.5            # 太陽天頂角超過此值的像元畫成黑色 (同 VIIRS)
# CREFL 沒有 FY-3 成像儀的官方係數：O3 / H2O 吸收借用 MODIS 相近波段 (CREFL 1.7.1，同 satima MODIS-HDF4)，
# Rayleigh 光學厚度以 Hansen & Travis (1974) 公式由等效中心波長計算。  紅 ~ MODIS B1、綠 ~ B4、藍 ~ B3
TC_GAS = {'R': (-5.60723, 0.820175, 0.0715289), 'G': (0., 0., 0.089691), 'B': (0., 0., 0.00743232)}
TC_UO3, TC_UH2O = 0.319, 2.93


def crefl(refl, color, wl_um, sza, vza, phi, height=None):
    """太陽天頂角已校正的 TOA 反射率 -> CREFL 校正後的反射率 (Rayleigh + O3 / H2O)。"""
    ah2o, bh2o, ao3 = TC_GAS[color]
    taur = tc.rayleigh_tau(wl_um)
    if height is not None:
        taur = taur * np.exp(-np.clip(height, 0., None) / tc.SCALE_HEIGHT)
    with np.errstate(invalid='ignore', divide='ignore'):
        mus = np.cos(np.deg2rad(sza))
        mus = np.where(mus > 0, mus, np.nan)
        muv = np.cos(np.deg2rad(vza))
        air_mass = np.minimum(1. / mus + 1. / muv, tc.MAXAIRMASS)
        t_o3 = tc.o3_transmittance(air_mass, TC_UO3, ao3)
        t_h2o = tc.h2o_transmittance_modis(air_mass, TC_UH2O, ah2o, bh2o)
        out = tc.crefl_core(refl, mus, muv, phi, taur, t_o3, t_h2o)
    return np.clip(out, tc.REFLMIN, tc.REFLMAX)


class TrueColor:
    """真實色：輸入 R / G / B 的定標反射率 (NSMC 的 Ref*，= π L / E0，未除 cos 太陽天頂角)。
    General = 只做太陽天頂角校正 (含日地距離)，線性；CIRA / Worldview = 再做 CREFL 與非線性增強。"""

    def __init__(self, refl, wls, des, angles, height=None):
        """refl：{'R','G','B': 陣列}；wls：等效中心波長 (µm)；des：日地距離比；angles：(sza, saz, vza, vaz)。"""
        sza, saz, vza, vaz = angles
        self.sza = sza
        corr = tc.sunz_corr(sza, 88., 95.) * des ** 2
        self.sunz = {c: refl[c] * corr for c in 'RGB'}
        self.args = (wls, sza, vza, saz - vaz, height)
        self._crefl = None

    def corrected(self):
        if self._crefl is None:
            wls, sza, vza, phi, height = self.args
            self._crefl = {c: crefl(self.sunz[c], c, wls[c], sza, vza, phi, height) for c in 'RGB'}
        return self._crefl

    def rgb(self, style):
        if style == 1:
            src = self.sunz
            out = np.dstack([np.clip(np.nan_to_num(src[c], nan=0.), 0., 1.) for c in 'RGB'])
        else:
            src = self.corrected()
            lin = np.dstack([np.clip(src[c], 0., None) for c in 'RGB'])
            out = tc.cira_stretch(lin) if style == 2 else tc.worldview_stretch(lin)
            out = np.nan_to_num(np.asarray(out, dtype=np.float64), nan=0.)
        with np.errstate(invalid='ignore'):
            out[~(self.sza < TC_NIGHT_SZA)] = 0.
        return np.clip(out, 0., 1.).astype(np.float32)
