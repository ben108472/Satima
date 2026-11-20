"""MSG SEVIRI (Meteosat-8~11) 共用：頻道、定標、導航 — SEVIRI-HRIT / SEVIRI-NAT / SEVIRI-NC 讀取器共用
(檔名沒有 '-'，不會被當成讀取器載入)

依據 EUMETSAT「MSG Level 1.5 Image Data Format Description」(EUM/MSG/ICD/105)、
「Radiometric Calibration of MSG SEVIRI Level 1.5 Image Data in Equivalent Spectral Blackbody Radiance」
(EUM/OPS-MSG/TEN/03/0064)、「The Conversion from Effective Radiances to Equivalent Brightness Temperatures」
(EUM/MET/TEN/11/0569)；常數與慣例與 satpy (seviri_l1b_*) 相同。

頻道 1~12：VIS006 VIS008 IR_016 IR_039 WV_062 WV_073 IR_087 IR_097 IR_108 IR_120 IR_134 HRV (選單 BAND-01~12)。
影像方向 (HRIT、native)：第 0 列 = 最南、第 0 行 = 最東。行號 L 由南、行號 C 由東起算 (由 1 起)：
  可見光 / 紅外 (3 km)：x = (1856 - C + 東西偏移) x 格距、y = (L - 1856 + 南北偏移) x 格距
  HRV (1 km)：中心 5566 (= 11136/2 - 2)；全圓盤的 HRV 分上下兩個視窗，各有自己的東界。
  2017/12 以前處理的資料 (EarthModel TypeOfEarthModel = 1) 影像向西北偏 1.5 km：
  導航東西偏移 +0.5、南北 -0.5 (可見光 / 紅外像素；HRV 為 1.5 像素)；TypeOfEarthModel = 2 不偏移。
  地球：a = 赤道半徑、b = (北極 + 南極半徑)/2 (檔頭，6378.169 / 6356.5838 km)；衛星高度 35785.831 km。
定標：計數值 -> 輻射 L = 計數值 x 斜率 + 截距 (計數值 0 為 NaN，負值取 0)，mW m-2 sr-1 (cm-1)-1。
  [seviri] calibration：nominal = 檔頭的 CalSlope / CalOffset (預設)；
    gsics = 紅外用 MPEF 的 GSICS 係數 (斜率 = GSICSCalCoeff、截距 = GSICSOffsetCount x GSICSCalCoeff；沒有時用 nominal)；
    meirink = 可見光 (VIS006、VIS008、IR_016) 用 Meirink 2023 的斜率 (A + B x 2000 年起的天數/1000)/1000。
  反射率 = π x L x d^2 / F (d = 日地距離 (AU)，F = 頻道太陽輻照度)，0~1。
  亮溫：T_L15 = C2 VC / ln(C1 VC^3 / L + 1)；有效輻射 (PlannedChanProcessing = 2)：T = (T_L15 - β)/α，
    光譜輻射 (= 1，2008/05 以後)：T = A T_L15^2 + B T_L15 + C。
  壞線：line_validity 為 2 或 3 (缺 / 壞資料) 且輻射、幾何品質皆為 4 (不可用) 的行為 NaN。
"""
import datetime
import math

import numpy as np

from satima.geo.nav import GeosNav, Nav
from satima.geo.reader import Channel

SATELLITES = {321: 'Meteosat-8', 322: 'Meteosat-9', 323: 'Meteosat-10', 324: 'Meteosat-11'}
SATID = {'MSG1': 321, 'MSG2': 322, 'MSG3': 323, 'MSG4': 324, 'MET08': 321, 'MET09': 322, 'MET10': 323, 'MET11': 324}
NAMES = ('VIS006', 'VIS008', 'IR_016', 'IR_039', 'WV_062', 'WV_073', 'IR_087', 'IR_097', 'IR_108', 'IR_120',
         'IR_134', 'HRV')
NUMBER = {n: i + 1 for i, n in enumerate(NAMES)}
# (中心波長 µm、種類、功能代稱)
_BANDS = {1: (0.635, 'vis', 'Visible'), 2: (0.81, 'nir', 'Veggie'), 3: (1.64, 'nir', 'Snow/Ice'),
          4: (3.92, 'swir', 'Shortwave Window'), 5: (6.25, 'wv', 'Upper-level Water Vapor'),
          6: (7.35, 'wv', 'Lower-level Water Vapor'), 7: (8.70, 'ir', 'Cloud-top Phase'), 8: (9.66, 'ir', 'Ozone'),
          9: (10.8, 'ir', 'LW Window'), 10: (12.0, 'ir', 'Dirty LW Window'), 11: (13.4, 'ir', 'CO2 LW'),
          12: (0.75, 'vis', 'HRV (High Resolution Visible)')}
REFLECTIVE = (1, 2, 3, 12)
H = 35785831.0                         # 衛星離地高度 (m)
CENTER_VISIR, CENTER_HRV = 1856, 5566
C1, C2 = 1.19104273e-5, 1.43877523

# 太陽輻照度 F (mW m-2 (cm-1)-1)：VIS006、VIS008、IR_016、HRV
_F = {321: (65.2296, 73.0127, 62.3715, 78.7599), 322: (65.2065, 73.1869, 61.9923, 79.0113),
      323: (65.5148, 73.1807, 62.0208, 78.9416), 324: (65.2656, 73.1692, 61.9416, 79.0035)}
# 紅外 (頻道 4~11)：中心波數 VC、α、β
_IR = {
    321: ((2567.33, 0.9956, 3.41), (1598.103, 0.9962, 2.218), (1362.081, 0.9991, 0.478), (1149.069, 0.9996, 0.179),
          (1034.343, 0.9999, 0.06), (930.647, 0.9983, 0.625), (839.66, 0.9988, 0.397), (752.387, 0.9981, 0.578)),
    322: ((2568.832, 0.9954, 3.438), (1600.548, 0.9963, 2.185), (1360.330, 0.9991, 0.47), (1148.620, 0.9996, 0.179),
          (1035.289, 0.9999, 0.056), (931.7, 0.9983, 0.64), (836.445, 0.9988, 0.408), (751.792, 0.9981, 0.561)),
    323: ((2547.771, 0.9915, 2.9002), (1595.621, 0.9960, 2.0337), (1360.337, 0.9991, 0.4340),
          (1148.130, 0.9996, 0.1714), (1034.715, 0.9999, 0.0527), (929.842, 0.9983, 0.6084),
          (838.659, 0.9988, 0.3882), (750.653, 0.9982, 0.5390)),
    324: ((2555.280, 0.9916, 2.9438), (1596.080, 0.9959, 2.0780), (1361.748, 0.9990, 0.4929),
          (1147.433, 0.9996, 0.1731), (1034.851, 0.9998, 0.0597), (931.122, 0.9983, 0.6256),
          (839.113, 0.9988, 0.4002), (748.585, 0.9981, 0.5635))}
# 光譜輻射 -> 亮溫的二次式 [A, B, C] (頻道 4~11)
_BTFIT = ((0.0, 1.011751900, -3.550400), (0.00001805700, 1.000255533, -1.790930),
          (0.00000231818, 1.000668281, -0.456166), (-0.00002332000, 1.011803400, -1.507390),
          (-0.00002055330, 1.009370670, -1.030600), (-0.00007392770, 1.032889800, -3.296740),
          (-0.00007009840, 1.031314600, -3.181090), (-0.00007293450, 1.030424800, -2.645950))
# Meirink 2023 可見光斜率 (A, B)：VIS006、VIS008、IR_016
_MEIRINK = {321: ((24.346, 0.3739), (30.989, 0.3111), (22.869, 0.0065)),
            322: ((21.026, 0.2556), (26.875, 0.1835), (21.394, 0.0498)),
            323: ((19.829, 0.5856), (25.284, 0.6787), (23.066, -0.0286)),
            324: ((20.515, 0.3600), (25.803, 0.4844), (22.354, -0.0187))}


def make_channel(n):
    """頻道 n -> Channel (選單、標註、檔名 BAND-01~12；HRV 只在選單的功能代稱標示 "HRV")。"""
    wl, kind, alias = _BANDS[n]
    return Channel(NAMES[n - 1], wl, 1.0 if n == 12 else 3.0, kind, band='BAND-%02d' % n, alias=alias)


def cds(days, msec, usec=0):
    """CCSDS 日數 (1958/01/01 起) + 毫秒 -> datetime；0 日 0 毫秒 (填充) 為 None。"""
    days, msec = int(days), int(msec)
    if days == 0 and msec == 0:
        return None
    return datetime.datetime(1958, 1, 1) + datetime.timedelta(days=days, milliseconds=msec, microseconds=int(usec))


def nominal_time(start, reduced_scan=False):
    """名目時間：觀測開始取整到 15 分鐘 (快速掃描 RSS 為 5 分鐘)。"""
    step = 5 if reduced_scan else 15
    t = start.replace(second=0, microsecond=0)
    return t - datetime.timedelta(minutes=t.minute % step)


def unpack10(raw):
    """10 位元 (每 5 bytes 4 個值，高位元在前) -> uint16。"""
    b = np.frombuffer(raw, np.uint8)
    n = len(b) // 5
    g = b[:n * 5].reshape(n, 5).astype(np.uint16)
    out = np.empty((n, 4), np.uint16)
    out[:, 0] = (g[:, 0] << 2) | (g[:, 1] >> 6)
    out[:, 1] = ((g[:, 1] & 0x3F) << 4) | (g[:, 2] >> 4)
    out[:, 2] = ((g[:, 2] & 0x0F) << 6) | (g[:, 3] >> 2)
    out[:, 3] = ((g[:, 3] & 0x03) << 8) | g[:, 4]
    return out.ravel()


def sun_earth_distance(t):
    """日地距離 (AU) (同 pyorbital.astronomy.sun_earth_distance_correction)。"""
    jd2000 = (t - datetime.datetime(2000, 1, 1, 12)).total_seconds() / 86400.
    return 1. - 0.0167 * math.cos(2. * math.pi * (jd2000 - 3.) / 365.25636)


def bad_lines(validity, rquality, gquality):
    v, r, g = (np.asarray(a).astype(np.int64) for a in (validity, rquality, gquality))
    return (v >= 2) & (v <= 3) & (r == 4) & (g == 4)


class Calibration:
    """一個頻道的定標：計數值 -> 反射率 (0~1) 或亮溫 (°C)。"""

    def __init__(self, satid, channel, slope, offset, rad_type=2, gsics=(0., 0.), mode='nominal', time=None):
        self.n = channel
        self.satid = satid if satid in _F else 324
        mode = str(mode).lower()
        self.mode = 'nominal'
        gain, off = float(slope), float(offset)
        if mode == 'gsics' and 4 <= channel <= 11 and gsics[0] and gsics[1]:
            gain, off = float(gsics[0]), float(gsics[1]) * float(gsics[0])
            self.mode = 'gsics'
        elif mode == 'meirink' and channel in (1, 2, 3) and time is not None:
            a, b = _MEIRINK[self.satid][channel - 1]
            days = (time - datetime.datetime(2000, 1, 1)).total_seconds() / 86400.
            gain = (a + b * days / 1000.) / 1000.
            self.mode = 'meirink'
        self.gain, self.offset = gain, off
        self.rad_type = int(rad_type)
        self.d2 = sun_earth_distance(time) ** 2 if time is not None else 1.

    def radiance(self, counts):
        c = np.asarray(counts, dtype=np.float64)
        return np.where(c > 0, np.maximum(c * self.gain + self.offset, 0.), np.nan)

    def __call__(self, counts):
        L = self.radiance(counts)
        n = self.n
        if n in REFLECTIVE:
            f = _F[self.satid][3 if n == 12 else n - 1]
            return (np.pi * L * self.d2 / f).astype(np.float32)
        vc, alpha, beta = _IR[self.satid][n - 4]
        with np.errstate(divide='ignore', invalid='ignore'):
            t = C2 * vc / np.log(C1 * vc ** 3 / L + 1.)
        if self.rad_type == 1:
            a, b, c = _BTFIT[n - 4]
            t = a * t * t + b * t + c
        else:
            t = (t - beta) / alpha
        return np.where(L > 0, t - 273.15, np.nan).astype(np.float32)


def geos_nav(shape, first_line, first_col, line_step, col_step, sub_lon, a, b, earth_model=2, hrv=False):
    """影像 (第 0 列 = 南端第 first_line 行、第 0 行 = 東端第 first_col 行) 的 GeosNav。
    line_step / col_step：掃描角格距 (弧度)；a、b：地球半軸 (m)。"""
    center = CENTER_HRV if hrv else CENTER_VISIR
    shift = (1.5 if hrv else 0.5) if int(earth_model) == 1 else 0.
    x0 = (center - first_col + shift) * col_step
    y0 = (first_line - center - shift) * line_step
    return GeosNav(shape, x0, -col_step, y0, line_step, sub_lon, H, a, b, sweep='y')


class SplitNav(Nav):
    """HRV 全圓盤的上下視窗：列 [starts[k], starts[k+1]) 用 navs[k] (各視窗南北相同、東界不同)。"""

    def __init__(self, navs, starts, shape):
        self.navs, self.starts = list(navs), np.asarray(starts, dtype=np.float64)
        self.shape = (int(shape[0]), int(shape[1]))

    def key(self):
        return ('split',) + tuple(n.key() for n in self.navs) + tuple(self.starts.tolist())

    def _which(self, row):
        return np.clip(np.searchsorted(self.starts, np.floor(np.asarray(row) + 0.5), side='right') - 1,
                       0, len(self.navs) - 1)

    def to_lonlat(self, row, col):
        row, col = np.broadcast_arrays(np.asarray(row, dtype=np.float64), np.asarray(col, dtype=np.float64))
        k = self._which(row)
        lon, lat = np.full(row.shape, np.nan), np.full(row.shape, np.nan)
        for i, nav in enumerate(self.navs):
            m = k == i
            if m.any():
                lon[m], lat[m] = nav.to_lonlat(row[m], col[m])
        return lon, lat

    def to_rowcol(self, lon, lat):
        lon, lat = np.broadcast_arrays(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        r, c = self.navs[0].to_rowcol(lon, lat)
        r, c = np.array(r, dtype=np.float64), np.array(c, dtype=np.float64)
        ok = np.isfinite(r)
        k = np.where(ok, self._which(np.where(ok, r, 0.)), 0)
        for i, nav in enumerate(self.navs[1:], 1):
            m = k == i
            if m.any():
                c[m] = nav.to_rowcol(lon[m], lat[m])[1]
        return r, c
