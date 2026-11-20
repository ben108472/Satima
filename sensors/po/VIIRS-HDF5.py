"""
VIIRS-HDF5 — S-NPP / NOAA-20 (JPSS-1) / NOAA-21 (JPSS-2) VIIRS SDR (HDF5)
原始腳本：SatimaV3NPPJPSS-AMAZS3-loop.py

所需檔案 (同一目錄)
  SVxnn_<sat>_d<日期>_t<開始>_e<結束>_b..._c..._ops.h5   例：SVI05、SVM15、SVDNB
  地理定位：GITCO (I-band)、GMTCO (M-band)、GDNBO (DNB)
  * NOAA CLASS 把定位與資料包在同一個檔的合併檔 (檔名前段列出所含的產品，例如
    GMTCO-SVM03-SVM04-SVM05_j01_d..._ops.h5、GDNBO-SVDNB_...) 也可讀，分開的檔與合併檔可以混用；
    幾個片段併成一個檔的聚合檔照檔名的起訖時間當成一個較長的片段。
  * 定位優先用地形校正版 (GMTCO / GITCO)；找不到時改用未做地形校正的 GMODO / GIMGO (畫面上會註明)。
    海面上兩者沒有差別；高山在掃描帶邊緣會偏移 (約為地形高度 x tan(視角))。
  * DNB：GDNBO 有地形校正的經緯度 (Latitude_TC / Longitude_TC) 時優先使用。
  * 一個片段 (granule) 約 85 秒：同衛星、時間相接的片段自動接成一筆 (資料列表的備註欄顯示檔數與缺檔數)，
    範圍確定後只讀與地圖相交的片段 (前後各多帶一個)。
  * 有資料但找不到定位檔的頻道不列入選單，資料列表會標出缺哪一種定位。
  (2026/11 以 S-NPP、NOAA-20、NOAA-21 的真實檔案驗證過 I5、I1、M5、DNB、TrueColor-M / -I 與 Latitude_TC；合併檔以真實資料組出來的檔
   測過。CLASS 的原始合併檔與聚合檔、GMODO / GIMGO、M 頻道亮溫 (M12~M16) 還沒有真實檔案可測。)

產品
  單頻道  反射率 M1–M11、I1–I3 / 亮溫 M12–M16、I4–I5 / DNB (4 種縮放)
  真實色  TrueColor-M (M5/M4/M3, 750 m) / TrueColor-I (I1 比值銳化, 375 m)

與原腳本的差異
  - 頻道名稱：原寫法 'M10'.replace('0', '') 會變成 'M1'，改為只去掉前導 0
  - 反射率頻道補上 M8–M11、I3 (原本會被當成亮溫而讀不到資料)；所有亮溫頻道都有色表選單
  - 單頻道的資料列與經緯度列一起去除填充列，避免兩者列數不一致
  - 真實色：CREFL 加上 GMTCO/GITCO 的地形高度修正與空氣質量上限 (同 MODIS)，新增 CIRA 樣式；
    CREFL 共用部分移到 satima/truecolor.py
  - DNB 縮放與真實色樣式可選 0.ALL 一次輸出全部
  - 中心時間：取最接近輸入中心點的像元，用地理定位檔逐掃描帶的 MidTime (IET，已扣閏秒)
    (原本以檔名起訖時間依裁切範圍的中間列線性推算，去除填充列或選不連續的 granule 時會偏移)
"""
import datetime
import os
import re

import h5py
import numpy as np
import pandas as pd
from matplotlib import colors

from satima.core import segments as seg
from satima.core.i18n import tr
from satima.po import SatReader, Segment, ui, swath, obstime, truecolor as tc
from satima.core.colormaps import ask_colormaps

DNB_THRESHOLD = 4e-8

# 頻道分類 (VIIRS SDR：M1–M11、I1–I3 為反射率，M12–M16、I4–I5 為亮溫)
REFL_BANDS = tuple('M%d' % i for i in range(1, 12)) + ('I1', 'I2', 'I3')
BT_BANDS = tuple('M%d' % i for i in range(12, 17)) + ('I4', 'I5')
SAT_NAMES = {'npp': 'S-NPP', 'j01': 'JPSS-1', 'j02': 'JPSS-2', 'j03': 'JPSS-3'}

# 檔名：<產品[-產品…]>_<衛星>_d<日期>_t<起 HHMMSSs>_e<訖 HHMMSSs>_b<軌道>_c<產製時間>_<來源>_ops.h5
# 產品：SVMnn / SVInn / SVDNB (SDR) 與定位 GMTCO / GMODO (M)、GITCO / GIMGO (I)、GDNBO (DNB)；CLASS 的合併檔有多個
_NAME = re.compile(r'^(?P<prods>[A-Z0-9]+(?:-[A-Z0-9]+)*)_(?P<sat>[a-z0-9]{3})_d(?P<date>[0-9]{8})_t(?P<t>[0-9]{7})'
                   r'_e(?P<e>[0-9]{7})_b(?P<orbit>[0-9]{5})_c(?P<created>[0-9]+)_(?P<src>[A-Za-z0-9]+)_ops[.]h5$')
# 定位產品 -> HDF5 群組；各類頻道依序嘗試 (地形校正版優先)
GEO_GROUPS = {'GMTCO': 'VIIRS-MOD-GEO-TC_All', 'GMODO': 'VIIRS-MOD-GEO_All', 'GITCO': 'VIIRS-IMG-GEO-TC_All',
              'GIMGO': 'VIIRS-IMG-GEO_All', 'GDNBO': 'VIIRS-DNB-GEO_All'}
GEO_ORDER = {'M': ('GMTCO', 'GMODO'), 'I': ('GITCO', 'GIMGO'), 'D': ('GDNBO',)}
_NOTICED = set()            # 已經提醒過「改用未做地形校正的定位」的種類


def channel_name(code):
    """檔名頻道代碼 -> 頻道名稱：'I05' -> 'I5'、'M10' -> 'M10'、'DNB' -> 'DNB'。"""
    if code == 'DNB':
        return code
    return code[0] + str(int(code[1:]))


# ===== ERF Dynamic DNB (CIRA/SSEC, 移植自 satpy ERFDNB / Polar2Grid dynamic_dnb) =====
# 依太陽天頂角、月亮天頂角與月亮照明比例動態決定縮放範圍; 輻射值單位為 SDR 原生的 W cm-2 sr-1
DNB_SAT_PCT = 0.005    # 飽和修正: 允許超過上限的像元比例
DNB_SAT_STEP = 1.1     # 飽和修正: 上限每次放大倍率
DNB_PHASE_CORRECTION = True   # True: 月光項依相位角計算; False: satpy 原版 (依照明比例)
DNB_MODES = {1: 'Original (fixed threshold)', 2: 'ERF Dynamic DNB',
             3: 'ERF Dynamic DNB + saturation correction', 4: 'Lunar reflectance (Miller & Turner 2009)'}


def lunar_phase_factor(phase):
    # Allen 月亮星等公式: 相位角 phase 時月光比滿月暗 dm 星等
    # 換算成 ERF 月光項 (夜間 10^(-2*factor) = 相對滿月的亮度); 上限 0.7 為 ERF 原本的新月 (氣輝) 設定
    p = abs(phase)
    dm = 0.026 * p + 4e-9 * p ** 4
    return min(0.2 * dm, 0.7)


def erf_dnb(rad, sza, lza, moon_frac, saturation=False, moon_phase=None):
    from scipy.special import erf
    rad = np.where(rad < -900., np.nan, rad)
    if moon_phase is not None:
        moon_factor1 = lunar_phase_factor(moon_phase)
    else:
        moon_factor1 = 0.7 * (1.0 - moon_frac)
    moon_factor2 = 0.0022 * lza
    erf_portion = 1 + erf((sza - 95.0) / (5.0 * np.sqrt(2.0)))
    max_val = np.power(10, -1.7 - (2.65 + moon_factor1 + moon_factor2) * erf_portion)
    min_val = np.power(10, -4.0 - (2.95 + moon_factor2) * erf_portion)
    if saturation:
        valid = np.isfinite(rad) & np.isfinite(max_val)
        nvalid = max(valid.sum(), 1)
        while np.sum(rad[valid] > max_val[valid]) / nvalid > DNB_SAT_PCT:
            max_val = max_val * DNB_SAT_STEP
    inner = (rad - min_val) / (max_val - min_val)
    inner = np.where(inner < 0, 0., inner)
    return np.clip(np.sqrt(inner), 0., 1.)


# ===== DNB 月光反射率 (Miller & Turner 2009, 移植自 DNB.py) =====
DNB_REFL_MAX = 1.0         # 顯示範圍 0..DNB_REFL_MAX (線性)
# DNB 等效頻寬 (um) = 官方 RSR (峰值正規化為 1, VIIRS SDR ATBD) 對波長積分, 16 個偵測器平均;
# SDR 的 DNB Radiance 為頻段積分值 (W/cm2/sr), 除以此值才是與 MT2009 表格一致的頻段平均光譜輻射.
# RSR 來源: pyspectral rsr_viirs_*.h5 (NOAA/NASA DAWG 發布版本)
DNB_BANDWIDTH_UM = {'npp': 0.3204, 'j01': 0.3331, 'j02': 0.3266}
DNB_BANDWIDTH_DEFAULT = 0.3266
NIGHT_SZA_MIN = 108.0      # 只轉換太陽天頂角大於此值 (天文夜) 的像元
LUNAR_ZENITH_MAX = 85.0    # 月亮高度低於 5 度的像元遮蔽 (反射率不可靠)

# MT2009 DNB 光譜響應加權的大氣頂月光照度 (mW m-2 um-1), 平均日地/地月距離, 相位角 0-180 度 (0=滿月)
MT2009_PHASE_DEG = np.arange(181, dtype=float)
MT2009_IRRAD_MW_M2_UM = np.array([
    4.11896, 4.02183, 3.92698, 3.83438, 3.74399, 3.65569, 3.56613, 3.47868,
    3.39348, 3.31033, 3.22924, 3.15167, 3.07599, 3.00212, 2.93009, 2.85973,
    2.7911, 2.7241, 2.65873, 2.59492, 2.53266, 2.47284, 2.41448, 2.35748,
    2.30179, 2.24746, 2.19441, 2.14263, 2.09204, 2.0427, 1.99445, 1.94665,
    1.9, 1.85443, 1.81001, 1.76663, 1.7243, 1.68297, 1.64264, 1.6033,
    1.56489, 1.52914, 1.49426, 1.46014, 1.42685, 1.39432, 1.36253, 1.33143,
    1.30106, 1.27137, 1.2424, 1.21369, 1.1856, 1.15817, 1.1314, 1.10524,
    1.07969, 1.05476, 1.03037, 1.0065, 0.983266, 0.958765, 0.934901, 0.91163,
    0.888928, 0.866796, 0.845213, 0.824174, 0.803657, 0.78365, 0.764138, 0.743059,
    0.722556, 0.702626, 0.683243, 0.664396, 0.646067, 0.628245, 0.610916, 0.594064,
    0.577678, 0.559427, 0.541756, 0.524635, 0.508063, 0.492009, 0.47647, 0.461418,
    0.446844, 0.432723, 0.419053, 0.403895, 0.389277, 0.375193, 0.361616, 0.348532,
    0.335923, 0.323766, 0.312051, 0.300761, 0.289879, 0.277654, 0.265944, 0.254728,
    0.243985, 0.233696, 0.223839, 0.214399, 0.205358, 0.196696, 0.1884, 0.179173,
    0.170397, 0.162053, 0.154113, 0.146566, 0.13939, 0.132561, 0.126068, 0.119897,
    0.114023, 0.107488, 0.101332, 0.0955243, 0.0900534, 0.0848946, 0.080031, 0.0754463,
    0.0711242, 0.0670496, 0.0632087, 0.0593831, 0.0557891, 0.0524129, 0.0492405, 0.0462609,
    0.0434608, 0.0408308, 0.0383591, 0.0360379, 0.0338565, 0.0315612, 0.0294209, 0.027426,
    0.0255659, 0.0238325, 0.0222165, 0.0207097, 0.0193054, 0.0179962, 0.0167755, 0.0155703,
    0.014451, 0.0134123, 0.0124481, 0.0115535, 0.0107233, 0.0099526, 0.00923708, 0.00857314,
    0.00795684, 0.00737225, 0.00683057, 0.00632876, 0.00586376, 0.00543289, 0.00503373, 0.0046639,
    0.00432122, 0.0040037, 0.00370952, 0.00344246, 0.00319467, 0.00296473, 0.00275132, 0.0025533,
    0.00236958, 0.0021991, 0.00204086, 0.00189405, 0.00175778
])
AU_KM = 149597870.7
MEAN_ES_KM = 149598022.6071
MEAN_EM_KM = 384400.0
R_EARTH_KM = 6378.140


def moon_geometry(dt):
    # 日地距離 km, 地月距離 km, 月亮相位角 (度, 0=滿月); 需要 ephem
    import ephem
    import math
    d = ephem.Date(dt)
    sun, moon = ephem.Sun(d), ephem.Moon(d)
    rs, rm = sun.earth_distance * AU_KM, moon.earth_distance * AU_KM
    e = abs(float(moon.elong))
    dsm = math.sqrt(rs ** 2 + rm ** 2 - 2 * rs * rm * math.cos(e))
    cosi = (rm ** 2 + dsm ** 2 - rs ** 2) / (2 * rm * dsm)
    return rs, rm, math.degrees(math.acos(max(-1.0, min(1.0, cosi))))


def lunar_irradiance_dnb(dt, phase_angle=None):
    # DNB 波段大氣頂月光照度 (W m-2 um-1) 與所用相位角; phase_angle 可用 GDNBO 的 MoonPhaseAngle
    import math
    rs, rm, ph_eph = moon_geometry(dt)
    ph = ph_eph if phase_angle is None else float(phase_angle)
    e1au = float(np.interp(ph, MT2009_PHASE_DEG, MT2009_IRRAD_MW_M2_UM))
    c = math.cos(math.radians(ph))
    t1 = MEAN_ES_KM ** 2 + MEAN_EM_KM ** 2 + 2 * MEAN_EM_KM * MEAN_ES_KM * c
    t2 = rs ** 2 + rm ** 2 + 2 * rm * rs * c
    t3 = ((MEAN_EM_KM - R_EARTH_KM) / (rm - R_EARTH_KM)) ** 2
    return e1au * (t1 / t2) * t3 * 1e-3, ph


def lunar_reflectance(rad, sza, lza, e_moon, bandwidth=DNB_BANDWIDTH_DEFAULT):
    # 月光反射率 = pi*L / (E_moon * cos(月亮天頂角)); rad 為 SDR 的 W/cm2/sr
    rad = np.where(rad < -900., np.nan, rad)
    radiance = rad * 1e4 / bandwidth
    with np.errstate(invalid='ignore', divide='ignore'):
        refl = np.pi * radiance / (e_moon * np.cos(np.deg2rad(lza)))
    refl[(sza <= NIGHT_SZA_MIN) | (lza >= LUNAR_ZENITH_MAX) | np.isnan(lza) | np.isnan(sza)] = np.nan
    if not np.isfinite(refl).any():
        print(tr('注意：沒有有效的月光反射率像元 (月亮在地平線下，或不是夜間的資料)。',
                 'WARNING: no valid lunar-reflectance pixels (moon below horizon or scene not at night).'))
    return refl


# ===== True Color (可見光真實色) =====
# M5=紅(0.672um) M4=綠(0.555um) M3=藍(0.488um)；I1(0.640um)僅用於解析度增強
# 大氣校正為 NASA CREFL (Corrected Reflectance) 演算法, 即 WorldView 影像所用,
# 移植自 NASA crefl.c / satpy crefl_utils.py, VIIRS 官方係數
# 注意：VIIRS SDR 的 Reflectance 已除以 cos(太陽天頂角) (satpy viirs_sdr 亦標記為 sunz_corrected)，不可再除
TC_COEFF = {           # taur0,   aH2O,        bH2O,     aO3
    'M5': (0.04350, 0.000406601, 0.812659, 0.0433461),
    'M4': (0.09740, 1.78644e-05, 0.867785, 0.0853012),
    'M3': (0.16176, 0.0,         1.0,      0.0178299),
    'I1': (0.05373, 0.00222253,  0.791204, 0.0813531),
}
TC_UO3 = 0.285    # 臭氧總量 (atm-cm)
TC_UH2O = 2.93    # 可降水量 (g/cm2)
TC_MAXSOLZ = 86.5
TC_MAXAIRMASS = tc.MAXAIRMASS   # 空氣質量上限 (同 MODIS / METimage)
TC_USE_HEIGHT = True            # 以 GMTCO/GITCO 的 Height 修正 Rayleigh 光學厚度 (同 crefl 的 DEM 修正)
TC_SHARPEN_MAX = 5.0            # I1 比值銳化的比值上限 (原設定；satpy 為 1.5)
# NASA Rapid Response / WorldView 增強 (Gumley, Descloitres & Schmaltz 2010, Table 3):
# 反射率 0.0-1.1 線性映射至 0-255 後套非線性增強表 (表格見 satima/truecolor.py)
TC_STRETCH_MAXREFL = 1.1
TC_STRETCH_IN = np.array(tc.WV_XP)
TC_STRETCH_OUT = np.array(tc.WV_FP)
TC_NIGHT_SZA = 88.5
TC_STYLES = {1: 'Simple', 2: 'WorldView', 3: 'CIRA'}
TC_STYLE_DESC = {1: 'Simple (CREFL correction + sRGB)',
                 2: 'NASA WorldView (CREFL correction + NASA enhancement)',
                 3: 'CIRA (CREFL correction + CIRA log stretch)'}

# 共用版本 (satima/truecolor.py)，保留原名稱
csalbr = tc.csalbr
chand = tc.chand
srgb_encode = tc.srgb_encode
worldview_stretch = tc.worldview_stretch


def crefl_correct(refl, ch, sz, saz, vz, vaz, height=None):
    # NASA CREFL 大氣校正: refl 需已除以 cos(太陽天頂角); VIIRS SDR 的 Reflectance 本身即是如此, 不可再除
    # height (m): 地形高度, Rayleigh 光學厚度乘上 exp(-height/8000)
    taur0, ah2o, bh2o, ao3 = TC_COEFF[ch]
    mus = np.cos(np.deg2rad(np.clip(sz, 0., TC_MAXSOLZ)))
    muv = np.cos(np.deg2rad(np.clip(vz, 0., TC_MAXSOLZ)))
    if height is not None:
        taur = taur0 * np.exp(-np.clip(height, 0., None) / tc.SCALE_HEIGHT)
    else:
        taur = taur0
    air_mass = np.minimum(1.0 / mus + 1.0 / muv, TC_MAXAIRMASS)
    t_o3 = tc.o3_transmittance(air_mass, TC_UO3, ao3)
    t_h2o = tc.h2o_transmittance_viirs(air_mass, TC_UH2O, ah2o, bh2o)
    corr_refl = tc.crefl_core(refl, mus, muv, saz - vaz, taur, t_o3, t_h2o)
    return np.clip(corr_refl, tc.REFLMIN, tc.REFLMAX)


# ===== 單一檔案 =====
def _name_time(date, hms):
    """檔名的日期與 HHMMSSs (最後一位是十分之一秒) -> datetime。"""
    return (datetime.datetime.strptime(date + hms[:6], '%Y%m%d%H%M%S')
            + datetime.timedelta(seconds=int(hms[6]) / 10.))


class ViirsPack:
    """一個實體檔案：檔名解析的結果 (一般的檔只有一個產品；CLASS 的合併檔有多個)。"""

    def __init__(self, fullname):
        m = _NAME.match(os.path.basename(fullname))
        if not m:
            raise ValueError(tr('不認得的 VIIRS 檔名', 'unexpected VIIRS file name'))
        self.fullname = fullname
        self.sat = m['sat'].lower()
        self.satellite = SAT_NAMES.get(self.sat, self.sat.upper())
        self.time = _name_time(m['date'], m['t'])
        self.end = _name_time(m['date'], m['e'])
        if self.end < self.time:            # 檔名的日期是起始時間的日期：跨日
            self.end += datetime.timedelta(days=1)
        self.created = m['created']
        self.source = m['src']
        self.products = m['prods'].split('-')
        self.sdr = [p for p in self.products if p.startswith('SV')]
        self.geo = [p for p in self.products if p in GEO_GROUPS]


class ViirsSlot:
    """一個片段 (同衛星、同起始時間)：各頻道的檔與各種定位檔。"""

    def __init__(self, pack):
        self.sat, self.satellite, self.time, self.end = pack.sat, pack.satellite, pack.time, pack.end
        self.bands = {}             # 頻道名稱 (M15、I5、DNB) -> ViirsFile
        self.geo = {}               # 定位產品 (GMTCO …) -> 檔名
        self._newest = {}

    def add(self, pack):
        """登記一個檔所含的產品；同一產品有多個檔 (重新下載、不同產製時間) 時留最新的。"""
        self.end = max(self.end, pack.end)
        for code in pack.sdr + pack.geo:
            if code in self._newest and self._newest[code] >= pack.created:
                continue
            self._newest[code] = pack.created
            if code in GEO_GROUPS:
                self.geo[code] = pack.fullname
            else:
                try:
                    ch = channel_name(code[2:5])
                except ValueError:
                    continue
                self.bands[ch] = ViirsFile(pack.fullname, ch, self)

    def usable(self):
        """有定位檔可用的頻道 (M 需要 GMTCO / GMODO，I 需要 GITCO / GIMGO，DNB 需要 GDNBO)。"""
        return {ch for ch in self.bands
                if any(code in self.geo for code in GEO_ORDER['D' if ch == 'DNB' else ch[0]])}

    def lacking(self):
        """有資料但缺定位檔的頻道種類 -> 需要的定位產品，例如 {'I': 'GITCO / GIMGO'}。"""
        out = {}
        for ch in set(self.bands) - self.usable():
            kind = 'D' if ch == 'DNB' else ch[0]
            out[kind] = ' / '.join(GEO_ORDER[kind])
        return out


class ViirsFile:
    """某個片段的一個頻道 (資料在 fullname；定位檔由所屬的片段 slot 查)。"""

    def __init__(self, fullname, channel, slot):
        self.fullname = fullname
        self.channel = channel
        self.slot = slot
        self.sat = slot.sat
        self.satellite = slot.satellite
        self.time = slot.time
        if self.channel == 'DNB':
            self._var = 'Radiance'
        else:
            self._var = 'Reflectance' if self.channel in REFL_BANDS else 'BrightnessTemperature'

    @property
    def end(self):
        return self.slot.end

    def _geo(self):
        """(定位檔, HDF5 群組名稱)：地形校正版優先，找不到時用未做地形校正的版本。"""
        kind = 'D' if self.channel == 'DNB' else self.channel[0]
        for n, code in enumerate(GEO_ORDER[kind]):
            path = self.slot.geo.get(code)
            if path:
                if n and code not in _NOTICED:
                    _NOTICED.add(code)
                    print(tr('注意：找不到地形校正的定位檔 %s，改用未做地形校正的 %s (高山在掃描帶邊緣會偏移)。',
                             'NOTE: terrain-corrected geolocation %s not found; using %s without terrain correction '
                             '(mountains shift near the swath edge).') % (GEO_ORDER[kind][0], code))
                return path, GEO_GROUPS[code]
        raise FileNotFoundError(tr('找不到定位檔：%s %s %s (%s；須放在同一個資料夾，或使用含定位的合併檔)',
                                   'Geolocation file not found: %s %s %s (%s; put it in the same folder, or use a '
                                   'packaged file that contains the geolocation)')
                                % (self.satellite, self.time.strftime('%Y/%m/%d %H:%M:%S'), self.channel,
                                   ' / '.join(GEO_ORDER[kind])))

    @property
    def _geoname(self):
        return self._geo()[1]

    def __lt__(self, other):
        return self.time < other.time

    def extract(self, drop_fill_rows=True):
        """反射率 (0–1.5) / 亮溫 (°C) / DNB 輻射 (W cm-2 sr-1)；填充值設為 NaN。
        drop_fill_rows=True 時去除整列皆為填充值的資料列 (原作法)。"""
        if not hasattr(self, '_sdr'):
            self._sdr = h5py.File(self.fullname, 'r')
        grp = self._sdr['All_Data']['VIIRS-%s-SDR_All' % (self.channel)]
        raw = grp[self._var][:]
        fill_value = grp[self._var].fillvalue

        print('Before remove NAN Line (data):', raw.shape, fill_value)
        temp = raw.copy().astype(float)
        temp[raw == fill_value] = np.nan
        if drop_fill_rows:
            temp = temp[~np.isnan(temp).all(axis=1)]
            print('After remove NAN Line (data):', temp.shape)
        if self.channel == 'DNB':
            return temp

        factors_name = '%sFactors' % (self._var)
        if factors_name in grp:
            factors = grp[factors_name][:]
            k, b = factors[0], factors[1]
        else:                       # 部分頻道 (如 M13) 為浮點數，沒有縮放係數
            k, b = 1., 0.
        if self.channel in REFL_BANDS:
            return np.ma.masked_greater(temp * k + b, 1.5)
        return np.ma.masked_greater(temp * k + b - 273.15, 100)

    def geofile(self):
        return self._geo()[0]

    def footprint(self, step=64):
        """粗略的涵蓋範圍：每隔 step 個像元取一點的 (緯度, 經度)；讀不到時 None。"""
        try:
            path, group = self._geo()
            with h5py.File(path, 'r') as geo:
                grp = geo['All_Data'][group]
                lat = np.asarray(grp['Latitude'][::step, ::step], dtype=np.float64)
                lon = np.asarray(grp['Longitude'][::step, ::step], dtype=np.float64)
        except Exception:
            return None
        ok = (np.abs(lat) <= 90.) & (np.abs(lon) <= 360.)
        return (lat[ok], lon[ok]) if ok.any() else None

    @staticmethod
    def _read_geo(grp, name):
        arr = grp[name][:]
        if arr.dtype.kind != 'f':
            arr = arr.astype(np.float32)
        arr[arr == grp[name].fillvalue] = np.nan
        arr[arr < -999.] = np.nan       # VIIRS 浮點填充值為 -999.x 系列
        return arr

    def geocoord(self, angles=False, drop_fill_rows=True):
        """lats, lons, sz (+ saz, vz, vaz)；DNB 另設定 self.lza / moonillum / moonphase。"""
        with h5py.File(self.geofile(), 'r') as geo:
            grp = geo['All_Data'][self._geoname]
            # DNB：較新的 GDNBO 另有地形校正的經緯度，有就用 (與 M / I 頻道的 GMTCO / GITCO 一致)
            tc = self.channel == 'DNB' and 'Latitude_TC' in grp and 'Longitude_TC' in grp
            lats = self._read_geo(grp, 'Latitude_TC' if tc else 'Latitude')
            lons = self._read_geo(grp, 'Longitude_TC' if tc else 'Longitude')
            sz = self._read_geo(grp, 'SolarZenithAngle')
            print('Before remove NAN Line (lats/lons):', lats.shape, lons.shape)
            drop = (lambda a: a[~np.isnan(a).all(axis=1)]) if drop_fill_rows else (lambda a: a)
            lats, lons, sz = drop(lats), drop(lons), drop(sz)
            if drop_fill_rows:
                print('After remove NAN Line (lats/lons):', lats.shape, lons.shape, sz.shape)
            if self.channel == 'DNB':
                self.lza = drop(self._read_geo(grp, 'LunarZenithAngle'))
                self.moonillum = float(grp['MoonIllumFraction'][0])
                self.moonphase = float(grp['MoonPhaseAngle'][0])
            if angles:
                extras = [drop(self._read_geo(grp, name))
                          for name in ('SolarAzimuthAngle', 'SatelliteZenithAngle', 'SatelliteAzimuthAngle')]
                return lats, lons, sz, extras[0], extras[1], extras[2]
        return lats, lons, sz

    def height(self):
        """地形高度 (m)；無效值設 0。地理定位檔沒有 Height 時回傳 None。"""
        with h5py.File(self.geofile(), 'r') as geo:
            grp = geo['All_Data'][self._geoname]
            if 'Height' not in grp:
                return None
            h = self._read_geo(grp, 'Height')
        return np.nan_to_num(np.clip(h, 0., None), nan=0.)

    def load(self):
        """單頻道：資料與經緯度一起去除填充列 (原本兩者分別去除，列數可能不一致)。"""
        data = self.extract(drop_fill_rows=False)
        lats, lons, sz = self.geocoord(drop_fill_rows=False)
        if np.ma.getdata(data).shape[0] != lats.shape[0]:
            raise ValueError(tr('%s：資料的列數 %d 與定位的列數 %d 不同', '%s: data rows %d != geolocation rows %d')
                             % (self.fullname, data.shape[0], lats.shape[0]))
        keep = ~np.isnan(np.ma.getdata(data)).all(axis=1) & ~np.isnan(lats).all(axis=1)
        print(tr('有效的列：%d / %d', 'Valid lines: %d / %d') % (keep.sum(), keep.size))
        if self.channel == 'DNB':
            self.lza = self.lza[keep]
        rowtime = self.row_times(lats.shape[0])
        return data[keep], lats[keep], lons[keep], sz[keep], rowtime[keep]

    @staticmethod
    def fill_gaps(data):
        """masked / NaN 以沿軌前後值補上 (pandas ffill / bfill)。"""
        data = np.ma.filled(data, np.nan)
        df = pd.DataFrame(data).ffill().bfill()
        return df.values

    def plot(self, canvas, lats, lons, data, sz, variant):
        """單頻道繪圖；variant：DNB 為縮放模式，亮溫為色表。回傳色表/縮放名稱。"""
        print(np.nanmax(lats), np.nanmax(lons), np.nanmin(lats), np.nanmin(lons))
        if data.size == 0:
            return ''
        if self.channel in REFL_BANDS:
            data = np.ma.masked_greater(data, 1.5)
        else:
            data = np.ma.masked_greater(data, 100)
        if self.channel == 'DNB':
            print(np.nanmax(data), np.nanmin(data))
            dnbmode = variant
            if dnbmode == 4:
                canvas.mesh(lons, lats, np.ma.masked_invalid(data), cmap='gray',
                            norm=colors.Normalize(vmin=0, vmax=DNB_REFL_MAX, clip=True))
                return 'LunarRefl'
            if dnbmode > 1:
                canvas.mesh(lons, lats, data, cmap='gray', vmin=0., vmax=1.)
                return 'ERF-SAT' if dnbmode == 3 else 'ERF'
            canvas.mesh(lons, lats, data * 100, cmap='gray',
                        norm=colors.Normalize(vmin=0, vmax=DNB_THRESHOLD, clip=True))
            return ''
        if self.channel in REFL_BANDS:
            canvas.mesh(lons, lats, data, cmap='gray', vmin=0., vmax=1.)
            return ''
        return canvas.bt_mesh(lons, lats, data, variant)

    def row_times(self, nrows):
        """每一列的觀測時間 (UTC)：地理定位檔逐掃描帶的 MidTime (IET，TAI 微秒)，
        同一掃描帶 (M / DNB 16 列、I 32 列) 同一時間；缺掃描 (填充值為負數) 為 NaT。"""
        with h5py.File(self.geofile(), 'r') as geo:
            mid = geo['All_Data'][self._geoname]['MidTime'][:]
        return obstime.per_scan(obstime.iet(np.ma.masked_less_equal(mid, 0)), nrows)


# ===== 讀取器 =====
class ViirsReader(SatReader):
    NAME = 'VIIRS-HDF5'
    SENSOR = 'VIIRS'
    DESCRIPTION = 'S-NPP / NOAA-20 / NOAA-21 VIIRS SDR (HDF5；NOAA CLASS 的合併檔也可讀)'
    DESCRIPTION_EN = 'S-NPP / NOAA-20 / NOAA-21 VIIRS SDR (HDF5; NOAA CLASS packaged files too)'
    # SV*：SDR (含以 SDR 開頭的合併檔)；G*：定位檔與以定位開頭的合併檔。(舊的 'S*ops*h5' 會誤認 ATMS 的 SATMS_*。)
    PATTERNS = ('SV*_d*_t*_e*_b*_c*_ops.h5', 'GMTCO*_d*_t*_ops.h5', 'GMODO*_d*_t*_ops.h5', 'GITCO*_d*_t*_ops.h5',
                'GIMGO*_d*_t*_ops.h5', 'GDNBO*_d*_t*_ops.h5')
    MULTI_FILE = True
    SOURCE = 'NOAA'
    FOOTPRINT_MARGIN = 2.0      # 判斷片段是否與地圖相交時，地圖範圍外加的度數

    def make_file(self, name):
        """一個實體檔案 (ViirsPack)；頻道與定位的對應在 segments() 整理。"""
        return ViirsPack(name)

    # ---- 片段與接續 ----------------------------------------------------------
    def segments(self, packs):
        """把檔案整理成片段 (同衛星、同起始時間的各頻道與定位)，再把時間相接的片段接成一段。"""
        slots = {}
        for pack in packs:
            key = (pack.sat, pack.time)
            if key not in slots:
                slots[key] = ViirsSlot(pack)
            slots[key].add(pack)
        by_sat = {}
        for (sat, _), slot in slots.items():
            if slot.bands:                  # 只有定位、沒有 SDR 的片段不列
                by_sat.setdefault(sat, []).append(slot)
        out = []
        for items in by_sat.values():
            for members, missing in seg.chain(items, lambda x: x.time, lambda x: x.end):
                union = set().union(*(m.usable() for m in members))
                lacking = {}
                for m in members:
                    lacking.update(m.lacking())
                # 列表的產品欄：可用的頻道數；有頻道缺定位檔時標出需要哪一種 (例如「缺 GITCO」)
                product = '%d bands' % len(union)
                if lacking:
                    product += tr('，缺 %s', ', no %s') % ' '.join(v.split(' / ')[0] for v in lacking.values())
                out.append(Segment(members[0].satellite, members[0].time, max(m.end for m in members),
                                   files=members, slots=len(members), missing=missing, product=product,
                                   short=True, extra={'lacking': lacking}))
        return out

    def begin(self, segment):
        self.slots = list(segment.files)
        self.lacking = dict(segment.extra.get('lacking', {}))
        self.files = []
        if getattr(self, 'product', None) is not None:      # 第 2 段起：沿用第一段選的產品
            self._resolve()

    @staticmethod
    def _order(ch):
        return ('MID'.index(ch[0]), int(ch[1:]) if ch[1:].isdigit() else 0)

    def _bands(self):
        """(整段每個片段都有的頻道, 只有部分片段有的頻道)，依 M、I、DNB 排序。"""
        union = set().union(*(sl.usable() for sl in self.slots))
        common = set(union).intersection(*(sl.usable() for sl in self.slots))
        return sorted(common, key=self._order), sorted(union - common, key=self._order)

    def _in_range(self, slots, band):
        """只留與地圖範圍相交的片段，前後各多帶一個 (銳化、掃描交界在接縫處需要鄰居)。判斷不出來時全留。"""
        if len(slots) < 2 or self.georange is None:
            return slots
        latmin, latmax, lonmin, lonmax = self.georange
        m = self.FOOTPRINT_MARGIN
        hit = []
        for k, sl in enumerate(slots):
            f = sl.bands.get(band)
            fp = f.footprint() if f is not None else None
            if fp is None:
                return slots
            lat, lon = fp[0], swath.wrap_lon(fp[1], self.center[1])
            if ((lat > latmin - m) & (lat < latmax + m) & (lon > lonmin - m) & (lon < lonmax + m)).any():
                hit.append(k)
        if not hit:
            return slots
        keep = slots[max(0, hit[0] - 1):hit[-1] + 2]
        if len(keep) < len(slots):
            print(tr('只讀與地圖範圍相交的片段：%d / %d 個 (%s ~ %s)', 'Reading only the granules over the map: '
                     '%d of %d (%s ~ %s)') % (len(keep), len(slots), keep[0].time.strftime('%H:%M:%S'),
                                             keep[-1].end.strftime('%H:%M:%S')))
        return keep

    def _resolve(self):
        """依所選的產品，決定這一段要讀的檔 (self.files；真實色另有 self.tc_files)。"""
        if self.tcmode:
            chs = ['M3', 'M4', 'M5'] + (['I1'] if self.tcmode == 2 else [])
            slots = self._in_range(self.slots, 'M5')
            self.tc_files = dict()
            for ch in chs:
                missing = [sl for sl in slots if ch not in sl.usable()]
                if missing:
                    raise IOError(tr('真實色缺少 %s：%s %s 的片段', 'TrueColor is missing %s: granule %s %s')
                                  % (ch, missing[0].satellite, missing[0].time.strftime('%Y/%m/%d %H:%M:%S')))
                self.tc_files[ch] = [sl.bands[ch] for sl in slots]
            self.files = self.tc_files['M5']
            return
        slots = [sl for sl in self.slots if self.channel in sl.usable()]
        if len(slots) < len(self.slots):
            print(tr('注意：%d 個片段沒有 %s，略過。', 'NOTE: %d granule(s) have no %s; skipped.')
                  % (len(self.slots) - len(slots), self.channel))
        self.files = [sl.bands[self.channel] for sl in self._in_range(slots, self.channel)]

    # ---- 選單 --------------------------------------------------------------
    def select_band(self):
        """產品選單 (單選)：單頻道 (整段都有的頻道)，以及 M3 / M4 / M5 (加 I1) 都有時的真實色。"""
        common, partial = self._bands()
        if partial:
            print(tr('注意：以下頻道只有部分片段有：%s', 'NOTE: these channels exist in only some granules: %s')
                  % ' '.join(partial))
        bands = common + partial
        if getattr(self, 'lacking', None):
            print(tr('注意：有些頻道的資料沒有對應的定位檔，不列入選單 (需要 %s；須放在同一個資料夾)。',
                     'NOTE: some channels have no matching geolocation file and are not offered (need %s; it must '
                     'be in the same folder).') % tr('、', ', ').join(self.lacking.values()))
        if not bands:
            raise SystemExit(tr('這一段沒有可用的頻道。', 'No usable channel in this segment.'))
        options = [(ch, 0) for ch in bands]
        if all(b in common for b in ('M3', 'M4', 'M5')):
            options.append(('TrueColor-M', 1))
            if 'I1' in common:
                options.append(('TrueColor-I', 2))
        text = {'TrueColor-M': 'TrueColor-M (M-band only, 750m)', 'TrueColor-I': 'TrueColor-I (I1 enhanced, 375m)'}
        groups = (('Reflectance', [k for k, (n, t) in enumerate(options, 1) if not t and n in REFL_BANDS]),
                  ('Brightness temperature', [k for k, (n, t) in enumerate(options, 1) if not t and n in BT_BANDS]),
                  ('Day/Night Band', [k for k, (n, t) in enumerate(options, 1) if n == 'DNB']),
                  ('RGB', [k for k, (n, t) in enumerate(options, 1) if t]))
        print('Products:')
        ui.print_columns([(title, ['%2d.%s' % (k, text.get(options[k - 1][0], options[k - 1][0])) for k in ks])
                          for title, ks in groups if ks])
        # 只有一個產品時不詢問 (同色表選單的做法)
        k = 1 if len(options) == 1 else ui.ask_int(tr('輸入產品編號：', 'Input the product: '),
                                                   range(1, len(options) + 1))
        self.product, self.tcmode = options[k - 1]
        self.channel = None if self.tcmode else self.product
        self._resolve()

    def select_variants(self):
        if self.tcmode:
            return ui.ask_menu('Style:', TC_STYLE_DESC, tr('選擇樣式 (可複選，0 = 全部)：', 'Select style: '))
        ch = self.channel
        if ch == 'DNB':
            return ui.ask_menu('DNB scaling:', DNB_MODES, tr('選擇 DNB 縮放 (可複選，0 = 全部)：',
                                                             'Select DNB scaling: '))
        if ch in BT_BANDS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    # ---- 資料 --------------------------------------------------------------
    def prepare(self):
        if not self.files:
            raise ValueError(tr('這一段沒有 %s 的資料', 'this segment has no %s data') % self.product)
        if self.tcmode:
            self.prepare_tc()
            timeindex = self.slice_tc()
            if timeindex is None:
                return False
            self.time = self.center_time(self.rowtime)
            return True

        loaded = [f.load() for f in self.files]
        self.data = self.concat([d[0] for d in loaded])
        lats = np.concatenate([d[1] for d in loaded], axis=0)
        lons = np.concatenate([d[2] for d in loaded], axis=0)
        self.sz_full = np.concatenate([d[3] for d in loaded], axis=0)
        rowtime = np.concatenate([d[4] for d in loaded])
        if len(self.files) > 1:
            print(tr('已合併 %d 個片段：', 'Merged %d granules:') % len(self.files), lats.shape)

        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        self.time = self.center_time(rowtime)
        self.sz = self.crop(self.sz_full)
        ch = self.files[0].channel
        if ch == 'DNB':
            # DNB 縮放在整條軌道上計算 (與 satpy 相同)，之後才裁切範圍 (見 dnb_scaled)
            self.lza_full = np.concatenate([f.lza for f in self.files], axis=0)
            self._dnb_cache = dict()
        else:
            self.data = self.crop(self.data)
            if ch in BT_BANDS:
                self.report_eyetemp(self.data)
        return True

    def slice(self, georange, lats, lons):
        """裁出包住 georange 的最小矩形；回傳 timeindex [yi, yj, 總列數]，範圍內無資料回傳 None。"""
        return self.slice_swath(lats, lons)

    def dnb_scaled(self, dnbmode):
        """DNB 縮放 (原 prepare 中的 DNB 分支)；結果依模式快取。"""
        if dnbmode in self._dnb_cache:
            return self._dnb_cache[dnbmode]
        f0 = self.files[0]
        if dnbmode == 4:
            phase = np.mean([f.moonphase for f in self.files])
            e_moon, phase = lunar_irradiance_dnb(f0.time, phase)
            bw = DNB_BANDWIDTH_UM.get(f0.sat, DNB_BANDWIDTH_DEFAULT)
            print(tr('月亮相位角 %.2f 度，MT2009 DNB 月光照度 %.4e W/m2/um，DNB 頻寬 %.4f um',
                     'Lunar phase angle %.2f deg, MT2009 DNB lunar irradiance %.4e W/m2/um, DNB bandwidth %.4f um')
                  % (phase, e_moon, bw))
            out = lunar_reflectance(self.data, self.sz_full, self.lza_full, e_moon, bw)
        elif dnbmode > 1:
            moon = np.mean([f.moonillum for f in self.files]) * 0.01
            phase = np.mean([f.moonphase for f in self.files]) if DNB_PHASE_CORRECTION else None
            print(tr('DNB ERF：月亮照明比例 %.1f%%，相位角 %s，月亮天頂角 %.1f-%.1f',
                     'DNB ERF: moon illumination %.1f%%, phase %s, LZA %.1f-%.1f')
                  % (moon * 100, tr('不使用', 'off') if phase is None else '%.1f' % phase,
                     np.nanmin(self.lza_full), np.nanmax(self.lza_full)))
            out = erf_dnb(self.data, self.sz_full, self.lza_full, moon, saturation=(dnbmode == 3), moon_phase=phase)
        else:
            out = self.data
        out = self.crop(out)
        self._dnb_cache[dnbmode] = out
        return out

    def plot(self, canvas, variant):
        if self.tcmode:
            return self.plot_truecolor(canvas, variant)
        f0 = self.files[0]
        data = self.dnb_scaled(variant) if f0.channel == 'DNB' else self.data
        cm = f0.plot(canvas, self.lats, self.lons, data, self.sz, variant)
        return dict(band=f0.channel, colorbar=cm, label=f0.channel)

    # ---- 真實色 ------------------------------------------------------------
    def prepare_tc(self):
        # 先完整讀入 (不各自去除填充列)，最後再以 M-band 經緯度一起去除填充列，
        # I-band (2 倍解析度) 同步去除對應的 2 列，維持 M/I 列號對應
        self.tc_data = dict()
        for ch in self.tc_files:
            arrs = [f.extract(drop_fill_rows=False) for f in self.tc_files[ch]]
            self.tc_data[ch] = self.concat(arrs)
        geos = [f.geocoord(angles=True, drop_fill_rows=False) for f in self.tc_files['M3']]
        self.lats = np.vstack([g[0] for g in geos])
        self.lons = np.vstack([g[1] for g in geos])
        self.sz = np.vstack([g[2] for g in geos])
        self.saz = np.vstack([g[3] for g in geos])
        self.vz = np.vstack([g[4] for g in geos])
        self.vaz = np.vstack([g[5] for g in geos])
        self.rowtime = np.concatenate([f.row_times(g[0].shape[0]) for f, g in zip(self.tc_files['M3'], geos)])
        self.height = None
        if TC_USE_HEIGHT:
            heights = [f.height() for f in self.tc_files['M3']]
            if all(h is not None for h in heights):
                self.height = np.vstack(heights)
            else:
                print(tr('注意：GMTCO 裡找不到 Height (地形高度)，CREFL 以海平面計算。',
                         'WARNING: Height not found in GMTCO; CREFL uses sea level.'))
        if self.tcmode == 2:
            igeos = [f.geocoord(drop_fill_rows=False) for f in self.tc_files['I1']]
            self.ilats = np.vstack([g[0] for g in igeos])
            self.ilons = np.vstack([g[1] for g in igeos])

        keep = ~np.isnan(self.lats).all(axis=1)
        if not keep.all():
            print(tr('去除 %d 列填充列 (共 %d 列)', 'Remove %d fill lines (of %d)') % ((~keep).sum(), keep.size))
            for ch in ('M3', 'M4', 'M5'):
                self.tc_data[ch] = self.tc_data[ch][keep]
            for name in ('lats', 'lons', 'sz', 'saz', 'vz', 'vaz', 'rowtime'):
                setattr(self, name, getattr(self, name)[keep])
            if self.height is not None:
                self.height = self.height[keep]
            if self.tcmode == 2:
                ikeep = np.repeat(keep, 2)
                if self.tc_data['I1'].shape[0] == ikeep.size:
                    self.tc_data['I1'] = self.tc_data['I1'][ikeep]
                    self.ilats, self.ilons = self.ilats[ikeep], self.ilons[ikeep]
                else:
                    print(tr('注意：I 頻道的列數 (%d) 不是 M 頻道列數 (%d) 的 2 倍',
                             'WARNING: I-band rows (%d) != 2 x M-band rows (%d)')
                          % (self.tc_data['I1'].shape[0], keep.size))
        self._tc_linear = None

    def slice_tc(self):
        lons = swath.wrap_lon(self.lons, self.center[1])
        box = swath.swath_box(self.georange, self.lats, lons)
        if box is None:
            return None
        yi, yj, xi, xj = box
        xx = self.lats.shape[0]
        self.box = box
        for ch in ('M3', 'M4', 'M5'):
            self.tc_data[ch] = self.tc_data[ch][yi:yj, xi:xj]
        self.lats = self.fill_gaps(np.ma.masked_outside(self.lats[yi:yj, xi:xj], -90., 90.))
        self.lons = self.fill_gaps(np.ma.masked_invalid(lons[yi:yj, xi:xj]))
        for name in ('sz', 'saz', 'vz', 'vaz'):
            # 零星的角度填充值以鄰近值補上 (否則 CREFL 為 NaN，畫成黑點)
            setattr(self, name, self.fill_gaps(np.ma.masked_invalid(getattr(self, name)[yi:yj, xi:xj])))
        if self.height is not None:
            self.height = self.height[yi:yj, xi:xj]
        if self.tcmode == 2:
            ilons = swath.wrap_lon(self.ilons, self.center[1])
            self.tc_data['I1'] = self.tc_data['I1'][2 * yi:2 * yj, 2 * xi:2 * xj]
            self.ilats = self.fill_gaps(np.ma.masked_outside(self.ilats[2 * yi:2 * yj, 2 * xi:2 * xj], -90., 90.))
            self.ilons = self.fill_gaps(np.ma.masked_invalid(ilons[2 * yi:2 * yj, 2 * xi:2 * xj]))
        return [yi, yj, xx]

    fill_gaps = staticmethod(ViirsFile.fill_gaps)

    def truecolor_linear(self):
        """CREFL 校正後的線性反射率 RGB 與夜間判斷用的太陽天頂角 (各樣式共用，只算一次)。"""
        if self._tc_linear is not None:
            return self._tc_linear

        def clean(arr):
            # 補掉 bowtie deletion 與 fill value 造成的 NaN 條紋
            return self.fill_gaps(np.ma.masked_invalid(np.ma.filled(arr.astype(float), np.nan)))

        chans = dict()
        for ch in ('M3', 'M4', 'M5'):
            refl = crefl_correct(clean(self.tc_data[ch]), ch, self.sz, self.saz, self.vz, self.vaz, self.height)
            chans[ch] = np.clip(refl, 0., None)

        if self.tcmode == 2:
            up = lambda a: np.repeat(np.repeat(a, 2, axis=0), 2, axis=1)
            i1 = clean(self.tc_data['I1'])
            ny = min(i1.shape[0], 2 * self.sz.shape[0], self.ilats.shape[0])
            nx = min(i1.shape[1], 2 * self.sz.shape[1], self.ilats.shape[1])
            fit = lambda a: a[:ny, :nx]
            i1 = fit(i1)
            self.ilats, self.ilons = fit(self.ilats), fit(self.ilons)
            szu, sazu = fit(up(self.sz)), fit(up(self.saz))
            vzu, vazu = fit(up(self.vz)), fit(up(self.vaz))
            hu = fit(up(self.height)) if self.height is not None else None
            i1 = np.clip(crefl_correct(i1, 'I1', szu, sazu, vzu, vazu, hu), 0., None)
            # ratio pan-sharpening: I1 與 M5 同為紅光波段
            red = fit(up(chans['M5']))
            ratio = np.clip(i1 / np.clip(red, 1e-4, None), 0., TC_SHARPEN_MAX)
            r = i1
            g = fit(up(chans['M4'])) * ratio
            b = fit(up(chans['M3'])) * ratio
            szmask = szu
        else:
            r, g, b = chans['M5'], chans['M4'], chans['M3']
            szmask = self.sz
        self._tc_linear = (np.dstack([r, g, b]), szmask)
        return self._tc_linear

    def truecolor_rgb(self, style):
        rgb, szmask = self.truecolor_linear()
        if style == 2:
            out = worldview_stretch(rgb)
        elif style == 3:
            out = tc.cira_stretch(rgb)
        else:
            out = srgb_encode(np.clip(rgb / TC_STRETCH_MAXREFL, 0., 1.))
        out = np.array(out, dtype=float)
        out[np.isnan(out)] = 0.
        with np.errstate(invalid='ignore'):
            out[szmask >= TC_NIGHT_SZA] = 0.
        return np.clip(out, 0., 1.)

    def plot_truecolor(self, canvas, style):
        rgb = self.truecolor_rgb(style)
        if self.tcmode == 2:
            lats, lons = self.ilats, self.ilons
        else:
            lats, lons = self.lats, self.lons
        canvas.rgb(lons, lats, rgb)
        label = 'TrueColor-I' if self.tcmode == 2 else 'TrueColor-M'
        return dict(band=label, colorbar=TC_STYLES[style], label=label)


READER = ViirsReader
