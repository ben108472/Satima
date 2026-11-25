"""
MVISR-L1B — FY-1C / FY-1D MVISR (多通道可見光紅外掃描輻射計) HRPT 1B (NSMC 定長紀錄格式)
原始腳本：SatimaV3-FY1D-NSMC.py
格式說明：國家衛星氣象中心《風雲一號 (02) 批極軌氣象衛星數據格式》§3.1 FY-1C/D HRPT AVHRR 1B 數據格式

所需檔案
  FY1D_AVHRR_HRPT_L1_ORB_MLT_NUL_<YYYYMMDD>_<HHMM>_1100M_PJ.L1B   (FY1C_… 亦可；FY-1C 沒有樣本，未實測)

檔案結構 (每筆紀錄 28400 bytes，big-endian)
  第 1 筆 TBM 頭紀錄、第 2 筆資料頭紀錄 (衛星代碼 113 = FY-1C、114 = FY-1D，起訖時間、軌道根數)、
  第 3 筆起每條掃描線一筆 (每秒 6 條)：時間、10 個頻道的定標係數、51 個定位點的太陽 / 衛星天頂角、方位角、經緯度、
  HRPT 遙測、2048 像元 × 10 頻道的 10-bit 計數值 (每 4 bytes 放 3 個，像元內依頻道 1~10 排列)

頻道
  1.CH1 0.63  2.CH2 0.865  6.CH6 1.61  7.CH7 0.455  8.CH8 0.505  9.CH9 0.555  10.CH10 0.93   反射率
  3.CH3 3.75  4.CH4 10.8   5.CH5 12.0                                                       亮溫
  0.真實色 R = CH1 (0.58–0.68)、G = CH9 (0.53–0.58)、B = CH7 (0.43–0.48)，比照 SGLI / MODIS / OLCI 提供
    1.RGB   TOA 反射率線性 0–1 (只套用頻道間增益，不做大氣校正)
    2.CORR  頻道間增益 + 太陽天頂角校正 + Rayleigh / 臭氧校正 (CREFL 公式，satima/po/truecolor.py) + CIRA 增強
    3.WV    校正同 CORR，NASA Worldview 增強

定標 (原腳本的問題與修正)
  - 定標係數是 Fortran 的 I*4 (2,10) 陣列：記憶體中為 斜率1, 截距1, 斜率2, 截距2, … 交錯存放
    (格式說明的「(1,10) 存放斜率、(2,10) 存放截距」是 Fortran 的下標，不是存放順序)。
    原腳本把前 10 個當斜率、後 10 個當截距，配對全錯 (例如 CH4 的「斜率」其實是 CH2 的截距)。
  - 斜率 / 2^30、截距 / 2^22 之後，意義與 NOAA pre-KLM 1b 相同：
      可見光 / 近紅外：反照率 (%) = 斜率 × 計數 + 截距   (零點 = 太空計數 12~29；整軌固定，應為發射前係數)
      紅外：輻射 mW/(m2 sr cm-1) = 斜率 × 計數 + 截距   (斜率為負；每條掃描線依星上黑體更新)
  - 反射率 = 反照率 / 100，未除 cos(太陽天頂角) (同 AVHRR-NAT)。
    原腳本把計數換成的值當成輻射，再除以假設的太陽輻照度與 cos(太陽天頂角)。
  - 亮溫：T* = c2·ν / ln(1 + c1·ν³ / L)，T = (T* − A) / B。ν、A、B 取自 RTTOV 係數檔
    rtcoef_fy1_4_mvisr.dat (FY-1D) / rtcoef_fy1_3_mvisr.dat (FY-1C) 的 FILTER_FUNCTIONS (濾波函數由 NSMC 提供)。
    原腳本以 1 / 波長 (10.8、12.0 µm) 當中心波數；FY-1D CH5 的實際中心波數是 850.4 cm-1 (12.0 µm 為 833.3)，
    冷雲頂的 CH5 會偏冷約 3 K。
  - 檢驗 (5 個 FY-1D 檔，2003~2006)：亮溫 190~302 K；分裂窗差 CH4 − CH5 在暖洋面 (CH4 > 290 K) 中位數 1.5~2.3 K，
    冷雲頂 (CH4 < 215 K) 中位數約 1 K，都在正常範圍。中心波數改用 1/波長時冷雲頂約 3.5 K (不合理)。
  - 未做非線性修正：格式說明與檔案內都沒有相關係數。

真實色 (感測器衰減與大氣校正)
  - 2003~2006 年檔案的可見光係數完全相同 (發射前的值)，沒有反映衰減；2011 年的檔案係數已更新 (CH1 斜率 0.0893 -> 0.1292)。
    厚冰雲 (CH4 < −30°C) 的 TOA 反射率在可見光應接近中性 (以 CREFL 的正向模式估算：雲頂上方有 200~400 hPa 的
    Rayleigh 時 CH7/CH1 仍在 1.00 ± 1%)，實測 TOA 的 CH7/CH1：2003/09 0.95、2004/10 0.94、2005/07 0.84、
    2005/09~10 三軌 0.82~0.83、2006/04 0.77、2006/10 0.72、2011/05 (新係數) 0.85；CH9/CH1 1.00 -> 0.92~0.94。
    直接合成會偏黃，且年份越晚越嚴重。
  - 頻道間增益 (只用於真實色；單頻道仍輸出檔案的原始定標)：以整軌的厚冰雲 (CH4 < −30°C、CH1/cos(θs) > 0.5、
    θs < 75°、衛星天頂角 < 50°，於定位點取樣) 的 TOA 比值中位數，令 CH9、CH7 與 CH1 一致。同期的三軌彼此相差 < 2%。
    樣本不足 VIS_GAIN_MIN_SAMPLES 點時 (例如晨昏)，FY-1D 改用同一組定標係數的實測值依日期內插 (範圍外取端點)；
    FY-1C 沒有參考值，不修正。CH1 本身的絕對衰減無法由此得知，不修正。
  - 已知的限制 (CORR / WV)：
      厚雲校正後 B/R 約 0.9 (略帶米黃)：CREFL 以整層大氣的 Rayleigh 校正高雲，定標正確的 MODIS 也一樣。
      太陽反光在太陽天頂角大時偏卡其色：CREFL 假設朗伯地表，對鏡面反光的藍光補得不夠。
      2003~2006 年的晴空海面藍光偏低 (校正後 B < G，海面偏橄欖色)：海面的藍光對 CH7 的定標與 τR 都很敏感，
      CH7 只有標稱的波段範圍、衰減也可能不是單純的倍數。τR 改為 0.18 時這三個場景的海面 B/G 都約 1.1，
      但 2011 年 (新係數) 用目前的 τR 海面已經偏藍 (B/G 1.4)，所以不以調整 τR 處理。
  - Rayleigh 光學厚度與臭氧吸收係數：沒有可見光頻道的光譜響應函數，以格式說明的波段範圍 (方波) 並以太陽光譜
    (5778 K 黑體) 加權平均；Rayleigh 用 Hansen & Travis (1974)，臭氧用 Chappuis 帶截面的近似表。
    同樣方法重算 MODIS B1/B3/B4、VIIRS M3/M4/M5/I1 的 CREFL 係數：τR 誤差 < 4%，臭氧多在 6% 以內 (B3 10%)。
    臭氧 0.319 cm-atm (同 MODIS CREFL 氣候值)，海平面 (無地形修正)，不做水氣修正。
  - 幾何：太陽天頂角、衛星天頂角、相對方位角取自定位點，與定位相同以地心角線性內插 (衛星天頂角以有號值內插、
    方位角以 cos / sin 內插)。已與 pyorbital 的天文計算比對：太陽天頂角相差 < 0.02°；
    相對方位角 = 太陽方位角 − 衛星方位角 (0° 為背向散射，同 CREFL 的約定)。

定位
  - 51 個定位點位於第 25, 65, …, 2025 像元 (1 起算；同 NOAA LAC 1b)。由衛星天頂角反推掃描角驗證：第 1 點約 −54.0°，
    間隔 2.16° = 40 像元 × 0.0541°；以留一法比較，「51 點平均分布在 1~2048 像元」(原腳本的假設) 的誤差明顯較大。
    原腳本的假設使掃描帶邊緣偏約 24 個像元。
  - 沿掃描內插以「地心角」為座標 (掃描角 θ 對應的地心角 γ = asin((R+h)/R·sinθ) − θ；地面位置對 γ 近似線性)，
    在地心直角座標線性內插，兩端各 24 個像元外插。留一法：內部約 0.4 km，外插 40 個像元到邊緣約 1~4 km
    (定位點本身量化為 1/128° ≈ 0.9 km)。以像元為座標的三次 Lagrange 外插誤差約 4~15 km。
  - 經緯度是有號整數 (I*2)。原腳本讀成無號整數，南半球 / 西半球的掃描線會被當成無效而略過。

其他
  - 每條掃描線最後一個 4-byte 字的 3 個槽依序為 CH9、CH10、0 (實際資料；格式說明寫成 0、CH9、CH10)。
  - 各頻道邊緣有幾個像元是配準時補上的 0 或固定值 (CH5/9/10 第 1~2、CH7 第 1~4、CH3/CH8 最後 2 個像元)，設為無效。
  - 資料頭的「掃描線數」與檔案的實際筆數不符，掃描線數依檔案大小計算。品質旗標在樣本中全為 0，不使用。
  - 時間為 UTC (以定位點的太陽天頂角對照)。
"""
import datetime
import os
import re

import numpy as np

from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr
from satima.po import SatReader, obstime, swath, truecolor as tc, ui

RECORD = 28400              # 每筆紀錄的長度 (bytes)
NPIX = 2048                 # 每條掃描線的像元數
NCH = 10                    # 頻道數
NWORD = 6827                # 每條掃描線的計數值所占的 4-byte 字數 (2048 × 10 個 10-bit 值，每字 3 個)
TIE_COLS = 24 + 40 * np.arange(51)      # 定位點所在的像元 (0 起算)

# ---- 紀錄格式 (格式說明 §3.1；位元組位置為 1 起算) ---------------------------------
DATA_HEADER = np.dtype([            # 第 2 筆紀錄的前 288 bytes
    ('satellite', 'u1'),            # 1      衛星代碼 113 = FY-1C、114 = FY-1D
    ('data_type', 'u1'),            # 2
    ('start_year', '>i2'),          # 3-4    年 (2003 存成 3)
    ('start_day', '>i2'),           # 5-6    年積日
    ('start_ms', '>i4'),            # 7-10   當日毫秒
    ('scan_lines', '>u2'),          # 11-12  (與實際筆數不符，不使用)
    ('end_year', '>i2'),            # 13-14
    ('end_day', '>i2'),             # 15-16
    ('end_ms', '>i4'),              # 17-20
    ('spare1', 'V178'),             # 21-198
    ('orbit', '>i2'),               # 199-200 軌道序號
    ('epoch', '>i2', (6,)),         # 201-212 曆元 年 月 日 時 分 秒×100
    ('semi_major_axis', '>i4'),     # 213-216 ×1000 km
    ('eccentricity', '>i4'),        # 217-220 ×1e8
    ('inclination', '>i4'),         # 221-224 ×1e6 度
    ('raan', '>i4'),                # 225-228
    ('perigee', '>i4'),             # 229-232
    ('mean_anomaly', '>i4'),        # 233-236
    ('period', '>i4'),              # 237-240
    ('orbit_count', '>i2'),         # 241-242
    ('ascending', '>i2'),           # 243-244 升降軌標識
    ('attitude', '>i4', (3,)),      # 245-256 姿態角 ×1e6 度
    ('corners', '>i4', (8,)),       # 257-288 四個角點的經緯度 ×1e4 度
])

SCAN_LINE = np.dtype([              # 第 3 筆紀錄起，每條掃描線一筆
    ('line', '>u2'),                # 1-2     掃描線序號
    ('year', '>i2'),                # 3-4
    ('day', '>i2'),                 # 5-6
    ('ms', '>i4'),                  # 7-10    當日毫秒 (UTC)
    ('quality', '>u2'),             # 11-12
    ('spare1', 'V4'),               # 13-16
    ('cal', '>i4', (NCH, 2)),       # 17-96   Fortran (2,10)：[頻道, 0] = 斜率 ×2^30、[頻道, 1] = 截距 ×2^22
    ('sun_zenith', '>i2', (51,)),   # 97-198  ×128 度
    ('sat_zenith', '>i2', (51,)),   # 199-300
    ('rel_azimuth', '>i2', (51,)),  # 301-402
    ('latlon', '>i2', (51, 2)),     # 403-606 Fortran (2,51)：[點, 0] = 緯度、[點, 1] = 經度，×128 度
    ('spare2', 'V2'),               # 607-608
    ('telemetry', 'V300'),          # 609-908 HRPT 幀遙測
    ('spare3', 'V92'),              # 909-1000
    ('counts', '>u4', (NWORD,)),    # 1001-28308
    ('spare4', 'V92'),              # 28309-28400
])
assert DATA_HEADER.itemsize == 288 and SCAN_LINE.itemsize == RECORD

SATELLITES = {113: 'FY-1C', 114: 'FY-1D'}
NAME_PREFIX = {'FY1C': 'FY-1C', 'FY1D': 'FY-1D'}

BAND_MENU = ('0.TrueColor(CH1/CH9/CH7) 1.CH1(0.63) 2.CH2(0.865) 3.CH3(3.75) 4.CH4(10.8) 5.CH5(12.0) '
             '6.CH6(1.61) 7.CH7(0.455) 8.CH8(0.505) 9.CH9(0.555) 10.CH10(0.93)')
BT_BANDS = (3, 4, 5)

# ---- 真實色 ----------------------------------------------------------------
TRUE_COLOR = 0                      # 頻道選單的真實色編號
TC_BANDS = (1, 9, 7)                # R, G, B
# Rayleigh 光學厚度 τR (海平面)、臭氧吸收係數 aO3 (cm-atm-1)：波段範圍 (方波) 以太陽光譜加權平均 (見檔頭說明)
TC_COEFF = {1: (0.05769, 0.09120),  # CH1 0.58–0.68 µm
            9: (0.09453, 0.09113),  # CH9 0.53–0.58 µm
            7: (0.21328, 0.00537)}  # CH7 0.43–0.48 µm
TC_UO3 = 0.319                      # 臭氧氣候值 (cm-atm，同 MODIS CREFL)
TC_MODES = {1: 'RGB', 2: 'CORR', 3: 'WV'}
TC_MODE_DESC = {1: 'RGB (raw composite)', 2: 'CORR (sunz + Rayleigh/O3 + CIRA stretch)',
                3: 'WV (NASA Worldview style)'}
TC_DESC = {'RGB': 'RGB', 'CORR': 'CORR', 'WV': 'WORLDVIEW'}

# 頻道間增益 (相對於 CH1)：厚冰雲的取樣條件與最少點數
VIS_GAIN_BT_MAX = -30.              # CH4 亮溫上限 (°C)
VIS_GAIN_NREFL_MIN = 0.5            # CH1 / cos(太陽天頂角) 下限
VIS_GAIN_SZA_MAX = 75.
VIS_GAIN_VZA_MAX = 50.
VIS_GAIN_MIN_SAMPLES = 100
# 樣本不足時的備援 (FY-1D)：實測的厚冰雲 TOA 比值 (日期, CH1 斜率的原始整數值, CH9/CH1, CH7/CH1)。
# NSMC 在 2006/10 ~ 2011/05 之間更新過可見光係數 (CH1 斜率 0.0893 -> 0.1292)，只用同一組係數的值依日期內插 (範圍外取端點)
VIS_RATIO_FY1D = (('2003-09-01', 95885144, 1.004, 0.953), ('2004-10-07', 95885144, 0.989, 0.936),
                  ('2005-07-16', 95885144, 0.974, 0.838), ('2005-09-30', 95885144, 0.946, 0.831),
                  ('2006-04-28', 95885144, 0.920, 0.772), ('2006-10-28', 95885144, 0.939, 0.720),
                  ('2011-05-26', 138727440, 0.933, 0.850))

# 亮溫換算 (RTTOV 係數檔 FUNDAMENTAL_CONSTANTS、FILTER_FUNCTIONS)：中心波數 ν (cm-1)、頻帶修正 T* = A + B·T
C1 = 1.191042972E-05        # mW/(m2 sr cm-4)
C2 = 1.438776878            # cm K
PLANCK = {
    'FY-1C': {3: (2691.459674, 3.638547839, 0.9946753160),      # rtcoef_fy1_3_mvisr.dat
              4: (912.0737561, 0.8457617504, 0.9969959516),
              5: (829.8681502, 0.7261620749, 0.9972656285)},
    'FY-1D': {3: (2699.586211, 3.389442697, 0.9952749930),      # rtcoef_fy1_4_mvisr.dat
              4: (925.5527421, 0.7551803441, 0.9973879923),
              5: (850.4203747, 0.7725137077, 0.9971325074)},
}

# 配準時補上的固定值 (計數 0 另外一律視為無效)：頻道 -> 像元 (0 起算)
EDGE_FILL = {3: (2046, 2047), 8: (2046, 2047)}

# 掃描幾何 (沿掃描內插用)：每像元的掃描角、衛星高度 (格式說明：掃描 ±55.4°、高度 863 km)
SCAN_STEP_DEG = 0.0541
ALTITUDE_KM = 863.
EARTH_KM = 6371.
MARGIN_DEG = 3.             # 讀檔前依定位點預選掃描線時，範圍外多帶的緯度 (經度依緯度放大)


def _year(y):
    y = int(y)
    return y + 2000 if y < 50 else (y + 1900 if y < 100 else y)


def _header_time(year, day, ms):
    """資料頭的時間 (年、年積日、當日毫秒) -> datetime；無效時 None。"""
    try:
        if not (1 <= int(day) <= 366 and 0 <= int(ms) < 86_400_000):
            return None
        return (datetime.datetime(_year(year), 1, 1) + datetime.timedelta(days=int(day) - 1, milliseconds=int(ms)))
    except (ValueError, OverflowError):
        return None


def earth_angle(cols):
    """像元 (0 起算) -> 地心角 γ (弧度)：掃描角 θ = (像元 − 1023.5) × 0.0541°，γ = asin((R+h)/R·sinθ) − θ。"""
    th = np.deg2rad((np.asarray(cols, dtype=np.float64) - (NPIX - 1) / 2.) * SCAN_STEP_DEG)
    return np.arcsin((EARTH_KM + ALTITUDE_KM) / EARTH_KM * np.sin(th)) - th


_G_TIE = earth_angle(TIE_COLS)
_G_PIX = earth_angle(np.arange(NPIX))
_SEG = np.clip(np.searchsorted(_G_TIE, _G_PIX, side='right') - 1, 0, TIE_COLS.size - 2)
_W = (_G_PIX - _G_TIE[_SEG]) / (_G_TIE[_SEG + 1] - _G_TIE[_SEG])     # 兩端 < 0 或 > 1 = 外插


def interp_tiepoints(lats, lons, chunk=1024):
    """定位點 (列, 51) -> 全解析度 (列, 2048) 的 (lats, lons)，float32；經度 -180~180。
    以地心角為座標，在地心直角座標線性內插 (兩端外插)。無效的定位點整列為 NaN。"""
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)
    out_lat = np.full((lats.shape[0], NPIX), np.nan, dtype=np.float32)
    out_lon = np.full_like(out_lat, np.nan)
    for r0 in range(0, lats.shape[0], chunk):
        la, lo = np.deg2rad(lats[r0:r0 + chunk]), np.deg2rad(lons[r0:r0 + chunk])
        xyz = [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
        x, y, z = (v[:, _SEG] * (1. - _W) + v[:, _SEG + 1] * _W for v in xyz)
        out_lat[r0:r0 + chunk] = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
        out_lon[r0:r0 + chunk] = np.rad2deg(np.arctan2(y, x))
    return out_lat, out_lon


def interp_tiefield(field):
    """定位點的純量場 (列, 51) -> (列, 2048) float32；與經緯度相同以地心角為座標線性內插 (兩端外插)。"""
    field = np.asarray(field, dtype=np.float64)
    return (field[:, _SEG] * (1. - _W) + field[:, _SEG + 1] * _W).astype(np.float32)


def unpack_counts(words, band, cols=None):
    """(列, 6827) 的 4-byte 字 -> 某頻道 (1~10) 的 (列, 像元) 10-bit 計數值；cols = 像元 (預設全部 2048 個)。
    計數值依 像元 0 的 CH1~CH10、像元 1 的 CH1~CH10 … 連續排列，每字 3 個 (最前面 2 bits 為 0)。"""
    cols = np.arange(NPIX) if cols is None else np.asarray(cols)
    i = cols * NCH + (band - 1)
    w = np.asarray(words)[:, i // 3].astype(np.uint32)
    return ((w >> ((2 - i % 3) * 10).astype(np.uint32)) & 0x3FF).astype(np.uint16)


def make_truecolor(stack, mode):
    """真實色 stack (列, 行, 6) = R、G、B 的 TOA 反射率 (已乘頻道間增益)、太陽天頂角、衛星天頂角、相對方位角
    -> 0–1 RGB (夜間 NaN)。
    RGB : TOA 反射率線性 0–1
    CORR: 太陽天頂角校正 + Rayleigh / 臭氧校正 (CREFL 公式) + CIRA 增強
    WV  : 校正同 CORR，NASA Worldview 增強"""
    toa = stack[..., :3].astype(np.float64)
    sza = stack[..., 3].astype(np.float64)
    night = ~(sza < 90.)
    if mode == 'RGB':
        rgb = np.clip(toa, 0., 1.)
    else:
        refl = toa * tc.sunz_corr(sza)[..., None]
        with np.errstate(invalid='ignore', divide='ignore', over='ignore'):
            mus = np.cos(np.deg2rad(sza))
            mus = np.where(mus > 0, mus, np.nan)
            muv = np.cos(np.deg2rad(stack[..., 4].astype(np.float64)))
            phi = stack[..., 5].astype(np.float64)
            air_mass = np.minimum(1. / mus + 1. / muv, tc.MAXAIRMASS)
            corr = np.empty_like(refl)
            for k, band in enumerate(TC_BANDS):
                taur, ao3 = TC_COEFF[band]
                t_o3 = tc.o3_transmittance(air_mass, TC_UO3, ao3)
                corr[..., k] = tc.crefl_core(refl[..., k], mus, muv, phi, taur, t_o3)
        corr = np.clip(corr, tc.REFLMIN, tc.REFLMAX)
        rgb = tc.cira_stretch(corr) if mode == 'CORR' else tc.worldview_stretch(corr)
    rgb = np.array(rgb, dtype=np.float32)
    rgb[night] = np.nan
    return rgb


def brightness_temperature(radiance, nu, a, b):
    """輻射 mW/(m2 sr cm-1) -> 亮溫 (K)：T* = c2·ν / ln(1 + c1·ν³/L)，T = (T* − A) / B；L <= 0 為 NaN。"""
    with np.errstate(invalid='ignore', divide='ignore'):
        rad = np.where(radiance > 0, radiance, np.nan)
        tstar = C2 * nu / np.log1p(C1 * nu ** 3 / rad)
    return (tstar - a) / b


# ===== 單一檔案 =====
class MvisrFile:

    def __init__(self, fullname):
        self.fullname = fullname
        size = os.path.getsize(fullname)
        if size % RECORD or size < 3 * RECORD:
            raise ValueError(tr('檔案大小 %d 不是 %d 的整數倍 (或少於 3 筆紀錄)',
                                'file size %d is not a multiple of %d (or fewer than 3 records)') % (size, RECORD))
        self.nlines = size // RECORD - 2
        with open(fullname, 'rb') as f:
            f.seek(RECORD)
            hdr = np.frombuffer(f.read(DATA_HEADER.itemsize), dtype=DATA_HEADER)[0]
        base = os.path.basename(fullname)
        self.satellite = SATELLITES.get(int(hdr['satellite'])) or NAME_PREFIX.get(base[:4].upper(), 'FY-1')
        self.time = _header_time(hdr['start_year'], hdr['start_day'], hdr['start_ms'])
        self.end = _header_time(hdr['end_year'], hdr['end_day'], hdr['end_ms'])
        if self.time is None:               # 資料頭的時間無效：用檔名的 <YYYYMMDD>_<HHMM>
            m = re.search(r'_(\d{8})_(\d{4})_', base)
            if not m:
                raise ValueError(tr('無法取得觀測時間', 'cannot determine the observation time'))
            self.time = datetime.datetime.strptime(m.group(1) + m.group(2), '%Y%m%d%H%M')
            self.end = None
        if self.end is not None and self.end < self.time:
            self.end = None
        self.orbit = int(hdr['orbit'])
        self._nav = None
        self._gains = None

    def __lt__(self, other):
        return self.time < other.time

    def _records(self):
        return np.memmap(self.fullname, dtype=SCAN_LINE, mode='r', offset=2 * RECORD, shape=(self.nlines,))

    def navigation(self):
        """全部掃描線的定位點經緯度 (列, 51)、觀測時間 (UTC datetime64)、有效列；讀一次後留著。"""
        if self._nav is None:
            rec = self._records()
            ll = np.array(rec['latlon'], dtype=np.float64) / 128.
            year = np.array(rec['year'], dtype=np.int64)
            day = np.array(rec['day'], dtype=np.int64)
            ms = np.array(rec['ms'], dtype=np.int64)
            del rec
            lat, lon = ll[..., 0], ll[..., 1]
            good = ((np.abs(lat) <= 90.).all(axis=1) & (np.abs(lon) <= 180.).all(axis=1)
                    & ~((lat == 0.) & (lon == 0.)).all(axis=1))
            year = np.where(year < 50, year + 2000, np.where(year < 100, year + 1900, year))
            tgood = (year > 1900) & (day >= 1) & (day <= 366) & (ms >= 0) & (ms < 86_400_000)
            jan1 = (np.where(tgood, year, 2000) - 1970).astype('datetime64[Y]').astype(obstime.US)
            times = jan1 + ((day - 1) * 86_400_000 + ms).astype('timedelta64[ms]')
            times[~tgood] = obstime.NAT
            lat[~good], lon[~good] = np.nan, np.nan
            self._nav = dict(lat=lat, lon=lon, time=times, good=good)
            bad = int((~good).sum())
            if bad:
                print(tr('%s：%d 條掃描線的定位無效', '%s: %d scanlines with invalid navigation') % (self.fullname, bad))
        return self._nav

    def rows_in(self, georange, clon):
        """與地圖範圍可能相交的掃描線 (起, 訖)；沒有時 None。依定位點判斷，範圍外多帶 MARGIN_DEG。"""
        nav = self.navigation()
        latmin, latmax, lonmin, lonmax = georange
        lat = nav['lat']
        lon = swath.wrap_lon(nav['lon'], clon)
        mlon = MARGIN_DEG / max(np.cos(np.deg2rad(min(max(abs(latmin), abs(latmax)) + MARGIN_DEG, 80.))), 0.1)
        with np.errstate(invalid='ignore'):
            inside = ((lat > latmin - MARGIN_DEG) & (lat < latmax + MARGIN_DEG)
                      & (lon > lonmin - mlon) & (lon < lonmax + mlon))
        rows = np.flatnonzero(inside.any(axis=1))
        if rows.size == 0:
            return None
        return max(int(rows[0]) - 2, 0), min(int(rows[-1]) + 3, self.nlines)

    def extract(self, band, rows=None, cols=None):
        """某頻道的反射率 (0~1) 或亮溫 (°C)，(列, 像元) float32；無效為 NaN。
        rows = (起, 訖)，預設全部；cols = 像元 (預設全部 2048 個)。"""
        r0, r1 = rows or (0, self.nlines)
        cols = np.arange(NPIX) if cols is None else np.asarray(cols)
        rec = self._records()[r0:r1]
        counts = unpack_counts(rec['counts'], band, cols)
        cal = np.array(rec['cal'][:, band - 1, :], dtype=np.float64)
        del rec
        slope = cal[:, 0] / 2. ** 30
        intercept = cal[:, 1] / 2. ** 22
        value = counts * slope[:, None] + intercept[:, None]

        bad = counts == 0
        bad[:, np.isin(cols, EDGE_FILL.get(band, ()))] = True
        bad |= (slope == 0.)[:, None]
        bad |= ~self.navigation()['good'][r0:r1, None]

        if band in BT_BANDS:
            nu, a, b = PLANCK.get(self.satellite, PLANCK['FY-1D'])[band]
            out = brightness_temperature(value, nu, a, b) - 273.15
        else:
            out = value / 100.
        out[bad] = np.nan
        return out.astype(np.float32)

    def geocoord(self, rows=None):
        r0, r1 = rows or (0, self.nlines)
        nav = self.navigation()
        return interp_tiepoints(nav['lat'][r0:r1], nav['lon'][r0:r1])

    def row_times(self, rows=None):
        """每條掃描線的觀測時間 (UTC datetime64[us]；無效為 NaT)。"""
        r0, r1 = rows or (0, self.nlines)
        return self.navigation()['time'][r0:r1]

    def tie_angles(self, rows=None):
        """定位點的 (太陽天頂角, 衛星天頂角, 相對方位角) (列, 51)，度。"""
        r0, r1 = rows or (0, self.nlines)
        rec = self._records()[r0:r1]
        out = tuple(np.array(rec[k], dtype=np.float64) / 128. for k in ('sun_zenith', 'sat_zenith', 'rel_azimuth'))
        del rec
        return out

    def angles(self, rows=None):
        """全解析度的 (太陽天頂角, 衛星天頂角, 相對方位角)，(列, 2048) float32，度。
        衛星天頂角在星下點兩側以正負號區分後內插再取絕對值；相對方位角在星下點跳 180°，以 cos / sin 內插。"""
        sza, vza, raa = self.tie_angles(rows)
        side = np.where(np.arange(TIE_COLS.size) < np.argmin(vza, axis=1)[:, None], -1., 1.)
        vza = np.abs(interp_tiefield(vza * side))
        r = np.deg2rad(raa)
        raa = np.rad2deg(np.arctan2(interp_tiefield(np.sin(r)), interp_tiefield(np.cos(r))))
        bad = ~self.navigation()['good'][slice(*(rows or (0, self.nlines)))]
        out = [interp_tiefield(sza), vza.astype(np.float32), raa.astype(np.float32)]
        for a in out:
            a[bad] = np.nan
        return out

    def vis_gains(self):
        """真實色的頻道間增益 {頻道: 增益} (CH1 = 1)：令整軌厚冰雲的 TOA 反射率 CH9、CH7 與 CH1 一致。
        取樣在定位點；樣本不足時 FY-1D 用 VIS_RATIO_FY1D 依日期內插，其他衛星不修正。讀一次後留著。"""
        if self._gains is not None:
            return self._gains
        sza, vza, _ = self.tie_angles()
        refl = {b: self.extract(b, cols=TIE_COLS).astype(np.float64) for b in TC_BANDS}
        bt = self.extract(4, cols=TIE_COLS)
        with np.errstate(invalid='ignore', divide='ignore'):
            nrefl = refl[1] / np.cos(np.deg2rad(sza))
            sel = ((bt < VIS_GAIN_BT_MAX) & (nrefl > VIS_GAIN_NREFL_MIN) & (sza < VIS_GAIN_SZA_MAX)
                   & (vza < VIS_GAIN_VZA_MAX) & np.isfinite(refl[9]) & np.isfinite(refl[7]) & (refl[9] > 0) & (refl[7] > 0))
        n = int(sel.sum())
        rec = self._records()
        slopes = np.array(rec['cal'][:, 0, 0])
        del rec
        slope1 = int(np.median(slopes[slopes != 0])) if np.any(slopes != 0) else 0
        table = [r for r in VIS_RATIO_FY1D if r[1] == slope1] if self.satellite == 'FY-1D' else []
        if n >= VIS_GAIN_MIN_SAMPLES:
            ratio = {b: float(np.median(refl[b][sel] / refl[1][sel])) for b in (9, 7)}
            how = tr('整軌厚冰雲 %d 點', 'thick ice cloud in this pass, %d samples') % n
        elif table:
            days = np.array([np.datetime64(r[0], 'D').astype(np.int64) for r in table], dtype=float)
            t = np.datetime64(self.time.date(), 'D').astype(np.int64)
            ratio = {9: float(np.interp(t, days, [r[2] for r in table])),
                     7: float(np.interp(t, days, [r[3] for r in table]))}
            how = tr('厚冰雲只有 %d 點，改用同一組定標係數的實測值依日期內插',
                     'only %d thick ice cloud samples; using measurements with the same calibration set, '
                     'interpolated by date') % n
        else:
            ratio = {9: 1., 7: 1.}
            how = tr('厚冰雲只有 %d 點，沒有參考值，不修正',
                     'only %d thick ice cloud samples; no reference values, not corrected') % n
        self._gains = {1: 1., 9: 1. / ratio[9], 7: 1. / ratio[7]}
        print(tr('頻道間增益 (%s)：CH9 x%.3f  CH7 x%.3f', 'Inter-channel gains (%s): CH9 x%.3f  CH7 x%.3f')
              % (how, self._gains[9], self._gains[7]))
        return self._gains

    def truecolor_stack(self, rows=None):
        """(列, 2048, 6)：R、G、B 的 TOA 反射率 (已乘頻道間增益)、太陽天頂角、衛星天頂角、相對方位角。"""
        gains = self.vis_gains()
        layers = [self.extract(b, rows) * np.float32(gains[b]) for b in TC_BANDS]
        layers += self.angles(rows)
        return np.stack(layers, axis=-1)


# ===== 讀取器 =====
class MvisrReader(SatReader):
    NAME = 'MVISR-L1B'
    SENSOR = 'MVISR'
    DESCRIPTION = 'FY-1C/1D MVISR HRPT 1B (NSMC 格式，.L1B)'
    DESCRIPTION_EN = 'FY-1C/1D MVISR HRPT 1B (NSMC format, .L1B)'
    SOURCE = 'NSMC'
    PATTERNS = ('FY1?_AVHRR_HRPT_L1_*.L1B',)

    def make_file(self, name):
        return MvisrFile(name)

    def listing(self, f):
        return (f.satellite, '10 bands', f.time)

    def select_band(self):
        for f in self.files:
            print(f.fullname)
        print(BAND_MENU)
        self.band = ui.ask_int(tr('輸入頻道編號：', 'Input the band: '), range(0, NCH + 1))

    @property
    def is_tc(self):
        return self.band == TRUE_COLOR

    def select_variants(self):
        if self.is_tc:
            return [TC_MODES[k] for k in ui.ask_menu('', TC_MODE_DESC, tr('選擇真實色模式 (可複選，0 = 全部)：',
                                                                          'Select true color mode: '))]
        if self.band in BT_BANDS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare(self):
        datas, latss, lonss, times = [], [], [], []
        for f in self.files:
            rows = f.rows_in(self.georange, self.center[1])
            if rows is None:
                continue
            datas.append(f.truecolor_stack(rows) if self.is_tc else f.extract(self.band, rows))
            lats, lons = f.geocoord(rows)
            latss.append(lats)
            lonss.append(lons)
            times.append(f.row_times(rows))
        if not datas:
            return False
        lats, lons = self.concat(latss), self.concat(lonss)
        if self.slice_swath(lats, lons) is None:
            return False
        self.data = self.crop(self.concat(datas))
        self.time = self.center_time(np.concatenate(times))
        if self.time is None:
            self.time = self.files[0].time
        if self.band in BT_BANDS:
            self.report_eyetemp(self.data)
        return True

    def plot(self, canvas, variant):
        if self.data.size == 0:
            return None
        if self.is_tc:
            canvas.rgb(self.lons, self.lats, make_truecolor(self.data, variant))
            return dict(band='TrueColor', colorbar=TC_DESC[variant], label='TRUE COLOR %s' % TC_DESC[variant])
        if variant is None:
            canvas.mesh(self.lons, self.lats, self.data, cmap='gray', vmin=0., vmax=1.)
            cm = ''
        else:
            cm = canvas.bt_mesh(self.lons, self.lats, self.data, variant)
        return dict(band='B%02d' % self.band, colorbar=cm, label='BAND%02d' % self.band)


READER = MvisrReader
