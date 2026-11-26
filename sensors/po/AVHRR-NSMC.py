"""
AVHRR-NSMC — NOAA-15~19 AVHRR/3 HRPT 1B (國家衛星氣象中心 NSMC 格式，.L1B)
原始腳本：SatimaV3-NOAA-NSMC_accnumpy.py (另有 SatimaV3-NOAA-NSMC.py、SatimaV3-NOAA-NSMC_numpy.py)
格式說明：國家衛星氣象中心《NOAA 系列極軌氣象衛星數據格式》§3.2 NOAA_K/L/M/N 衛星 1B 數據格式
  https://img.nsmc.org.cn/PORTAL/NSMC/DATASERVICE/DataFormat/Other/noaa_1b_format.pdf

所需檔案
  NA<衛星>_AVHRR_HRPT_L1_ORB_MLT_NUL_<YYYYMMDD>_<HHMM>_1100M_PJ.L1B   (NA16、NA17、NA18 已實測；NA15、NA19 未實測)

檔案結構 (每筆 22016 bytes，big-endian；與 NOAA KLM 1b 的 16-bit 格式相近，但部分欄位與 KLM 不同)
  第 1 筆表頭、其後每條掃描線一筆：時間、品質旗標、可見光 / 紅外定標係數、51 個定位點的角度與經緯度、
  HRPT 遙測、2048 像元 × 5 通道的 16-bit 計數值 (像元內依通道 1~5 排列)。
  表頭的紀錄數與實際筆數可能差幾筆，筆數依檔案大小計算。

頻道
  1.CH1 0.63  2.CH2 0.86  3.CH3A 1.61   反射率 (反照率 / 100，未除 cos(太陽天頂角)，同 AVHRR-NAT)
  4.CH3B 3.74 5.CH4 10.8  6.CH5 12.0    亮溫
  CH3A / CH3B 共用通道 3 的資料位置，依掃描線切換 (NOAA-17 白天 3A；NOAA-16、18 一直是 3B)

定標 (實際檔案與 NSMC 格式說明對照)
  - 可見光 (CH1、CH2、CH3A)：反照率 (%) = 斜率1 × C + 截距1 (C <= 交點)，否則 斜率2 × C + 截距2。
    縮放為斜率 ×10^10、截距 ×10^7 (NSMC 格式說明；NOAA KLM 為 ×10^7、×10^6)。
    交點欄位在 2008 年以前的檔案是 0，改由兩條直線的交點計算 (約 498 count)。
    NOAA-17 CH3A 的斜率2 (約 0.222) 超過 int32 能存的 0.2147，檔案中為 2147483647，改由交點處連續的條件反推。
  - 紅外 (CH3B、CH4、CH5)：線性輻射 = 係數3 + 係數2 × C + 係數1 × C²，各 ×10^6 (NSMC 的順序與 KLM 相反，
    樣本中係數1 一律為 0)。以典型太空計數 (約 990) 代入，CH4 / CH5 約為 −5.2 / −2.1 (NOAA-18 的 Ns = −5.53 / −2.22)、
    CH3B 約為 0，表示線性定標已依 NOAA KLM 的程序含太空輻射 Ns，所以 CH4、CH5 再做非線性修正
    N = b0 + (1 + b1)·N_lin + b2·N_lin² (NOAA KLM User's Guide 附錄 D 的係數)。
  - 亮溫：T* = c2·ν / ln(1 + c1·ν³ / N)，T = (T* − A) / B (NOAA KLM 的約定；ν、A、B 與 RTTOV
    rtcoef_noaa_18_avhrr.dat 的 FILTER_FUNCTIONS 相同，RTTOV 的定義為 T* = A + B·T)。
    表頭的輻射換算常數不可靠 (NOAA-18 的檔填的是 NOAA-16 的太陽輻照度與 NOAA-17 的中心波數)，不使用。
  - 通道 3 的 3A / 3B：掃描線 bit 域 bit 0 的意義在不同檔案相反 (2008 年以前照 NSMC 格式說明 0 = 3A，2010 年以後照 KLM
    0 = 3B)，同一檔內也有不一致 (NOAA-17 2005/10/01 前 95 條)，不使用。改由 3B 定標係數的線性項判斷：
    負 = 3B (熱紅外)，正 = 3A (此時該欄位存的是 3A 的反照率係數)；以 11 條掃描線的中位數濾掉單獨跳動的掃描線。
    已以計數分布驗證 (3A 的暗處約 40 count；3B 在夜間 770~990 count)。
  - 定標係數的品質：個別掃描線的係數會跳動 (CH4 在 600 count 的輻射跳 1~4 mW/(m2 sr cm-1)，約 1~4 K)，
    也有係數為 2147483647 或符號錯誤的掃描線。紅外係數與 51 條掃描線的移動中位數相差超過門檻時改用中位數；
    無效的掃描線以前後的有效值補上。整檔都無效 (例如 NOAA-17 2010/10/28，AVHRR 故障期間) 時該通道為無效。

資料受損 (接收或儀器造成的掃描線錯位、位元錯誤)
  - 品質旗標沒有標出來。NOAA-16 2005/09/28 整軌、2005/07/16 北緯 32 度以北約 23% 的像元是錯位的亂碼，其他檔案 < 1%。
  - 整條錯位：CH4 − CH5 不合理 (正常約 −3 ~ +12 K，SPLIT_RANGE = −5 ~ +15 K) 的像元超過 10%，或 CH4 計數與上下兩條
    掃描線中較接近者的差的中位數超過 15 count (兩個紅外通道一起偏移時 CH4 − CH5 看不出來；正常的檔 < 8 count)
    -> 整條遮蔽 (所有頻道)。錯位的掃描線不合理的比例大多在 20~100%，單純位元錯誤的掃描線在 10% 以下。
  - 單點位元錯誤 (例如周圍都是 610 count 的計數中夾一個 111)：只遮蔽 CH4 / CH5 中 CH4 − CH5 不合理的像元
    (單點雜訊會被當成颱風眼的最高亮溫)。紅外不用空間濾波，以免把只有一兩個像元的颱風眼當成雜訊。
  - 可見光 (CH1、CH2、CH3A) 的單條線段雜訊：與上下兩條都差超過 0.1 (同方向)、而上下兩條相差不到 0.05 的像元，
    以上下的平均取代。CH3B 不處理 (火點常只有一兩個像元)。
  - 連續被遮蔽的掃描線不超過 2 條時，以上下兩條內插補上 (所有頻道；內插值介於上下之間，不影響颱風眼的最高亮溫)，
    否則留空。

定位
  - 51 個定位點位於第 25、65、…、2025 像元 (1 起算；NSMC 格式說明明載)。沿掃描內插以「地心角」為座標
    (掃描角 θ = (像元 − 1023.5) × 0.0541° 對應的地心角 γ = asin((R+h)/R·sinθ) − θ)，在地心直角座標線性內插，
    兩端各 24 個像元外插 (方法同 MVISR-L1B；FY-1D 的留一法驗證內部約 0.4 km、邊緣約 1~4 km)。

與原腳本的差異 (原腳本的錯誤)
  - 亮溫的頻帶修正寫成 T = A + B·T*，NOAA KLM 為 T = (T* − A) / B。以 NOAA-18 為例，CH4 在 300 K 偏高 0.04 K、
    200 K 偏高 0.3 K；CH5 0.0~0.2 K；CH3B 在 300 K 偏高約 1.6 K
  - 可見光把反照率 (%) 當成輻射再乘 100π / F，再在繪圖時除以 100：反射率被放大 π/F × 100 倍 (NOAA-18 CH1 約 2.4 倍)，
    反照率 41% 以上就飽和。反照率本身就是 πL/F，反射率 = 反照率 / 100
  - 可見光的高低增益交點寫死為 500 count (2010 年以後的檔案有交點欄位，NOAA-18 為 501)
  - 沒有通道 3 (選單有 3，但選了會因為變數未定義而中斷)
  - 51 個定位點平均分布在 0~2047 像元 (應為第 25~2025 像元)，掃描帶邊緣約偏 24 個像元；且使用 scipy 1.14 已移除的 interp2d
  - QC 以第 2 / 98 百分位剔除係數：只套用在 CH1、CH4，且每個檔固定有 4% 的掃描線被改掉 (不論是否異常)
  - 掃描線時間都用第一條掃描線的日期 (跨日的軌道會錯)；中心時間以裁切範圍的中間列推算
  - 範圍內沒有資料時 slice 回傳 3 個值而呼叫端要 4 個 (會中斷)；裁切範圍少一列一行；eval(input())
"""
import datetime
import os
import re

import numpy as np
from scipy.ndimage import median_filter

from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr
from satima.po import SatReader, obstime, swath, ui

RECORD = 22016              # 每筆紀錄的長度 (bytes)
NPIX = 2048
TIE_COLS = 24 + 40 * np.arange(51)      # 定位點所在的像元 (0 起算)
INT32_MAX = 2 ** 31 - 1

# ---- 紀錄格式 (NSMC 格式說明 §3.2；位元組位置為 1 起算) -------------------------------
SCAN_LINE = np.dtype([
    ('line', '>u2'),                # 1-2      掃描線數
    ('year', '>u2'),                # 3-4
    ('day', '>u2'),                 # 5-6
    ('clock_drift', '>i2'),         # 7-8      (毫秒)
    ('ms', '>u4'),                  # 9-12     當日毫秒 (UTC)
    ('bitfield', '>u2'),            # 13-14    bit 15：0 = 升軌、1 = 降軌；bit 0：通道 3 選擇 (不可靠，不使用)
    ('spare1', 'V10'),              # 15-24
    ('quality', '>u4'),             # 25-28
    ('line_quality', '>u4'),        # 29-32
    ('ir_cal_quality', '>u2', (3,)),    # 33-38
    ('frame_sync_errors', '>u2'),   # 39-40
    ('spare2', 'V8'),               # 41-48
    ('vis_cal', '>i4', (3, 3, 5)),  # 49-228   [CH1/CH2/CH3A, 業務/測試/發射前, 斜率1 截距1 斜率2 截距2 交點]
    ('ir_cal', '>i4', (3, 2, 3)),   # 229-300  [CH3B/CH4/CH5, 業務/發射前, 係數1 係數2 係數3]
    ('spare3', 'V12'),              # 301-312
    ('nav_status', '>u4'),          # 313-316
    ('euler_time', '>u4'),          # 317-320
    ('euler', '>i2', (3,)),         # 321-326  ×10^3 度
    ('altitude', '>u2'),            # 327-328  ×10 km (樣本中為 0)
    ('angles', '>i2', (51, 3)),     # 329-634  太陽天頂角、衛星天頂角、相對方位角 ×10^2 度
    ('spare4', 'V6'),               # 635-640
    ('latlon', '>i4', (51, 2)),     # 641-1048 緯度、經度 ×10^4 度
    ('spare5', 'V8'),               # 1049-1056
    ('telemetry', '>u2', (104,)),   # 1057-1264 HRPT 副幀遙測
    ('counts', '>u2', (NPIX, 5)),   # 1265-21744
    ('spare6', 'V272'),             # 21745-22016 TIP 副幀頭、CPU 遙測
])
assert SCAN_LINE.itemsize == RECORD

BAND_MENU = '1.CH1(0.63) 2.CH2(0.86) 3.CH3A(1.61) 4.CH3B(3.74) 5.CH4(10.8) 6.CH5(12.0)'
BAND_LABEL = {1: '01', 2: '02', 3: '03A', 4: '03B', 5: '04', 6: '05'}     # 檔名 / 圖說的頻道名稱 (同 AVHRR-NAT)
VIS_BANDS = {1: 0, 2: 1, 3: 2}      # 頻道 -> vis_cal 的通道索引
IR_BANDS = {4: 0, 5: 1, 6: 2}       # 頻道 -> ir_cal 的通道索引
COUNT_INDEX = {1: 0, 2: 1, 3: 2, 4: 2, 5: 3, 6: 4}
BT_BANDS = (4, 5, 6)

# 衛星高度 (km；NSMC 格式說明，NOAA-18 / 19 為標稱值)：只用於沿掃描內插的地心角
ALTITUDE_KM = {'NOAA-15': 808., 'NOAA-16': 850., 'NOAA-17': 811., 'NOAA-18': 854., 'NOAA-19': 870.}
SCAN_STEP_DEG = 0.0541
EARTH_KM = 6371.
MARGIN_DEG = 3.             # 讀檔前依定位點預選掃描線時，範圍外多帶的緯度

# 亮溫換算 (NOAA KLM User's Guide 附錄 D)：通道 -> (中心波數 ν, A, B)；T = (T* − A) / B
C1 = 1.1910427E-5           # mW/(m2 sr cm-4)
C2 = 1.4387752              # cm K
PLANCK = {
    'NOAA-15': {4: (2695.9743, 1.621256, 0.998015), 5: (925.4075, 0.337810, 0.998719), 6: (839.8979, 0.304558, 0.999024)},
    'NOAA-16': {4: (2700.1148, 1.592459, 0.998147), 5: (917.2289, 0.332380, 0.998522), 6: (838.1255, 0.674623, 0.998363)},
    'NOAA-17': {4: (2669.3554, 1.702380, 0.997378), 5: (926.2947, 0.271683, 0.998794), 6: (839.8246, 0.309180, 0.999012)},
    'NOAA-18': {4: (2659.7952, 1.698704, 0.996960), 5: (928.1460, 0.436645, 0.998607), 6: (833.2532, 0.253179, 0.999057)},
    'NOAA-19': {4: (2670.0, 1.67396, 0.997364), 5: (928.9, 0.53959, 0.998534), 6: (831.9, 0.36064, 0.998913)},
}
# 非線性修正 (同上)：通道 -> (b0, b1, b2)；N = b0 + (1 + b1)·N_lin + b2·N_lin² (CH3B 不修正)
NONLINEAR = {
    'NOAA-15': {5: (4.76, -0.0932, 0.0004524), 6: (3.83, -0.0659, 0.0002811)},
    'NOAA-16': {5: (2.96, -0.05411, 0.00024532), 6: (2.25, -0.03665, 0.00014854)},
    'NOAA-17': {5: (8.22, -0.15795, 0.00075579), 6: (4.31, -0.07318, 0.00030976)},
    'NOAA-18': {5: (5.82, -0.11069, 0.00052337), 6: (2.67, -0.04360, 0.00017715)},
    'NOAA-19': {5: (5.70, -0.11187, 0.00054668), 6: (3.58, -0.05991, 0.00024985)},
}

# 紅外定標係數的檢查：移動中位數的掃描線數、與中位數的輻射差門檻 (mW/(m2 sr cm-1)，在 400 / 800 count 比較)
CAL_WINDOW = 51
CAL_TOLERANCE = {4: 0.005, 5: 0.2, 6: 0.2}
MODE_WINDOW = 11            # 3A / 3B 判斷的中位數濾波 (掃描線數)

# 資料受損的檢查 (見檔頭說明)：CH4 − CH5 (K) 超出 SPLIT_RANGE 視為不合理 (正常約 −3 ~ +12 K)。
# 整條遮蔽 (所有頻道)：不合理的像元超過 BAD_LINE_FRACTION，或 CH4 計數與上下兩條掃描線中較接近者的差的中位數
# 超過 LINE_JUMP_COUNTS (正常的檔 < 8 count)。其餘掃描線只遮蔽 CH4 / CH5 不合理的像元。
SPLIT_RANGE = (-5., 15.)
BAD_LINE_FRACTION = 0.10
LINE_JUMP_COUNTS = 15.
MAX_GAP_FILL = 2            # 連續被遮蔽的掃描線不超過幾條時，以上下兩條內插補上 (所有頻道)
STREAK_REFL = 0.1           # 可見光的單條線段雜訊：與上下兩條都差超過此值、而上下兩條相差不到一半時，以上下的平均取代


def _year_day_ms(year, day, ms):
    try:
        if not (1 <= int(day) <= 366 and 0 <= int(ms) < 86_400_000 and int(year) > 1970):
            return None
        return datetime.datetime(int(year), 1, 1) + datetime.timedelta(days=int(day) - 1, milliseconds=int(ms))
    except (ValueError, OverflowError):
        return None


def earth_angle(cols, altitude_km):
    """像元 (0 起算) -> 地心角 γ (弧度)：掃描角 θ = (像元 − 1023.5) × 0.0541°，γ = asin((R+h)/R·sinθ) − θ。"""
    th = np.deg2rad((np.asarray(cols, dtype=np.float64) - (NPIX - 1) / 2.) * SCAN_STEP_DEG)
    return np.arcsin((EARTH_KM + altitude_km) / EARTH_KM * np.sin(th)) - th


def interp_tiepoints(lats, lons, altitude_km, chunk=1024):
    """定位點 (列, 51) -> 全解析度 (列, 2048) 的 (lats, lons)，float32；經度 -180~180。
    以地心角為座標，在地心直角座標線性內插 (兩端外插)。無效的定位點整列為 NaN。"""
    g_tie = earth_angle(TIE_COLS, altitude_km)
    g_pix = earth_angle(np.arange(NPIX), altitude_km)
    seg = np.clip(np.searchsorted(g_tie, g_pix, side='right') - 1, 0, TIE_COLS.size - 2)
    w = (g_pix - g_tie[seg]) / (g_tie[seg + 1] - g_tie[seg])
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)
    out_lat = np.full((lats.shape[0], NPIX), np.nan, dtype=np.float32)
    out_lon = np.full_like(out_lat, np.nan)
    for r0 in range(0, lats.shape[0], chunk):
        la, lo = np.deg2rad(lats[r0:r0 + chunk]), np.deg2rad(lons[r0:r0 + chunk])
        xyz = [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
        x, y, z = (v[:, seg] * (1. - w) + v[:, seg + 1] * w for v in xyz)
        out_lat[r0:r0 + chunk] = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
        out_lon[r0:r0 + chunk] = np.rad2deg(np.arctan2(y, x))
    return out_lat, out_lon


def _fill_from_neighbours(values, valid):
    """無效的掃描線以前後最近的有效值線性內插 (兩端取端點)；沒有任何有效值時原樣回傳。"""
    if valid.all() or not valid.any():
        return values
    idx = np.arange(values.shape[0])
    out = values.copy()
    for k in range(values.shape[1]):
        out[~valid, k] = np.interp(idx[~valid], idx[valid], values[valid, k])
    return out


def vis_albedo(counts, cal):
    """可見光計數 -> 反照率 (%)。cal (列, 5) = 原始整數 [斜率1 ×10^10, 截距1 ×10^7, 斜率2, 截距2, 交點]。
    交點為 0 時由兩條直線的交點計算；斜率2 為 int32 上限 (溢位) 時由交點處連續反推。無效的掃描線為 NaN。"""
    cal = np.asarray(cal, dtype=np.float64)
    s1, i1 = cal[:, 0] / 1e10, cal[:, 1] / 1e7
    s2, i2 = cal[:, 2] / 1e10, cal[:, 3] / 1e7
    x = cal[:, 4].copy()
    overflow = cal[:, 2] >= INT32_MAX
    with np.errstate(invalid='ignore', divide='ignore'):
        s2 = np.where(overflow & (x > 0), (s1 * x + i1 - i2) / np.where(x > 0, x, np.nan), s2)
        x = np.where(x > 0, x, (i2 - i1) / (s1 - s2))
    ok = (s1 > 0) & (s2 > s1) & np.isfinite(x) & (x > 0) & (x < 1023)
    c = counts.astype(np.float64)
    albedo = np.where(c <= x[:, None], s1[:, None] * c + i1[:, None], s2[:, None] * c + i2[:, None])
    albedo[~ok] = np.nan
    return albedo


def fill_short_gaps(arr, lines, max_gap=MAX_GAP_FILL):
    """被遮蔽的掃描線 (lines) 連續不超過 max_gap 條、且上下都有掃描線時，逐像元以上下兩條線性內插 (就地修改)。
    內插值介於上下兩條之間，不會產生比鄰近值更極端的數值。"""
    idx = np.flatnonzero(lines)
    if idx.size == 0:
        return arr
    runs = np.split(idx, np.flatnonzero(np.diff(idx) > 1) + 1)
    for run in runs:
        a, b = run[0] - 1, run[-1] + 1
        if run.size > max_gap or a < 0 or b >= arr.shape[0]:
            continue
        w = (run - a) / float(b - a)
        arr[run] = arr[a][None, :] * (1. - w[:, None]) + arr[b][None, :] * w[:, None]
    return arr


def repair_streaks(arr, threshold=STREAK_REFL):
    """可見光的單條線段雜訊：與上下兩條都差超過 threshold (同方向)、而上下兩條相差不到一半時，以上下的平均取代 (就地修改)。"""
    if arr.shape[0] < 3:
        return arr
    up, mid, dn = arr[:-2], arr[1:-1], arr[2:]
    with np.errstate(invalid='ignore'):
        a, b = mid - up, mid - dn
        streak = (np.abs(a) > threshold) & (np.abs(b) > threshold) & (a * b > 0) & (np.abs(up - dn) < threshold / 2.)
    mid[streak] = 0.5 * (up[streak] + dn[streak])
    return arr


def brightness_temperature(radiance, nu, a, b):
    """輻射 mW/(m2 sr cm-1) -> 亮溫 (K)：T* = c2·ν / ln(1 + c1·ν³/N)，T = (T* − A) / B；N <= 0 為 NaN。"""
    with np.errstate(invalid='ignore', divide='ignore'):
        rad = np.where(radiance > 0, radiance, np.nan)
        tstar = C2 * nu / np.log1p(C1 * nu ** 3 / rad)
    return (tstar - a) / b


# ===== 單一檔案 =====
class AvhrrNsmcFile:

    def __init__(self, fullname):
        self.fullname = fullname
        base = os.path.basename(fullname)
        m = re.match(r'NA(\d\d)_', base, re.IGNORECASE)
        self.satellite = 'NOAA-%s' % m.group(1) if m else 'NOAA'
        self.channel = 'IR/VIS'
        size = os.path.getsize(fullname)
        if size % RECORD or size < 2 * RECORD:
            raise ValueError(tr('檔案大小 %d 不是 %d 的整數倍', 'file size %d is not a multiple of %d') % (size, RECORD))
        self.nlines = size // RECORD - 1
        with open(fullname, 'rb') as f:
            hdr = f.read(104)
        be = lambda off, dt: np.frombuffer(hdr[off:off + np.dtype(dt).itemsize], dt)[0]
        self.time = _year_day_ms(be(84, '>u2'), be(86, '>u2'), be(88, '>u4'))
        self.end = _year_day_ms(be(96, '>u2'), be(98, '>u2'), be(100, '>u4'))
        if self.time is None:               # 表頭的時間無效：用檔名的 <YYYYMMDD>_<HHMM>
            m = re.search(r'_(\d{8})_(\d{4})_', base)
            if not m:
                raise ValueError(tr('無法取得觀測時間', 'cannot determine the observation time'))
            self.time = datetime.datetime.strptime(m.group(1) + m.group(2), '%Y%m%d%H%M')
            self.end = None
        if self.end is not None and self.end < self.time:
            self.end = None
        self._nav = None
        self._cal = {}

    def __lt__(self, other):
        return self.time < other.time

    def _records(self):
        return np.memmap(self.fullname, dtype=SCAN_LINE, mode='r', offset=RECORD, shape=(self.nlines,))

    def navigation(self):
        """全部掃描線的定位點經緯度 (列, 51)、觀測時間 (UTC datetime64)、有效列；讀一次後留著。"""
        if self._nav is None:
            rec = self._records()
            ll = np.array(rec['latlon'], dtype=np.float64) / 1e4
            year = np.array(rec['year'], dtype=np.int64)
            day = np.array(rec['day'], dtype=np.int64)
            ms = np.array(rec['ms'], dtype=np.int64)
            del rec
            lat, lon = ll[..., 0], ll[..., 1]
            good = ((np.abs(lat) <= 90.).all(axis=1) & (np.abs(lon) <= 180.).all(axis=1)
                    & ~((lat == 0.) & (lon == 0.)).all(axis=1))
            tgood = (year > 1970) & (day >= 1) & (day <= 366) & (ms >= 0) & (ms < 86_400_000)
            jan1 = (np.where(tgood, year, 2000) - 1970).astype('datetime64[Y]').astype(obstime.US)
            times = jan1 + ((day - 1) * 86_400_000 + ms).astype('timedelta64[ms]')
            times[~tgood] = obstime.NAT
            good &= tgood
            lat[~good], lon[~good] = np.nan, np.nan
            self._nav = dict(lat=lat, lon=lon, time=times, good=good)
            bad = int((~good).sum())
            if bad:
                print(tr('%s：%d 條掃描線的定位或時間無效', '%s: %d scanlines with invalid navigation or time')
                      % (self.fullname, bad))
        return self._nav

    def ch3_mode(self):
        """每條掃描線的通道 3 模式：True = 3A、False = 3B (由 3B 定標係數的線性項正負判斷，見檔頭說明)。"""
        if 'mode' not in self._cal:
            rec = self._records()
            slope = np.array(rec['ir_cal'][:, 0, 0, 1], dtype=np.float64)
            del rec
            is3a = np.where(slope > 0, 1., np.where(slope < 0, 0., np.nan))
            ok = np.isfinite(is3a)
            if ok.any():
                is3a = np.interp(np.arange(is3a.size), np.flatnonzero(ok), is3a[ok])
                is3a = median_filter(is3a, size=MODE_WINDOW, mode='nearest') > 0.5
            else:
                is3a = np.zeros(is3a.size, dtype=bool)
            self._cal['mode'] = is3a
        return self._cal['mode']

    def ir_coefficients(self, band):
        """紅外通道的線性定標係數 (列, 3) = [二次項, 一次項, 常數項] (已除 10^6)，及有效列。
        無效 (int32 上限、一次項 >= 0) 的掃描線以前後的有效值補上；與移動中位數相差超過門檻的掃描線改用中位數。"""
        if band in self._cal:
            return self._cal[band]
        rec = self._records()
        raw = np.array(rec['ir_cal'][:, IR_BANDS[band], 0, :], dtype=np.float64)
        del rec
        coef = raw / 1e6
        valid = (np.abs(raw) < INT32_MAX).all(axis=1) & (coef[:, 1] < 0)
        if band == 4:                       # 3A 的掃描線沒有 3B 的資料
            valid &= ~self.ch3_mode()
        usable = valid.copy()
        replaced = 0
        if valid.sum() >= 2:
            filled = _fill_from_neighbours(coef, valid)
            med = np.stack([median_filter(filled[:, k], size=CAL_WINDOW, mode='nearest') for k in range(3)], axis=1)
            ref = np.array([400., 800.])
            rad = lambda c: c[:, 2:3] + c[:, 1:2] * ref + c[:, 0:1] * ref ** 2
            spike = valid & (np.abs(rad(coef) - rad(med)) > CAL_TOLERANCE[band]).any(axis=1)
            replaced = int(spike.sum() + (~valid).sum())
            coef = np.where((spike | ~valid)[:, None], med, coef)
            usable = np.ones(coef.shape[0], dtype=bool)
            if band == 4:
                usable = ~self.ch3_mode()
        else:
            usable[:] = False
        if replaced and usable.any():
            print(tr('%s CH%s：%d 條掃描線的定標係數異常，改用前後 %d 條的中位數',
                     '%s CH%s: %d scanlines with abnormal calibration coefficients replaced by the %d-line median')
                  % (os.path.basename(self.fullname), BAND_LABEL[band], replaced, CAL_WINDOW))
        if band == 4 and self.ch3_mode().all():
            print(tr('%s：這一軌的通道 3 全部是 3A，沒有 3B', '%s: channel 3 is 3A for the whole pass; no 3B')
                  % os.path.basename(self.fullname))
        elif not usable.any():
            print(tr('%s CH%s：沒有有效的定標係數 (儀器異常？)', '%s CH%s: no valid calibration coefficients (instrument problem?)')
                  % (os.path.basename(self.fullname), BAND_LABEL[band]))
        self._cal[band] = (coef, usable)
        return self._cal[band]

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

    def _calibrate(self, band, r0, r1):
        """定標 (不含資料受損的檢查)：反射率 (0~1) 或亮溫 (°C)，(列, 2048) float64；無效為 NaN。"""
        rec = self._records()[r0:r1]
        counts = np.array(rec['counts'][:, :, COUNT_INDEX[band]])
        vis_cal = np.array(rec['vis_cal'][:, VIS_BANDS[band], 0, :]) if band in VIS_BANDS else None
        del rec
        bad = (counts == 0) | ~self.navigation()['good'][r0:r1, None]
        if band in VIS_BANDS:
            out = vis_albedo(counts, vis_cal) / 100.
            if band == 3:
                bad |= ~self.ch3_mode()[r0:r1, None]
        else:
            coef, usable = self.ir_coefficients(band)
            coef = coef[r0:r1]
            c = counts.astype(np.float64)
            rad = coef[:, 2:3] + coef[:, 1:2] * c + coef[:, 0:1] * c ** 2
            nlc = NONLINEAR.get(self.satellite, {}).get(band)
            if nlc:
                b0, b1, b2 = nlc
                rad = b0 + (1. + b1) * rad + b2 * rad ** 2
            nu, a, b = PLANCK.get(self.satellite, PLANCK['NOAA-18'])[band]
            out = brightness_temperature(rad, nu, a, b) - 273.15
            bad |= ~usable[r0:r1, None]
        out[bad] = np.nan
        return out

    def damage(self, r0, r1, split=None):
        """資料受損的檢查 (見 SPLIT_RANGE、LINE_JUMP_COUNTS)：回傳 (整條受損的掃描線, CH4 − CH5 不合理的像元)。
        split = 已算好的 (CH4, CH5) 亮溫，可省略。CH4 或 CH5 沒有有效定標時只做相鄰掃描線的檢查。"""
        t4, t5 = split if split is not None else (self._calibrate(5, r0, r1), self._calibrate(6, r0, r1))
        with np.errstate(invalid='ignore'):
            d = t4 - t5
            odd = (d < SPLIT_RANGE[0]) | (d > SPLIT_RANGE[1])
            n = np.isfinite(d).sum(axis=1)
            line = (n > 0) & (odd.sum(axis=1) > BAD_LINE_FRACTION * np.maximum(n, 1))
        # 相鄰掃描線：整條錯位時兩個紅外通道一起偏移，CH4 − CH5 看不出來
        a0, a1 = max(r0 - 1, 0), min(r1 + 1, self.nlines)
        rec = self._records()[a0:a1]
        c = np.array(rec['counts'][:, :, COUNT_INDEX[5]], dtype=np.int32)
        del rec
        if c.shape[0] >= 2:
            step = np.abs(np.diff(c, axis=0))
            big = np.full((1, c.shape[1]), 1 << 20, dtype=np.int64)
            prev = np.vstack([big, step])
            nxt = np.vstack([step, big])
            jump = np.median(np.minimum(prev, nxt), axis=1)[r0 - a0:r0 - a0 + (r1 - r0)]
            line |= jump > LINE_JUMP_COUNTS
        return line, odd

    def extract(self, band, rows=None):
        """反射率 (0~1) 或亮溫 (°C)，(列, 2048) float32；無效為 NaN。rows = (起, 訖)，預設全部。
        資料受損的掃描線整條為 NaN；CH4 / CH5 另外遮蔽分裂窗差不合理的像元。"""
        r0, r1 = rows or (0, self.nlines)
        out = self._calibrate(band, r0, r1)
        if band in (5, 6):
            other = self._calibrate(11 - band, r0, r1)
            line, odd = self.damage(r0, r1, (out, other) if band == 5 else (other, out))
            out[odd] = np.nan
        else:
            line, _ = self.damage(r0, r1)
        out[line] = np.nan
        if band in VIS_BANDS:
            repair_streaks(out)
        fill_short_gaps(out, line)
        if line.any():
            print(tr('%s：%d / %d 條掃描線的資料受損 (CH4 − CH5 不合理或與相鄰掃描線不連續)，已遮蔽',
                     '%s: %d / %d scanlines damaged (implausible CH4 − CH5 or discontinuous with the adjacent lines); masked')
                  % (os.path.basename(self.fullname), int(line.sum()), line.size))
        return out.astype(np.float32)

    def geocoord(self, rows=None):
        r0, r1 = rows or (0, self.nlines)
        nav = self.navigation()
        return interp_tiepoints(nav['lat'][r0:r1], nav['lon'][r0:r1], ALTITUDE_KM.get(self.satellite, 850.))

    def row_times(self, rows=None):
        """每條掃描線的觀測時間 (UTC datetime64[us]；無效為 NaT)。"""
        r0, r1 = rows or (0, self.nlines)
        return self.navigation()['time'][r0:r1]


# ===== 讀取器 =====
class AvhrrNsmcReader(SatReader):
    NAME = 'AVHRR-NSMC'
    SENSOR = 'AVHRR3'
    DESCRIPTION = 'NOAA-15~19 AVHRR/3 HRPT 1B (NSMC 格式，.L1B)'
    DESCRIPTION_EN = 'NOAA-15~19 AVHRR/3 HRPT 1B (NSMC format, .L1B)'
    SOURCE = 'NSMC'
    PATTERNS = ('NA??_AVHRR_HRPT_L1_*.L1B',)

    def make_file(self, name):
        return AvhrrNsmcFile(name)

    def listing(self, f):
        return (f.satellite, f.channel, f.time)

    def select_band(self):
        for f in self.files:
            print(f.fullname)
        print(BAND_MENU)
        self.band = ui.ask_int(tr('輸入頻道編號：', 'Input the band: '), range(1, 7))
        if self.band in (3, 4):
            n3a = sum(int(f.ch3_mode().sum()) for f in self.files)
            n = sum(f.nlines for f in self.files)
            print(tr('通道 3：3A %d 條、3B %d 條掃描線', 'Channel 3: %d 3A and %d 3B scanlines') % (n3a, n - n3a))

    def select_variants(self):
        if self.band in BT_BANDS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare(self):
        datas, latss, lonss, times = [], [], [], []
        for f in self.files:
            rows = f.rows_in(self.georange, self.center[1])
            if rows is None:
                continue
            datas.append(f.extract(self.band, rows))
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
        if variant is None:
            canvas.mesh(self.lons, self.lats, self.data, cmap='gray', vmin=0., vmax=1.)
            cm = ''
        else:
            cm = canvas.bt_mesh(self.lons, self.lats, self.data, variant)
        channel = BAND_LABEL[self.band]
        return dict(band='BAND-%s' % channel, colorbar=cm, label='BAND-%s' % channel)


READER = AvhrrNsmcReader
