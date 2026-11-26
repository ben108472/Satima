"""NOAA AVHRR Level 1b (POD / KLM 格式) 的檔案結構、衛星辨識與定標常數 (AVHRR-L1B 讀取器的輔助模組)

格式說明 (NCEI：https://www.ncei.noaa.gov/pub/data/cdo/documentation/podguides/)：
  POD Guide (TIROS-N ~ NOAA-14)：§2.0 資料集表頭、TBM 檔頭；§3.1 GAC、§3.2 LAC/HRPT；附錄 K、L (1992、1994 以前的版本)
  NOAA KLM User's Guide (NOAA-15 ~ 19、MetOp)：§8.3.1 AVHRR 1b；§7.1 定標；附錄 D 各衛星參數
位元組位置以下都是 0 起算 (指南為 1 起算)。所有欄位 big-endian。

檔案結構
  [前置檔頭] + 表頭記錄 + 資料記錄 (每條掃描線一筆，與表頭等長)
  - 前置檔頭：0、122 (TBM，ASCII，第 31~74 byte 為資料集名稱)、512 (ARS)。依檔案大小與表頭欄位判斷。
  - POD：GAC 3220 bytes (409 像元)，兩筆邏輯記錄一個實體記錄：表頭占前兩筆 (第二筆是殘留的某條掃描線)，
         資料從第 3 筆起；LAC/HRPT 14800 bytes (2048 像元)，資料從第 2 筆起。
  - KLM：GAC 4608、LAC/HRPT/FRAC 15872 bytes (10-bit)；資料從 表頭第 15-16 byte 的表頭筆數 之後起。
  - 感測資料：每個 32-bit 字放 3 個 10-bit 值 (bits 29-20、19-10、9-0)，像元內依頻道 1~5 排列 (BIP)。
衛星：資料集名稱 (NSS.<類型>.<衛星代碼>.D<yyddd>…) 優先，其次表頭的衛星 ID (POD 與 KLM 的編號不同；
  千葉大學把 NOAA-18 存成 POD 格式，ID 仍是 KLM 的 7 = NOAA-18，在 POD 表裡是 NOAA-9，所以依年份判斷)。
"""
import codecs
import datetime
import os
import re

import numpy as np

# ---- 衛星 ----------------------------------------------------------------------------------
# 資料集名稱的衛星代碼 -> (衛星, AVHRR 世代 1/2/3)。世代 1 (4 頻道，頻道 5 = 頻道 4 的複本)、2 (5 頻道)、3 (3A/3B)
SATELLITES = {
    'TN': ('TIROS-N', 1), 'NA': ('NOAA-6', 1), 'NC': ('NOAA-7', 2), 'NE': ('NOAA-8', 1), 'NF': ('NOAA-9', 2),
    'NG': ('NOAA-10', 1), 'NH': ('NOAA-11', 2), 'ND': ('NOAA-12', 2), 'NI': ('NOAA-13', 2), 'NJ': ('NOAA-14', 2),
    'NK': ('NOAA-15', 3), 'NL': ('NOAA-16', 3), 'NM': ('NOAA-17', 3), 'NN': ('NOAA-18', 3), 'NP': ('NOAA-19', 3),
    'M2': ('METOP-A', 3), 'M1': ('METOP-B', 3), 'M3': ('METOP-C', 3),
}
GENERATION = {name: gen for name, gen in SATELLITES.values()}
# 表頭的衛星 ID：POD (指南 Table 2.0.4-3；1 = TIROS-N 與 NOAA-11，NOAA-13 的 ID 指南前後不一) 與 KLM
POD_IDS = {1: ('TIROS-N', 'NOAA-11'), 2: ('NOAA-6', 'NOAA-13'), 3: ('NOAA-14', 'NOAA-13'), 4: ('NOAA-7',),
           5: ('NOAA-12',), 6: ('NOAA-8',), 7: ('NOAA-9',), 8: ('NOAA-10',)}
KLM_IDS = {4: 'NOAA-15', 2: 'NOAA-16', 6: 'NOAA-17', 7: 'NOAA-18', 8: 'NOAA-19', 12: 'METOP-A', 11: 'METOP-B',
           13: 'METOP-C'}
# 運作年份 (只用於 ID 判斷衛星)
SERVICE = {'TIROS-N': (1978, 1981), 'NOAA-6': (1979, 1987), 'NOAA-7': (1981, 1986), 'NOAA-8': (1983, 1985),
           'NOAA-9': (1984, 1998), 'NOAA-10': (1986, 2001), 'NOAA-11': (1988, 2004), 'NOAA-12': (1991, 2007),
           'NOAA-13': (1993, 1993), 'NOAA-14': (1994, 2007), 'NOAA-15': (1998, 2025), 'NOAA-16': (2000, 2014),
           'NOAA-17': (2002, 2013), 'NOAA-18': (2005, 2025), 'NOAA-19': (2009, 2025), 'METOP-A': (2006, 2021),
           'METOP-B': (2012, 2040), 'METOP-C': (2018, 2040)}
# 平均高度 (km)：只用於沿掃描內插的地心角
ALTITUDE_KM = {'TIROS-N': 850., 'NOAA-6': 833., 'NOAA-7': 847., 'NOAA-8': 813., 'NOAA-9': 850., 'NOAA-10': 833.,
               'NOAA-11': 841., 'NOAA-12': 804., 'NOAA-13': 850., 'NOAA-14': 845., 'NOAA-15': 808., 'NOAA-16': 850.,
               'NOAA-17': 811., 'NOAA-18': 854., 'NOAA-19': 870., 'METOP-A': 817., 'METOP-B': 817., 'METOP-C': 817.}
DATA_TYPES = {1: 'LAC', 2: 'GAC', 3: 'HRPT', 13: 'FRAC'}
NAME_TYPES = {'LHRR': 'LAC', 'GHRR': 'GAC', 'HRPT': 'HRPT', 'FRAC': 'FRAC'}

# ---- 記錄結構 (0 起算) -------------------------------------------------------------------------
POD = 'POD'
KLM = 'KLM'
RECORD = {(POD, 'GAC'): 3220, (POD, 'LAC'): 14800, (KLM, 'GAC'): 4608, (KLM, 'LAC'): 15872}
NPIX = {'GAC': 409, 'LAC': 2048}
HEADER_BYTES = 4096                 # 偵測時讀取的開頭長度


def _words(npix):
    return -(-npix * 5 // 3)


def pod_record(kind):
    """POD 資料記錄 (POD Guide Table 3.1.2.1-1 / 3.2.2.1-1)。kind = 'GAC' / 'LAC' (HRPT 同 LAC)。"""
    n = NPIX[kind]
    fields = [('scan', '>u2', 0), ('yearday', '>u2', 2), ('ms', '>u4', 4), ('quality', '>u4', 8),
              ('cal', '>i4', 12, (5, 2)),           # [頻道 1~5, 斜率 ×2^30 / 截距 ×2^22]
              ('npoints', 'u1', 52), ('sza', 'u1', 53, (51,)),          # 太陽天頂角 ×2
              ('latlon', '>i2', 104, (51, 2)),      # 緯度、經度 ×128
              ('telemetry', '>u4', 308, (35,)),     # HRPT minor frame word 1~103 (10-bit，每字 3 個)
              ('video', '>u4', 448, (_words(n),))]
    return _dtype(fields, RECORD[(POD, kind)])


def klm_record(kind):
    """KLM 資料記錄 (KLM Guide Table 8.3.1.3.3.1-1 / 8.3.1.4.3.1-1；GAC 的 1~1264 byte 與 LAC 相同)。"""
    n = NPIX[kind]
    fields = [('scan', '>u2', 0), ('year', '>u2', 2), ('day', '>u2', 4), ('drift', '>i2', 6), ('ms', '>u4', 8),
              ('bitfield', '>u2', 12),              # bit 15：1 = 南行；bits 1-0：頻道 3 (0 = 3B、1 = 3A、2 = 切換中)
              ('quality', '>u4', 24), ('line_quality', '>u4', 28),
              ('vis_cal', '>i4', 48, (3, 3, 5)),    # [CH1/CH2/CH3A, 業務/測試/發射前, 斜率1 截距1 斜率2 截距2 交點]
              ('ir_cal', '>i4', 228, (3, 2, 3)),    # [CH3B/CH4/CH5, 業務/測試, a0 a1 a2]
              ('angles', '>i2', 328, (51, 3)),      # 太陽天頂角、衛星天頂角、相對方位角 ×100
              ('latlon', '>i4', 640, (51, 2)),      # 緯度、經度 ×1e4
              ('prt', '>u2', 1090, (3,)), ('ict', '>u2', 1100, (30,)), ('space', '>u2', 1160, (50,)),
              ('video', '>u4', 1264, (_words(n),))]
    return _dtype(fields, RECORD[(KLM, kind)])


def _dtype(fields, size):
    names, formats, offsets = [], [], []
    for f in fields:
        names.append(f[0])
        formats.append((f[1], f[3]) if len(f) > 3 else f[1])
        offsets.append(f[2])
    return np.dtype({'names': names, 'formats': formats, 'offsets': offsets, 'itemsize': size})


def unpack10(words, nvalues):
    """32-bit 字 (列, n) -> 10-bit 值 (列, nvalues)，uint16。"""
    w = np.asarray(words, dtype=np.uint32)
    out = np.empty(w.shape[:-1] + (w.shape[-1] * 3,), dtype=np.uint16)
    out[..., 0::3] = (w >> 20) & 0x3FF
    out[..., 1::3] = (w >> 10) & 0x3FF
    out[..., 2::3] = w & 0x3FF
    return out[..., :nvalues]


# ---- 偵測 ----------------------------------------------------------------------------------
def _pod_time(b):
    yd = int.from_bytes(b[0:2], 'big')
    ms = int.from_bytes(b[2:6], 'big') & 0x7FFFFFF
    return yd >> 9, yd & 0x1FF, ms


def full_year(y):
    return y + (1900 if y >= 70 else 2000) if y < 100 else y


def _datetime(year, day, ms):
    try:
        if not (1 <= day <= 366 and 0 <= ms < 86_400_000 and 1970 < year < 2100):
            return None
        return datetime.datetime(year, 1, 1) + datetime.timedelta(days=day - 1, milliseconds=ms)
    except (ValueError, OverflowError):
        return None


_NAME = re.compile(r'NSS\.(GHRR|LHRR|HRPT|FRAC)\.([A-Z]\w)\.D(\d{5})\.S(\d{4})\.E(\d{4})\.B(\d{7})\.(\w\w)')


def _names(raw):
    """一段 bytes 中出現的資料集名稱 (ASCII 或 EBCDIC)。"""
    out = []
    for text in (raw.decode('latin1'), codecs.decode(raw, 'cp500')):
        out += [m for m in _NAME.finditer(text)]
    return out


class Structure:
    """檔案結構：layout (POD / KLM)、kind (GAC / LAC)、dtype ('GAC'、'LAC'、'HRPT'、'FRAC')、prefix、record、first (第一筆資料
    的位置)、nlines、satellite、generation、version (KLM 表頭版本)、time / end (表頭的起訖時間)、station (前置檔頭或
    B0000000 = 接收站自行產生，例如千葉大學)、header (表頭記錄 bytes)。"""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def sniff(path, head=None):
    """依檔案內容判斷是否為 NOAA AVHRR 1b；不是時回傳 None (不丟例外，用於掃描資料夾)。"""
    try:
        size = os.path.getsize(path)
        if size < 3 * 3220:
            return None
        if head is None:
            with open(path, 'rb') as f:
                head = f.read(HEADER_BYTES)
    except OSError:
        return None
    for pre in (0, 122, 512):
        h = head[pre:pre + 200]
        if len(h) < 120:
            continue
        s = _try_klm(path, size, head, pre) or _try_pod(path, size, head, pre)
        if s is not None:
            return s
    return None


def _try_klm(path, size, head, pre):
    h = head[pre:]
    ver = int.from_bytes(h[4:6], 'big')
    nhdr = int.from_bytes(h[14:16], 'big')
    sc = int.from_bytes(h[72:74], 'big')
    dt = int.from_bytes(h[76:78], 'big')
    if not (1 <= ver <= 10 and sc in KLM_IDS and dt in DATA_TYPES and 1 <= nhdr <= 4):
        return None
    kind = 'GAC' if dt == 2 else 'LAC'
    rec = RECORD[(KLM, kind)]
    if (size - pre) % rec:
        return None
    t0 = _datetime(int.from_bytes(h[84:86], 'big'), int.from_bytes(h[86:88], 'big'), int.from_bytes(h[88:92], 'big'))
    if t0 is None:
        return None
    t1 = _datetime(int.from_bytes(h[96:98], 'big'), int.from_bytes(h[98:100], 'big'),
                   int.from_bytes(h[100:104], 'big'))
    nlines = (size - pre) // rec - nhdr
    if nlines < 1:
        return None
    names = _names(head[:pre] + h[22:64])
    sat, station = _satellite(names, KLM_IDS.get(sc), t0, path)
    if pre == 512 and not station:
        station = 'B0000000' in head[:pre].decode('latin1')
    return Structure(layout=KLM, kind=kind, dtype=DATA_TYPES[dt], prefix=pre, record=rec, first=pre + nhdr * rec,
                     nlines=nlines, satellite=sat, generation=GENERATION.get(sat, 3), version=ver, time=t0,
                     end=t1 if t1 and t1 >= t0 else None, station=station, header=bytes(h[:rec][:1024]))


def _try_pod(path, size, head, pre):
    h = head[pre:]
    sc, dt = h[0], h[1] >> 4
    if dt not in (1, 2, 3):
        return None
    kind = 'GAC' if dt == 2 else 'LAC'
    rec = RECORD[(POD, kind)]
    if (size - pre) % rec:
        return None
    y, d, ms = _pod_time(h[2:8])
    year = int.from_bytes(h[38:40], 'big')        # 1998/12 以後的表頭有四位數年份；更早的版本這裡可能是殘留資料
    if not (1970 < year < 2100 and year % 100 == y):
        year = full_year(y)
    t0 = _datetime(year, d, ms)
    if t0 is None:
        return None
    y1, d1, ms1 = _pod_time(h[10:16])
    t1 = _datetime(year + (1 if d1 < d else 0), d1, ms1)
    nskip = 2 if kind == 'GAC' else 1
    nlines = (size - pre) // rec - nskip
    if nlines < 1:
        return None
    names = _names(head[:pre] + h[40:84])
    if not names:
        # 沒有資料集名稱：第一筆資料的時間應接近表頭的起始時間 (避免把其他格式誤認)
        try:
            with open(path, 'rb') as f:
                f.seek(pre + nskip * rec)
                r = f.read(8)
            y2, d2, ms2 = _pod_time(r[2:8])
            if d2 != d or abs(ms2 - ms) > 600_000:
                return None
        except OSError:
            return None
    guess = POD_IDS.get(sc, ())
    sat, station = _satellite(names, None, t0, path, guess + ((KLM_IDS[sc],) if sc in KLM_IDS else ()))
    if pre == 122:
        station = True
    if sat is None:
        return None
    return Structure(layout=POD, kind=kind, dtype=DATA_TYPES[dt], prefix=pre, record=rec, first=pre + nskip * rec,
                     nlines=nlines, satellite=sat, generation=GENERATION.get(sat, 2), version=0, time=t0,
                     end=t1 if t1 and t1 >= t0 else None, station=station, header=bytes(h[:1024]))


def _satellite(names, by_id, t0, path, candidates=()):
    """衛星：檔內的資料集名稱 > 表頭 ID (依年份在候選中挑) > 檔名。回傳 (衛星, 是否為接收站產品)。"""
    station = False
    for m in names:
        if m.group(6) == '0000000':
            station = True
        if m.group(2) in SATELLITES:
            return SATELLITES[m.group(2)][0], station
    if by_id:
        return by_id, station
    year = t0.year
    for sat in candidates:
        a, b = SERVICE[sat]
        if a <= year <= b:
            return sat, station
    m = _NAME.search(os.path.basename(path))
    if m and m.group(2) in SATELLITES:
        return SATELLITES[m.group(2)][0], station or m.group(6) == '0000000'
    return (candidates[0] if candidates else None), station


# ---- 定標常數 ------------------------------------------------------------------------------
# Planck (KLM Guide Eq. 7.1.2.4-8/9)：T* = c2·ν / ln(1 + c1·ν³/N)，T = (T* − A) / B
C1 = 1.1910427E-5           # mW/(m2 sr cm-4)
C2 = 1.4387752              # cm K
# 頻道 (3 = 3B / POD 的頻道 3、4、5) -> (ν, A, B)
#   NOAA-15~19、MetOp-A/B：KLM Guide 附錄 D (與 1b 表頭的 −A/B、1/B 一致)；MetOp-C：指南沒有，用 pygac 1.8 (PATMOS-x)
#   TIROS-N ~ NOAA-14：POD Guide 只有依溫度分段的中心波數 (沒有頻帶修正)，用 pygac 的單一 ν、A、B
#   (與指南分段值換算的亮溫差：頻道 4/5 ≤ 0.03 K、頻道 3 −0.3~+0.6 K)；NOAA-13 用指南 270–310 K 段、A = 0、B = 1
PLANCK = {
    'TIROS-N': {3: (2655.7409, 1.645107, 0.997915), 4: (913.05397, 0.530593, 0.998568)},
    'NOAA-6': {3: (2671.5433, 1.762406, 0.997563), 4: (913.46088, 0.503276, 0.998643)},
    'NOAA-7': {3: (2684.5233, 1.943141, 0.997083), 4: (928.23757, 0.527340, 0.998598),
               5: (841.52137, 0.405093, 0.998822)},
    'NOAA-8': {3: (2651.3776, 1.772111, 0.997580), 4: (915.3033, 0.499508, 0.998656)},
    'NOAA-9': {3: (2690.0451, 1.877825, 0.997111), 4: (930.5023, 0.510840, 0.998645), 5: (845.75, 0.387780, 0.998880)},
    'NOAA-10': {3: (2672.6164, 1.793970, 0.997374), 4: (910.49626, 0.456510, 0.998774)},
    'NOAA-11': {3: (2680.05, 1.733160, 0.996657), 4: (927.462, 0.320810, 0.998788), 5: (840.746, 0.048620, 0.999336)},
    'NOAA-12': {3: (2651.7708, 1.899556, 0.996999), 4: (922.36261, 0.632961, 0.998295),
                5: (838.02678, 0.410373, 0.998800)},
    'NOAA-13': {3: (2643.153, 0., 1.), 4: (924.9732, 0., 1.), 5: (836.7651, 0., 1.)},
    'NOAA-14': {3: (2654.25, 1.878120, 0.996176), 4: (928.349, 0.307940, 0.998559), 5: (833.04, -0.022159, 0.999462)},
    'NOAA-15': {3: (2695.9743, 1.621256, 0.998015), 4: (925.4075, 0.337810, 0.998719),
                5: (839.8979, 0.304558, 0.999024)},
    'NOAA-16': {3: (2700.1148, 1.592459, 0.998147), 4: (917.2289, 0.332380, 0.998522),
                5: (838.1255, 0.674623, 0.998363)},
    'NOAA-17': {3: (2669.3554, 1.702380, 0.997378), 4: (926.2947, 0.271683, 0.998794),
                5: (839.8246, 0.309180, 0.999012)},
    'NOAA-18': {3: (2659.7952, 1.698704, 0.996960), 4: (928.1460, 0.436645, 0.998607),
                5: (833.2532, 0.253179, 0.999057)},
    'NOAA-19': {3: (2670.0, 1.67396, 0.997364), 4: (928.9, 0.53959, 0.998534), 5: (831.9, 0.36064, 0.998913)},
    'METOP-A': {3: (2687.0, 2.06699, 0.996577), 4: (927.2, 0.55126, 0.998509), 5: (837.7, 0.34716, 0.998947)},
    'METOP-B': {3: (2684.32, 1.763611, 0.997018), 4: (933.63, 0.504183, 0.998638), 5: (839.62, 0.381279, 0.998610)},
    'METOP-C': {3: (2707.6457, 1.782461, 0.997638), 4: (931.89092, 0.564729, 0.998492),
                5: (832.69445, 0.391622, 0.998851)},
}
for _sat in ('TIROS-N', 'NOAA-6', 'NOAA-8', 'NOAA-10'):        # 4 頻道儀器：頻道 5 是頻道 4 的複本
    PLANCK[_sat][5] = PLANCK[_sat][4]

# 非線性修正 N = b0 + (1 + b1)·N_lin + b2·N_lin²，N_lin 以太空輻射 Ns 為零點的線性定標：頻道 -> (Ns, b0, b1, b2)
#   NOAA-15~19、MetOp-A/B：KLM 附錄 D (MetOp-B 頻道 4 的 b1 印成 −0.0096771，應為 −0.096771)；MetOp-C：pygac
#   TIROS-N ~ NOAA-12：Walton et al. 1998 (pygac；用指南的 Ns = 0 修正表模擬比較，差 0.1~0.4 K)；
#   NOAA-13 / 14：POD Guide Table 1.4.9-4 / 1.4.10-3 (RAD = A·R + B·R² + D，即 b0 = D、b1 = A − 1、b2 = B)
#   沒有列出的頻道 (頻道 3) 不做非線性修正、Ns = 0
NONLINEAR = {
    'TIROS-N': {4: (-8.13, 6.13, -0.131942, 0.000673193)},
    'NOAA-6': {4: (-3.26, 2.24, -0.03964, 0.00016925)},
    'NOAA-7': {4: (-5.16, 5.25, -0.10217, 0.0004819), 5: (-4.28, 3.93, -0.06317, 0.0002425)},
    'NOAA-8': {4: (-3.26, 2.24, -0.03964, 0.00016925)},
    'NOAA-9': {4: (-5.53, 5.24, -0.1136, 0.0006033), 5: (-3.06, 2.42, -0.0469, 0.0002198)},
    'NOAA-10': {4: (-7.29, 5.76, -0.1157, 0.0005882)},
    'NOAA-11': {4: (-8.055, 7.21, -0.1588, 0.0008739), 5: (-3.51, 2.92, -0.054, 0.0002504)},
    'NOAA-12': {4: (-5.51, 5.11, -0.1107, 0.0005968), 5: (-2.51, 1.91, -0.037, 0.0001775)},
    'NOAA-13': {4: (-5.31, 5.01, -0.08841, 0.0003820), 5: (-3.28, 3.24, -0.05216, 0.0002057)},
    'NOAA-14': {3: (0.0069, -0.0031, 0.00359, 0.), 4: (-4.05, 3.72, -0.07622, 0.0003822),
                5: (-2.29, 2.00, -0.03806, 0.0001742)},
    'NOAA-15': {4: (-4.50, 4.76, -0.0932, 0.0004524), 5: (-3.61, 3.83, -0.0659, 0.0002811)},
    'NOAA-16': {4: (-2.467, 2.96, -0.05411, 0.00024532), 5: (-2.009, 2.25, -0.03665, 0.00014854)},
    'NOAA-17': {4: (-8.55, 8.22, -0.15795, 0.00075579), 5: (-3.97, 4.31, -0.07318, 0.00030976)},
    'NOAA-18': {4: (-5.53, 5.82, -0.11069, 0.00052337), 5: (-2.22, 2.67, -0.04360, 0.00017715)},
    'NOAA-19': {4: (-5.49, 5.70, -0.11187, 0.00054668), 5: (-3.39, 3.58, -0.05991, 0.00024985)},
    'METOP-A': {4: (-4.98, 5.44, -0.10152, 0.00046964), 5: (-3.40, 3.84, -0.06249, 0.00025239)},
    'METOP-B': {4: (-4.75, 4.85, -0.096771, 0.00048091), 5: (-4.39, 4.36, -0.0766350, 0.00033524)},
    'METOP-C': {4: (-6.27, 6.58, -0.13203, 0.00065922), 5: (-2.55, 3.23, -0.05692, 0.00024963)},
}
for _sat in ('TIROS-N', 'NOAA-6', 'NOAA-8', 'NOAA-10'):
    NONLINEAR[_sat][5] = NONLINEAR[_sat][4]

# PRT 計數 -> 溫度 T = d0 + d1·C + d2·C² (+ d3·C³ + d4·C⁴) (K)，四顆 PRT (KLM 附錄 D；MetOp-C：pygac)。
# 用於由星上黑體自行定標：KLM 衛星的資料缺定標係數 (例如存成 POD 格式) 時 (KLM 格式的檔用表頭的係數)，以及原始框 (AVHRR-L0)。
# POD 衛星 (TIROS-N ~ NOAA-14)：指南沒有，用 pygac 的值；以 CLASS 1b 驗證 (NOAA-9/10/11/14)，換算的黑體溫度與 1b 係數
# 在黑體計數處的溫度相差 ±0.01 K 以內 (中位數)。
_POD_PRT_A = [(276.659, 0.051275, 1.363e-6)] * 4
_POD_PRT_B = [(276.597, 0.051275, 1.363e-6)] * 4
PRT = {
    'TIROS-N': _POD_PRT_A, 'NOAA-6': _POD_PRT_A, 'NOAA-8': _POD_PRT_A, 'NOAA-10': _POD_PRT_A,
    'NOAA-7': [(277.099, 0.05048, 2.823e-6), (276.734, 0.05069, 2.493e-6), (276.876, 0.05148, 1.04e-6),
               (276.16, 0.05128, 1.414e-6)],
    'NOAA-9': [(277.018, 0.05128, 0.), (276.750, 0.05128, 0.), (276.862, 0.05128, 0.), (276.546, 0.05128, 0.)],
    'NOAA-11': _POD_PRT_B, 'NOAA-12': _POD_PRT_B, 'NOAA-13': _POD_PRT_B, 'NOAA-14': _POD_PRT_B,
    'NOAA-15': [(276.60157, 0.051045, 1.36328e-6), (276.62531, 0.050909, 1.47266e-6),
                (276.67413, 0.050907, 1.47656e-6), (276.59258, 0.050966, 1.47656e-6)],
    'NOAA-16': [(276.355, 5.562e-2, -1.590e-5, 2.486e-8, -1.199e-11), (276.142, 5.605e-2, -1.707e-5, 2.595e-8, -1.224e-11),
                (275.996, 5.486e-2, -1.223e-5, 1.862e-8, -0.853e-11), (276.132, 5.494e-2, -1.344e-5, 2.112e-8, -1.001e-11)],
    'NOAA-17': [(276.628, 0.05098, 1.371e-6), (276.538, 0.05098, 1.371e-6), (276.761, 0.05097, 1.369e-6),
                (276.660, 0.05100, 1.348e-6)],
    'NOAA-18': [(276.601, 0.05090, 1.657e-6), (276.683, 0.05101, 1.482e-6), (276.565, 0.05117, 1.313e-6),
                (276.615, 0.05103, 1.484e-6)],
    'NOAA-19': [(276.6067, 0.051111, 1.405783e-6), (276.6119, 0.051090, 1.496037e-6),
                (276.6311, 0.051033, 1.496990e-6), (276.6268, 0.051058, 1.493110e-6)],
    'METOP-A': [(276.6194, 0.050919, 1.470892e-6), (276.6511, 0.050892, 1.489000e-6),
                (276.6597, 0.050845, 1.520646e-6), (276.3685, 0.050992, 1.482390e-6)],
    'METOP-B': [(276.5853, 0.05093323, 1.543330e-6), (276.5335, 0.05103343, 1.497510e-6),
                (276.5721, 0.05109724, 1.429280e-6), (276.5720, 0.05102045, 1.508410e-6)],
    'METOP-C': [(276.5862, 0.051051, 1.474208e-6), (276.6136, 0.051029, 1.472138e-6),
                (276.5975, 0.051065, 1.469268e-6), (276.4595, 0.050990, 1.506223e-6)],
}

# KLM 衛星的發射前可見光係數 (反照率 %；KLM 附錄 D)：頻道 (1、2、3 = 3A) -> (斜率1, 截距1, 斜率2, 截距2, 交點)。
# 只用於 KLM 衛星的資料存成 POD 格式 (單一斜率的欄位只放得下低增益段) 時。
VIS_PRELAUNCH = {
    'NOAA-15': {1: (0.0568, -2.1874, 0.1633, -54.9928, 496), 2: (0.0596, -2.4096, 0.1629, -55.2436, 511),
                3: (0.0275, -1.0684, 0.1846, -78.1691, 491)},
    'NOAA-16': {1: (0.0523, -2.016, 0.1528, -51.91, 497.5), 2: (0.0513, -1.943, 0.1510, -51.77, 500.3),
                3: (0.0262, -1.01, 0.1920, -84.2, 501.8)},
    'NOAA-17': {1: (0.0555, -2.2193, 0.1627, -55.9635, 497.53), 2: (0.0543, -2.1227, 0.1621, -56.2160, 500.32),
                3: (0.0265, -1.1153, 0.1860, -81.2520, 498.66)},
    'NOAA-18': {1: (0.05400, -2.130, 0.16100, -55.37, 500.54), 2: (0.05290, -2.084, 0.15873, -55.31, 500.40),
                3: (0.02618, -0.982, 0.18492, -80.88, 500.56)},
    'NOAA-19': {1: (0.055091, -2.1415, 0.16253, -55.863, 496.43), 2: (0.054892, -2.1288, 0.16325, -56.445, 500.37),
                3: (0.027174, -1.0881, 0.18798, -81.491, 496.11)},
    'METOP-A': {1: (0.053714, -2.171861, 0.158723, -54.782424, 500.54), 2: (0.0545174, -2.167, 0.16192, -55.91, 500.40),
                3: (0.026446, -1.087, 0.18369, -80.01, 500.56)},
    'METOP-B': {1: (0.05176, -2.039, 0.15335, -52.899, 500.54), 2: (0.05678, -2.291, 0.16727, -57.812, 500.40)},
}


# POD 衛星的發射前可見光係數 (反照率 % = S·C + I；POD Guide 表 3.3.2-1)：頻道 1、2 -> (S, I)。
# 只用於原始框 (AVHRR-L0)。沒有在軌衰減修正 (NOAA-14 自 1996/11 起 1b 逐月更新，原始框沒有這些係數)。
VIS_PRELAUNCH_POD = {
    'TIROS-N': {1: (0.1071, -3.9), 2: (0.1051, -3.5)}, 'NOAA-6': {1: (0.1071, -4.1136), 2: (0.1058, -3.4539)},
    'NOAA-7': {1: (0.1068, -3.4400), 2: (0.1069, -3.488)}, 'NOAA-8': {1: (0.1060, -4.1619), 2: (0.1060, -4.1492)},
    'NOAA-9': {1: (0.1063, -3.8464), 2: (0.1075, -3.8770)}, 'NOAA-10': {1: (0.1059, -3.5279), 2: (0.1061, -3.4766)},
    'NOAA-11': {1: (0.0906, -3.730), 2: (0.0900, -3.390)}, 'NOAA-12': {1: (0.1042, -4.4491), 2: (0.1014, -3.9925)},
    'NOAA-13': {1: (0.1076, -3.9747), 2: (0.1035, -3.8280)}, 'NOAA-14': {1: (0.1081, -3.8648), 2: (0.1090, -3.6749)},
}


# 在軌可見光係數：NOAA STAR (F. Wu) 以利比亞沙漠做的替代定標，每月一組 (AAPP 附的檔，放在 conversion_table/AVHRR/)
#   <衛星>_AVHRR_Libya_ch<1|2|3a>.txt：每行 Update、Data、Center 三個日期 ... Slope_lo Int_lo Slope_hi Int_hi (行末可能有註解)
#   NOAA-16 2003/05~2013/11、NOAA-17 2004/02~2010/10、NOAA-18 2005/10~2017/08、NOAA-19 2009/05~2017/08；NOAA-15 以前沒有。
STAR_NAMES = {'NOAA-16': 'N16', 'NOAA-17': 'N17', 'NOAA-18': 'N18', 'NOAA-19': 'N19'}
STAR_MAX_DAYS = 400             # 觀測與最近一組的 Center 相差超過此天數時不用 (改用發射前係數，或依下面的內插)
_STAR = {}
_DATE3 = re.compile(r'^\s*(\d{2})/(\d{2})/(\d{4})\s+\d{2}/\d{2}/\d{4}\s+(\d{2})/(\d{2})/(\d{4})\s')
_FLOAT = re.compile(r'^-?\d+\.\d+$')


def _star_table(table_dir, satellite, ch):
    """STAR 係數檔 -> (發射日, [(Center 日期, (s1, i1, s2, i2))])；沒有檔時 None。"""
    key = (table_dir, satellite, ch)
    if key in _STAR:
        return _STAR[key]
    out = None
    tag = STAR_NAMES.get(satellite)
    path = os.path.join(table_dir or '', 'AVHRR', '%s_AVHRR_Libya_ch%s.txt' % (tag, '3a' if ch == 3 else ch)) if tag else ''
    if path and os.path.isfile(path):
        launch, rows = None, []
        pre = VIS_PRELAUNCH.get(satellite, {}).get(ch)
        for line in open(path, encoding='latin1'):
            m = _DATE3.match(line)
            if not m:
                d = re.match(r'^\s*(\d{2})/(\d{2})/(\d{4})\s', line)
                if d and launch is None and not rows:
                    launch = datetime.date(int(d.group(3)), int(d.group(1)), int(d.group(2)))
                continue
            nums = [float(t) for t in line.split()[3:] if _FLOAT.match(t)]
            if len(nums) < 4:
                continue
            s1, i1, s2, i2 = nums[-4:]
            ok = 0 < s1 < s2 and i2 < i1 < 0 and 300 < (i1 - i2) / (s2 - s1) < 700
            if ok and pre:
                ok = 0.6 < s1 / pre[0] < 1.8               # 排除明顯錯誤的一組 (例如 NOAA-19 的 3A)
            if ok:
                center = datetime.date(int(m.group(6)), int(m.group(4)), int(m.group(5)))
                rows.append((center, (s1, i1, s2, i2)))
        if rows:
            out = (launch, sorted(rows))
    _STAR[key] = out
    return out


def star_vis(table_dir, satellite, ch, when):
    """觀測時間 when 的在軌可見光係數 (s1, i1, s2, i2, 交點) 與說明；沒有可用的時 None。
    取 Center 最接近的一組；早於第一組超過 120 天時，在發射日的發射前係數與第一組之間線性內插。"""
    tab = _star_table(table_dir, satellite, ch)
    if tab is None:
        return None
    launch, rows = tab
    day = when.date() if isinstance(when, datetime.datetime) else when
    center, coef = min(rows, key=lambda r: abs((r[0] - day).days))
    pre = VIS_PRELAUNCH.get(satellite, {}).get(ch)
    first = rows[0][0]
    if day < first and (first - day).days > 120 and launch and pre and launch <= day:
        w = (day - launch).days / max(1, (first - launch).days)
        coef = tuple(p + w * (c - p) for p, c in zip(pre[:4], rows[0][1]))
        note = 'STAR %s~%s' % (launch, first)
    elif abs((center - day).days) > STAR_MAX_DAYS:
        return None
    else:
        note = 'STAR %s' % center
    s1, i1, s2, i2 = coef
    return (s1, i1, s2, i2, (i1 - i2) / (s2 - s1)), note


def planck_bt(radiance, nu, a, b):
    """輻射 mW/(m2 sr cm-1) -> 亮溫 (K)；N <= 0 為 NaN。"""
    with np.errstate(invalid='ignore', divide='ignore'):
        rad = np.where(radiance > 0, radiance, np.nan)
        tstar = C2 * nu / np.log1p(C1 * nu ** 3 / rad)
    return (tstar - a) / b


def planck_radiance(temp, nu, a, b):
    """亮溫 (K) -> 輻射。"""
    with np.errstate(over='ignore'):
        return C1 * nu ** 3 / np.expm1(C2 * nu / (a + b * np.asarray(temp, dtype=np.float64)))
