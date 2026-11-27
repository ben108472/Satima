"""
AVHRR-L1B — TIROS-N ~ NOAA-19、MetOp AVHRR 的 NOAA Level 1b (NOAA CLASS 的 GAC / LAC / HRPT / FRAC、千葉大學 HRPT)
格式說明：NOAA POD Guide (TIROS-N ~ NOAA-14)、NOAA KLM User's Guide (NOAA-15 以後)；結構與常數見 avhrr_l1b_noaa.py

所需檔案
  NSS.GHRR / LHRR / HRPT / FRAC.<衛星>.D<yyddd>.S<hhmm>.E<hhmm>.B<nnnnnnn>.<站>   (CLASS)
  不依檔名辨識：資料夾中的每個檔依內容判斷 (POD / KLM、GAC / LAC、前置檔頭)，改過檔名也能讀。
  千葉大學的檔 (處理區塊編號 B0000000) 前面多 122 bytes (TBM，POD 格式) 或 512 bytes (KLM 格式)，
  而且 NOAA-18 存成 POD 格式 (表頭的衛星 ID 7 是 KLM 的編號；衛星以檔內的資料集名稱為準)。
  也可以是 .gz / .bz2 / .xz 壓縮檔 (例如千葉的 n1297121720.l1b.bz2)：列表時只解壓開頭判斷，畫圖時整個解壓進記憶體
  (不在磁碟留下解壓檔；.bz2 解壓 70 MB 約 2~3 秒)。.bz2 事先不知道解壓後的大小，掃描線數先用表頭的記錄數，
  讀檔時以實際長度為準。同一個檔的壓縮與未壓縮版本都在時只列未壓縮的。

頻道
  AVHRR/3 (NOAA-15 以後)：1.CH1 0.63  2.CH2 0.86  3.CH3A 1.61  4.CH3B 3.74  5.CH4 10.8  6.CH5 12.0
    CH3A / CH3B 共用通道 3，依掃描線切換 (KLM 掃描線 bit 域 bits 1-0；POD 格式依頻道 3 定標斜率的正負)
  AVHRR/2 (NOAA-7、9、11、12、13、14)：1.CH1 2.CH2 3.CH3 3.7 4.CH4 10.8 5.CH5 12.0
  AVHRR/1 (TIROS-N、NOAA-6、8、10)：CH1~CH4 (頻道 5 是頻道 4 的複本，不列出)
  可見光為反射率 = 反照率 / 100 (未除 cos(太陽天頂角)，同其他 AVHRR 讀取器)；紅外為亮溫

定標
  可見光
    - KLM：記錄中的業務 (operational) 係數，雙增益：反照率 = 斜率1·C + 截距1 (C <= 交點)，否則斜率2、截距2；
      斜率 /1e7、截距 /1e6 (NOAA 隨時間更新，已含衰減修正)。無效的掃描線用前後的值；整檔都無效時用發射前係數。
    - POD：反照率 = 斜率·C + 截距 (斜率 /2^30、截距 /2^22)；CLASS 1b 的係數多半是發射前值 (NOAA-14 1996/11 起逐月更新)。
    - KLM 衛星存成 POD 格式 (千葉 NOAA-18)：單一斜率的欄位只放得下低增益段，改用 NOAA STAR 的在軌係數
      (conversion_table/AVHRR，利比亞沙漠替代定標，每月一組；NOAA-16~19，取自 AAPP)，沒有時用附錄 D 的發射前雙增益係數。
  紅外 (亮溫：T* = c2·ν / ln(1 + c1·ν³/N)，T = (T* − A) / B；ν、A、B 見 avhrr_l1b_noaa.PLANCK)
    - KLM：N = a0 + a1·C + a2·C² (記錄中的業務係數，已含太空輻射 Ns 與非線性修正，不再另外修正)；
      a0、a1 /1e6，a2 在格式版本 <= 2 為 /1e6、>= 3 (NOAA-N 格式) 的頻道 4/5 為 /1e7。
      接收站產品 (千葉) 或係數整軌固定不變 (千葉 NOAA-17 的紅外係數與可見光都是固定值，不是逐線的星上定標)
      時，改用下面的自行定標 (以 CLASS 的 NOAA-15/18/19 驗證：與業務係數相差 ±0.05 K 以內，NOAA-15 1998 年 0.1 K)。
    - POD：1b 的係數只是太空點與星上黑體 (ICT) 兩點的線性定標，太空點的輻射依衛星與年代不同 (0 或負的 Ns；
      NOAA-9 在 1986~87 年改過幾次)，且非線性修正另外以表格提供。為了各衛星一致，改以 KLM 的模型重算：
        x = (C − C太空) / (C黑體 − C太空)              (C太空、C黑體：遙測的太空與黑體計數，取 51 條掃描線的中位數)
        N黑體 = 斜率·C黑體 + 截距                     (1b 的線性定標通過黑體點，等於 NOAA 由 PRT 算出的黑體輻射)
        N_lin = Ns + (N黑體 − Ns)·x，N = b0 + (1 + b1)·N_lin + b2·N_lin²   (Walton et al. 1998 的 Ns、b；NOAA-13/14 為指南值)
      NOAA-14 與指南的作法 (二次式) 完全相同；其他衛星與指南的 Ns = 0 修正表相差約 0.1~0.4 K。遙測無效時退回 1b 的線性。
      (POD 指南沒有 PRT 係數；pygac 的 POD PRT 係數換算的黑體溫度與 CLASS 1b 的黑體點相差 ±0.01 K 以內
      (NOAA-9/10/11/14)，用於沒有 1b 係數的原始框 AVHRR-L0。)
    - KLM 衛星存成 POD 格式 (千葉 NOAA-18)：檔中的紅外斜率不對 (頻道 4 只有正常值的約 1/3.5)，
      改由遙測的 PRT 計數 (附錄 D 的 PRT 係數) 求黑體溫度，依 KLM Guide §7.1.2.4 自行定標 (同 Ns、b 的非線性修正)。
      KLM 格式的檔若整檔沒有有效的紅外係數，也用同樣的方式 (PRT 係數取表頭)。
      PRT 每 5 條掃描線有一條參考讀數：KLM 指南寫 0，實際資料常是 1~5 (POD 衛星都不是 0)，三個讀數的中位數 < 20 即視為參考；
      接收站的資料常有個別讀數的位元錯誤，每條取三個讀數的中位數，四顆 PRT 各自去掉離群值、平滑後再平均。
定位
  - 51 個定位點：LAC/HRPT 在第 25、65、…、2025 像元，GAC 在第 5、13、…、405 像元 (1 起算)；GAC 像元 n 為
    LAC 樣本 5n−2 ~ 5n+1 的平均 (POD Guide)。沿掃描以地心角內插 (同 AVHRR-NSMC)。
  - 定位點無效、或與前後掃描線不連續 (位元錯誤) 的掃描線整條為無資料。
時間
  - 殘差「時間碼 − 記錄序 × 掃描週期」(LAC/HRPT 1/6 秒、GAC 0.5 秒) 在連續段內固定、資料缺口處跳一階；與它的
    15 條中位數相差超過 2 個掃描週期的時間碼視為受損，改用推算值 (千葉的檔在接收開頭常有受損的時間碼；CLASS 的檔也有
    個別位元錯誤，例如日數錯成 19)。掃描線序號在缺口處不可靠 (CLASS 的檔也會跳號或重複)，不使用。
資料受損
  - 計數值 0 為無資料。接收不良造成的條紋 (接收站資料的開頭、結尾) 與位元錯誤的雜點不另外處理，照原值畫出。
不同來源的合併 (merge_sources)
  - 同衛星的 LAC、HRPT、FRAC 解析度相同 (merge_key)，在列表中一起選取就合併成一張；GAC 另成一組。
  - 同一次過境的多個來源 (CLASS 的 LAC 與接收站的 HRPT、CLASS 與千葉…)：掃描線依時間對齊並以頻道 4 計數驗證
    (同一筆資料幾乎逐位元相同；時間差整數條掃描線時依計數校正)；每條掃描線選無資料與尖峰最少的來源
    (以整體最好的來源為主)，選中的掃描線裡無資料或尖峰的像元再以其他來源同一位置的好值補上。
  - 只有一個來源時不做任何挑選或修補。
"""
import lzma
import os
import re
import warnings
import zlib

import numpy as np
from scipy.ndimage import median_filter

import avhrr_hrpt as hrpt
import avhrr_l1b_noaa as fmt
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr
from satima.po import SatReader, obstime, swath, ui

SCAN_STEP_DEG = 0.0541      # LAC 取樣間隔 (度)
EARTH_KM = 6371.
MARGIN_DEG = 3.             # 讀檔前依定位點預選掃描線時，範圍外多帶的緯度
TELEMETRY_WINDOW = 51       # 太空 / 黑體計數、黑體溫度的平滑 (掃描線數)
MODE_WINDOW = 11            # 3A / 3B 判斷的中位數濾波 (掃描線數)
TIME_TOLERANCE = 2.        # 時間碼與推算值相差超過幾個掃描週期時，視為受損
TIME_WINDOW = 15            # 推算時間的中位數濾波 (掃描線數)
NAV_JUMP_DEG = 0.5          # 定位點與前後掃描線的中位數相差超過此值時，視為受損
PRT_REFERENCE = 20          # PRT 三個讀數的中位數小於此值的掃描線為參考讀數 (每 5 條一條)
PRT_SPREAD = 5              # 三個讀數中至少兩個與中位數相差在此值 (計數) 以內才採用
PRT_OUTLIER_K = 1.          # 同一顆 PRT 的溫度與前後 PRT_OUTLIER_WINDOW 個讀數的中位數相差超過此值 (K) 時捨去
PRT_OUTLIER_WINDOW = 11
# 多個來源合併 (見 merge_sources)
SPIKE_COUNTS = 50           # 計數與左右相鄰都相差超過此值 (同方向) 的像元視為尖峰
ALIGN_RANGE = 30            # 依計數對齊時，掃描線位移的搜尋範圍
SWITCH_MARGIN = 0.002       # 其他來源的分數 (無資料與尖峰的比例) 好過主要來源此值以上才換
SKIP_EXT = {'.png', '.jpg', '.jpeg', '.gif', '.txt', '.pdf', '.py', '.pyc', '.toml', '.nc', '.h5', '.he5', '.hdf',
            '.zip', '.7z', '.tar', '.tif', '.tiff', '.xml', '.json', '.md', '.csv', '.log', '.ini'}
# 壓縮檔 (hrpt.COMPRESSED：.gz / .bz2 / .xz / .lzma) 依內層的副檔名略過：SKIP_EXT，以及 HSD 的 .DAT.bz2 (數量多，不必逐一解壓)
SKIP_INNER_EXT = SKIP_EXT | {'.dat'}
COMPRESSED_HEAD = 2 ** 17   # 壓縮檔辨識時解壓的開頭長度 (含最長的前置檔頭、表頭與第一筆資料)
# USGS EarthExplorer 的檔名 = Entity ID：AO + 衛星編號 + mmddyy + hhmmss (例：ao14101797050456_012257.l1b)
USGS_ENTITY = re.compile(r'^ao\d{14}(?!\d)', re.I)

# 頻道 (選單編號，同 AVHRR-NSMC)：1 CH1、2 CH2、3 CH3A、4 CH3B (AVHRR/1、2 的頻道 3)、5 CH4、6 CH5
COUNT_INDEX = {1: 0, 2: 1, 3: 2, 4: 2, 5: 3, 6: 4}
IR_CHANNEL = {4: 3, 5: 4, 6: 5}     # 選單頻道 -> AVHRR 頻道
VIS_CHANNEL = {1: 1, 2: 2, 3: 3}
BT_BANDS = (4, 5, 6)
MENUS = {
    3: ('1.CH1(0.63) 2.CH2(0.86) 3.CH3A(1.61) 4.CH3B(3.74) 5.CH4(10.8) 6.CH5(12.0)', {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6},
        {1: '01', 2: '02', 3: '03A', 4: '03B', 5: '04', 6: '05'}),
    2: ('1.CH1(0.63) 2.CH2(0.86) 3.CH3(3.7) 4.CH4(10.8) 5.CH5(12.0)', {1: 1, 2: 2, 3: 4, 4: 5, 5: 6},
        {1: '01', 2: '02', 4: '03', 5: '04', 6: '05'}),
    1: ('1.CH1(0.63) 2.CH2(0.86) 3.CH3(3.7) 4.CH4(10.8)', {1: 1, 2: 2, 3: 4, 4: 5}, {1: '01', 2: '02', 4: '03', 5: '04'}),
}
SENSORS = {1: 'AVHRR', 2: 'AVHRR2', 3: 'AVHRR3'}

_SNIFF = {}


def is_compressed(path):
    return os.path.splitext(path)[1].lower() in hrpt.COMPRESSED


def sniff(path):
    """fmt.sniff 的快取 (資料夾掃描與建立檔案物件各做一次)。"""
    try:
        st = os.stat(path)
    except OSError:
        return None
    key = (os.path.abspath(path), st.st_size, st.st_mtime)
    if key not in _SNIFF:
        _SNIFF[key] = _sniff_compressed(path) if is_compressed(path) else fmt.sniff(path)
    return _SNIFF[key]


def _sniff_compressed(path):
    """壓縮檔：解壓開頭 COMPRESSED_HEAD bytes 判斷。解壓後的大小：整個檔都在開頭裡時即為其長度，
    .gz / .xz 讀容器記錄的大小，.bz2 不知道 (fmt.UNKNOWN_SIZE)。"""
    try:
        with hrpt.open_stream(path) as f:
            head = f.read(COMPRESSED_HEAD)
    except (OSError, EOFError, lzma.LZMAError, zlib.error, ValueError):
        return None
    if len(head) < COMPRESSED_HEAD:
        size = len(head)
    else:
        size = hrpt.raw_size(path)
        if size is None:
            size = fmt.UNKNOWN_SIZE
    return fmt.sniff(path, head=head, size=size)


# ===== 共用計算 =====
def earth_angle(cols, altitude_km):
    """LAC 像元座標 (0 起算，可為小數) -> 地心角 γ (弧度)：掃描角 θ = (像元 − 1023.5) × 0.0541°。"""
    th = np.deg2rad((np.asarray(cols, dtype=np.float64) - 1023.5) * SCAN_STEP_DEG)
    return np.arcsin((EARTH_KM + altitude_km) / EARTH_KM * np.sin(th)) - th


def interp_tiepoints(lats, lons, pix_cols, tie_cols, altitude_km, chunk=1024):
    """定位點 (列, 51) -> (列, 像元數) 的 (lats, lons)，float32：以地心角為座標在地心直角座標線性內插 (兩端外插)。"""
    g_tie = earth_angle(tie_cols, altitude_km)
    g_pix = earth_angle(pix_cols, altitude_km)
    seg = np.clip(np.searchsorted(g_tie, g_pix, side='right') - 1, 0, len(tie_cols) - 2)
    w = (g_pix - g_tie[seg]) / (g_tie[seg + 1] - g_tie[seg])
    lats = np.asarray(lats, dtype=np.float64)
    lons = np.asarray(lons, dtype=np.float64)
    out_lat = np.full((lats.shape[0], len(pix_cols)), np.nan, dtype=np.float32)
    out_lon = np.full_like(out_lat, np.nan)
    for r0 in range(0, lats.shape[0], chunk):
        la, lo = np.deg2rad(lats[r0:r0 + chunk]), np.deg2rad(lons[r0:r0 + chunk])
        xyz = [np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)]
        x, y, z = (v[:, seg] * (1. - w) + v[:, seg + 1] * w for v in xyz)
        out_lat[r0:r0 + chunk] = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
        out_lon[r0:r0 + chunk] = np.rad2deg(np.arctan2(y, x))
    return out_lat, out_lon


def fill_rows(values, valid):
    """無效的掃描線以前後的有效值線性內插 (兩端取端點)；values (列,) 或 (列, k)。沒有任何有效值時原樣回傳。"""
    values = np.asarray(values, dtype=np.float64)
    if valid.all() or not valid.any():
        return values
    idx = np.arange(values.shape[0])
    out = values.copy()
    if out.ndim == 1:
        out[~valid] = np.interp(idx[~valid], idx[valid], values[valid])
        return out
    for k in range(values.shape[1]):
        out[~valid, k] = np.interp(idx[~valid], idx[valid], values[valid, k])
    return out


def line_mean(samples, tolerance=15.):
    """每條掃描線 10 個太空 / 黑體樣本的平均：去掉 0、1023 與離中位數超過 tolerance 的樣本；有效樣本少於 3 個為 NaN。"""
    s = np.asarray(samples, dtype=np.float64).copy()
    s[(s <= 0) | (s >= 1023)] = np.nan
    with np.errstate(invalid='ignore'), warnings.catch_warnings():
        warnings.simplefilter('ignore', RuntimeWarning)
        med = np.nanmedian(s, axis=1)
        s[np.abs(s - med[:, None]) > tolerance] = np.nan
        n = np.isfinite(s).sum(axis=1)
        out = np.nanmean(s, axis=1)
    out[n < 3] = np.nan
    return out


def smooth_rows(v, window=TELEMETRY_WINDOW):
    """沿掃描線的穩健平滑：NaN 先以前後內插，再做中位數濾波；全為 NaN 時原樣回傳。"""
    ok = np.isfinite(v)
    if not ok.any():
        return v
    return median_filter(fill_rows(v, ok), size=window, mode='nearest')


def running_mean(v, window=TELEMETRY_WINDOW):
    """忽略 NaN 的移動平均 (視窗內沒有有效值時以前後內插)。"""
    ok = np.isfinite(v)
    if not ok.any():
        return v
    c = np.r_[0., np.cumsum(np.where(ok, v, 0.))]
    n = np.r_[0, np.cumsum(ok)]
    h = window // 2
    i = np.arange(v.size)
    a, b = np.clip(i - h, 0, v.size), np.clip(i + h + 1, 0, v.size)
    with np.errstate(invalid='ignore', divide='ignore'):
        m = (c[b] - c[a]) / (n[b] - n[a])
    return fill_rows(m, np.isfinite(m))


def spikes(counts):
    """沿掃描的尖峰 (位元錯誤)：與左右相鄰的計數都相差超過 SPIKE_COUNTS 且同方向。只用於多個來源之間挑選。"""
    c = np.asarray(counts, dtype=np.int32)
    out = np.zeros(c.shape, dtype=bool)
    if c.shape[-1] >= 3:
        a = c[:, 1:-1] - c[:, :-2]
        b = c[:, 1:-1] - c[:, 2:]
        out[:, 1:-1] = (np.abs(a) > SPIKE_COUNTS) & (np.abs(b) > SPIKE_COUNTS) & (a * b > 0)
    return out


def _identity(a, b):
    """兩組計數 (列, 像元) 逐像元相同的比例 (兩者皆非 0 的像元)。"""
    ok = (a > 0) & (b > 0)
    n = ok.sum()
    return float((ok & (a == b)).sum()) / n if n else 0.


def merge_sources(parts, period):
    """同一次過境的多個來源 (例如 CLASS 的 LAC 與接收站的 HRPT) 合併成一條掃描帶。
    parts：每個檔一筆 dict(name, data, lats, lons, times, counts, align, good)：data / lats / lons / counts (列, 像元)、
           times (列,) datetime64、align = 頻道 4 計數 (對齊用)、good = 定位有效的列；period：掃描週期 (毫秒)。
    1. 掃描線依時間對齊到 period 的格點；重疊的來源再以頻道 4 計數驗證 (同一筆資料應幾乎逐位元相同)，
       時間差了整數條掃描線時 (±ALIGN_RANGE 內) 依計數校正。
    2. 每條掃描線選用品質最好的來源：分數 = 無資料 (計數 0) 與尖峰的比例；以重疊段整體分數最好的來源為主，
       其他來源好 SWITCH_MARGIN 以上才換 (避免定位在來源之間來回跳)。
    3. 選中的掃描線裡無資料或尖峰的像元，以其他來源同一像元的好值補上。
    回傳 (data, lats, lons, times, 說明文字列表)。"""
    t0 = min(p['times'].min() for p in parts)
    for p in parts:
        p['slot'] = np.round((p['times'] - t0).astype('timedelta64[us]').astype(np.int64) / (period * 1000.)
                             ).astype(np.int64)
        p['spike'] = spikes(p['counts'])
        bad = (p['counts'] == 0).mean(axis=1) + p['spike'].mean(axis=1)
        p['score'] = np.where(p['good'] & np.isfinite(p['data']).any(axis=1), bad, np.inf)
    notes = []
    # 1. 對齊：每個來源與先前 (已對齊) 重疊最多的來源比對頻道 4 計數
    for i, p in enumerate(parts[1:], 1):
        q = max(parts[:i], key=lambda x: np.intersect1d(p['slot'], x['slot']).size)
        if np.intersect1d(p['slot'], q['slot']).size < 5 and not _near(p['slot'], q['slot']):
            continue
        qi = {int(s): k for k, s in enumerate(q['slot'])}
        sample = np.arange(0, p['slot'].size, max(1, p['slot'].size // 200))

        def ident(shift):
            pairs = [(int(k), qi[int(p['slot'][k]) + shift]) for k in sample if int(p['slot'][k]) + shift in qi]
            if len(pairs) < 5:
                return 0.
            return _identity(p['align'][[a for a, _ in pairs], ::4], q['align'][[b for _, b in pairs], ::4])

        i0 = ident(0)
        if i0 >= 0.5:
            continue
        scores = {s: ident(s) for s in range(-ALIGN_RANGE, ALIGN_RANGE + 1) if s}
        shift = max(scores, key=scores.get)
        if scores[shift] >= 0.5 and scores[shift] > i0 + 0.2:
            p['slot'] = p['slot'] + shift
            notes.append(tr('%s 的時間與 %s 差 %d 條掃描線 (%.2f 秒)，依計數對齊',
                            '%s: time offset from %s of %d scanlines (%.2f s); aligned by the counts')
                         % (p['name'], q['name'], shift, shift * period / 1000.))
    # 2. 每個格點的候選；整體品質排序 (重疊段的平均分數)
    slots = np.unique(np.concatenate([p['slot'] for p in parts]))
    index = {int(s): k for k, s in enumerate(slots)}
    cands = [[] for _ in slots]
    for i, p in enumerate(parts):
        for k, s in enumerate(p['slot']):
            cands[index[int(s)]].append((i, k))
    overall = []
    for p in parts:
        shared = np.array([len(cands[index[int(s)]]) > 1 for s in p['slot']])
        sc = p['score'][shared] if shared.any() else p['score']
        sc = sc[np.isfinite(sc)]
        overall.append(float(sc.mean()) if sc.size else np.inf)
    order = {int(i): n for n, i in enumerate(np.argsort(overall, kind='stable'))}
    primary = min(order, key=order.get)
    npix = parts[0]['data'].shape[1]
    out = np.full((slots.size, npix), np.nan, dtype=np.float32)
    lat = np.full((slots.size, npix), np.nan, dtype=np.float32)
    lon = np.full_like(lat, np.nan)
    times = np.empty(slots.size, dtype=parts[0]['times'].dtype)
    used = np.zeros(len(parts), dtype=np.int64)
    filled = 0
    for k, cs in enumerate(cands):
        cs = sorted(cs, key=lambda c: (parts[c[0]]['score'][c[1]], order[c[0]]))
        best = parts[cs[0][0]]['score'][cs[0][1]]
        choice = next((c for c in cs if c[0] == primary and parts[primary]['score'][c[1]] <= best + SWITCH_MARGIN),
                      cs[0])
        i, r = choice
        p = parts[i]
        times[k] = p['times'][r]
        if not np.isfinite(p['score'][r]):
            continue
        used[i] += 1
        out[k], lat[k], lon[k] = p['data'][r], p['lats'][r], p['lons'][r]
        need = ~np.isfinite(out[k]) | p['spike'][r]
        for j, rr in cs:
            if j == i or not need.any():
                continue
            q = parts[j]
            take = need & np.isfinite(q['data'][rr]) & ~q['spike'][rr]
            out[k][take] = q['data'][rr][take]
            need &= ~take
            filled += int(take.sum())
    keep = np.isfinite(lat).any(axis=1)
    sep = tr('、', ', ')
    notes.append(tr('合併 %d 個來源，採用的掃描線：%s；另以其他來源補上 %d 個無資料 / 尖峰像元',
                    'Merged %d sources; scanlines used: %s; %d missing / spike pixels filled from the other sources')
                 % (len(parts), sep.join('%s %d' % (p['name'], n) for p, n in zip(parts, used)), filled))
    return out[keep], lat[keep], lon[keep], times[keep], notes


def _near(a, b, gap=ALIGN_RANGE):
    """兩組格點是否在 gap 條掃描線內相鄰 (可能因時間偏移而看似不重疊)。"""
    return a.size and b.size and (min(a.max(), b.max()) + gap >= max(a.min(), b.min()))


def dual_gain_albedo(counts, s1, i1, s2, i2, x):
    """雙增益的反照率 (%)：C <= 交點用 (s1, i1)，否則 (s2, i2)；係數為 (列,) 陣列。"""
    c = counts.astype(np.float64)
    return np.where(c <= x[:, None], s1[:, None] * c + i1[:, None], s2[:, None] * c + i2[:, None])


# ===== 單一檔案 =====
class AvhrrL1bFile:
    compressed = False          # 壓縮檔：資料整個解壓進記憶體 (_load)
    _pending = False            # 壓縮檔解壓之前，nlines 只是表頭記錄的條數 (第一次用到 nlines 時解壓並以實際長度為準)

    def __init__(self, fullname):
        s = sniff(fullname)
        if s is None:
            raise ValueError(tr('不是 NOAA AVHRR 1b 檔', 'not a NOAA AVHRR 1b file'))
        self.fullname = fullname
        self.compressed = self._pending = is_compressed(fullname)
        self._buf = None
        self.s = s
        self.satellite = s.satellite
        self.generation = s.generation
        self.time = s.time
        self.end = s.end
        self.channel = s.dtype                  # 列表的產品欄：GAC / LAC / HRPT / FRAC
        self.station = s.station
        self.source = self._source(fullname, s)  # 列表的資料來源欄
        self.nlines = s.nlines
        self.npix = fmt.NPIX[s.kind]
        self.dtype = (fmt.klm_record if s.layout == fmt.KLM else fmt.pod_record)(s.kind)
        self._nav = None
        self._cache = {}
        self.label = os.path.basename(fullname)
        self.table_dir = None               # 校正表目錄 (STAR 可見光係數；讀取器建立檔案時設定)
        self.cfg = None                     # 讀取器的設定 (接收站 1b 以 TLE 重新定位時用；讀取器建立檔案時設定)

    def __lt__(self, other):
        return self.time < other.time

    @staticmethod
    def _source(fullname, s):
        """列表的資料來源 (只用於顯示)：接收站產品 (處理區塊 B0000000 等，樣本為千葉大學) 為 Chiba；
        USGS EarthExplorer 的檔 (檔名為 Entity ID，或有 122 bytes TBM 且資料集名稱的站碼為 SF = Sioux Falls) 為 USGS；其他為 CLASS。"""
        if s.station:
            return 'Chiba'
        name = getattr(s, 'name', None) or ''
        if USGS_ENTITY.match(os.path.basename(fullname)) or (s.prefix == 122 and name.upper().endswith('.SF')):
            return 'USGS'
        return 'CLASS'

    @property
    def klm(self):
        return self.s.layout == fmt.KLM

    @property
    def nlines(self):
        if self._pending:
            self._load()
        return self._nlines

    @nlines.setter
    def nlines(self, value):
        self._nlines = value

    def _load(self):
        """壓縮檔：整個解壓進記憶體 (不在磁碟留下解壓檔)；掃描線數以解壓後的實際長度為準。"""
        if self._buf is None:
            buf = hrpt.read_bytes(self.fullname)
            rec, first = self.s.record, self.s.first
            n = (len(buf) - first) // rec
            if n < 1:
                raise ValueError(tr('%s：解壓後不到一條掃描線', '%s: less than one scanline after decompression')
                                 % self.label)
            if (len(buf) - self.s.prefix) % rec:
                print(tr('%s：解壓後的長度不是記錄長度 (%d bytes) 的整數倍，略過最後不完整的記錄',
                         '%s: the decompressed length is not a multiple of the record length (%d bytes); the incomplete '
                         'last record is skipped') % (self.label, rec))
            self._buf = buf
            self._nlines = n
            self._pending = False
        return self._buf

    def release(self):
        """釋放解壓後的資料 (壓縮檔)；定位、遙測等已算好的結果留著，之後還要讀時再解壓。"""
        self._buf = None

    def _records(self):
        if self.compressed:
            buf = self._load()
            return np.frombuffer(buf, dtype=self.dtype, count=self._nlines, offset=self.s.first)
        return np.memmap(self.fullname, dtype=self.dtype, mode='r', offset=self.s.first, shape=(self.nlines,))

    def _field(self, name, r0=0, r1=None):
        rec = self._records()
        out = np.array(rec[name][r0:r1])
        del rec
        return out

    # ---- 時間與定位 ------------------------------------------------------------------------
    def navigation(self):
        """定位點經緯度 (列, 51)、觀測時間 (UTC datetime64[us])、有效列；讀一次後留著。"""
        if self._nav is not None:
            return self._nav
        rec = self._records()
        if self.klm:
            ll = np.array(rec['latlon'], dtype=np.float64) / 1e4
            year = np.array(rec['year'], dtype=np.int64)
            day = np.array(rec['day'], dtype=np.int64)
            no_nav_bit = 27
        else:
            ll = np.array(rec['latlon'], dtype=np.float64) / 128.
            yd = np.array(rec['yearday'], dtype=np.int64)
            year = np.array([fmt.full_year(int(y)) for y in (yd >> 9)], dtype=np.int64) if yd.size else yd
            day = yd & 0x1FF
            no_nav_bit = 26
        ms = np.array(rec['ms'], dtype=np.int64) & (0xFFFFFFFF if self.klm else 0x7FFFFFF)
        scan = np.array(rec['scan'], dtype=np.int64)
        quality = np.array(rec['quality'], dtype=np.int64)
        del rec
        lat, lon = ll[..., 0], ll[..., 1]
        good = ((np.abs(lat) <= 90.).all(axis=1) & (np.abs(lon) <= 180.).all(axis=1)
                & ~((lat == 0.) & (lon == 0.)).all(axis=1) & ((quality >> no_nav_bit) & 1 == 0))
        good &= self._continuous(lat, lon, good)
        times, fixed = self._times(year, day, ms, scan)
        renav = self._renavigate(times)
        if renav is not None:
            lat, lon = renav
            good = np.isfinite(lat).all(axis=1)
        lat[~good], lon[~good] = np.nan, np.nan
        self._nav = dict(lat=lat, lon=lon, time=times, good=good)
        bad = int((~good).sum())
        if bad or fixed:
            print(tr('%s：%d 條掃描線的定位無效、%d 條的時間碼受損 (改用推算值)',
                     '%s: %d scanlines with invalid navigation, %d with damaged time codes (replaced by the expected time)')
                  % (self.label, bad, fixed))
        return self._nav

    def _renavigate(self, times):
        """接收站產品 (千葉等) 的 LAC/HRPT：以 TLE 重新定位 (同 AVHRR-L0：大地天底 + 各衛星的定位參數)，
        取代檔內的定位點；[tle] station_l1b = false、GAC 或沒有 TLE 時 None (用檔內的定位)。
        千葉檔內的定位與海岸線相差 3.5~6 km (NOAA-12)、掃描邊緣到 8~11 km (NOAA-19 沒有修正 yaw)。"""
        cfg = self.cfg
        if not self.station or self.s.kind == 'GAC' or cfg is None or not cfg.get('tle') \
                or not cfg.tle.get('station_l1b', True):
            return None
        try:
            lat, lon = hrpt.navigate(cfg, self.satellite, self.time, times, 24. + 40. * np.arange(51), self.label)
        except ValueError as err:
            print(tr('%s (改用檔內的定位)', '%s (using the navigation in the file)') % err)
            return None
        print(tr('%s：接收站產品，以 TLE 重新定位 ([tle] station_l1b)', '%s: station product; navigated from the TLE '
                 '([tle] station_l1b)') % self.label)
        return lat, lon

    @staticmethod
    def _continuous(lat, lon, good):
        """中間定位點與前後 MODE_WINDOW 條掃描線的中位數相差 > NAV_JUMP_DEG 的掃描線視為受損 (位元錯誤)。"""
        ok = good.copy()
        if ok.sum() < MODE_WINDOW:
            return ok
        lon_c = lon[:, 25].copy()
        lon_c[good] = np.rad2deg(np.unwrap(np.deg2rad(lon_c[good])))
        for v in (lat[:, 25], lon_c):
            ref = median_filter(fill_rows(v, good), size=MODE_WINDOW, mode='nearest')
            ok &= np.abs(v - ref) <= NAV_JUMP_DEG
        return ok

    def _times(self, year, day, ms, scan):
        """每條掃描線的時間 (datetime64[us])；受損的時間碼改用推算值。回傳 (時間, 改掉的條數)。
        殘差 r = 時間 − 記錄序 × 掃描週期 在連續段內為常數、資料缺口處跳一階；r 的 TIME_WINDOW 條中位數濾波保留階梯，
        與它相差超過 TIME_TOLERANCE 個掃描週期的時間碼視為受損 (含檔頭、檔尾)。掃描線序號在缺口處不可靠 (CLASS 也會跳號、重複)，不使用。"""
        t0 = np.datetime64(self.time, 'ms')
        good = (year > 1970) & (day >= 1) & (day <= 366) & (ms >= 0) & (ms < 86_400_000)
        jan1 = (np.where(good, year, 2000) - 1970).astype('datetime64[Y]').astype('datetime64[ms]')
        t = jan1 + ((day - 1) * 86_400_000 + ms).astype('timedelta64[ms]')
        rel = (t - t0).astype(np.int64).astype(np.float64)
        rel[~good] = np.nan
        period = 500. if self.s.kind == 'GAC' else 1000. / 6.
        steps = np.arange(rel.size) * period
        r = rel - steps                     # 連續段內為常數，資料缺口處往上跳一階
        ok = np.isfinite(r)
        if not ok.any():
            return (t0 + np.round(steps).astype('timedelta64[ms]')).astype(obstime.US), rel.size
        tol = TIME_TOLERANCE * period
        bad = ~ok
        for _ in range(2):                  # 第二輪以換掉明顯壞值後的殘差為參考
            ref = median_filter(fill_rows(r, ~bad), size=TIME_WINDOW, mode='mirror')
            bad = ~(np.abs(r - ref) <= tol)
            if bad.all():
                break
        rel = np.where(bad, steps + ref, rel)
        times = (t0 + np.round(rel).astype('timedelta64[ms]')).astype(obstime.US)
        return times, int(bad.sum())

    def rows_in(self, georange, clon):
        """與地圖範圍可能相交的掃描線 (起, 訖)；沒有時 None。"""
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

    def geocoord(self, rows=None):
        r0, r1 = rows or (0, self.nlines)
        nav = self.navigation()
        if self.s.kind == 'GAC':
            pix = 5. * np.arange(self.npix) + 3.5            # GAC 像元 n 的中心在 LAC 座標 5n + 3.5 (0 起算)
            tie = 5. * (4 + 8 * np.arange(51)) + 3.5
        else:
            pix = np.arange(self.npix, dtype=np.float64)
            tie = 24. + 40. * np.arange(51)
        return interp_tiepoints(nav['lat'][r0:r1], nav['lon'][r0:r1], pix, tie,
                                fmt.ALTITUDE_KM.get(self.satellite, 850.))

    def row_times(self, rows=None):
        r0, r1 = rows or (0, self.nlines)
        return self.navigation()['time'][r0:r1]

    # ---- 計數值與遙測 ----------------------------------------------------------------------
    def counts(self, band, r0, r1):
        w = self._field('video', r0, r1)
        vals = fmt.unpack10(w, self.npix * 5).reshape(-1, self.npix, 5)
        return vals[:, :, COUNT_INDEX[band]]

    def telemetry(self):
        """遙測 (全部掃描線)：prt (列, 3)、ict {頻道 3/4/5: (列, 10)}、space {頻道 1~5: (列, 10)}。"""
        if 'tm' not in self._cache:
            if self.klm:
                prt, ict, space = self._field('prt'), self._field('ict'), self._field('space')
            else:
                t = fmt.unpack10(self._field('telemetry'), 103)
                prt, ict, space = t[:, 17:20], t[:, 22:52], t[:, 52:102]     # word 18~20、23~52、53~102 (1 起算)
            self._cache['tm'] = dict(prt=prt.astype(np.int64),
                                     ict={c: ict[:, k::3].astype(np.int64) for k, c in enumerate((3, 4, 5))},
                                     space={c: space[:, c - 1::5].astype(np.int64) for c in range(1, 6)})
        return self._cache['tm']

    def ch3_mode(self):
        """每條掃描線的通道 3：True = 3A、False = 3B；AVHRR/1、2 全為 False。另回傳切換中 (兩者皆無效) 的掃描線。"""
        if 'mode' not in self._cache:
            n = self.nlines
            transition = np.zeros(n, dtype=bool)
            if self.generation < 3:
                is3a = np.zeros(n, dtype=bool)
            elif self.klm:
                sel = self._field('bitfield') & 3
                is3a, transition = sel == 1, sel == 2
            else:                           # KLM 衛星存成 POD 格式：頻道 3 的定標斜率 正 = 3A (反照率)、負 = 3B
                slope = self._field('cal')[:, 2, 0].astype(np.float64)
                v = np.where(slope > 0, 1., np.where(slope < 0, 0., np.nan))
                ok = np.isfinite(v)
                is3a = (median_filter(fill_rows(v, ok), size=MODE_WINDOW, mode='nearest') > 0.5) if ok.any() \
                    else np.zeros(n, dtype=bool)
            self._cache['mode'] = (is3a, transition)
        return self._cache['mode']

    # ---- 定標 ------------------------------------------------------------------------------
    def _vis_coefficients(self, ch):
        """可見光頻道 ch (1、2、3 = 3A) 每條掃描線的 (s1, i1, s2, i2, 交點)，或 POD 的 (s, i)。"""
        key = ('vis', ch)
        if key in self._cache:
            return self._cache[key]
        sat = self.satellite
        if self.klm:
            raw = self._field('vis_cal')[:, ch - 1].astype(np.float64)      # (列, 業務/測試/發射前, 5)
            out = None
            for k in (0, 2):
                c = raw[:, k] / np.array([1e7, 1e6, 1e7, 1e6, 1.])
                ok = (c[:, 0] > 0) & (c[:, 2] > c[:, 0]) & (c[:, 4] > 0) & (c[:, 4] < 1023)
                if ok.any():
                    c = fill_rows(c, ok)
                    out = tuple(c[:, j] for j in range(5))
                    if k:
                        print(tr('%s CH%d：沒有業務定標係數，改用記錄中的發射前係數', '%s CH%d: no operational coefficients; '
                                 'using the pre-launch coefficients in the records') % (self.label, ch))
                    break
            if out is None:
                out = self._vis_fallback(ch, tr('檔中沒有可見光係數', 'no visible coefficients in the file'))
        elif self.generation == 3:          # KLM 衛星存成 POD 格式
            out = self._vis_fallback(ch, tr('POD 格式只有單一斜率', 'the POD format holds a single slope'))
        else:
            cal = self._field('cal')[:, ch - 1].astype(np.float64)
            s, i = cal[:, 0] / 2 ** 30, cal[:, 1] / 2 ** 22
            ok = s > 0
            out = (fill_rows(s, ok), fill_rows(i, ok)) if ok.any() else None
        self._cache[key] = out
        return out

    def _vis_fallback(self, ch, why):
        """檔中沒有可用的可見光係數時：NOAA STAR 的在軌係數 (conversion_table/AVHRR，見 fmt.star_vis)，
        沒有時用附錄 D 的發射前雙增益係數。回傳 (s1, i1, s2, i2, 交點) 各為 (列,) 陣列，或 None。"""
        star = fmt.star_vis(self.table_dir, self.satellite, ch, self.time)
        if star is not None:
            coef, note = star
            print(tr('%s CH%d：%s，改用 NOAA STAR 的在軌係數 (%s)', '%s CH%d: %s; using the NOAA STAR in-orbit coefficients (%s)')
                  % (self.label, ch, why, note))
        else:
            coef = fmt.VIS_PRELAUNCH.get(self.satellite, {}).get(ch)
            if coef is None:
                return None
            print(tr('%s CH%d：%s，改用附錄 D 的發射前雙增益係數 (沒有在軌衰減修正)',
                     '%s CH%d: %s; using the pre-launch dual-gain coefficients of Appendix D (no degradation correction)')
                  % (self.label, ch, why))
        return tuple(np.full(self.nlines, v, dtype=np.float64) for v in coef)

    def _ir_operational(self, ch):
        """KLM 記錄的紅外業務係數 (列, 3) = [a0, a1, a2]；整檔無效時 None。"""
        raw = self._field('ir_cal')[:, ch - 3, 0].astype(np.float64)
        a2_scale = 1e7 if (self.s.version >= 3 and ch in (4, 5)) else 1e6
        c = raw / np.array([1e6, 1e6, a2_scale])
        ok = (c[:, 1] < 0) & (np.abs(raw) < 2 ** 31 - 1).all(axis=1)
        if ch == 3:
            ok &= ~self.ch3_mode()[0]
        if not ok.any():
            return None
        return fill_rows(c, ok)

    def _prt_coefficients(self):
        """四顆 PRT 的係數 (4, 6) (T = d0 + d1·C + … + d5·C⁵)；KLM 檔取表頭，其他取 fmt.PRT；沒有時 None。"""
        if self.klm:
            raw = np.frombuffer(self.s.header[200:248], '>i2').reshape(4, 6).astype(np.float64)
            c = raw / 10. ** np.array([2, 5, 8, 11, 14, 17])
            if np.all((c[:, 0] > 250) & (c[:, 0] < 300)):
                return c
        if self.satellite in fmt.PRT:
            return np.array([list(p) + [0.] * (6 - len(p)) for p in fmt.PRT[self.satellite]])
        return None

    def _blackbody(self):
        """星上黑體溫度 (列,) (K)：PRT 計數 (每第 5 條掃描線是參考讀數，之後依序為 PRT1~4) 換算成溫度。
        每條掃描線取三個讀數的中位數 (接收站資料常有個別讀數的位元錯誤)；四顆 PRT 各自去掉與前後中位數相差
        PRT_OUTLIER_K 以上的值、平滑 (TELEMETRY_WINDOW 條) 後再平均 (某顆 PRT 缺值時不會偏向其他幾顆)。"""
        if 'tbb' in self._cache:
            return self._cache['tbb']
        coefs = self._prt_coefficients()
        if coefs is None:
            self._cache['tbb'] = None
            return None
        prt = self.telemetry()['prt']
        med = np.median(prt, axis=1)
        agree = (np.abs(prt - med[:, None]) <= PRT_SPREAD).sum(axis=1) >= 2
        ref = med < PRT_REFERENCE
        val = np.where(~ref & agree, med, np.nan)
        n = prt.shape[0]
        if ref.any():
            last = np.maximum.accumulate(np.where(ref, np.arange(n), -10 ** 9))
            k = np.arange(n) - last                         # 1~4 = PRT1~4
            groups = [(k == j + 1, coefs[j]) for j in range(4)]
        else:                                               # 找不到參考掃描線：四顆 PRT 的平均係數 (KLM 相差 < 0.1 K)
            groups = [(np.ones(n, dtype=bool), coefs.mean(axis=0))]
        series = []
        for sel, c in groups:
            sel = sel & np.isfinite(val)
            if sel.sum() < 3:
                continue
            temp = np.full(n, np.nan)
            temp[sel] = np.polyval(c[::-1], val[sel])
            temp[(temp < 260.) | (temp > 320.)] = np.nan
            ok = np.flatnonzero(np.isfinite(temp))
            if ok.size >= 3:
                ref_t = median_filter(temp[ok], size=PRT_OUTLIER_WINDOW, mode='nearest')
                temp[ok[np.abs(temp[ok] - ref_t) > PRT_OUTLIER_K]] = np.nan
            if np.isfinite(temp).any():
                series.append(running_mean(temp))
        tbb = np.mean(series, axis=0) if series else None
        self._cache['tbb'] = tbb if tbb is not None and np.isfinite(tbb).any() else None
        return self._cache['tbb']

    def _space_ict(self, ch):
        """頻道 ch 每條掃描線 (平滑後) 的太空與黑體計數。"""
        key = ('sp', ch)
        if key not in self._cache:
            tm = self.telemetry()
            self._cache[key] = (smooth_rows(line_mean(tm['space'][ch])), smooth_rows(line_mean(tm['ict'][ch])))
        return self._cache[key]

    def _raw_possible(self, ch):
        """能否由星上黑體自行定標：有 PRT 黑體溫度與太空、黑體計數。"""
        sp, ict = self._space_ict(ch)
        return self._blackbody() is not None and np.isfinite(sp).any() and np.isfinite(ict).any()

    def _ir_radiance(self, ch, counts, r0, r1):
        """紅外頻道 ch (3、4、5) 的輻射 (列, 像元) mW/(m2 sr cm-1)；無法定標時 None。"""
        sat = self.satellite
        c = counts.astype(np.float64)
        ns, b0, b1, b2 = fmt.NONLINEAR.get(sat, {}).get(ch, (0., 0., 0., 0.))
        method = self._cache.get(('irmethod', ch))
        if self.klm:
            coef = self._cache.get(('ir', ch))
            if coef is None and method is None:
                coef = self._ir_operational(ch)
                self._cache[('ir', ch)] = coef
                # 係數整軌固定 (不是逐線的星上定標，例如千葉的檔) 或為接收站產品：能自行定標時以星上黑體定標
                fixed = coef is not None and np.ptp(coef[:, 1]) == 0.
                if (fixed or self.station) and self._raw_possible(ch):
                    method = 'raw'
                    print(tr('%s CH%d：%s，改由 PRT、黑體與太空計數自行定標 (KLM Guide §7.1.2.4)',
                             '%s CH%d: %s; calibrating from the PRT, blackbody and space counts (KLM Guide 7.1.2.4)')
                          % (self.label, ch, tr('檔中的紅外係數整軌固定', 'the IR coefficients are constant over the pass')
                             if fixed else tr('接收站產品', 'station product')))
                    self._cache[('rawmsg', ch)] = True
                else:
                    method = 'operational' if coef is not None else 'raw'
                self._cache[('irmethod', ch)] = method
            if method == 'operational':
                a = coef[r0:r1]
                return a[:, 0:1] + a[:, 1:2] * c + a[:, 2:3] * c * c
        elif self.generation < 3 and not self._station_raw(ch):
            sp, ict = self._space_ict(ch)
            cal = self._field('cal')[:, ch - 1].astype(np.float64)
            s, i = cal[:, 0] / 2 ** 30, cal[:, 1] / 2 ** 22
            ok = s < 0
            if not ok.any():
                return None
            s, i = fill_rows(s, ok), fill_rows(i, ok)
            if np.isfinite(sp).any() and np.isfinite(ict).any():
                if method is None:
                    self._cache[('irmethod', ch)] = 'walton'
                n_ict = (s * ict + i)[r0:r1, None]
                x = (c - sp[r0:r1, None]) / (ict - sp)[r0:r1, None]
                n_lin = ns + (n_ict - ns) * x
                return b0 + (1. + b1) * n_lin + b2 * n_lin * n_lin
            if method is None:
                self._cache[('irmethod', ch)] = 'linear'
                print(tr('%s CH%d：遙測中沒有太空 / 黑體計數，只用 1b 的線性定標 (不做非線性修正)',
                         '%s CH%d: no space / blackbody counts in the telemetry; using the 1b linear calibration only '
                         '(no nonlinearity correction)') % (self.label, ch))
            return s[r0:r1, None] * c + i[r0:r1, None]
        # 由星上黑體自行定標 (KLM 衛星存成 POD 格式，或 KLM 檔沒有有效的係數)
        if method == 'raw' and ('rawmsg', ch) not in self._cache:
            self._cache[('rawmsg', ch)] = True
            print(tr('%s CH%d：檔中沒有可用的紅外定標係數，由 PRT、黑體與太空計數自行定標 (KLM Guide §7.1.2.4)',
                     '%s CH%d: no usable IR calibration coefficients; calibrating from the PRT, blackbody and space counts '
                     '(KLM Guide 7.1.2.4)') % (self.label, ch))
        self._cache[('irmethod', ch)] = 'raw'
        return self._ir_raw(ch, c, r0, r1)

    def _station_raw(self, ch):
        """POD 格式的接收站產品 (千葉 NOAA-12/14 等)：能由 PRT 自行定標時改用它 (同 AVHRR-L0 與 AAPP 的 PRT 係數)。
        千葉 1b 係數隱含的黑體溫度比 PRT 換算的低 0.3~0.7 K (NOAA-12)；CLASS 的 1b 則與 PRT 一致 (±0.01 K)。"""
        if not self.station:
            return False
        if ('stationraw', ch) not in self._cache:
            ok = self._raw_possible(ch)
            self._cache[('stationraw', ch)] = ok
            if ok:
                self._cache[('rawmsg', ch)] = True
                print(tr('%s CH%d：接收站產品，改由 PRT、黑體與太空計數自行定標 (KLM Guide §7.1.2.4)',
                         '%s CH%d: station product; calibrating from the PRT, blackbody and space counts (KLM Guide 7.1.2.4)')
                      % (self.label, ch))
        return self._cache[('stationraw', ch)]

    def _ir_raw(self, ch, c, r0, r1):
        """由 PRT 黑體溫度、黑體與太空計數自行定標 (KLM Guide §7.1.2.4)：
        N_lin = Ns + (N黑體 − Ns)·(C太空 − C) / (C太空 − C黑體)，N = b0 + (1 + b1)·N_lin + b2·N_lin²。無法定標時 None。"""
        tbb = self._blackbody()
        sp, ict = self._space_ict(ch)
        if tbb is None or not (np.isfinite(sp).any() and np.isfinite(ict).any()):
            return None
        ns, b0, b1, b2 = fmt.NONLINEAR.get(self.satellite, {}).get(ch, (0., 0., 0., 0.))
        nu, a, b = fmt.PLANCK[self.satellite][ch]
        n_bb = fmt.planck_radiance(tbb, nu, a, b)[r0:r1, None]
        spr, ictr = sp[r0:r1, None], ict[r0:r1, None]
        n_lin = ns + (n_bb - ns) * (spr - c) / (spr - ictr)
        return b0 + (1. + b1) * n_lin + b2 * n_lin * n_lin

    def calibrate(self, band, r0, r1):
        """反射率 (0~1) 或亮溫 (°C)，(列, 像元) float32；無效為 NaN。"""
        counts = self.counts(band, r0, r1)
        bad = (counts == 0) | ~self.navigation()['good'][r0:r1, None]
        is3a, transition = self.ch3_mode()
        if band in (3, 4) and self.generation == 3:
            bad |= (~is3a if band == 3 else is3a)[r0:r1, None] | transition[r0:r1, None]
        if band in VIS_CHANNEL:
            coef = self._vis_coefficients(VIS_CHANNEL[band])
            if coef is None:
                return np.full(counts.shape, np.nan, dtype=np.float32)
            coef = tuple(v[r0:r1] for v in coef)
            if len(coef) == 5:
                albedo = dual_gain_albedo(counts, *coef)
            else:
                albedo = coef[0][:, None] * counts + coef[1][:, None]
            out = albedo / 100.
        else:
            ch = IR_CHANNEL[band]
            rad = self._ir_radiance(ch, counts, r0, r1)
            if rad is None:
                print(tr('%s CH%d：無法定標', '%s CH%d: cannot be calibrated') % (self.label, ch))
                return np.full(counts.shape, np.nan, dtype=np.float32)
            nu, a, b = fmt.PLANCK[self.satellite][ch]
            out = fmt.planck_bt(rad, nu, a, b) - 273.15
        out = np.asarray(out, dtype=np.float32)
        out[bad] = np.nan
        return out

    def extract(self, band, rows=None):
        r0, r1 = rows or (0, self.nlines)
        return self.calibrate(band, r0, r1)


# ===== 讀取器 =====
class AvhrrL1bReader(SatReader):
    NAME = 'AVHRR-L1B'
    SENSOR = 'AVHRR'
    DESCRIPTION = 'TIROS-N~NOAA-19、MetOp AVHRR NOAA 1b (CLASS GAC/LAC/HRPT/FRAC、千葉大學 HRPT，可為 .gz/.bz2/.xz；依內容辨識)'
    DESCRIPTION_EN = ('TIROS-N~NOAA-19, MetOp AVHRR NOAA 1b (CLASS GAC/LAC/HRPT/FRAC, Chiba Univ. HRPT, also .gz/.bz2/.xz; '
                      'detected by content)')
    SOURCE = 'CLASS'
    PATTERNS = ('*',)

    @classmethod
    def find(cls, folder=None):
        """資料夾中所有依內容判斷為 NOAA AVHRR 1b 的檔 (不看檔名)。"""
        root = folder or '.'
        names = []
        try:
            entries = os.listdir(root)
        except OSError:
            return []
        for name in entries:
            ext = os.path.splitext(name)[1].lower()
            if ext in hrpt.COMPRESSED:
                if hrpt._inner_ext(name) in SKIP_INNER_EXT or name.lower().endswith('.tar.' + ext[1:]):
                    continue
            elif ext in SKIP_EXT:
                continue
            path = os.path.join(root, name)
            if os.path.isfile(path) and sniff(path) is not None:
                names.append(name)
        # 同一個檔的壓縮與未壓縮版本 (X 與 X.bz2) 都在時只列未壓縮的
        have = set(names)
        names = [n for n in names if not (is_compressed(n) and os.path.splitext(n)[0] in have)]
        return sorted(names, key=cls.sort_key)

    def make_file(self, name):
        f = AvhrrL1bFile(name)
        f.table_dir = self.cfg.general.conversion_table_dir
        f.cfg = self.cfg
        return f

    def listing(self, f):
        return (f.satellite, f.channel, f.time)

    def merge_key(self, segment):
        """一次選到多段時：同衛星的 LAC、HRPT、FRAC (解析度相同) 可以合併，GAC 另成一組；
        同一次過境的不同來源由 merge_sources 逐條掃描線挑選較好的。"""
        return (segment.satellite, 'GAC' if segment.product == 'GAC' else 'LAC/HRPT/FRAC')

    def source_label(self, segment):
        # 各檔的 source (見 AvhrrL1bFile._source)；一段有多個來源時依 Chiba、USGS、CLASS 的順序取第一個
        sources = {getattr(f, 'source', 'CLASS') for f in segment.files}
        return next(s for s in ('Chiba', 'USGS', 'CLASS') if s in sources or s == 'CLASS')

    def begin(self, segment):
        super().begin(segment)
        self.generation = max(f.generation for f in self.files)
        self.SENSOR = SENSORS[self.generation]

    def select_band(self):
        for f in self.files:
            print(f.fullname)
        menu, mapping, self.labels = MENUS[self.generation]
        print(menu)
        self.band = mapping[ui.ask_int(tr('輸入頻道編號：', 'Input the band: '), sorted(mapping))]
        if self.generation == 3 and self.band in (3, 4):
            n3a = sum(int(f.ch3_mode()[0].sum()) for f in self.files)
            n = sum(f.nlines for f in self.files)
            print(tr('通道 3：3A %d 條、3B %d 條掃描線', 'Channel 3: %d 3A and %d 3B scanlines') % (n3a, n - n3a))

    def select_variants(self):
        if self.band in BT_BANDS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare(self):
        try:
            return self._prepare()
        finally:
            for f in self.files:            # 壓縮檔解壓進記憶體的資料 (畫圖只用 self.data)
                f.release()

    def _prepare(self):
        parts = []
        for f in self.files:
            rows = f.rows_in(self.georange, self.center[1])
            if rows is None:
                continue
            lats, lons = f.geocoord(rows)
            parts.append(dict(name=f.label, data=f.extract(self.band, rows), lats=lats, lons=lons, times=f.row_times(rows),
                              counts=f.counts(self.band, *rows), align=f.counts(5, *rows),
                              good=f.navigation()['good'][rows[0]:rows[1]], kind=f.s.kind))
        if not parts:
            return False
        if len({p['kind'] for p in parts}) > 1:
            raise ValueError(tr('GAC 與 LAC 不能接在一起', 'GAC and LAC cannot be joined'))
        if len(parts) == 1:
            p = parts[0]
            data, lats, lons, times = p['data'], p['lats'], p['lons'], p['times']
        else:
            period = 500. if parts[0]['kind'] == 'GAC' else 1000. / 6.
            data, lats, lons, times, notes = merge_sources(parts, period)
            for line in notes:
                print(line)
        if self.slice_swath(lats, lons) is None:
            return False
        self.data = self.crop(data)
        self.time = self.center_time(times)
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
        channel = self.labels[self.band]
        return dict(band='BAND-%s' % channel, colorbar=cm, label='BAND-%s' % channel)


READER = AvhrrL1bReader
