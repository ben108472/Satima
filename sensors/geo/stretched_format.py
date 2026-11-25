"""SMS-1/2、GOES-1~7 的 Mode A / AA / AAA 原始資料 (stretched VISSR / VAS) 的格式與解碼
— GOES-STRETCHED 讀取器用 (檔名沒有 '-'，不會被當成讀取器載入)

依據 NOAA 的格式文件：
  Mode A    VISSR Digital Archive Users' Guide (NCC 1978) 的文件欄位與 O/A、
            MODE-A Modified IR Documentation (1982/10/01) 的 128 字 IR 文件
  Mode AA   GOES Mode AA Output Format Rev. B / Rev. C (VAS，GOES-4~6)
  Mode AAA  Operational VAS Mode AAA Format SFP 002 v2.5A (1987/09，GOES-7)

NCEI 2024 重處理版的檔案 (goesNN.yyyy.ddd.hhmm.A.data、.AA-MSI.data、.AAA-DS.data…)：
  一筆記錄 = 一圈 (自旋約 0.6 秒) 去掉同步碼 (PN 序列、AAA 的區塊標頭與 FCS) 後的資料，位元由高到低連續排列，
  各區塊都從 byte 邊界開始。記錄大致等長，但偶有不完整的記錄 (例如 AA-DS 樣本中缺 VIS 第 1 段、之後的區塊整體前移的一筆，
  以及檔尾截斷的記錄)：有同名 .indx (NCEI 附的區塊索引) 時依它逐區塊定位，否則以 VIS 段碼定位完整的記錄 (RawFile)。
  同名的 .meta (JSON) 記錄處理步驟、導航來源 (navigation_source：native = 原始資料流的 O&A)、整份檔 O&A 的眾數
  (與 assemble() 的結果相同) 與外推到第 0 掃描的時刻 (= 導航的 pictim)；讀取器只把它放進 header()。

  Mode A    99244 bytes  位元 0~1151      IR 文件 128 字 x 9 位元 (1 個旗標位元 + 8 位元資料)
                         位元 1152~       IR 3818 個樣本 x 9 位元 (8 位元的已校正計數值；第 3819~3822 樣本為行尾標記)
                         byte 4444 起     VIS 8 段 x 11850 bytes
  Mode AA   109536 bytes 位元 0           IR1 (上偵測器)：文件 16 + 影像 3822 + EOLC 8 + 共同文件 2048 字 x 10 位元，再 4 位元補零
                         位元 58944       IR2 (下偵測器)，同上
                         byte 14736 起    VIS 8 段
  Mode AAA  109516 bytes 位元 0 / 58864   兩個 IR 區塊各 58864 位元：文件 16 + 影像 3822 + EOLC 8 + 共同文件 512
                                          + 格網 1024 + Mode A 文件 128 + 備用 (截斷)，10 位元字
                         byte 14716 起    VIS 8 段
  VIS 段：文件 512 字 (384 bytes) + 影像 15288 字 (11466 bytes)，6 位元；文件第 1~3 字為段碼 (ONE = 111110、ZERO = 000001)。
  文件中的 ONE / ZERO 旗標以 $FE / $01 表示。

VIS 線號 (由 1 起算，同 McIDAS AREA：掃描 N 的 8 條線為 (N-1) x 8 + 1 ~ + 8；s / N = 記錄的掃描數，k / j = 段 0~7)：
  AAA   (N-1) x 8 + 1 + j：每段 VIS 文件第 9~12 字記錄該線的掃描數 = N (表 9)；與 CLASS 的 SSEC AREA_VIS 逐線相同。
  A     (s-1) x 8 + 5 + k：段 0~3 = 掃描 s 的 V5~V8、段 4~7 = 掃描 s+1 的 V1~V4；與 CLASS 的 SSEC AREA_VIS 逐線相同
        (SMS-2 1979、GOES-6 1985；pregvar 以 GOES-1 也得到 AREA 線 = 8 (s-1) + 偵測器 + 4)。
  AA    N x 8 + 5 + j：同時次的 Mode A 掃描 s 第 k 段與 AA 掃描 s-1 第 k 段是同一圈 (GOES-6 1985 19:31、20:31，
        24 線區塊相關 0.98，80 個區塊全部對在同一線；兩者由不同的地面處理產生，像素值不完全相同)。
        (CLASS 的 SSEC AA AREA_VIS 是 8 (N-1) + 5 + j，比這裡北 8 線；但同組 AA AREA_IR 的線號與這裡相同、導航區塊在
        同一線號的位置與這裡相差 < 0.5 線，所以是 SSEC 的 AA 可見光偏了一個掃描，不採用。)
  依據：在這個編號下，GOES-6 (Mode A、AA) 與 SMS-2 的 IR 與 VIS 南北配準都在 2 線內 (AAA +2 線)，可見光海岸線殘差 +2、-6 線；
  若照 Rev. C 表 2 的字面 (段 0~3 屬於記錄的掃描 N) 排，配準差 9~11 線、殘差 -6、-14 線。
IR 線位置 (掃描步，一步 = 8 條 VIS 線)：
  A     記錄的掃描數 s。
  AA    小偵測器 (HgCdTe 小) 上 = s、下 = s+1；大偵測器 (HgCdTe 大、InSb) 各涵蓋 2 步，上 = s、s+1，下 = s+4、s+5。
        (以影像相關性量得，與 AAA 的 adjusted earth count 一致；各格式、各波段的位置與 CLASS 的 SSEC AREA_IR 逐像素相同。)
  AAA   IR 文件第 9~10 字的 adjusted earth count (已依上下偵測器位移修正) + 整檔的 (預測掃描數 - earth count) 中位數；
        大偵測器涵蓋 2 步；adjusted earth count 為 0 時改用 AA 的規則。
截短的區塊：NCEI 2024 年以後的處理 (.meta 的 processing_stages 有 zero-padded) 把截短的區塊補零到完整長度，
  並在 .meta 的 truncated_line_bits 記下 {.indx 行號 (由 1 起算): 有效位元數}；補零的部分 (ir_valid、vis_valid 之後) 為 NaN。
"""
import json
import mmap
import os

import numpy as np

# ---- 記錄配置 -----------------------------------------------------------------------------------
VIS_SECTOR = 11850              # VIS 一段的 bytes (文件 384 + 影像 11466)
VIS_DOC = 384
NVIS = 15288                    # VIS 每線樣本數 (= 4 x 3822)
NIR = 3822                      # IR 每線樣本數


class Layout:
    def __init__(self, mode, rec, vis0, ir_slots, ir_bits, ir_doc_words, ir_valid, ir_bytes, ir_types, vis_type0):
        self.mode, self.rec, self.vis0 = mode, rec, vis0
        self.ir_slots = ir_slots            # 完整記錄中各 IR 區塊的起點 (位元)
        self.ir_bits = ir_bits              # IR 字長 (9 或 10)
        self.ir_doc_words = ir_doc_words    # IR 文件字數 (影像緊接其後)
        self.ir_valid = ir_valid            # 有資料的 IR 樣本數
        self.ir_bytes = ir_bytes            # 一個 IR 區塊的 bytes
        self.ir_types = ir_types            # .indx 的區塊類型：各 IR 區塊
        self.vis_type0 = vis_type0          # .indx 的區塊類型：VIS 第 0 段 (第 k 段 = vis_type0 + k)


LAYOUTS = {
    'A': Layout('A', 99244, 4444, (0,), 9, 128, 3818, 4444, (0,), 1),
    'AA': Layout('AA', 109536, 14736, (0, 58944), 10, 16, NIR, 7368, (0, 1), 2),
    'AAA': Layout('AAA', 109516, 14716, (0, 58864), 10, 16, NIR, 7358, (2, 3), 4),     # AAA 的區塊號 2~11
}
COMMON_DOC = 16 + NIR + 8       # AA / AAA：共同文件在 IR 區塊內的字位置 (位元 38460)

ONE, ZERO = 0xFE, 0x01
UPPER, LOWER = 0x14A, 0x2B5     # IR 區塊段碼第 1 字：上 (IR1) / 下 (IR2) 偵測器

# 濾鏡輪位置碼 -> VAS 波段 (預測標頭的低 4 位元)
FILTER_BAND = {0b0000: 1, 0b1011: 2, 0b0100: 3, 0b0010: 4, 0b1010: 5, 0b0110: 6, 0b0101: 7, 0b1000: 8,
               0b0001: 9, 0b0111: 10, 0b0011: 11, 0b1001: 12}
SMALL, LARGE, INSB = 0, 1, 3    # 預測標頭的偵測器類型 (小 HgCdTe、大 HgCdTe、InSb)
_SMALL_BANDS = (3, 4, 5, 7, 8, 9, 10)

# VAS 波段：(中心波長 µm、種類、功能代稱)
VAS_BANDS = {1: (14.73, 'ir', 'CO2 LW'), 2: (14.48, 'ir', 'CO2 LW'), 3: (14.25, 'ir', 'CO2 LW'),
             4: (14.01, 'ir', 'CO2 LW'), 5: (13.33, 'ir', 'CO2 LW'), 6: (4.525, 'swir', 'CO2 SW'),
             7: (12.66, 'ir', 'Dirty LW Window'), 8: (11.17, 'ir', 'LW Window'), 9: (7.261, 'wv', 'Lower-level Water Vapor'),
             10: (6.725, 'wv', 'Water Vapor'), 11: (4.444, 'swir', 'CO2 SW'), 12: (3.945, 'swir', 'Shortwave Window')}
# 輻射 (mW/(m2 sr cm-1)) -> 亮溫 的係數 (舊版 goes_pregvar_tools.py 的一組 VAS 係數，GOES-4~7 共用)
FK1 = (.37407E+4, .39157E+4, .40871E+4, .43417E+4, .50296E+4, .12819E+6, .58519E+4, .84911E+4, .30936E+5, .38873E+5,
       .13611E+6, .19511E+6)
FK2 = (.97802E+3, .99304E+3, .10073E+4, .10278E+4, .10795E+4, .31767E+4, .11354E+4, .12853E+4, .19778E+4, .21343E+4,
       .32408E+4, .36542E+4)
TC1 = (-.89414E-3, .24306E-2, .24917E-2, .34902E-2, .39097E-2, .66916E-1, .63888E-2, .34408, .70558E-1, .78113,
       .58431E-1, .43968)
TC2 = (.99998, .99995, .99995, .99993, .99993, .99983, .99991, .99722, .99973, .99717, .99986, .99903)
# Mode AA：各波段的放大位數 F 與 DELTA SCALE FACTOR (低 4 位元) -> DELTAF
F_AA = (8, 8, 8, 8, 8, 4, 8, 8, 7, 7, 4, 2)
DELTAF = (0, -1, -2, -3, -4, -5, 0, 0, 0, 1, 2, 3, 4, 5, 0, 0)


def brightness_temperature(rad, band):
    """輻射 -> 亮溫 (°C)：T = (FK2 / ln(FK1/R + 1) - TC1) / TC2 (McIDAS 的形式)；R <= 0 為 NaN。"""
    i = band - 1
    with np.errstate(divide='ignore', invalid='ignore'):
        teff = FK2[i] / np.log(FK1[i] / rad + 1.)
    return np.where(rad > 0, (teff - TC1[i]) / TC2[i] - 273.15, np.nan)


def mode_a_temperature(count):
    """Mode A 的 8 位元已校正計數值 -> 溫度 (°C) (NCC 1978 附錄 II：0~175 為 330 - c/2，176~255 為 418 - c K)。"""
    c = np.asarray(count, dtype=np.float32)
    return np.where(c >= 176, 418. - c, 330. - c / 2.) - 273.15


def detector_index(det, upper, band):
    """偵測器 / 濾鏡組合 k (1~38，AAA 表 5C1)：大 HgCdTe 與 InSb 上 = 波段、下 = 12 + 波段；
    小 HgCdTe 上 = 25~31、下 = 32~38 (波段 3、4、5、7、8、9、10)。"""
    if det == SMALL:
        if band not in _SMALL_BANDS:
            return 0
        return (25 if upper else 32) + _SMALL_BANDS.index(band)
    return band if upper else 12 + band


# ---- 位元工具 -----------------------------------------------------------------------------------
def bits_at(buf, bitpos, nbits):
    """buf：uint8 (n, L) (每列後面至少多 2 bytes)；bitpos：各值起點的位元位置 -> (n, len(bitpos)) 的 nbits 位元值 (nbits <= 17)。"""
    bitpos = np.asarray(bitpos, dtype=np.int64)
    byte = bitpos >> 3
    shift = (24 - nbits - (bitpos & 7)).astype(np.uint32)
    b = buf.astype(np.uint32)
    v = (b[:, byte] << 16) | (b[:, byte + 1] << 8) | b[:, byte + 2]
    return (v >> shift) & np.uint32((1 << nbits) - 1)


def gather(data, offsets, start, nbytes):
    """各位置 (offsets；-1 = 沒有這個區塊) 從 start byte 起取 nbytes (+2 bytes 補零) -> uint8 (n, nbytes + 2)；
    沒有的區塊與超出檔尾的部分為 0。"""
    offsets = np.asarray(offsets, dtype=np.int64)
    idx = offsets[:, None] + start + np.arange(nbytes + 2)
    bad = (offsets[:, None] < 0) | (idx >= len(data))
    out = data[np.clip(idx, 0, len(data) - 1)]
    out[bad] = 0
    return out


def words(data, offsets, bitpos, count, nbits):
    """各記錄從位元 bitpos 起連續 count 個 nbits 位元的字 -> (n, count)。"""
    byte0 = bitpos >> 3
    nbytes = (bitpos - byte0 * 8 + count * nbits + 7) // 8
    buf = gather(data, offsets, byte0, nbytes)
    return bits_at(buf, (bitpos - byte0 * 8) + nbits * np.arange(count), nbits)


def bcd(x):
    x = np.asarray(x, dtype=np.int64)
    return (x >> 4) * 10 + (x & 15)


def _sector_code(k):
    v = 0
    for j in range(3):
        v = (v << 6) | (0b111110 if (k >> (2 - j)) & 1 else 0b000001)
    return v


_CODES = [_sector_code(k) for k in range(8)]


# ---- 檔案 -----------------------------------------------------------------------------------------
class RawFile:
    """一個 Mode A / AA / AAA 檔：各記錄的區塊位置與文件欄位 (meta)。
      offsets  各記錄的起點 (byte)
      ir_pos   (記錄數, IR 區塊數) 各 IR 區塊的起點；vis_pos (記錄數, 8) 各 VIS 段的起點；-1 = 這筆記錄沒有該區塊
      ir_valid、vis_valid  同形狀，各區塊的有效樣本數 (NCEI 補零的截短區塊之後為補零，見 _valid_samples)
    同名的 .indx (NCEI 附的區塊索引：每列 檔名、位置 (32 位元字)、位元偏移、區塊類型 [掃描數]) 存在時以它為準，
    可讀到不完整的記錄 (缺段的短記錄、檔尾截斷的記錄)；每個區塊仍以段碼驗證，不符的區塊不用。
    沒有 .indx 或索引與資料不符時，改以 VIS 段碼自行定位 (只用完整的記錄)。
    probe = (前 n 筆, 後 m 筆)：只取檔頭、檔尾的記錄 (資料列表用)。每筆記錄的文件分散在整個檔 (每筆約 100 KB)，
    全部讀取在傳統硬碟上冷啟動每檔要 3~5 秒；波段在檔頭就會全部出現 (MSI 前 14 筆、DS 前 57 筆)。"""

    def __init__(self, path, mode, probe=None):
        self.path, self.mode = path, mode
        self.layout = LAYOUTS[mode]
        with open(path, 'rb') as f:
            self._mm = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
        self.data = np.frombuffer(self._mm, np.uint8)
        self.index_file = None
        blocks = self._read_indx()
        if blocks is None:
            blocks = self._sync_index(probe)
        elif probe:
            n = len(blocks[0])
            keep = np.unique(np.r_[np.arange(min(probe[0], n)), np.arange(max(0, n - probe[1]), n)])
            blocks = tuple(b[keep] for b in blocks)
        self.offsets, self.ir_pos, self.vis_pos, ir_line, vis_line = blocks
        self.ir_valid, self.vis_valid = self._valid_samples(ir_line, vis_line)
        self.meta = parse(self.data, self.ir_pos, mode)

    def _valid_samples(self, ir_line, vis_line):
        """各 IR 區塊、VIS 段的有效樣本數 (ir_valid (記錄數, IR 區塊數)、vis_valid (記錄數, 8))：
        .meta 的 truncated_line_bits {.indx 行號: 有效位元數} 記錄的補零區塊只有前面的樣本是資料，其餘為完整長度。
        AA / AAA 的 IR1 截在共同文件 (掃描數、時間、校正) 之前時，這筆記錄的 IR 都不用。"""
        L = self.layout
        ir_n = np.full(ir_line.shape, NIR, np.int64)
        vis_n = np.full(vis_line.shape, NVIS, np.int64)
        try:
            with open(os.path.splitext(self.path)[0] + '.meta', encoding='utf-8') as f:
                trunc = {int(k): int(v) for k, v in json.load(f).get('truncated_line_bits', {}).items()}
        except (OSError, ValueError, AttributeError, TypeError):
            trunc = {}
        if not trunc:
            return ir_n, vis_n
        lines = np.array(sorted(trunc), np.int64)
        bits = np.array([trunc[k] for k in lines], np.int64)
        for line, n, full, valid in ((ir_line, ir_n, NIR, lambda b: b // L.ir_bits - L.ir_doc_words),
                                     (vis_line, vis_n, NVIS, lambda b: (b - VIS_DOC * 8) // 6)):
            i = np.minimum(np.searchsorted(lines, line), len(lines) - 1)
            hit = (line > 0) & (lines[i] == line)
            n[hit] = np.clip(valid(bits[i][hit]), 0, full)
        if L.mode != 'A':
            # IR1 截在共同文件 (第 COMMON_DOC 字起、parse 用到的前 384 字) 結束之前：文件欄位不可靠
            i = np.minimum(np.searchsorted(lines, ir_line[:, 0]), len(lines) - 1)
            cut = (ir_line[:, 0] > 0) & (lines[i] == ir_line[:, 0]) & (bits[i] < (COMMON_DOC + 384) * L.ir_bits)
            ir_n[cut] = 0
        return ir_n, vis_n

    # ---- 區塊位置 ------------------------------------------------------------------------------
    def _read_indx(self):
        path = os.path.splitext(self.path)[0] + '.indx'
        if not os.path.exists(path):
            return None
        L = self.layout
        starts, irs, viss = [], [], []
        irl, visl = [], []              # 各區塊在 .indx 的行號 (由 1 起算；.meta 的 truncated_line_bits 用)
        last = None                     # 目前記錄中最後一個區塊的類型；None = 不在記錄中
        try:
            with open(path) as f:
                for lno, line in enumerate(f, 1):
                    s = line.split()
                    if len(s) < 4:
                        continue
                    pos, t = int(s[1]) * 4 + int(s[2]) // 8, int(s[3])
                    if t == L.ir_types[0]:
                        starts.append(pos)
                        irs.append([pos] + [-1] * (len(L.ir_types) - 1))
                        viss.append([-1] * 8)
                        irl.append([lno] + [0] * (len(L.ir_types) - 1))
                        visl.append([0] * 8)
                        last = t
                    elif last is None or t <= last:
                        # 一筆記錄內區塊類型遞增；回頭 = 沒有 IR1 開頭的損壞記錄 (例如 AAA-MSI 掃描 1120)，略過到下一個 IR1
                        last = None
                    else:
                        if t in L.ir_types:
                            irs[-1][L.ir_types.index(t)] = pos
                            irl[-1][L.ir_types.index(t)] = lno
                        elif L.vis_type0 <= t < L.vis_type0 + 8:
                            viss[-1][t - L.vis_type0] = pos
                            visl[-1][t - L.vis_type0] = lno
                        last = t
        except (OSError, ValueError):
            return None
        if not starts:
            return None
        start = np.array(starts, np.int64)
        ir = np.array(irs, np.int64).reshape(len(starts), -1)
        vis = np.array(viss, np.int64).reshape(len(starts), 8)
        # 驗證：區塊完整在檔內；不在正常配置位置的區塊 (例如短記錄中前移的段) 還要段碼相符
        # (在正常位置的區塊即使標頭段碼有位元錯誤也保留，影像本身沒有問題)。不符的區塊標成 -1，太多不符時不用這個索引
        # 只讀被移動的區塊 (逐筆讀每一段的段碼，冷啟動時在傳統硬碟上每檔要 6~8 秒)
        size = len(self.data)
        vis_bad = (vis < 0) | (vis + VIS_SECTOR > size)
        chk = np.nonzero((vis != start[:, None] + L.vis0 + VIS_SECTOR * np.arange(8)[None, :]) & ~vis_bad)
        if chk[0].size:
            b = gather(self.data, vis[chk], 0, 1).astype(np.int64)
            vis_bad[chk] = ((b[:, 0] << 10) | (b[:, 1] << 2) | (b[:, 2] >> 6)) != np.array(_CODES)[chk[1]]
        ir_bad = (ir < 0) | (ir + L.ir_bytes > size)
        if L.mode != 'A':
            chk = np.nonzero((ir != start[:, None] + np.array(L.ir_slots, np.int64)[None, :] // 8) & ~ir_bad)
            if chk[0].size:
                ir_bad[chk] = ~np.isin(words(self.data, ir[chk], 0, 1, 10)[:, 0], (UPPER, LOWER))
        checked = (vis >= 0).sum() + (ir >= 0).sum()
        if checked == 0 or ((vis_bad & (vis >= 0)).sum() + (ir_bad & (ir >= 0)).sum()) > .1 * checked:
            return None
        vis[vis_bad] = -1
        ir[ir_bad] = -1
        keep = ir[:, 0] >= 0
        self.index_file = path
        irl = np.array(irl, np.int64).reshape(len(starts), -1)
        visl = np.array(visl, np.int64).reshape(len(starts), 8)
        return start[keep], ir[keep], vis[keep], irl[keep], visl[keep]

    def _sync_index(self, probe=None):
        L = self.layout
        if probe:
            start = np.unique(np.r_[self._index(probe[0]), self._last_records(probe[1])]).astype(np.int64)
        else:
            start = self._index()
        ir = start[:, None] + np.array(L.ir_slots, np.int64)[None, :] // 8
        vis = start[:, None] + L.vis0 + VIS_SECTOR * np.arange(8)[None, :]
        return start, ir, vis, np.zeros(ir.shape, np.int64), np.zeros(vis.shape, np.int64)     # 沒有 .indx 行號

    def _sector_ok(self, p, k):
        q = p + self.layout.vis0 + k * VIS_SECTOR
        if q + 3 > len(self.data):
            return False
        b0, b1, b2 = (int(x) for x in self.data[q:q + 3])
        return ((b0 << 10) | (b1 << 2) | (b2 >> 6)) == _CODES[k]

    def _ok(self, p, sectors=(0, 1, 6, 7)):
        return all(self._sector_ok(p, k) for k in sectors)

    def _mostly_ok(self, p):
        return sum(self._sector_ok(p, k) for k in range(8)) >= 6

    def _index(self, limit=None):
        """依 VIS 第 0、1 段的段碼確認每筆記錄；不符時，8 段中至少 6 段的段碼相符也算 (段碼的位元錯誤，
        例如 SMS-2 1979/243 12:30 掃描 1267、1670：SSEC AREA 有這兩圈，內容與此處讀到的相同)，
        否則往後找下一個第 0、1、6、7 段都相符的位置重新對齊 (短記錄、雜訊)。
        記錄末端 (第 6、7 段) 只在下一筆不在正常位置 (及最後一筆) 時確認：下一筆接在正常位置就表示這筆完整。
        (每筆只讀記錄開頭，冷啟動時在傳統硬碟上快一倍。) limit：只找前 limit 筆。"""
        size, rec, vis0 = len(self.data), self.layout.rec, self.layout.vis0
        offs, p = [], 0
        while p + rec <= size and (limit is None or len(offs) <= limit):
            if self._ok(p, (0, 1)) or (offs and p == offs[-1] + rec and self._mostly_ok(p)):
                offs.append(p)
                p += rec
                continue
            # 上一筆只確認過開頭，可能是短記錄 (下一筆在 p 之前就開始)：從上一筆之後找
            base = offs[-1] if offs and p == offs[-1] + rec else p
            q = self._mm.find(b'\x04\x10', base + vis0 + 1)
            while q != -1 and not (q - vis0 > base and q - vis0 + rec <= size and self._ok(q - vis0)):
                q = self._mm.find(b'\x04\x10', q + 1)
            if q == -1:
                break
            p = q - vis0
        out = [p for i, p in enumerate(offs)
               if (i + 1 < len(offs) and offs[i + 1] == p + rec) or self._ok(p, (6, 7))]
        return np.array(out[:limit], dtype=np.int64)

    def _last_records(self, n):
        """檔尾的最後 n 筆完整記錄 (由檔尾往前找段碼)。"""
        size, rec, vis0 = len(self.data), self.layout.rec, self.layout.vis0
        q = self._mm.rfind(b'\x04\x10')
        while q != -1 and not (q - vis0 >= 0 and q - vis0 + rec <= size and self._ok(q - vis0)):
            q = self._mm.rfind(b'\x04\x10', 0, q)
        if q == -1:
            return np.zeros(0, np.int64)
        out = [q - vis0]
        while len(out) < n and out[-1] - rec >= 0 and self._ok(out[-1] - rec):
            out.append(out[-1] - rec)
        return np.array(out[::-1], dtype=np.int64)

    # ---- 影像 ------------------------------------------------------------------------------
    def vis_lines(self, recs, sectors, c0, c1):
        """各 (記錄序號, 段) 的 VIS 計數值 (0~63)，行 c0~c1-1 -> uint8 (n, c1-c0)。"""
        a = c0 - c0 % 4
        pos = self.vis_pos[np.asarray(recs), np.asarray(sectors)]
        buf = gather(self.data, pos, VIS_DOC + a * 3 // 4, ((c1 - a) * 6 + 7) // 8)
        return bits_at(buf, 6 * np.arange(c1 - a), 6)[:, c0 - a:].astype(np.uint8)

    def ir_video(self, recs, slot=0):
        """各記錄 IR 區塊 (slot) 的影像字 -> (n, 3822) (Mode A 為 8 位元計數值，3819~3822 為 0)。"""
        L = self.layout
        bitpos = L.ir_doc_words * L.ir_bits                     # 文件之後 (Mode A 1152、AA / AAA 160 位元，皆在 byte 邊界)
        pos = self.ir_pos[np.asarray(recs), slot]
        out = np.zeros((len(pos), NIR), np.uint16)
        if L.mode == 'A':                       # 9 位元字 = 1 個旗標位元 + 8 位元計數值
            buf = gather(self.data, pos, bitpos >> 3, (L.ir_valid * 9 + 7) // 8 + 1)
            out[:, :L.ir_valid] = bits_at(buf, 1 + 9 * np.arange(L.ir_valid), 8)
        else:
            out[:] = words(self.data, pos, bitpos, NIR, 10)
        return out


# ---- 文件欄位 -------------------------------------------------------------------------------------
def _time(year, day, hour, minute, second, msec):
    """各記錄的時間 -> (YYDDD 整數陣列, 當日秒數陣列)。"""
    yyddd = (np.asarray(year) % 100) * 1000 + np.asarray(day)
    sec = np.asarray(hour) * 3600. + np.asarray(minute) * 60. + np.asarray(second) + np.asarray(msec) / 1000.
    return yyddd.astype(np.int64), sec


def _oa_words(doc):
    """文件第 99~124 字 (0 起算 98~123)：O&A 區塊號、minor frame index、6 個 32 位元整數 (有號)。"""
    blk, mfi = doc[:, 98].astype(np.int64), doc[:, 99].astype(np.int64)
    b = doc[:, 100:124].astype(np.int64).reshape(-1, 6, 4)
    v = (b[:, :, 0] << 24) | (b[:, :, 1] << 16) | (b[:, :, 2] << 8) | b[:, :, 3]
    v = np.where(v >= 1 << 31, v - (1 << 32), v)
    return blk, mfi, v


def repair_scan(scan, limit=20):
    """掃描數的位元錯誤：前後兩筆是連續掃描 (相差 2)、本筆卻偏離超過 limit 時，改為前一筆 + 1。
    (SMS-2 1979/243 12:30 掃描 670、940 的掃描數錯成 2673；SSEC AREA 有這兩圈，內容與這兩筆記錄相同。)
    只修明顯的錯誤，DS 駐留 (同一掃描數重複) 與回掃不受影響。"""
    s = np.array(scan, dtype=np.int64)
    if len(s) >= 3:
        prev, cur, nxt = s[:-2], s[1:-1].copy(), s[2:]
        fix = (nxt - prev == 2) & (np.abs(cur - prev) > limit) & (prev >= 1) & (nxt <= 1821)
        s[1:-1] = np.where(fix, prev + 1, cur)
    return s


def parse(data, ir_pos, mode):
    """各記錄的文件欄位 -> dict (陣列，長度 = 記錄數)。ir_pos：各記錄 IR 區塊的起點 (記錄數, 區塊數)。"""
    n = len(ir_pos)
    m = {}
    if mode == 'A':
        doc = bits_at(gather(data, ir_pos[:, 0], 0, 145), 1 + 9 * np.arange(128), 8).astype(np.int64)
        m['frame'] = doc[:, 3] == ONE
        m['ir_select'] = doc[:, 7]
        m['scan'] = repair_scan(((doc[:, 24] >> 4) & 3) * 1000 + (doc[:, 24] & 15) * 100 + bcd(doc[:, 25]))
        m['yyddd'], m['seconds'] = _time(bcd(doc[:, 26]) * 100 + bcd(doc[:, 27]), bcd(doc[:, 28]) * 100 + bcd(doc[:, 29]),
                                         bcd(doc[:, 30]), bcd(doc[:, 31]), bcd(doc[:, 32]), bcd(doc[:, 33]) * 10)
        m['blk'], m['mfi'], m['oa'] = _oa_words(doc)
        m['doc'] = doc
        m['submode'] = np.zeros(n, np.int64)
        return m
    L = LAYOUTS[mode]
    irdoc = [words(data, ir_pos[:, s], 0, 16, 10).astype(np.int64) for s in range(len(L.ir_slots))]
    common = (words(data, ir_pos[:, 0], COMMON_DOC * 10, 384, 10) & 0xFF).astype(np.int64)
    f1, f2, f3 = common[:, :128], common[:, 128:256], common[:, 256:384]
    m['scan'] = repair_scan(((f1[:, 24] >> 4) & (7 if mode == 'AA' else 3)) * 1000 + (f1[:, 24] & 15) * 100
                            + bcd(f1[:, 25]))
    m['yyddd'], m['seconds'] = _time(bcd(f1[:, 27]), (f1[:, 28] & 15) * 100 + bcd(f1[:, 29]), bcd(f1[:, 30]),
                                     bcd(f1[:, 31]), bcd(f1[:, 32]),
                                     (f1[:, 33] >> 4) * 100 + (f1[:, 33] & 15) * 10 + (f1[:, 28] >> 4))
    m['blk'], m['mfi'], m['oa'] = _oa_words(f1)
    m['submode'] = f2[:, 22]
    m['vasmode'] = f2[:, 23]            # 第二欄第 24 字 Mode：bit 9 (& 2) 1 = VAS、0 = VISSR 模式；bit 10 (& 1) 1 = MSI、0 = DS
    hdr_word, frame_word = (4, 2) if mode == 'AA' else (6, 3)
    m['frame'] = (irdoc[0][:, frame_word] == ONE) & (m['submode'] != 4)
    m['slot_frame'] = np.stack([d[:, frame_word] == ONE for d in irdoc], 1)
    nslot = len(L.ir_slots)
    m['sector'] = np.stack([d[:, 0] for d in irdoc], 1)            # $14A 上 / $2B5 下
    hdr = np.stack([d[:, hdr_word] & 0xFF for d in irdoc], 1)
    m['band'] = np.vectorize(lambda c: FILTER_BAND.get(int(c), 0))(hdr & 15) if n else np.zeros((0, nslot), np.int64)
    m['det'] = (hdr >> 6) & 3
    # 預測標頭 = $80 (文件所寫 VISSR 模式的標頭，偵測器類型 10 未定義) 時波段無法由濾鏡碼判斷，不使用。
    # (實際的 VISSR 模式樣本 GOES-6 1987/183 .AAA 的標頭是一般的 band 8 小偵測器 $38，見 GOES-STRETCHED 的 _Info)
    m['present'] = np.isin(m['sector'], (UPPER, LOWER)) & (m['band'] > 0) & (m['det'] != 2)
    if mode == 'AA':
        # 共同文件第 3 段 (校正)：第 60~63 字 Y_Z + Q (IR1、IR2)、第 64~65 字 DELTA SCALE FACTOR
        m['yzq'] = np.stack([(f3[:, 59] << 8) | f3[:, 60], (f3[:, 61] << 8) | f3[:, 62]], 1)
        m['dsf'] = np.stack([f3[:, 63], f3[:, 64]], 1)
        m['step'] = np.where(m['det'] == SMALL, m['scan'][:, None] + np.array([0, 1]),
                             m['scan'][:, None] + np.array([0, 4]))
    else:
        m['earth'] = bcd(f1[:, 38]) * 100 + bcd(f1[:, 39])
        adj = np.stack([bcd(d[:, 8]) * 100 + bcd(d[:, 9]) for d in irdoc], 1)
        # 掃描數 - earth count 在一個畫格內是常數 (= EQUAT - 836)：取整檔的中位數，不受個別記錄位元錯誤影響
        ok = (m['earth'] > 0) & (m['scan'] > 0)
        k = int(np.median((m['scan'] - m['earth'])[ok])) if ok.any() else 0
        aa_rule = m['scan'][:, None] + np.where(m['sector'] == LOWER, np.where(m['det'] == SMALL, 1, 4), 0)
        # adjusted earth count 為 0 (未填) 時改用 AA 的上下偵測器規則
        m['step'] = np.where(adj > 0, adj + k, aa_rule)
        m['ircal'] = f3[:, 62:112]                                   # 第 3 欄第 63~112 字：IR 校正參數 (20 圈輪替)
    m['nsteps'] = np.where(m['det'] == SMALL, 1, 2)
    m['common'] = common
    m['irdoc'] = irdoc
    return m


# ---- 軌道姿態 (O&A) -----------------------------------------------------------------------------
def assemble(blk, mfi, values, nwords=6):
    """各記錄的 (區塊號 1~10, MFI 1~2, 內容) -> {(區塊, MFI): 出現最多次的內容} (20 圈一組，重複出現時取眾數以排除雜訊)。"""
    count = {}
    for b, f, v in zip(blk, mfi, values):
        if 1 <= b <= 10 and f in (1, 2):
            key = (int(b), int(f))
            t = tuple(int(x) for x in v[:nwords])
            count.setdefault(key, {}).setdefault(t, 0)
            count[key][t] += 1
    return {k: max(c, key=c.get) for k, c in count.items()}


def oa_fields(blocks):
    """O&A 區塊 -> {名稱: 整數} (Mode AAA 表 3A；Mode A / AA 相同，block 8、9 的 Kepler 元素為 VIRGS 產生)。"""
    names = {1: ('DATE1', 'TIME1', None, None, None, None, 'X1', 'Y1', 'Z1', 'VX1', 'VY1', 'VZ1'),
             2: ('SPER', 'SPRA1', 'SPDC1', 'ZETA', 'RHO', 'ETA', 'GAMMA', 'NAMES', 'ID', 'SRA1', 'SDC1', 'GRA1'),
             3: ('EST', 'EET', 'FPER', 'TC', 'CE0', 'CE1', 'CE2', 'CE3', 'CR0', 'CR1', 'CR2', 'CR3'),
             7: (None,) * 11 + ('EPY',),
             8: ('EPH', 'SMA', 'ECC', 'INC', 'MA', 'AP', 'RAN', 'SBSCAN', 'SBSAMP', 'SBLAT', 'SBLONG', 'YAW'),
             9: ('SBSCANB', 'SBSAMPB', 'SBLATB', 'SBLONGB', 'YAWB', None, 'X2', 'Y2', 'Z2', 'VX2', 'VY2', 'VZ2'),
             10: (None, 'SPRA2', 'SPDC2', 'DELTAL', 'DELTAS', None, None, None, None, 'SRA2', 'SDC2', 'GRA2')}
    out = {}
    for (b, f), v in blocks.items():
        for i, x in enumerate(v):
            name = names.get(b, (None,) * 12)[(f - 1) * 6 + i]
            if name:
                out[name] = x
    return out


def nav_params(oa, yyddd, seconds, scan):
    """O&A -> McIDAS 'GOES' 導航參數 (mcidas_goes_nav.params() 的 dict)。yyddd、seconds、scan：參考記錄的日期、時間、掃描數。
    軌道用 VIRGS 的 Kepler 元素 (EPY、EPH、SMA、ECC、INC、MA、AP、RAN，同 McIDAS AREA 導航區塊)；
    自旋軸 = SPRA、SPDC 依 13 小時區間線性內插的反向 (赤經 + 180°、赤緯取負；沒有第二組時直接用 SPRA1、SPDC1)；
    姿態 pitch = -ZETA (line bias)、yaw = ETA (skew)、gamma = RHO (element bias，換成 VIS 像素數)。
    (同日 AREA 導航區塊的 pitch / yaw 與 ZETA / ETA 數值相同；GAMMA = sun pulse 角 = AREA 的 ISEANG。)
    沒有 Kepler 元素或自旋資料時 ValueError。"""
    from satima.core.i18n import tr
    from mcidas_goes_nav import params
    need = ('DATE1', 'TIME1', 'SPER', 'SPRA1', 'SPDC1', 'EPY', 'EPH', 'SMA', 'ECC', 'INC', 'MA', 'AP', 'RAN')
    miss = [k for k in need if k not in oa]
    if miss or oa['SMA'] <= 0 or oa['SPER'] <= 0:
        raise ValueError(tr('O&A 不完整或沒有軌道參數 (%s)', 'orbit & attitude data incomplete or without orbit elements (%s)')
                         % (', '.join(miss) or 'SMA = %d' % oa.get('SMA', 0)))
    s21 = 2. ** 21
    yyddd, scan = int(yyddd), int(scan)
    hours = seconds / 3600.
    spra1, spdc1 = oa['SPRA1'] / s21, oa['SPDC1'] / s21
    if oa['DATE1'] > 0 and (oa.get('SPRA2') or oa.get('SPDC2')):
        # 參考時刻相對 TIME1 (DATE1) 的歸一化時間 u (13 小時區間 -1~1)
        days = _daynum(yyddd) - _daynum(int(oa['DATE1']))
        u = 2. * ((hours + 24. * days) - oa['TIME1'] / 360000.) / 13. - 1.
        spra2, spdc2 = oa['SPRA2'] / s21, oa['SPDC2'] / s21
        if abs(spra2 - spra1) > 180.:
            spra2 += 360. if spra2 < spra1 else -360.
        spra = .5 * (spra2 + spra1 + (spra2 - spra1) * u)
        spdc = .5 * (spdc2 + spdc1 + (spdc2 - spdc1) * u)
    else:
        # 只有一組自旋軸 (GOES-1 1979/130 等 NCEI 移植的 SYSNAV 導航：區塊 1、9、10 全為 0)：不內插
        spra, spdc = spra1, spdc1
    a = np.zeros(128, dtype=np.int64)
    eph = int(oa['EPH'])
    a[1], a[3], a[4] = yyddd, 1, int(oa['EPY'])
    a[5] = (eph // 100) * 100 + int(round((eph % 100) / .6))      # McIDAS ETIMH：HHMM + 分鐘的百分位
    a[6:12] = [oa['SMA'], oa['ECC'], oa['INC'], oa['MA'], oa['AP'], oa['RAN']]
    a[12], a[13], a[14], a[15] = 895800, 1000000, 7285, oa['SPER']
    a[16], a[17], a[18], a[19] = 200000, 801821, 182230, 15288    # 20°、8 感測器 x 1821 掃描、18.375°、15288 像素
    p = params(a)
    p['declin'], p['rascen'] = -spdc, (spra + 180.) % 360.
    p['pitch'], p['yaw'], p['roll'] = -oa.get('ZETA', 0) / s21, oa.get('ETA', 0) / s21, 0.
    p['spinra'] = oa['SPER'] / 1000.                                # 毫秒
    p['tmpscl'] = p['spinra'] / 3600000.
    p['pictim'] = hours - scan * p['tmpscl']                        # 掃描 k 的時刻 = pictim + k 圈
    p['gamma'] = (oa.get('RHO', 0) / s21) / (p['degele'] / (p['eletot'] - 1.))
    p['beta'] = (1., p['pictim'], 1821., p['pictim'] + 1821 * p['tmpscl'])
    p['oa'] = dict(oa)
    return p


def _daynum(yyddd):
    y, d = yyddd // 1000 % 100, yyddd % 1000
    year = 1900 + y if y >= 50 else 2000 + y
    return (year - 1970) * 365 + (year - 1969) // 4 + d


# ---- Mode AA 可見光 --------------------------------------------------------------------------------
def vis_inverse(m):
    """Mode AA 的可見光是地面處理器 (VIP) 以可見光非線性係數換算後的值 Y = 63 Σ c_j (d/63)^j (j = 0~3)；
    係數在共同文件第 3 段第 77~110 字 (校正參數，20 圈輪替) 的資訊 byte 585~648 (8 個偵測器 x 4 個 16 位元有號數，
    c_j = 值 / 2^15)。回傳反算表 (64,)：Y -> 原始計數值 d (0~63，同 Mode A)；係數為恆等、不完整或不合理時 None。
    用 8 個偵測器係數的平均：偵測器與可見光線的對應無法由資料確認 (GOES-5 1982 與同時次 Mode A 比對，
    8 種排法的差異都在 0.15 個計數值內)。"""
    sets = assemble(m['blk'], m['mfi'], m['common'][:, 256 + 76:256 + 110], nwords=34)
    need = [(9, 2), (10, 1), (10, 2)]           # 資訊 byte 579~680 (第 18~20 圈)
    if not all(k in sets for k in need):
        return None
    info = np.concatenate([np.array(sets[k], np.int64) for k in need])
    b = info[584 - 578:648 - 578].reshape(8, 4, 2)
    w = (b[..., 0] << 8) | b[..., 1]
    c = (np.where(w >= 32768, w - 65536, w) / 32768.).mean(0)
    d = np.linspace(0., 63., 6301)
    y = 63. * sum(c[j] * (d / 63.) ** j for j in range(4))
    if np.abs(y - d).max() < .1 or not np.all(np.diff(y) > 0) or y[-1] < 8.:
        return None
    return np.interp(np.arange(64.), y, d)


# ---- Mode AAA IR 校正 ---------------------------------------------------------------------------
def aaa_calibration(blk, mfi, ircal):
    """第 3 欄第 63~112 字 (每圈 25 個 16 位元參數，20 圈輪替) -> (AB (38, 2), FAB (38,))；缺資料時 None。
    資訊 byte 編號 = 2 x (50 x (區塊 - 1) + (MFI - 1) x 25 + i) + 1 (表 5C)；AB 在 753~904、FAB 在 905~980。"""
    sets = assemble(blk, mfi, ircal, nwords=50)
    need = [(8, 2), (9, 1), (9, 2), (10, 1), (10, 2)]
    if not all(k in sets for k in need):
        return None
    info = np.zeros(1000, np.uint8)
    for (b, f), v in sets.items():
        start = 2 * (50 * (b - 1) + (f - 1) * 25)
        info[start:start + 50] = np.asarray(v, np.uint8)
    i16 = info.view('>i2').astype(np.int64)                     # 第 i 個 16 位元參數 = byte 2i+1、2i+2
    ab = i16[376:452].reshape(38, 2)                               # byte 753~904
    fab = i16[452:490]                                             # byte 905~980
    if not ab.any():
        return None
    return ab, fab
