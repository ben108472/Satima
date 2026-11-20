"""GMS-1 VISSR 磁帶映像 (FGGE 時期，1978~1979) — 例如 VS2303、VR0123、V10585

檔案：一卷磁帶 = 一個檔案 = 一個時次，可見光與紅外都在裡面；檔名為磁帶的卷號 (V 開頭共 6 個字，沒有副檔名)。
  是不是磁帶映像由檔案開頭的 VOL1 標籤判斷；衛星名稱與時間取自檔頭 (Mode Block 的觀測開始時刻)。
容器 (同隨磁帶提供的 vissr.c / gms_area.c)：每筆紀錄前有 4 bytes 大端序的長度，長度 0 = tape mark。
    VOL1、HDR1、HDR2 (各 80 bytes 的 EBCDIC 標籤)、tape mark、
    參數區塊 x 3 (各 12096 bytes)、影像紀錄 (紅外 7884、可見光 10152 bytes)、結束紀錄 (3024 bytes)、tape mark、
    EOF1、EOF2、tape mark x 2
格式：內容就是 JMA VISSR 格式 (GMS-VISSR 讀的 VS/IR*.Z??；JMA「FORMAT OF VISSR ARCHIVE DATA」)，只是每個 4-byte word
  存成一個 36-bit word，整筆紀錄是連續的位元串 (9 bytes = 2 個 word)：
    I*4   36 位元二補數整數
    R*4   高 27 位元為二補數尾數、低 9 位元為指數：值 = 尾數 x 2^(指數 - 255 - 27)
    R*8   2 個 word：尾數 63 位元 (第 1 個 word 加上第 2 個 word 的高 27 位元)、低 9 位元為指數：尾數 x 2^(指數 - 255 - 63)
    文字  每個 word 4 個字元，每個字元 9 位元 (EBCDIC)
  參數區塊 = 4 個參數段，每段 672 個 word，欄位與 JMA 格式 2688 bytes 的參數段一一對應 (定義在 sensors/vissr_common.py)：
    第 1 塊  Mode、S/DB、IR 校正、VIS 校正
    第 2 塊  座標轉換、姿態預報、軌道預報 1、軌道預報 2
    第 3 塊  位置對應 JMA 文件的 ⑧ (Mode Block 狀態表的 Data segment 8)、保留 x 2、⑨ β 角取樣；文件對 ⑧、⑨ 的內容
             都寫「Omitted」，沒有欄位定義，所以以原始 word 保留 (h['Data_segment_8'])。磁帶上只有 ⑧ 有內容：
             開頭 12 個 word 是各子表的起始位置，子表包含與 Mode Block 相同的運用模式 / 掃描行號上下限 / 赤道行號、
             感測器故障表、各種紀錄 (影像、黑體、註記…) 的行號範圍、分區的行號範圍與紅外 / 可見光的畫面常數，
             看起來是 S/DB 的畫面設定 (文件的紅外檔圖示把 ⑧ 寫成 DCD communication、可見光檔圖示寫成 S/DB operation)；
             定標與導航都用不到。⑨ β 角取樣全為 0 (各線的 β 角在 LCW 裡)。
  影像紀錄 = 一條掃描線：Line Control Word (16 個 word) + DOC + 計數值，全部以該波段的樣本位元數連續存放
  (位元數、LCW / DOC 的樣本數、像素數取自 Mode Block)：
    紅外    9 位元：LCW 64 + DOC 254 + 6688 個像素 (計數值 0~255)
    可見光  6 位元：LCW 96 + DOC 61 + 13376 個像素 (計數值 0~63)
  每一圈 (自旋) 依序是 1 筆紅外 (行號 N) 與同一圈的可見光 (行號 4N-3 ~ 4N)。
    Data ID 高 18 位元 = 影像區分、低 18 位元 = 資料區分 (1 = 紅外，2、4、8、16 = 可見光感測器 1~4)
    Line name：1 = 影像、4 = 黑體、8 = 測試、16 = 註記、32 = 灰階；只取影像 (1)，行號沿用 LCW 的行號。
  已知的磁帶 (1979/01、1979/10) 都用備援感測器組 (Mode Block 的 Sensor selection = 2)，其中可見光第 2 感測器故障
  (VIS redundant 2ch = 0)，磁帶上沒有它的紀錄：可見光每 4 條掃描線缺 1 條 (行號 4N-2)。
校正：同 GMS-VISSR — 紅外為檔頭的計數值->亮溫表 (256 項)；可見光為各感測器自己的計數值->反照率表 (64 項)，
  每條線依 LCW 的 Data ID 選用。
  [visible] fill_missing_lines = true (預設) 時，單獨缺的可見光掃描線 (上下兩條都有資料) 以上下兩條反照率的平均補上
  (McIDAS 的 GMS AREA 檔也是以上下兩條的平均補這條線)；false 則保留為無資料。
導航：同 GMS-VISSR — JMA MGIVSR (sensors/vissr_nav.py) 與地球邊緣修正，參數取自第 2 個參數區塊。
完整檔頭：reader.header(scene, 'IR') 回傳 {紀錄名稱: numpy 結構}，名稱與欄位同 GMS-VISSR (例如
    h['Coordinate_Conversion']['Longitude of SSP']、h['IR_BT_Table']、h['Line_Control_Word']['β angle'])，
  數值保留磁帶的精確度：整數為 int64、實數為 float64、文字轉成 ASCII；另有 h['Labels'] (磁帶標籤)。
"""
import os

import numpy as np

from satima.core.i18n import tr
from satima.geo.reader import Scene
from vissr_common import LINE_CONTROL_WORD, MODE, SEGMENTS, VissrBase, observation_time, parse_segment, satellite_name

_VOL1 = (80).to_bytes(4, 'big') + 'VOL1'.encode('cp037')         # 檔案開頭：80 bytes 的 VOL1 標籤
_LABEL = 80
_SEGMENT_WORDS = 672                                             # 一個參數段的 word 數
_SEGMENT_BYTES = _SEGMENT_WORDS * 9 // 2
_BLOCK = 4 * _SEGMENT_BYTES                                      # 參數區塊 (4 個參數段) = 12096 bytes
_LCW_BYTES = LINE_CONTROL_WORD.itemsize // 4 * 9 // 2            # LCW 16 個 word = 72 bytes
# 第 3 個參數區塊的 4 段 (JMA 文件的 ⑧、保留 x 2、⑨；文件沒有定義內容，全為 0 的不保留)
_SEGMENTS = SEGMENTS['GMS'] + ('*Data_segment_8', '*Reserved_segment_10', '*Reserved_segment_11', '*Beta_angle_sampling')
_EBCDIC = np.frombuffer(bytes(range(256)).decode('cp037').encode('latin-1', 'replace'), np.uint8)


# ==================================================================================================
#  36-bit word
# ==================================================================================================
def words36(buf):
    """bytes -> 36-bit word (uint64)；每 9 bytes 為 2 個 word，不足 9 bytes 的尾端捨去。"""
    b = np.frombuffer(buf, np.uint8)
    b = b[:len(b) // 9 * 9].reshape(-1, 9).astype(np.uint64)
    w = np.empty((len(b), 2), np.uint64)
    w[:, 0] = (b[:, 0] << 28) | (b[:, 1] << 20) | (b[:, 2] << 12) | (b[:, 3] << 4) | (b[:, 4] >> 4)
    w[:, 1] = ((b[:, 4] & 15) << 32) | (b[:, 5] << 24) | (b[:, 6] << 16) | (b[:, 7] << 8) | b[:, 8]
    return w.ravel()


def _int(w):
    """36 位元二補數整數 -> int64。"""
    return (w.astype(np.int64) ^ (1 << 35)) - (1 << 35)


def _real(w):
    """單精確度實數：尾數 27 位元 (二補數)、指數 9 位元。"""
    s = _int(w)
    return np.ldexp((s >> 9).astype(np.float64), ((s & 0x1FF) - 282).astype(np.int32))


def _double(hi, lo):
    """倍精確度實數 (2 個 word)：尾數 63 位元 (二補數)、指數 9 位元。"""
    lo = lo.astype(np.int64)
    return np.ldexp(_int(hi).astype(np.float64) * float(1 << 27) + (lo >> 9), ((lo & 0x1FF) - 318).astype(np.int32))


def _text(w):
    """(n, k) 個 word -> n 個長度 4k 的字串 (每字元 9 位元的 EBCDIC -> ASCII)。"""
    ch = np.stack([(w >> s) & 0xFF for s in (27, 18, 9, 0)], axis=-1).astype(np.intp)
    return np.ascontiguousarray(_EBCDIC[ch.reshape(len(w), -1)]).view('S%d' % (4 * w.shape[1]))[:, 0]


def _field(w, ftype):
    kind, size = ftype.kind, ftype.itemsize
    if kind == 'i' and size == 4:
        return _int(w[:, 0])
    if kind == 'f' and size == 4:
        return _real(w[:, 0])
    if kind == 'f' and size == 8:
        return _double(w[:, 0], w[:, 1])
    if kind == 'S' and size % 4 == 0:
        return _text(w)
    raise TypeError(tr('36-bit word 無法對應的欄位型別 %s', 'field type %s cannot be mapped to 36-bit words') % ftype)


def decode(words, dtype, count=1):
    """36-bit word -> count 筆 numpy 資料。dtype 是 JMA 格式 (bytes) 的定義，每 4 bytes 對應 1 個 word；
    回傳的欄位名稱相同、型別加寬：I*4 -> int64、R*4 / R*8 -> float64、文字 -> ASCII。"""
    dtype = np.dtype(dtype)
    n = dtype.itemsize // 4
    if len(words) < count * n:
        raise ValueError(tr('資料長度不足 (%d 個 word，需要 %d 個)', 'not enough data (%d words, %d needed)')
                         % (len(words), count * n))
    w = np.asarray(words[:count * n], np.uint64).reshape(count, n)
    if dtype.names is None:
        return _field(w, dtype)
    wide = {'i': 'i8', 'f': 'f8'}
    out = np.zeros(count, [(name, wide.get(dtype[name].kind, dtype[name].str)) for name in dtype.names])
    for name in dtype.names:
        ftype, offset = dtype.fields[name][:2]
        out[name] = _field(w[:, offset // 4:(offset + ftype.itemsize) // 4], ftype)
    return out


class Words36:
    """參數段的來源：36-bit word (vissr_common.parse_segment 用)。offset 為 JMA 格式的 bytes 位移 (4 bytes = 1 個 word)。"""

    def __init__(self, words):
        self.words = words

    def record(self, dtype, offset=0):
        return decode(self.words[offset // 4:], dtype)[0]

    def array(self, dtype, count, offset=0):
        return decode(self.words[offset // 4:], dtype, count)

    def raw(self):
        return self.words.copy()


# ---- 影像樣本 --------------------------------------------------------------------
def _unpack6(raw):
    """(n, 3k) bytes -> (n, 4k) 個 6 位元樣本。"""
    b = raw.reshape(raw.shape[0], -1, 3)
    out = np.empty(b.shape[:2] + (4,), np.uint8)
    out[..., 0] = b[..., 0] >> 2
    out[..., 1] = ((b[..., 0] & 3) << 4) | (b[..., 1] >> 4)
    out[..., 2] = ((b[..., 1] & 15) << 2) | (b[..., 2] >> 6)
    out[..., 3] = b[..., 2] & 63
    return out.reshape(raw.shape[0], -1)


def _unpack9(raw):
    """(n, 9k) bytes -> (n, 8k) 個 9 位元樣本。"""
    b = raw.reshape(raw.shape[0], -1, 9).astype(np.uint16)
    out = np.empty(b.shape[:2] + (8,), np.uint16)
    for k in range(8):
        out[..., k] = ((b[..., k] << (k + 1)) | (b[..., k + 1] >> (7 - k))) & 0x1FF
    return out.reshape(raw.shape[0], -1)


# 樣本位元數 -> (每組樣本數, 每組 bytes, 解開的函式)：一組 = 樣本與 byte 邊界對齊的最小單位
_SAMPLES = {6: (4, 3, _unpack6), 9: (8, 9, _unpack9)}


# ==================================================================================================
#  磁帶映像
# ==================================================================================================
def _length(f):
    """下一筆紀錄的長度 (0 = tape mark)；檔尾為 None。"""
    head = f.read(4)
    return int.from_bytes(head, 'big') if len(head) == 4 else None


def _labels(f):
    """讀到第一個 tape mark 為止的標籤 (VOL1、HDR1、HDR2…)。"""
    labels = []
    while True:
        n = _length(f)
        if n == 0:
            return labels
        if n != _LABEL:
            raise ValueError(tr('磁帶標籤不正確 (紀錄長度 %s)', 'invalid tape label (record length %s)') % n)
        labels.append(f.read(_LABEL).decode('cp037'))


def read_mode(path):
    """磁帶映像 -> Mode Block；不是 VISSR 磁帶映像時 None。只讀檔案開頭 (資料列表用)。"""
    with open(path, 'rb') as f:
        if f.read(len(_VOL1)) != _VOL1:
            return None
        f.seek(0)
        _labels(f)
        n = _length(f)
        if n is None or n < _SEGMENT_BYTES:
            raise ValueError(tr('標籤之後不是參數區塊 (紀錄長度 %s)',
                                'no parameter block after the labels (record length %s)') % n)
        return decode(words36(f.read(_SEGMENT_BYTES)), MODE['GMS'])[0]


def read_tape(path):
    """磁帶映像 -> dict(header=檔頭, offset / length=各影像紀錄的位置與長度 (bytes), lcw=各影像紀錄的 LCW)。
    只讀第一個資料檔 (到資料後的 tape mark 為止)；檔案被截斷時用到最後一筆完整的紀錄。"""
    size = os.path.getsize(path)
    params, offset, length, heads = [], [], [], []
    with open(path, 'rb') as f:
        labels = _labels(f)
        while True:
            n = _length(f)
            pos = f.tell()
            if not n or pos + n > size:
                break
            if n == _BLOCK and not offset:
                params.append(words36(f.read(n)))
            elif n >= _LCW_BYTES:
                offset.append(pos)
                length.append(n)
                heads.append(f.read(_LCW_BYTES))
            f.seek(pos + n)
    if len(params) < 2:
        raise ValueError(tr('找不到參數區塊 (只有 %d 個)', 'parameter blocks not found (only %d)') % len(params))
    h = {'Labels': labels}
    for j, name in enumerate(_SEGMENTS[:4 * len(params)]):
        words = params[j // 4][(j % 4) * _SEGMENT_WORDS:(j % 4 + 1) * _SEGMENT_WORDS]
        if j < len(SEGMENTS['GMS']) or words.any():
            parse_segment(name, Words36(words), 'GMS', h)
    lcw = decode(words36(b''.join(heads)), LINE_CONTROL_WORD, len(heads))
    return dict(header=h, offset=np.array(offset, np.int64), length=np.array(length, np.int64), lcw=lcw)


class GmsTapeReader(VissrBase):
    NAME = 'GMS-TAPE'
    DESCRIPTION = 'GMS-1 VISSR 磁帶映像 (36-bit word；檔名為卷號，例如 VS2303)'
    DESCRIPTION_EN = 'GMS-1 VISSR tape image (36-bit words; the file name is the volume serial, e.g. VS2303)'
    PATTERNS = ('V?????',)

    def __init__(self, cfg):
        super().__init__(cfg)
        self._tapes = {}
        self._info = {}

    # ---- 檔案 --------------------------------------------------------------
    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            path = os.path.join(folder, name)
            if not os.path.isfile(path):
                continue
            try:
                mode = read_mode(path)
            except (OSError, ValueError) as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                continue
            if mode is None:                # 剛好同名但不是磁帶映像
                continue
            # Sensor mode：1 = 只有可見光、2 = 只有紅外、3 = 兩者
            bands = {1: ('VIS',), 2: ('IR',)}.get(int(mode['Sensor mode']), ('VIS', 'IR'))
            out.append(Scene(satellite_name(mode['Satellite name'], 'GMS'), observation_time(mode), 'FD',
                             files={b: path for b in bands}, extra={'gen': 'GMS'}))
        return sorted(out, key=lambda s: (s.time, s.satellite))

    # ---- 檔頭與影像結構 ----------------------------------------------------
    def _load(self, scene, channel):
        path = scene.files[channel]
        if (path, channel) in self._info:
            return self._info[path, channel]
        if path not in self._tapes:
            self._tapes[path] = read_tape(path)
        tape = self._tapes[path]
        kind = 'VIS' if channel == 'VIS' else 'IR'
        mode = tape['header']['Mode_Block']
        bits, n_lcw, n_doc, npix = (int(mode['%s-%s' % (kind, item)]) for item in (
            'Bit length', 'LCW-pixel size', 'DOC-pixel size', 'Number of pixels'))
        if bits not in _SAMPLES or n_lcw * bits != _LCW_BYTES * 8 or npix <= 0:
            raise ValueError(tr('%s 的影像紀錄格式不支援 (每樣本 %d 位元、LCW %d、像素 %d)',
                                'unsupported image record format for %s (%d bits per sample, LCW %d, %d pixels)')
                             % (channel, bits, n_lcw, npix))
        per, nbytes, _ = _SAMPLES[bits]
        first = n_lcw + n_doc                           # 第 1 個像素是紀錄中的第幾個樣本 (由 0 起算)
        seg = tape['lcw']['Data ID'] & 0x3FFFF          # 資料區分：1 = 紅外，2、4、8、16 = 可見光感測器 1~4
        sel = np.flatnonzero((seg >= 2) & (seg <= 16) if kind == 'VIS' else seg == 1)
        lcw = tape['lcw'][sel]
        usable = (lcw['Line name'] == 1) & (tape['length'][sel] >= -(-(first + npix) // per) * nbytes)
        h = dict(tape['header'], Line_Control_Word=lcw)
        info = self._index(h, 'GMS', channel, usable)
        info.update(path=path, offset=tape['offset'][sel], bits=bits, first=first, npix=npix)
        if kind == 'IR':                                # 9 位元的樣本：超過 255 的計數值 (不應出現) 視為無資料
            info['lut'] = np.concatenate([info['lut'], np.full((1 << bits) - len(info['lut']), np.nan, np.float32)])
        self._info[path, channel] = info
        return info

    # ---- 讀取 --------------------------------------------------------------
    @staticmethod
    def _counts(mm, info, recs, c0, c1):
        """紀錄 recs 的第 c0~c1 個像素的計數值 (len(recs), c1 - c0)。"""
        per, nbytes, unpack = _SAMPLES[info['bits']]
        s0, s1 = info['first'] + c0, info['first'] + c1
        g0, g1 = s0 // per, -(-s1 // per)               # 涵蓋這些樣本的組
        raw = np.empty((len(recs), (g1 - g0) * nbytes), np.uint8)
        for i, o in enumerate(info['offset'][recs]):
            raw[i] = mm[o + g0 * nbytes:o + g1 * nbytes]
        return unpack(raw)[:, s0 - g0 * per:s1 - g0 * per]

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        info = self._load(scene, channel)
        nrow = len(info['rec_of_row'])
        vis = info['kind'] == 'VIS'
        # 補缺行需要上下相鄰的線：多讀一列
        margin = 1 if vis and self.cfg.visible.get('fill_missing_lines', True) else 0
        a0, a1 = max(r0 - margin, 0), min(r1 + margin, nrow)
        recs = info['rec_of_row'][a0:a1]
        have = recs >= 0
        data = np.full((max(a1 - a0, 0), c1 - c0), np.nan, dtype=np.float32)
        if have.any():
            rows = np.flatnonzero(have)
            mm = np.memmap(info['path'], np.uint8, 'r')
            try:
                for i in range(0, len(rows), 256):      # 分批解開，暫存的計數值不會太大
                    part = rows[i:i + 256]
                    counts = self._counts(mm, info, recs[part], c0, c1)
                    if vis:
                        data[part] = info['lut'][info['sensor'][recs[part]][:, None], counts]
                    else:
                        data[part] = info['lut'][counts]
            finally:
                del mm
            if margin and len(recs) > 2:
                gap = np.flatnonzero(~have[1:-1] & have[:-2] & have[2:]) + 1
                data[gap] = 0.5 * (data[gap - 1] + data[gap + 1])
        if 0 <= r0 and r1 <= nrow:
            return data[r0 - a0:r1 - a0]
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)       # 視窗超出影像
        lo, hi = max(r0, 0), min(r1, nrow)
        if hi > lo:
            out[lo - r0:hi - r0] = data[lo - a0:hi - a0]
        return out


READER = GmsTapeReader
