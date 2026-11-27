"""NOAA HRPT 原始框 (Level 0) 的結構、檔案辨識、TLE 與定位 (AVHRR-L0 讀取器的輔助模組)

HRPT 小框 (minor frame) 每條掃描線一框，11090 個 10-bit 字 (POD Guide 表 4.1.3-1；KLM Guide §8.3.1.3.3.1 的
「HRPT MINOR FRAME TELEMETRY」)。以下字的位置 1 起算 (陣列索引 = 字號 − 1)：
  1~6      框同步 0x284 0x16F 0x35C 0x19D 0x20F 0x095
  7        ID：bit 9 (最高位) AVHRR 同步；bits 8-7 小框編號 1~3；bits 6-3 衛星位址；bit 2 重新同步；
           KLM 另有 bit 1：1 = 正常 AVHRR 輸入 (0 = 測試用的偽隨機碼)、bit 0：1 = 頻道 3A、0 = 3B
  9~12     時間碼：字 9 的 bits 9-1 = 年中日數；字 10 的 bits 6-0、字 11、字 12 = 當日毫秒 (27 bits)
  13~17    斜坡定標；18~20 PRT (同一顆 PRT 的三個讀數，每條掃描線輪換；每第 5 條是參考讀數，接近 0)
  23~52    星上黑體 (ICT)：頻道 3、4、5 各 10 個，交錯排列
  53~102   太空：頻道 1~5 各 10 個，交錯排列
  104~623  TIP；751~10990 影像：2048 像元 × 頻道 1~5 (BIP)
沒有年份 (年份取路徑中的日期，或檔案的修改時間) 也沒有軌道資料：定位用 TLE (SGP4) 推算衛星位置，
再依 AVHRR 的掃描幾何 (大地天底、各衛星的步進角與姿態偏差、每像元 25 µs) 求視線與 WGS84 橢球的交點 (見 geolocate)。

檔案
  - 每個字存成 16 bits (高 6 bits 為 0)，小端或大端都可 (自動判斷)；每框 22180 bytes，框的位置依同步字找。
  - 可以是未壓縮的檔、.xz / .gz / .bz2 壓縮檔，或未壓縮的 .tar 封存檔中的這些檔 (不必先解開；路徑寫成 封存檔::成員)。
  - 開頭可以有 ANSI 磁帶標籤 (VOL1、HDR1、HDR2…，各 80 bytes 的 ASCII 記錄)：千葉大學 CEReS 的 L0 檔 (.l0、.gi)
    前面有 vol1 + hdr1 + hdr2 = 240 bytes，辨識時跳過。vol1 的第 13~16 byte 是整個檔 (含標籤) 的大小 (大端)，
    .bz2 無法事先得知解壓後的大小，用它推算框數與結束時間 (與框長對不上時不用)。
  - 衛星：ID 字的衛星位址 (AAPP 的對照表，同一位址依年份區分)，其次路徑中的 noaa18、n14 等字樣。
  - 年份：路徑中的日期 (見 infer_year)；千葉的檔名 20141007163202.n19.l0 (日期緊接時刻) 與 n1297121720.gi
    (n + 衛星編號 + yymmddhh) 也可以 (原本支援的寫法都先比對，對不上才看這兩種)。
"""
import bz2
import datetime
import gzip
import io
import lzma
import os
import re
import tarfile
import zipfile
import zlib

import numpy as np

import avhrr_l1b_noaa as fmt
from satima.core.i18n import tr

FRAME_WORDS = 11090
FRAME_BYTES = 2 * FRAME_WORDS
SYNC = np.array([0x284, 0x16F, 0x35C, 0x19D, 0x20F, 0x095], dtype=np.uint16)
SYNC_MIN = 4                    # 框同步 6 個字至少要對幾個才算有效的框
FIRST_FRAMES = 12               # 辨識時先讀的框數
HEAD_FRAMES = 40                # 好的框太少時最多讀到的框數
LINE_SECONDS = 1. / 6.
SAMPLE_SECONDS = 25e-6          # AVHRR 每個樣本 25 µs
SCAN_DEG = 55.37                # 第 1、2048 像元的掃描角 (KLM 附錄 J；只用於 scan_angles 沒有指定步進角時)
STEP_DEG = 0.0540865            # 步進角 (度)：鏡面每秒 6 轉、每秒 39936 個樣本 (AAPP satid.txt)
# 各衛星的定位參數 (roll, pitch, yaw 度，步進角 度)：以海岸線擬合 (陸海遮罩 vs 紅外影像，2003-02-04 與 2013-11-07
# 日本附近 14 次升降軌過境，見 docs/PO.md)，擬合後殘差 RMS 0.3~0.5 km；以 2013-08-11 的過境獨立驗證，多數區塊在 1 km 內。符號同 geolocate：roll > 0 視線往左
# (第 2048 像元那側)、yaw > 0 見 geolocate。與 NOAA 1b 表頭的固定姿態修正、AAPP satid.txt 的預設偏差同量級
# (NOAA-18 roll 0.12~0.15°、NOAA-19 yaw 0.25°、NOAA-16 步進角 0.053955°)。pitch 與時鐘誤差分不開，設為 0。
NAV_PARAMS = {
    'NOAA-12': (0.295, 0., 0.124, 0.054107),          # 千葉 1997-12-17 的 2 次降軌 (AAPP：roll −4.9 mrad，符號相反)
    'NOAA-14': (-0.076, 0., -0.056, 0.0540865),       # 2003 年的資料受損，只有 2 次降軌 16 個區塊；步進角未擬合
    'NOAA-15': (-0.054, 0., 0.180, 0.054092),
    'NOAA-16': (-0.034, 0., -0.130, 0.053960),
    'NOAA-17': (-0.026, 0., 0.047, 0.054080),
    'NOAA-18': (-0.181, 0., 0.156, 0.054059),
    'NOAA-19': (-0.034, 0., 0.326, 0.054075),
}
NAV_DEFAULT = (-0.03, 0., 0., STEP_DEG)    # 沒有擬合的衛星：各衛星共有的約 −0.03° roll (約半個樣本)
TAR_SEP = '::'
COMPRESSED = {'.xz': lzma.open, '.lzma': lzma.open, '.gz': gzip.open, '.bz2': bz2.open}
LABEL_BYTES = 80                # ANSI 磁帶標籤記錄的長度
LABEL_IDS = (b'VOL', b'HDR', b'UHL', b'UVL', b'EOF', b'EOV', b'UTL')
SKIP_EXT = {'.png', '.jpg', '.jpeg', '.gif', '.txt', '.pdf', '.py', '.pyc', '.toml', '.nc', '.h5', '.he5', '.hdf',
            '.zip', '.7z', '.tif', '.tiff', '.xml', '.json', '.md', '.csv', '.log', '.ini', '.tle', '.dat'}

# ID 字的衛星位址 (bits 6-3) -> 衛星 (AAPP satid.txt；NOAA-14~19 由樣本確認)。同一位址有不同年代的衛星，依年份區分
ADDRESSES = {2: ('TIROS-N',), 3: ('NOAA-6', 'NOAA-13', 'NOAA-16'), 7: ('NOAA-7', 'NOAA-15'), 13: ('NOAA-8', 'NOAA-18'),
             11: ('NOAA-9', 'NOAA-17'), 15: ('NOAA-10', 'NOAA-19'), 1: ('NOAA-11',), 9: ('NOAA-12',), 5: ('NOAA-14',)}


def satellite_from_address(addr, year):
    """衛星位址 + 年份 -> 衛星 (運作年份 fmt.SERVICE 涵蓋該年者；都不涵蓋時取最近的)；位址不明時 None。"""
    cands = ADDRESSES.get(addr, ())
    if not cands:
        return None
    def gap(sat):
        a, b = fmt.SERVICE[sat]
        return 0 if a <= year <= b else min(abs(year - a), abs(year - b))
    return min(cands, key=gap)
# NORAD 編號 (TLE 用)
NORAD = {'TIROS-N': 11060, 'NOAA-6': 11416, 'NOAA-7': 12553, 'NOAA-8': 13923, 'NOAA-9': 15427, 'NOAA-10': 16969,
         'NOAA-11': 19531, 'NOAA-12': 21263, 'NOAA-13': 22739, 'NOAA-14': 23455, 'NOAA-15': 25338, 'NOAA-16': 26536,
         'NOAA-17': 27453, 'NOAA-18': 28654, 'NOAA-19': 33591}


# ---- 檔案 (未壓縮 / 壓縮 / tar 成員) -----------------------------------------------------------
def split_member(path):
    """'封存檔.tar::成員' -> (封存檔, 成員)；一般的檔 -> (path, None)。"""
    if TAR_SEP in path:
        a, b = path.split(TAR_SEP, 1)
        return a, b
    return path, None


def _decoder(name):
    return COMPRESSED.get(os.path.splitext(name)[1].lower())


def open_stream(path):
    """解壓後的位元組串流 (需 close)。"""
    archive, member = split_member(path)
    if member is None:
        opener = _decoder(path)
        return opener(path, 'rb') if opener else open(path, 'rb')
    tar = tarfile.open(archive, 'r:*')
    raw = tar.extractfile(member)
    opener = _decoder(member)
    stream = opener(raw, 'rb') if opener else raw
    return _Closing(stream, tar)


class _Closing(io.RawIOBase):
    """關閉時一併關閉 tar 檔。"""

    def __init__(self, stream, tar):
        self.stream, self.tar = stream, tar

    def readable(self):
        return True

    def read(self, n=-1):
        return self.stream.read(n)

    def readinto(self, b):
        data = self.stream.read(len(b))
        b[:len(data)] = data
        return len(data)

    def close(self):
        try:
            self.stream.close()
            self.tar.close()
        finally:
            super().close()


def read_bytes(path, nbytes=None):
    with open_stream(path) as f:
        return f.read() if nbytes is None else f.read(nbytes)


def mtime(path):
    """修改時間 (tar 成員取成員記錄的時間)。"""
    archive, member = split_member(path)
    try:
        if member is None:
            return datetime.datetime.utcfromtimestamp(os.path.getmtime(path))
        with tarfile.open(archive, 'r:*') as tar:
            return datetime.datetime.utcfromtimestamp(tar.getmember(member).mtime)
    except (OSError, KeyError, tarfile.TarError, ValueError, OverflowError):
        return None


def raw_size(path):
    """解壓後的大小 (bytes)；不讀整個檔就無法得知時 None。.xz 讀檔尾的索引、.gz 讀檔尾的 ISIZE (mod 2^32)。"""
    archive, member = split_member(path)
    try:
        if member is None:
            name, size = path, os.path.getsize(path)
            with open(path, 'rb') as f:
                return _container_size(f, 0, size, name)
        with tarfile.open(archive, 'r:*') as tar:
            info = tar.getmember(member)
            if not tar.fileobj or not hasattr(info, 'offset_data'):
                return None
        with open(archive, 'rb') as f:
            return _container_size(f, info.offset_data, info.size, member)
    except (OSError, KeyError, tarfile.TarError, ValueError):
        return None


def _container_size(f, offset, size, name):
    ext = os.path.splitext(name)[1].lower()
    if ext not in COMPRESSED:
        return size
    if ext == '.gz' and size >= 18:
        f.seek(offset + size - 4)
        return int.from_bytes(f.read(4), 'little')
    if ext == '.xz':
        return _xz_size(f, offset, size)
    return None


def _varint(buf, i):
    v, shift = 0, 0
    while True:
        b = buf[i]
        v |= (b & 0x7F) << shift
        i += 1
        if not b & 0x80:
            return v, i
        shift += 7


def _xz_size(f, offset, size):
    """xz 檔 (單一串流) 的解壓後大小：檔尾 12 bytes 的 Backward Size 指到索引，索引記錄各區塊的解壓後大小。"""
    if size < 32:
        return None
    f.seek(offset + size - 12)
    foot = f.read(12)
    if foot[10:12] != b'YZ':
        return None
    index_size = (int.from_bytes(foot[4:8], 'little') + 1) * 4
    f.seek(offset + size - 12 - index_size)
    idx = f.read(index_size)
    if not idx or idx[0] != 0:
        return None
    n, i = _varint(idx, 1)
    total = 0
    for _ in range(n):
        _, i = _varint(idx, i)
        u, i = _varint(idx, i)
        total += u
    return total


# ---- 辨識 ---------------------------------------------------------------------------------
class Structure:
    """order ('<' / '>')、offset (第一框的位元組位置)、satellite、generation、time、end、nframes (None = 未知)、
    year、doy (開頭的年中日數)、year_from (年份的來源：'path' / 'mtime')、address (衛星位址)。"""

    def __init__(self, **kw):
        self.__dict__.update(kw)


def label_size(buf):
    """開頭的 ANSI 磁帶標籤 (VOL1、HDR1、HDR2…：前 3 個字元為標籤代號、第 4 個為數字，每筆 80 bytes) 的總長度；沒有時 0。"""
    k = 0
    while buf[k:k + 3].upper() in LABEL_IDS and buf[k + 3:k + 4].isdigit():
        k += LABEL_BYTES
    return k


def label_file_size(buf, offset):
    """千葉的 vol1 標籤第 13~16 byte 記錄的檔案大小 (大端)；與第一框的位置、框長對得上才採用，否則 None。"""
    if buf[:4].upper() != b'VOL1':
        return None
    size = int.from_bytes(buf[12:16], 'big')
    if size >= offset + FRAME_BYTES and (size - offset) % FRAME_BYTES == 0:
        return size
    return None


def find_frames(buf):
    """在一段位元組中找框：回傳 (字序, 第一框的位元組位置)；找不到 None。
    條件：字都 <= 1023 (16 bits 存 10-bit 值)，且至少兩個同步字 (前 4 個字) 的間隔是框長的整數倍。"""
    for order in ('<', '>'):
        for start in (0, 1):
            n = (len(buf) - start) // 2
            if n < FRAME_WORDS + 6:
                continue
            w = np.frombuffer(buf, dtype=order + 'u2', count=n, offset=start)
            if (w > 1023).mean() > 0.01:
                continue
            hit = np.flatnonzero((w[:-3] == SYNC[0]) & (w[1:-2] == SYNC[1]) & (w[2:-1] == SYNC[2]) & (w[3:] == SYNC[3]))
            if hit.size < 2:
                continue
            phase = hit % FRAME_WORDS
            vals, counts = np.unique(phase, return_counts=True)
            if counts.max() < 2:
                continue
            return order, start + 2 * int(vals[np.argmax(counts)])
    return None


def _ten_bit(buf):
    """開頭是否像 16 bits 存的 10-bit 字 (任一種字序，幾乎都 <= 1023)：先用來快速排除其他格式。"""
    for order in ('<', '>'):
        for start in (0, 1):
            w = np.frombuffer(buf, dtype=order + 'u2', count=(len(buf) - start) // 2, offset=start)
            if (w > 1023).mean() <= 0.01:
                return True
    return False


def frame_words(buf, order, offset=0):
    """位元組 -> (框數, 11090) uint16 (小端)。"""
    n = (len(buf) - offset) // FRAME_BYTES
    w = np.frombuffer(buf, dtype=order + 'u2', count=n * FRAME_WORDS, offset=offset).reshape(n, FRAME_WORDS)
    return w.astype(np.uint16) if order == '>' else w


def sync_count(w):
    return (w[:, :6] == SYNC).sum(axis=1)


def time_code(w):
    """(年中日數, 當日毫秒) (框,)。"""
    day = (w[:, 8].astype(np.int64) >> 1) & 0x1FF
    ms = ((w[:, 9].astype(np.int64) & 0x7F) << 20) | ((w[:, 10].astype(np.int64) & 0x3FF) << 10) | \
        (w[:, 11].astype(np.int64) & 0x3FF)
    return day, ms


_DATE = re.compile(r'(?<!\d)((?:19[6-9]|20[0-9])\d)[-_.]?(0[1-9]|1[0-2])[-_.]?(0[1-9]|[12]\d|3[01])(?!\d)')
_SAT_NAMES = [(re.compile(r'tiros[-_ ]?n', re.I), 'TIROS-N'),
              (re.compile(r'noaa[-_ ]?(\d{1,2})(?!\d)', re.I), None),
              (re.compile(r'(?:^|[^a-z0-9])n(\d{2})(?=[_.-])', re.I), None)]
# 千葉大學 CEReS 的檔名 (上面的寫法都對不上時才用，原本能解析的路徑結果不變)：
#   日期緊接時刻 yyyymmddhhmm[ss]，例如 20141007163202.n19.l0
#   n + 衛星編號 + yymmdd[hh]，例如 n1297121720.gi (NOAA-12，1997/12/17 20 時)
_DATE_TIME = re.compile(r'(?<!\d)((?:19[6-9]|20[0-9])\d)(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])'
                        r'(?:[01]\d|2[0-3])[0-5]\d(?:[0-5]\d)?(?!\d)')
_CHIBA_NAME = re.compile(r'(?:^|[^a-z0-9])n(\d{2})(\d{2})(0[1-9]|1[0-2])(0[1-9]|[12]\d|3[01])(?:[01]\d|2[0-3])?'
                         r'(?=[_.-]|$)', re.I)


def _path_parts(path):
    """路徑的各段，由檔名往上。"""
    return list(reversed(re.split(r'[\\/]|' + re.escape(TAR_SEP), path)))


def satellite_from_path(path):
    """路徑中的衛星名稱 (由檔名往上找)；原本的寫法都找不到時再找千葉的 nSSyymmddhh。"""
    parts = _path_parts(path)
    for part in parts:
        for rx, name in _SAT_NAMES:
            m = rx.search(part)
            if m:
                sat = name or 'NOAA-%d' % int(m.group(1))
                if sat in fmt.GENERATION:
                    return sat
    for part in parts:
        for m in _CHIBA_NAME.finditer(part):
            sat = 'NOAA-%d' % int(m.group(1))
            if sat in fmt.GENERATION:
                return sat
    return None


def _two_digit_year(yy):
    return 1900 + yy if yy >= 57 else 2000 + yy


def _dates_in_path(path):
    """路徑中的日期 (由檔名往上)：原本的寫法 (_DATE) 全部排在前面，千葉的兩種寫法排在最後。"""
    parts = _path_parts(path)
    out, extra = [], []
    for part in parts:
        for m in _DATE.finditer(part):
            try:
                out.append(datetime.date(int(m.group(1)), int(m.group(2)), int(m.group(3))))
            except ValueError:
                pass
    for part in parts:
        found = [(int(m.group(1)), m.group(2), m.group(3)) for m in _DATE_TIME.finditer(part)]
        found += [(_two_digit_year(int(m.group(2))), m.group(3), m.group(4)) for m in _CHIBA_NAME.finditer(part)]
        for y, mo, d in found:
            try:
                extra.append(datetime.date(y, int(mo), int(d)))
            except ValueError:
                pass
    return out + extra


def infer_year(path, doy):
    """年份：路徑中的日期與年中日數相差 2 天以內者 (跨年時取前後一年)，其次檔案的修改時間 (修改時間的年中日數
    比資料早時取前一年)。回傳 (年, 'path' / 'mtime') 或 (None, None)。"""
    for d in _dates_in_path(path):
        for year in (d.year, d.year - 1, d.year + 1):
            try:
                t = datetime.date(year, 1, 1) + datetime.timedelta(days=doy - 1)
            except (ValueError, OverflowError):
                continue
            if abs((t - d).days) <= 2:
                return year, 'path'
    t = mtime(path)
    if t is None:
        return None, None
    year = t.year
    if doy > t.timetuple().tm_yday + 2:
        year -= 1
    return year, 'mtime'


def sniff(path):
    """依內容判斷是否為 HRPT 原始框；不是時 None (不丟例外)。"""
    try:
        with open_stream(path) as f:
            head = f.read(8192)
            skip = label_size(head)                 # 千葉的 vol1/hdr1/hdr2 標籤 (沒有標籤時 0)
            if len(head) < 8192 or not _ten_bit(head[skip:]):
                return None
            # 先讀 FIRST_FRAMES 框；接收開頭的框常同步不良，好的框太少時再多讀到 HEAD_FRAMES 框
            for nframes in (FIRST_FRAMES, HEAD_FRAMES):
                head += f.read((nframes + 1) * FRAME_BYTES + skip - len(head))
                found = find_frames(head[skip:]) if len(head) - skip >= 2 * FRAME_BYTES else None
                if found is not None:
                    found = (found[0], found[1] + skip)
                    w = frame_words(head, *found)
                    good = sync_count(w) >= 5
                    if good.sum() >= 5:
                        break
    except (OSError, EOFError, lzma.LZMAError, tarfile.TarError, zlib.error, ValueError):
        return None
    if found is None or good.sum() < 2:
        return None
    order, offset = found
    w = w[good]
    rows = np.flatnonzero(good)
    day, ms = time_code(w)
    ok = (day >= 1) & (day <= 366) & (ms < 86_400_000)
    if ok.sum() < 2:
        return None
    vals, counts = np.unique(day[ok], return_counts=True)
    doy = int(vals[np.argmax(counts)])
    year, year_from = infer_year(path, doy)
    if year is None:
        return None
    address = (w[:, 6] >> 3) & 0xF
    vals, counts = np.unique(address, return_counts=True)
    addr = int(vals[np.argmax(counts)])
    sat = satellite_from_address(addr, year) if counts.max() >= 0.6 * address.size else None
    named = satellite_from_path(path)
    if named and (sat is None or named in ADDRESSES.get(addr, ())):
        sat = named
    if sat is None or sat not in NORAD:
        return None
    sel = ok & (day == doy)
    # 起始時間：各框「時間 − 框序 × 掃描週期」的中位數 (開頭的時間碼常有位元錯誤)
    t_rel = ms[sel] - (rows[sel] * LINE_SECONDS * 1000.)
    ms0 = float(np.median(t_rel))
    time = datetime.datetime(year, 1, 1) + datetime.timedelta(days=doy - 1, milliseconds=ms0)
    size = raw_size(path)
    if size is None and skip:                       # .bz2 等不知道解壓後大小的檔：用千葉 vol1 標籤記錄的大小
        size = label_file_size(head, offset)
    nframes = (size - offset) // FRAME_BYTES if size else None
    end = time + datetime.timedelta(seconds=nframes * LINE_SECONDS) if nframes else None
    return Structure(order=order, offset=offset, satellite=sat, generation=fmt.GENERATION.get(sat, 3), time=time,
                     end=end, nframes=nframes, year=year, doy=doy, year_from=year_from, address=addr)


def candidates(root, depth=2):
    """資料夾中可能是原始框的檔 (相對路徑)：一般檔、壓縮檔、.tar 封存檔中的成員；往下找 depth 層子資料夾。"""
    out = []
    root = root or '.'

    def walk(folder, rel, level):
        try:
            entries = sorted(os.listdir(folder))
        except OSError:
            return
        for name in entries:
            full = os.path.join(folder, name)
            relname = os.path.join(rel, name) if rel else name
            if os.path.isdir(full):
                if level < depth:
                    walk(full, relname, level + 1)
                continue
            ext = os.path.splitext(name)[1].lower()
            low = name.lower()
            if ext == '.tar':                   # 未壓縮的 tar 才讀成員 (壓縮的 tar 要整檔解壓才能列出成員)
                out.extend(_tar_members(full, relname))
            elif low.endswith(('.tar.gz', '.tgz', '.tar.xz', '.tar.bz2')):
                continue
            elif ext in COMPRESSED:
                if _inner_ext(name) not in SKIP_EXT:      # 例如 HSD 的 .DAT.bz2 不必解壓檢查
                    out.append(relname)
            elif ext not in SKIP_EXT and _size_ok(full):
                out.append(relname)

    walk(root, '', 0)
    return out


def _inner_ext(name):
    """壓縮檔內層的副檔名 (x.DAT.bz2 -> .dat)。"""
    return os.path.splitext(os.path.splitext(name)[0])[1].lower()


def _size_ok(path):
    try:
        return os.path.getsize(path) >= 3 * FRAME_BYTES
    except OSError:
        return False


def _tar_members(full, relname):
    try:
        with tarfile.open(full, 'r:*') as tar:
            out = []
            for m in tar.getmembers():
                ext = os.path.splitext(m.name)[1].lower()
                if not m.isfile() or ext in SKIP_EXT:
                    continue
                if ext in COMPRESSED and _inner_ext(m.name) not in SKIP_EXT or m.size >= 3 * FRAME_BYTES:
                    out.append(relname + TAR_SEP + m.name)
            return out
    except (OSError, tarfile.TarError):
        return []


# ---- TLE ----------------------------------------------------------------------------------
_TLE_CACHE = {}


def _tle_epoch(line1):
    yy = int(line1[18:20])
    doy = float(line1[20:32])
    year = 1900 + yy if yy >= 57 else 2000 + yy
    return datetime.datetime(year, 1, 1) + datetime.timedelta(days=doy - 1)


def _parse_tle_text(text, table):
    lines = [ln.rstrip() for ln in text.splitlines()]
    for a, b in zip(lines, lines[1:]):
        if len(a) >= 64 and len(b) >= 63 and a.startswith('1 ') and b.startswith('2 ') and a[2:7] == b[2:7]:
            try:
                num = int(a[2:7])
                epoch = _tle_epoch(a)
                float(b[52:63])
            except ValueError:
                continue
            table.setdefault(num, {})[epoch] = (a[:69], b[:69])


def load_tles(folders):
    """資料夾 (遞迴) 中所有 TLE；見 load_tle_files。"""
    files = []
    for folder in folders:
        if folder and os.path.isdir(folder):
            for dirpath, _, names in os.walk(folder):
                files.extend(os.path.join(dirpath, n) for n in sorted(names))
    return load_tle_files(files)


def load_tle_files(files):
    """TLE 檔 (任何副檔名的文字檔，或 .zip 中的文字檔；2 行或 3 行格式，可含多顆衛星、多個時間)
    -> {NORAD 編號: [(epoch, line1, line2)] (依 epoch 排序)}。依檔案的大小與修改時間快取。"""
    stamp = []
    for full in files:
        try:
            st = os.stat(full)
        except OSError:
            continue
        if st.st_size <= 200 * 2 ** 20:
            stamp.append((os.path.abspath(full), st.st_size, st.st_mtime))
    stamp = tuple(stamp)
    if stamp in _TLE_CACHE:
        return _TLE_CACHE[stamp]
    table = {}
    for full, _, _ in stamp:
        try:
            if full.lower().endswith('.zip'):
                with zipfile.ZipFile(full) as z:
                    for info in z.infolist():
                        if info.file_size < 200 * 2 ** 20:
                            _parse_tle_text(z.read(info).decode('latin1'), table)
                continue
            with open(full, 'rb') as f:
                head = f.read(4096)
                if b'\x00' in head or not re.search(rb'(^|\n)1 [ \d]{5}', head):
                    continue
                _parse_tle_text((head + f.read()).decode('latin1'), table)
        except (OSError, zipfile.BadZipFile):
            continue
    out = {num: sorted((e, a, b) for e, (a, b) in d.items()) for num, d in table.items()}
    _TLE_CACHE[stamp] = out
    return out


def _plausible(line2):
    """太陽同步極軌衛星的 TLE (傾角 95~105 度、每天 12~15.5 圈)；排除歷史檔中編號標錯的資料。"""
    try:
        return 95. <= float(line2[8:16]) <= 105. and 12. <= float(line2[52:63]) <= 15.5
    except ValueError:
        return False


TLE_EXT = {'.tle', '.txt', '.2le', '.3le'}     # 資料目錄中當作 TLE 讀的副檔名 ([tle] dir 則讀所有檔)
TLE_SPREAD_KM = 20.                             # 最接近的幾組 TLE 推算的位置相差超過此值時警告


def tle_table(cfg):
    """[tle] dir (遞迴，相對路徑以 Satima.py 所在目錄為準) 與目前目錄 (資料目錄；只讀 TLE_EXT 的檔) 中的 TLE。"""
    folder = str(cfg.tle.dir or '').strip()
    if folder and not os.path.isabs(folder):
        folder = os.path.join(cfg.get('_home', '') or '', folder)
    tables = [load_tles([folder])] if folder else []
    here = os.getcwd()
    try:
        names = [n for n in os.listdir(here) if os.path.splitext(n)[1].lower() in TLE_EXT]
    except OSError:
        names = []
    if names:
        tables.append(load_tle_files([os.path.join(here, n) for n in names]))
    out = {}
    for t in tables:
        for num, items in t.items():
            out.setdefault(num, []).extend(items)
    return {num: sorted(set(items)) for num, items in out.items()}


def navigate(cfg, satellite, when, times, cols, label):
    """以 TLE 定位：取 epoch 最接近 when 的 TLE (印出 epoch；超過 [tle] max_age 天、或最接近的幾組互相矛盾時警告)，
    套用衛星的定位參數 (NAV_PARAMS) 計算 cols 的經緯度 (列, 像元)。沒有該衛星的 TLE 時 ValueError。"""
    near = nearest_tles(tle_table(cfg), satellite, when)
    if not near:
        raise ValueError(tr('%s：找不到 %s (NORAD %d) 的 TLE，無法定位。請把 TLE 檔放到 %s 或資料目錄 '
                            '(python tools/fetch_spacetrack_tle.py 可從 Space-Track 下載)',
                            '%s: no TLE for %s (NORAD %d); cannot navigate. Put TLE files in %s or the data folder '
                            '(python tools/fetch_spacetrack_tle.py downloads them from Space-Track)')
                         % (label, satellite, NORAD[satellite], cfg.tle.dir))
    epoch, line1, line2 = near[0]
    age = (when - epoch).total_seconds() / 86400.
    if age >= 0:
        msg = tr('%s：TLE epoch %s (觀測前 %.1f 天)', '%s: TLE epoch %s (%.1f days before the observation)')
    else:
        msg = tr('%s：TLE epoch %s (觀測後 %.1f 天)', '%s: TLE epoch %s (%.1f days after the observation)')
    print(msg % (label, epoch.strftime('%Y-%m-%d %H:%M'), abs(age)))
    if abs(age) > float(cfg.tle.max_age):
        print(tr('注意：TLE 與觀測相差 %.1f 天 (超過 [tle] max_age = %g)，定位誤差可能達數公里以上',
                 'NOTE: the TLE is %.1f days from the observation (over [tle] max_age = %g); navigation errors may '
                 'reach several km or more') % (abs(age), float(cfg.tle.max_age)))
    spread = tle_spread(near, when)
    if spread > TLE_SPREAD_KM:
        print(tr('注意：最接近的 %d 組 TLE (epoch %s) 推算的衛星位置相差達 %.0f km，TLE 檔中可能有錯誤的根數，定位可能不準',
                 'NOTE: the %d nearest TLEs (epochs %s) disagree by up to %.0f km in satellite position; the TLE files '
                 'may contain bad elements and the navigation may be off')
              % (len(near), ', '.join(e.strftime('%m-%d %H:%M') for e, _, _ in near), spread))
    roll, pitch, yaw, step = NAV_PARAMS.get(satellite, NAV_DEFAULT)
    return geolocate(line1, line2, times, cols, attitude=(roll, pitch, yaw), step_deg=step)


def pick_tle(tles, satellite, when):
    """最接近 when 的 TLE (epoch, line1, line2)；沒有該衛星的 TLE 時 None。"""
    items = nearest_tles(tles, satellite, when, 1)
    return items[0] if items else None


def nearest_tles(tles, satellite, when, n=3):
    """epoch 最接近 when 的 n 組 TLE (依接近程度排序)。"""
    items = [x for x in tles.get(NORAD.get(satellite), ()) if _plausible(x[2])]
    return sorted(items, key=lambda x: abs((x[0] - when).total_seconds()))[:n]


def tle_spread(items, when):
    """幾組 TLE 推算到 when 的衛星位置，與第一組 (最接近的) 的最大距離 (km)；用來發現歷史檔中互相矛盾的 TLE。"""
    if len(items) < 2:
        return 0.
    t = np.array([np.datetime64(when, 'us')])
    ref = propagate(items[0][1], items[0][2], t)[0][0]
    return max(float(np.linalg.norm(propagate(a, b, t)[0][0] - ref)) for _, a, b in items[1:])


# ---- 軌道與定位 -----------------------------------------------------------------------------
WGS84_A = 6378.137
WGS84_F = 1. / 298.257223563
WGS84_B = WGS84_A * (1. - WGS84_F)
WGS84_E2 = WGS84_F * (2. - WGS84_F)


def propagate(line1, line2, times):
    """SGP4：times (datetime64) -> TEME 位置 (n, 3) km、速度 (n, 3) km/s。有 sgp4 套件時用它，否則用 pyorbital。"""
    t = np.asarray(times, dtype='datetime64[us]')
    try:
        from sgp4.api import Satrec
    except ImportError:
        Satrec = None
    if Satrec is not None:
        sat = Satrec.twoline2rv(line1, line2)
        us = (t - np.datetime64('2000-01-01T12:00:00', 'us')).astype(np.int64)
        days = us // 86_400_000_000
        jd = 2451545.0 + days.astype(np.float64)
        fr = (us - days * 86_400_000_000) / 86_400e6
        err, r, v = sat.sgp4_array(jd, fr)
        r, v = np.asarray(r, dtype=np.float64), np.asarray(v, dtype=np.float64)
        r[err != 0] = np.nan
        return r, v
    from pyorbital.orbital import Orbital
    orb = Orbital('NOAA', line1=line1, line2=line2)
    pos, vel = orb.get_position(t, normalize=False)
    return np.asarray(pos, dtype=np.float64).T, np.asarray(vel, dtype=np.float64).T


def gmst(times):
    """格林威治平恆星時 (弧度)，IAU 1982 (同 SGP4 的 TEME；UT1 以 UTC 代替)。"""
    t = np.asarray(times, dtype='datetime64[us]')
    d = (t - np.datetime64('2000-01-01T12:00:00', 'us')).astype(np.int64) / 86_400e6
    tu = d / 36525.
    sec = 67310.54841 + (876600. * 3600. + 8640184.812866) * tu + 0.093104 * tu ** 2 - 6.2e-6 * tu ** 3
    return np.deg2rad((sec % 86400.) / 240.)


def scan_angles(cols, step_deg=None):
    """像元 (0 起算) -> 掃描角 (弧度)；負 = 第 1 像元那一側 (前進方向的右邊)。"""
    step = SCAN_DEG / 1023.5 if step_deg is None else step_deg
    return np.deg2rad((np.asarray(cols, dtype=np.float64) - 1023.5) * step)


def geolocate(line1, line2, times, cols, attitude=(0., 0., 0.), nadir='geodetic', step_deg=None):
    """每條掃描線 (時間 times，datetime64) 在像元 cols (0 起算) 的大地緯度、經度 (度)，(列, 像元) float64。
    天底：nadir = 'geodetic' (橢球法線，KLM Guide §2.2 與 AAPP 的作法) 或 'geocentric' (地心方向)。以海岸線檢驗
    (2003 年 NOAA-15/16/17 升降軌)，大地天底的南北誤差約 −1 km，地心天底約 −3.5 km；2013 年以後 NOAA 業務定位也是大地天底
    (1998 年的 CLASS 1b 則較接近地心天底)。
    掃描面垂直於慣性速度；掃描角 = (像元 − 1023.5) × 步進角 (step_deg，預設 SCAN_DEG / 1023.5)；
    attitude = (roll, pitch, yaw) 度：roll > 0 視線往左 (第 2048 像元那側)，pitch > 0 往前，yaw 為繞天底的旋轉。
    每個像元的時間 = 掃描線時間 + 像元 × 25 µs (衛星位置以速度線性外推)。看不到地球的方向為 NaN。"""
    t = np.asarray(times, dtype='datetime64[us]')
    cols = np.asarray(cols, dtype=np.float64)
    r0, v0 = propagate(line1, line2, t)
    dt = cols * SAMPLE_SECONDS                                      # (k,)
    r = r0[:, None, :] + v0[:, None, :] * dt[None, :, None]         # (n, k, 3)
    v = np.broadcast_to(v0[:, None, :], r.shape)
    if nadir == 'geodetic':                                         # 橢球法線
        x, y, z = r[..., 0], r[..., 1], r[..., 2]
        p = np.hypot(x, y)
        lat = np.arctan2(z, p * (1. - WGS84_E2))
        for _ in range(5):
            n_r = WGS84_A / np.sqrt(1. - WGS84_E2 * np.sin(lat) ** 2)
            h = p / np.cos(lat) - n_r
            lat = np.arctan2(z, p * (1. - WGS84_E2 * n_r / (n_r + h)))
        lam = np.arctan2(y, x)
        nadir = -np.stack([np.cos(lat) * np.cos(lam), np.cos(lat) * np.sin(lam), np.sin(lat)], axis=-1)
    else:
        nadir = -r / np.linalg.norm(r, axis=-1, keepdims=True)      # 地心天底
    fwd = v - (v * nadir).sum(axis=-1, keepdims=True) * nadir
    fwd /= np.linalg.norm(fwd, axis=-1, keepdims=True)
    right = np.cross(nadir, fwd)
    roll, pitch, yaw = np.deg2rad(attitude)
    th = scan_angles(cols, step_deg)[None, :] + roll
    # 繞前進軸轉 θ (θ > 0 往左)，再依 pitch (往前)、yaw (繞天底) 修正
    look = nadir * np.cos(th)[..., None] - right * np.sin(th)[..., None]
    if pitch:
        look = look * np.cos(pitch) + fwd * np.sin(pitch)
    if yaw:
        c, s = np.cos(yaw), np.sin(yaw)
        f2 = fwd * c + right * s
        r2 = right * c - fwd * s
        a = (look * fwd).sum(-1, keepdims=True)
        b = (look * right).sum(-1, keepdims=True)
        d = (look * nadir).sum(-1, keepdims=True)
        look = a * f2 + b * r2 + d * nadir
    # 與橢球的交點
    scale = np.array([1. / WGS84_A, 1. / WGS84_A, 1. / WGS84_B])
    rs, ls = r * scale, look * scale
    qa = (ls * ls).sum(-1)
    qb = (rs * ls).sum(-1)
    qc = (rs * rs).sum(-1) - 1.
    disc = qb * qb - qa * qc
    with np.errstate(invalid='ignore'):
        s = (-qb - np.sqrt(disc)) / qa
    s[disc < 0] = np.nan
    pt = r + look * s[..., None]
    # TEME -> 地固 (繞 z 軸轉 GMST)
    tt = t[:, None] + (dt * 1e6).astype('timedelta64[us]')[None, :]
    lon = np.arctan2(pt[..., 1], pt[..., 0]) - gmst(tt)
    glat = np.arctan2(pt[..., 2], (1. - WGS84_E2) * np.hypot(pt[..., 0], pt[..., 1]))
    lon = (np.rad2deg(lon) + 180.) % 360. - 180.
    return np.rad2deg(glat), lon
