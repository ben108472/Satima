"""GOES-8~15 Imager — GVAR 原始資料 (NOAA CLASS 的 goes13.2015.296.001519、中國 NSMC 的 GOS9_GVR_MLT_NUL_*.dat)

格式 (GOES I-P GVAR，DRL 504-02 Section 3；結構與欄位名稱同舊版 GVAR_format.py，見 sensors/gvar_format.py)：
  資料由 GVAR 區塊組成，每個掃描 (可見光 8 條線) 為 Block 0 (文件，8 位元字) + Block 1~2 (紅外紀錄) +
  Block 3~10 (可見光 8 個偵測器，Block 3 最北)，10 位元字。每個區塊前有 30 bytes 的區塊頭，重複 3 次。
  Block 1~10 的每筆紀錄 = Line Documentation (16 字：頻道、偵測器、side、掃描序號 RISCT、像素數…) + 像素。
  NSMC 檔每個區塊固定佔 26624 bytes (開頭 2 bytes 對齊)；CLASS 檔則依區塊頭的字數依序排列
  (區塊頭 720 bits + 資料 (字數 - 2) x 字長 + CRC 16 bits，補到整數 bytes)。
  CLASS 另外提供的 <檔名>.indx (每個區塊一行：檔名、32 位元字位移、位元位移、區塊位元數、Block ID，Block 0 另有
  INSLN；資料損壞的區塊在最後標 ?) 不需要：依區塊頭找到的位置與 .indx 相同 (GOES-8 1998/10/26 四個檔逐塊比對)。
容錯 (原始資料常有位元錯誤、區塊或參數缺失，但影像本身多半有效)：
  區塊頭以 3 份逐 byte 多數決；長度不合理時往後找下一個區塊頭 (CLASS) 或跳到下一個固定位置 (NSMC)。
  每筆紀錄依自己 Line Documentation 的掃描序號 RISCT 定出行號，並與 Block 0 的 INSLN、RISCT 交叉驗證；
  導航參數、影像起點等取所有 Block 0 的眾數，個別 Block 0 損壞不影響。
  不以 CRC 捨棄區塊：CLASS 檔約 0.2~1.3% 的區塊 CRC 錯誤 (CRC-16/CCITT，初值 0xFFFF、結果取反)，但絕大多數
  只有少數像素的位元錯誤，其餘像素正確，照常讀出。
  地面站失鎖時填入的偽隨機資料 (一條線從某處起整段是 0~1023 的亂數，CLASS .indx 標 ? 的區塊即是) 不是影像；
  [gvar] mask_fill = true 時標為無資料 (預設 false：照原值畫出)。判定：每 32 像素中至少 24 個落在上下兩列範圍
  ±40 計數值之外，且沿線相鄰像素差平均 >= 100 計數值；真實影像 (含 CRC 錯誤區塊的零星位元錯誤) 最多 23 個、
  平均差 < 50，不受影響 (GOES-8/12/13/15 共 11 個檔驗證)。
波段 (BAND-01~06 = GVAR 頻道 1~6)：GOES-8~11 有 1~5、GOES-12~15 有 1、2、3、4、6；
  可見光 1 km、紅外 4 km (GOES-8~11 的 BAND-03 與 GOES-12/13 的 BAND-06 為 8 km)。
  紅外每條線涵蓋 4 (或 8) 條可見光線、每個像素 4 個可見光像素；行列號以可見光為單位換算後導航。
校正：紅外依 conversion_table/GOES 的計數值 -> 亮溫表 (衛星、side、頻道、偵測器)；
      可見光依 Block 0 各偵測器的係數 (GOES-14/15 沒有，用 NOAA 公布值)。
導航：GOES I-M O&A 模型 (sensors/goes_imager.py，同 McIDAS nvxgvar)。
完整檔頭：reader.header(scene) 回傳 {'Block_0': 所有掃描的 Block 0 (numpy 結構)、'Block_0_b': 後段 (I-N 或 O-P 格式)、
  'records': 所有紀錄的 Line Documentation 與位置、'nav': 導航參數}。
"""
import collections
import datetime
import os
import re

import numpy as np

import gvar_format as fmt
from goes_imager import (IR_LINE_CENTER, IR_PIXEL_CENTER, GoesNav, SATELLITE, bcd_time, fit_line_time, gould, ir_lut,
                         make_channel, params_from_block0, sector_label, unpack10, vis_coefficients, vis_scale)
from satima.core.i18n import tr
from satima.geo.reader import GeoReader, Scene

_CLASS =re.compile(r'^goes(\d\d)\.(\d{4})\.(\d{3})\.(\d{6})$', re.I)
_NSMC = re.compile(r'^GOS(\d+)_GVR_\w*?_(\d{8})_(\d{4})', re.I)
_SLOT = 26624                       # NSMC 每個區塊固定佔的 bytes
_HEADER = 90                        # 區塊頭 30 bytes x 3
_B0 = 8040                          # Block 0 資料 (8 位元字) 的 bytes
_REC = np.dtype([('block', 'u1'), ('offset', 'i8'), ('word', 'i4'), ('scan', 'i4'), ('SPCID', 'u2'), ('SPSID', 'u2'),
                 ('LSIDE', 'u2'), ('LIDET', 'u2'), ('LICHA', 'u2'), ('RISCT', 'u4'), ('LPIXLS', 'u4'), ('LWORDS', 'u4'),
                 ('LZCOR', 'u2'), ('LLAG', 'u2'), ('slot', 'u1')])
# 填充資料 (地面站失鎖時的亂數) 的判定：_FILL_W 像素中至少 _FILL_OUT 個落在上下兩列範圍 ±_FILL_TOL 之外，
# 且沿線相鄰像素差的平均 >= _FILL_JUMP (計數值)
_FILL_W, _FILL_OUT, _FILL_TOL, _FILL_JUMP = 32, 24, 40, 100


def _mode(values):
    values = [v for v in values if v is not None]
    return collections.Counter(values).most_common(1)[0][0] if values else None


def _vote(raw):
    """區塊頭 3 份逐 byte 多數決。"""
    a = np.frombuffer(raw, np.uint8).reshape(3, 30)
    out = a[0].copy()
    same = a[1] == a[2]
    out[same] = a[1][same]
    return out


def _ldoc(w):
    return (int(w[0]), int(w[1]), int(w[2]), int(w[3]), int(w[4]), int(w[5]) << 10 | int(w[6]),
            int(w[9]) << 10 | int(w[10]), int(w[11]) << 10 | int(w[12]), int(w[13]), int(w[14]))


class GvarFile:
    """一個 GVAR 檔案的索引：所有 Block 0 與各頻道每條線的紀錄位置。"""

    def __init__(self, path, slot=None):
        self.path = path
        self.size = os.path.getsize(path)
        with open(path, 'rb') as f:
            self._scan(f, slot)

    @staticmethod
    def _header_ok(h):
        bid, ws, wc = int(h[0]), int(h[1]), int(h[2]) << 8 | int(h[3])
        if bid == 240:
            bid = 0
        if bid > 11 or ws not in (6, 8, 10) or wc < 2:
            return None
        n = -(-(wc - 2) * ws // 8) + 2         # 10 位元字數 x 10 / 8 不是整數時無條件進位 (舊版奇數字數 +1)
        if (bid == 0 and ws == 8 and n < _B0) or n > 65536:
            return None
        return bid, ws, n, int(h[6])

    @classmethod
    def _find_header(cls, buf):
        """buf 中第一個區塊頭的位置：3 份中至少 2 份相同、多數決結果合理。
        不要求第 1 份完好 (常有位元錯誤)；取最早的位置，避免把第 2、3 份誤當成起點而錯開 30 bytes。"""
        a = np.frombuffer(buf, np.uint8)
        if a.size < _HEADER:
            return None
        e = np.r_[0, np.cumsum(a[:-30] == a[30:])]
        run = e[30:] - e[:-30]                       # run[k] == 30：buf[k:k+30] 與 buf[k+30:k+60] 相同
        n = a.size - _HEADER + 1
        same = np.zeros(n, bool)
        same |= run[:n] == 30                        # 第 1、2 份相同
        same |= run[30:30 + n] == 30                 # 第 2、3 份相同
        for k in np.flatnonzero(same):
            if cls._header_ok(_vote(buf[k:k + _HEADER])) is not None:
                return int(k)
        return None

    def _resync(self, f, start, chunk=1 << 20):
        """從 start 往後找下一個區塊頭。"""
        while start + _HEADER <= self.size:
            f.seek(start)
            k = self._find_header(f.read(chunk + _HEADER))
            if k is not None:
                return start + k
            start += chunk
        return None

    def _first_offset(self, f):
        """NSMC：固定槽位的起點 (第一個可辨識的區塊頭，不一定是 Block 0) 對槽位長度取餘數。"""
        f.seek(0)
        k = self._find_header(f.read(4 * _SLOT + _HEADER))
        return 0 if k is None else k % _SLOT

    def _scan(self, f, slot):
        """逐區塊建立索引。不檢查 CRC 與 Data Valid 旗標 (常有錯誤但影像正常)，只看紀錄本身是否合理。"""
        blocks0, recs = [], []
        bad = 0
        pos = self._first_offset(f) if slot else 0
        last_data = 0
        scan = -1
        while pos + _HEADER <= self.size:
            f.seek(pos)
            raw = f.read(_HEADER)
            info = self._header_ok(_vote(raw))
            if info is None:
                bad += 1
                if slot:
                    pos += slot
                    continue
                # 上一個區塊頭的字數可能有誤：從上一個區塊的資料起點往後找下一個區塊頭
                nxt = self._resync(f, last_data + 20 if last_data else pos + 1)
                if nxt is None or nxt <= last_data:
                    nxt = self._resync(f, pos + 1)
                if nxt is None:
                    break
                pos = nxt
                continue
            bid, ws, n, _ = info
            data_pos = last_data = pos + _HEADER
            if bid == 0 and ws == 8:
                b = f.read(_B0)
                if len(b) == _B0:
                    blocks0.append(b)
                    scan = len(blocks0) - 1
            elif 1 <= bid <= 10 and scan >= 0:
                if bid <= 2:
                    raw = f.read(max(n - 2, 0))
                    w = unpack10(raw)
                    k = j = 0
                    while k + 16 <= len(w) + 16 and j < 4:
                        if k + 16 > len(w):
                            break
                        d = _ldoc(w[k:k + 16])
                        lwords, npix = d[7], d[6]
                        if not (0 < npix <= 8000 and 16 + npix <= lwords <= 16 + npix + 64):
                            break
                        if k + lwords > len(w):          # 區塊頭的字數偏少：多讀一些
                            f.seek(data_pos)
                            w = unpack10(f.read((k + lwords) * 10 // 8 + 5))
                            if k + lwords > len(w):
                                break
                        recs.append((bid, data_pos, k + 16, scan) + d + (j,))
                        k += lwords
                        j += 1
                else:
                    d = _ldoc(unpack10(f.read(20)))
                    if 0 < d[6] <= 30680 and d[7] >= 16 + d[6]:
                        recs.append((bid, data_pos, 16, scan) + d + (0,))
            pos = pos + slot if slot else data_pos + n
        if not blocks0:
            raise ValueError(tr('找不到 GVAR Block 0', 'GVAR Block 0 not found'))
        self.bad_headers = bad
        self.block0 = np.frombuffer(b''.join(b[:fmt._Block_0.itemsize] for b in blocks0), fmt._Block_0)
        spcid = _mode(int(x) for x in self.block0['SPCID'] if 8 <= x <= 15) or 8
        self.satellite = SATELLITE[spcid]
        tail = fmt._Block_0_O if spcid >= 14 else fmt._Block_0_I_N
        self.block0b = np.frombuffer(b''.join(b[fmt._Block_0.itemsize:fmt._Block_0.itemsize + tail.itemsize]
                                              for b in blocks0), tail)
        self.records = np.array(recs, dtype=_REC) if recs else np.zeros(0, _REC)
        self._repair()
        self._lines()

    def _repair(self):
        """Line Documentation 個別欄位的位元錯誤：同一區塊、同一位置的紀錄在整個檔案中頻道與偵測器固定，
        以眾數修正 (可見光 Block 3~10 一定是頻道 1)；side 不是 0 / 1023 時用眾數。"""
        r = self.records
        if not len(r):
            return
        self.repaired = 0
        for bid in range(1, 11):
            for j in range(4):
                sel = np.flatnonzero((r['block'] == bid) & (r['slot'] == j))
                if not sel.size:
                    continue
                cha = 1 if bid >= 3 else _mode(r['LICHA'][sel].tolist())
                det = _mode(r['LIDET'][sel].tolist())
                wrong = (r['LICHA'][sel] != cha) | (r['LIDET'][sel] != det)
                self.repaired += int(wrong.sum())
                r['LICHA'][sel] = cha
                r['LIDET'][sel] = det
        side = _mode([int(x) for x in r['LSIDE'] if x in (0, 1023)]) or 0
        r['LSIDE'][(r['LSIDE'] != 0) & (r['LSIDE'] != 1023)] = side

    def _lines(self):
        """每筆紀錄的掃描起始行 (可見光行號)：Line Documentation 的 RISCT 與 Block 0 的 INSLN、RISCT 交叉驗證。"""
        b0, r = self.block0, self.records
        insln, risct = b0['INSLN'].astype(np.int64), b0['RISCT'].astype(np.int64)
        ok = (risct >= 1) & (risct <= 1974) & (insln >= 1) & (insln <= 15787)
        self.first_line = int(_mode((insln - 8 * (risct - 1))[ok].tolist()) or insln[0])
        self.west_pixel = int(_mode([int(x) for x in b0['IWFPX'] if 0 < x < 30681]) or 1)
        base_rec = self.first_line + 8 * (r['RISCT'].astype(np.int64) - 1)
        base_b0 = insln[r['scan']]
        base_b0r = self.first_line + 8 * (risct[r['scan']] - 1)
        rec_ok = (r['RISCT'] >= 1) & (r['RISCT'] <= 1974)
        base = np.where(rec_ok & ((base_rec == base_b0) | (base_rec == base_b0r)), base_rec,
                        np.where(base_b0 == base_b0r, base_b0, np.where(rec_ok, base_rec, -1)))
        self.base_line = base

    def line_time(self):
        """各掃描的時間 (Block 0 的 TCHED) 對掃描序號 RISCT 做穩健的直線擬合 (排除損壞的時間)
        -> (第一條可見光線 first_line 的時間, 每條可見光線的分鐘數)；資料不足時 None。"""
        r = self.block0['RISCT'].astype(np.int64)
        ok = (r >= 1) & (r <= 1974)
        lines = np.r_[self.first_line, self.first_line + 8 * (r[ok] - 1)]
        times = [None] + [bcd_time(x) for x in self.block0['TCHED'][ok]]
        return fit_line_time(lines, times)

    # ---- 各頻道的影像結構 -----------------------------------------------------
    def channel_layout(self, number):
        """頻道 -> (紀錄索引 (每列一筆，沒有為 -1), 每列可見光行數, 列數, 行數)。"""
        r = self.records
        sel = np.flatnonzero((r['LICHA'] == number) & (self.base_line >= 0))
        if sel.size == 0:
            return None
        if number == 1:
            lines = self.base_line[sel] + (r['block'][sel].astype(np.int64) - 3)
            step, rows_per_scan = 1, 8
        else:
            # 同一掃描中同頻道的紀錄數 = 偵測器數 (2 -> 4 km、1 -> 8 km)；偵測器依 LIDET 由北到南
            per_scan = collections.Counter(collections.Counter(zip(r['scan'][sel], self.base_line[sel])).values())
            ndet = per_scan.most_common(1)[0][0]
            step, rows_per_scan = 8 // ndet, ndet
            order = np.zeros(sel.size, np.int64)
            # 同一 base_line 內依 LIDET 排序取得偵測器順序
            idx = np.lexsort((r['LIDET'][sel], self.base_line[sel]))
            prev, k = None, 0
            for i in idx:
                b = self.base_line[sel][i]
                k = k + 1 if b == prev else 0
                prev = b
                order[i] = k
            flipped = bool(_mode([int(x) >> (31 - 16) & 1 for x in self.block0['ISCAN']]))
            if flipped:
                order = rows_per_scan - 1 - order
            lines = self.base_line[sel] + step * order
        row = (lines - self.first_line) // step
        good = (row >= 0) & (row < 20000)
        sel, row = sel[good], row[good]
        nrows = int(row.max()) + 1
        rec_of_row = np.full(nrows, -1, np.int64)
        rec_of_row[row[::-1]] = sel[::-1]           # 重複的列取第一筆
        ncols = int(np.percentile(r['LPIXLS'][sel], 99))
        return rec_of_row, step, nrows, ncols

    def read_words(self, rec, c0, c1):
        """一筆紀錄的像素 c0..c1 (10 位元計數值)。"""
        r = self.records[rec]
        c1 = min(c1, int(r['LPIXLS']))
        if c1 <= c0:
            return np.zeros(0, np.uint16)
        w0, w1 = int(r['word']) + c0, int(r['word']) + c1
        g0, g1 = w0 // 4, -(-w1 // 4)
        with open(self.path, 'rb') as f:
            f.seek(int(r['offset']) + 5 * g0)
            words = unpack10(f.read(5 * (g1 - g0)))
        return words[w0 - 4 * g0:w1 - 4 * g0]

    def fill_mask(self, rec_of_row, row, c0, words):
        """第 row 列 (像素 c0 起的計數值 words) 中地面站填充資料的位置 (bool)；沒有時 None。
        先看沿線跳動 (不必多讀)，夠大才讀上下兩列比較。"""
        n, w = words.size, _FILL_W
        if n < w:
            return None
        x = words.astype(np.int32)
        cs = np.r_[0, np.cumsum(np.abs(np.diff(x)))]
        jump = (cs[w - 1:] - cs[:n - w + 1]) >= _FILL_JUMP * (w - 1)      # 視窗 [k, k+w) 內 w-1 個相鄰差
        if not jump.any():
            return None
        lo = np.full(n, np.iinfo(np.int32).max, np.int32)
        hi = np.full(n, np.iinfo(np.int32).min, np.int32)
        for j in (row - 1, row + 1):
            if 0 <= j < len(rec_of_row) and rec_of_row[j] >= 0:
                ref = self.read_words(int(rec_of_row[j]), c0, c0 + n).astype(np.int32)
                lo[:ref.size] = np.minimum(lo[:ref.size], ref)
                hi[:ref.size] = np.maximum(hi[:ref.size], ref)
        known = lo <= hi
        if not known.any():
            return None
        out = known & ((x < lo - _FILL_TOL) | (x > hi + _FILL_TOL))
        co = np.r_[0, np.cumsum(out)]
        bad = jump & (co[w:] - co[:n - w + 1] >= _FILL_OUT)
        if not bad.any():
            return None
        mask = np.convolve(bad.astype(np.int32), np.ones(w, np.int32))[:n] > 0
        # 填充段的開頭 / 結尾常是較稀疏的亂數：往兩側延伸到最後一個異常像素 (中間最多隔 3 個正常像素)
        edges = np.flatnonzero(np.diff(np.r_[0, mask.view(np.int8), 0]))
        for s, e in zip(edges[::2], edges[1::2]):
            k, last = s - 1, s
            while k >= 0 and last - k <= 4:
                last = k if out[k] else last
                k -= 1
            k, end = e, e
            while k < n and k - end < 4:
                end = k + 1 if out[k] else end
                k += 1
            mask[last:end] = True
        return mask


def peek(path, slot):
    """只讀檔案開頭的幾個 Block 0 (資料列表用，不必建立整個檔案的索引) -> (Block 0, 後段) 的眾數紀錄。"""
    size = os.path.getsize(path)
    found = []
    with open(path, 'rb') as f:
        pos = GvarFile._first_offset(GvarFile.__new__(GvarFile), f) if slot else 0
        for _ in range(80):
            if pos + _HEADER > size or len(found) >= 7:
                break
            f.seek(pos)
            info = GvarFile._header_ok(_vote(f.read(_HEADER)))
            if info is None:
                if not slot:
                    break
                pos += slot
                continue
            bid, ws, n, _ = info
            if bid == 0 and ws == 8:
                b = f.read(_B0)
                if len(b) == _B0:
                    found.append(b)
            pos = pos + slot if slot else pos + _HEADER + n
    if not found:
        return None
    b = collections.Counter(x[294:1622] for x in found).most_common(1)[0][0]
    raw = next(x for x in found if x[294:1622] == b)
    b0 = np.frombuffer(raw[:fmt._Block_0.itemsize], fmt._Block_0)[0].copy()
    # 影像框 (INFLN、ISFLN、IWFPX、IEFPX) 也取眾數，避免第一個 Block 0 損壞
    frames = [tuple(int(np.frombuffer(x[:fmt._Block_0.itemsize], fmt._Block_0)[0][k])
                    for k in ('INFLN', 'ISFLN', 'IWFPX', 'IEFPX')) for x in found]
    for k, v in zip(('INFLN', 'ISFLN', 'IWFPX', 'IEFPX'), _mode(frames)):
        b0[k] = v
    tail = fmt._Block_0_O if int(b0['SPCID']) >= 14 else fmt._Block_0_I_N
    return b0, np.frombuffer(raw[fmt._Block_0.itemsize:fmt._Block_0.itemsize + tail.itemsize], tail)[0]


class GoesGvarReader(GeoReader):
    NAME = 'GOES-GVAR'
    SENSOR = 'IMAGER'
    DESCRIPTION = 'GOES-8~15 Imager GVAR (NOAA CLASS goes??.yyyy.ddd.hhmmss / NSMC GOS*_GVR_*.dat)'
    PATTERNS = ('goes??.????.???.??????', 'GOS*_GVR_*')
    CHANNELS = {'CH%d' % n: make_channel(n, 'GOES-8' if n != 6 else 'GOES-13', 4) for n in range(1, 7)}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._files = {}
        self._layouts = {}
        self._fill_seen = set()

    # ---- 檔案 --------------------------------------------------------------
    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m, n = _CLASS.match(name), _NSMC.match(name)
            if m:
                sat = 'GOES-%d' % int(m.group(1))
                t = datetime.datetime.strptime(m.group(2) + m.group(3) + m.group(4), '%Y%j%H%M%S')
                slot = None
            elif n:
                sat = 'GOES-%d' % int(n.group(1))
                t = datetime.datetime.strptime(n.group(2) + n.group(3), '%Y%m%d%H%M')
                slot = _SLOT
            else:
                continue
            path = os.path.join(folder, name)
            try:
                first = peek(path, slot)
            except OSError as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                continue
            if first is None:
                print(tr('  [略過] %s：找不到 GVAR Block 0', '  [skip] %s: GVAR Block 0 not found') % name)
                continue
            b0, b0b = first
            # 影像框 (可見光行、像素)：涵蓋全圓盤 (> 9000 條線) 為 FD，否則為區域 (地圖範圍自動帶入)
            frame = (int(b0['INFLN']), int(b0['ISFLN']), int(b0['IWFPX']), int(b0['IEFPX']))
            full = frame[1] - frame[0] > 9000 or frame[3] - frame[2] > 18000
            s = Scene(sat, t, 'FD' if full else 'SECTOR', files={'path': path},
                      extra={'slot': slot, 'frame': frame, 'first': first},
                      format='NSMC GVAR' if slot else 'GVAR')
            if not full:
                try:
                    s.extra['georange'] = self._image_range(s)
                    s.region, s.auto_range = sector_label(s.extra['georange'])
                except (SystemExit, ValueError):
                    pass
            out.append(s)
        return sorted(out, key=lambda s: (s.time, s.satellite))

    def _image_range(self, scene):
        """區域影像的經緯度範圍：用第一個 Block 0 的影像框與導航參數計算 (不必讀整個檔案)。"""
        n0, n1, w0, w1 = scene.extra['frame']
        p = params_from_block0(*scene.extra['first'])
        nav = GoesNav((n1 - n0 + 1, w1 - w0 + 1), p, n0, w0)
        ny, nx = nav.shape
        rows, cols = np.meshgrid(np.linspace(0, ny - 1, 41), np.linspace(0, nx - 1, 41), indexing='ij')
        lon, lat = nav.to_lonlat(rows, cols)
        ok = np.isfinite(lon)
        if not ok.any():
            raise SystemExit(tr('無法計算 %s 的影像範圍', 'Cannot determine the image range of %s') % scene.area)
        lon = lon[ok]
        if lon.max() - lon.min() > 180:
            lon = lon % 360.
        return float(lat[ok].min()), float(lat[ok].max()), float(lon.min()), float(lon.max())

    def _gvar(self, scene):
        path = scene.files['path']
        if path not in self._files:
            g = GvarFile(path, scene.extra.get('slot'))
            if g.bad_headers or getattr(g, 'repaired', 0):
                print(tr('  [GVAR] %s：略過 %d 個無法辨識的區塊頭、修正 %d 筆 Line Documentation',
                         '  [GVAR] %s: skipped %d unrecognized block headers, repaired %d Line Documentation records') % (
                    os.path.basename(path), g.bad_headers, getattr(g, 'repaired', 0)))
            self._files[path] = g
        return self._files[path]

    def _layout(self, scene, number):
        key = (scene.files['path'], number)
        if key not in self._layouts:
            self._layouts[key] = self._gvar(scene).channel_layout(number)
        if self._layouts[key] is None:
            raise ValueError(tr('%s 沒有頻道 %d 的資料', '%s has no data for channel %d')
                             % (os.path.basename(scene.files['path']), number))
        return self._layouts[key]

    def channels(self, scene):
        """各衛星的頻道 (GOES-8~11：1~5；GOES-12~15：1、2、3、4、6)；8 km 頻道：GOES-8~11 的 3、GOES-12/13 的 6。"""
        n = int(scene.satellite.split('-')[1])
        numbers = (1, 2, 3, 4, 5) if n <= 11 else (1, 2, 3, 4, 6)
        coarse = {3} if n <= 11 else ({6} if n <= 13 else set())
        return [make_channel(c, scene.satellite, 8 if c in coarse else 4) for c in numbers]

    def header(self, scene, channel=None):
        g = self._gvar(scene)
        return {'Block_0': g.block0, 'Block_0_b': g.block0b, 'records': g.records, 'nav': self._params(scene)}

    # ---- 導航 --------------------------------------------------------------
    def _params(self, scene):
        if 'nav' not in scene.extra:
            g = self._gvar(scene)
            # 各 Block 0 的 O&A 參數 (bytes 294~1622：參考經度 REFLO ~ 俯仰失準級數) 取眾數，個別損壞的 Block 0 不影響
            raw = [bytes(x) for x in np.frombuffer(g.block0.tobytes(), np.dtype((np.void, fmt._Block_0.itemsize)))]
            oa = slice(294, 1622)
            best = collections.Counter(r[oa] for r in raw if r[322:330] != b'\x00' * 8).most_common(1)
            k = next(i for i, r in enumerate(raw) if best and r[oa] == best[0][0]) if best else 0
            p = params_from_block0(g.block0[k], g.block0b[k])
            p['nad'] = _mode([tuple(int(x[n]) for n in ('IOFNC', 'IOFEC', 'IOFNI', 'IOFEI')) for x in g.block0b])
            p['iscan_mode'] = _mode([int(x) for x in g.block0['ISCAN']])
            p['imc'] = bool(p['iscan_mode'] >> (31 - 8) & 1)
            p['flip'] = -1. if p['iscan_mode'] >> (31 - 16) & 1 else 1.
            p['time'] = _mode([bcd_time(t) for t in g.block0['TCHED']][:50]) or p['time']
            scene.extra['nav'] = p
        return scene.extra['nav']

    def nav(self, scene, channel):
        number = int(channel[2:])
        g = self._gvar(scene)
        rec_of_row, step, nrows, ncols = self._layout(scene, number)
        p = self._params(scene)
        if number == 1:
            first_line, first_pixel, pstep = g.first_line, g.west_pixel, 1
        else:
            first_line = g.first_line + IR_LINE_CENTER * step / 4.
            first_pixel, pstep = g.west_pixel + IR_PIXEL_CENTER, 4
        if 'line_time' not in scene.extra:
            lt = g.line_time()
            # line_time 以 first_line 為起點；IMC 關閉時導航依各線的掃描時間計算
            scene.extra['line_time'] = lt
        lt = scene.extra['line_time']
        if lt is not None:
            lt = (lt[0] + datetime.timedelta(minutes=lt[1] * (first_line - g.first_line)), lt[1])
        return GoesNav((nrows, ncols), p, first_line, first_pixel, step, pstep,
                       key=(scene.files['path'], number), engine=str(self.cfg.navigation.engine).lower(), line_time=lt)

    def satellite_position(self, scene):
        st = self.nav(scene, 'CH4' if 'CH4' in [c.name for c in self.channels(scene)] else 'CH1').state
        return st['sublon'], st['height']

    # ---- 讀取與校正 --------------------------------------------------------
    def _calibrators(self, scene, number):
        key = ('cal', scene.files['path'], number)
        if key in self._layouts:
            return self._layouts[key]
        g = self._gvar(scene)
        if number == 1:
            b = g.block0b
            if 'IVCRB' in b.dtype.names:
                k = len(b) // 2
                coef = vis_coefficients(g.satellite, gould(np.frombuffer(bytes(b['IVCRB'][k]), '>u4')),
                                        gould(np.frombuffer(bytes(b['IVCR1'][k]), '>u4')),
                                        gould(np.frombuffer(bytes(b['IVCR2'][k]), '>u4')), gould(int(b['IVRAL'][k])))
            else:
                coef = vis_coefficients(g.satellite, [0] * 8, [0] * 8, [0] * 8, 0)
            bias, g1, g2, alb = coef
            cal = ('vis', (bias, g1, g2, alb * vis_scale(self.cfg, g.satellite, scene.time)))
        else:
            cal = ('ir', {})
        self._layouts[key] = cal
        return cal

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        number = int(channel[2:])
        g = self._gvar(scene)
        rec_of_row = self._layout(scene, number)[0]
        kind, cal = self._calibrators(scene, number)
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        table_dir = self.cfg.general.conversion_table_dir
        mask_fill = bool(self.cfg.get('gvar', {}).get('mask_fill', False))
        pad = _FILL_W if mask_fill else 0           # 填充資料的判定多讀兩側各一個視窗，視窗邊緣也判得準
        a0 = max(c0 - pad, 0)
        for i, rec in enumerate(rec_of_row[r0:r1]):
            if rec < 0:
                continue
            wide = g.read_words(int(rec), a0, c1 + pad)
            counts = wide[c0 - a0:c1 - a0].astype(np.int64)
            if counts.size == 0:
                continue
            fill = g.fill_mask(rec_of_row, r0 + i, a0, wide) if mask_fill else None
            if fill is not None:
                fill = fill[c0 - a0:c0 - a0 + counts.size]
                self._report_fill(scene, number, r0 + i, int(fill.sum()))
            r = g.records[rec]
            if kind == 'vis':
                bias, g1, g2, alb = cal
                d = max(0, min(7, int(r['LIDET']) - 1))
                vals = (bias[d] + g1[d] * counts + g2[d] * counts * counts) * alb
                out[i, :counts.size] = np.where(counts > 0, vals, np.nan)      # 計數值 0 = 填充
            else:
                side = 2 if int(r['LSIDE']) == 1023 else 1
                det = 2 if int(r['LIDET']) % 2 == 0 else 1
                key = (side, det)
                if key not in cal:
                    cal[key] = ir_lut(table_dir, g.satellite, side, number, det)
                lut = cal[key]
                if lut is not None:
                    out[i, :counts.size] = lut[np.clip(counts, 0, 1023)]
            if fill is not None:
                out[i, :counts.size][fill] = np.nan
        return out

    def _report_fill(self, scene, number, row, npix):
        """填充資料每個 (檔案、頻道、列) 只報告一次。"""
        key = (scene.files['path'], number, row)
        if key in self._fill_seen or not npix:
            return
        self._fill_seen.add(key)
        print(tr('  [GVAR] %s 頻道 %d 第 %d 列：%d 個像素是地面站填充的亂數，標為無資料',
                 '  [GVAR] %s channel %d row %d: %d pixels are ground-station fill (random data), set to no data') % (
            os.path.basename(scene.files['path']), number, row + 1, npix))


READER = GoesGvarReader
