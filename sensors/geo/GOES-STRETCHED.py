"""GVAR 以前的 GOES 與 SMS (SMS-1/2、GOES-1~7) 的原始 stretched 資料：Mode A (VISSR)、Mode AA / AAA (VAS)
檔名 (NCEI 2024 重處理版，NOAA CLASS)：goesNN.YYYY.DDD.HHMM.A.data、.AA-MSI.data、.AA-DS.data、.AAA-MSI.data、
  .AAA-DS.data、smsNN.YYYY.DDD.HHMM.A.data。一個檔 = 一筆資料 (VIS 與 IR 在同一個檔)。
格式 (記錄配置、文件欄位、線號) 見 sensors/stretched_format.py。
波段：
  Mode A     VIS (6 位元，15288 像素) 與 IR (8 位元已校正計數值，3818 個有效樣本)，每圈一條 IR 線。
  Mode AA    VIS 與 VAS 的 BAND-01~12：每圈的濾鏡輪位置決定波段 (MSI 例如 8、7、8、7、8、10、8、10 輪替；
  Mode AAA   DS (dwell sounding) 在部分掃描步駐留、依序取 12 個波段)，上下兩個偵測器各一條線。
             同一波段、同一步有多圈 (DS 駐留) 時輻射取平均 (各圈用自己的校正係數；SSEC AREA 則平均計數值)。
             VAS 的 VISSR 模式 (.AAA，沒有 -MSI / -DS) 只有 BAND-08，只用上偵測器 (IR1)。
校正：
  VIS        反射率 = (計數值 / 63)^2 (同 GOES-PREGVAR：平方根編碼、無絕對校正)。Mode AA 先以檔內的可見光非線性係數
             反算地面處理器的換算 (GOES-5 1982 斜率約 0.5；係數為恆等時不變；見 stretched_format.vis_inverse)。
  Mode A IR  NCC 1978 附錄 II：計數值 < 176 為 330 - c/2，否則 418 - c (K)。
  Mode AA    輻射 = (32 x V - YZQ) x 2^(F + DELTAF - 15)；YZQ = 共同文件第 3 段的 Y_Z + Q、DELTAF 由 DELTA SCALE FACTOR 查表，
             逐圈、上下偵測器各一組 (同 McIDAS AREA 的 YSUBZ / RAWDELTAF；AREA 的 P = 32 x V)。
  Mode AAA   輻射 = (AB_k(1) x V - AB_k(0)) / 2^(15 - FAB_k)；AB、FAB 由第 3 欄第 63~112 字 (20 圈輪替) 組成，
             k = 偵測器 / 濾鏡組合 (1~38)。
  輻射 -> 亮溫 同 GOES-PREGVAR 的 VAS 係數。
導航：McIDAS 'GOES' 導航 (mcidas_goes_nav.McGoesNav)，參數取自每圈文件中輪替的 O&A (20 圈一組)：
  VIRGS 的 Kepler 軌道、自旋軸、ZETA / ETA 姿態 (見 stretched_format.nav_params)。沒有 Kepler 元素的檔案 (例如 1977/11 以前的
  舊 O&A 格式) 列出原因後略過。IR 的位置照資料本身的掃描編號、比照 McIDAS 慣例 (同 GOES-PREGVAR：IR 線中心 = 該步第一條
  VIS 線 + 3.5)，不另外修正配準。
同名的 .indx (區塊索引)、.meta (JSON) 不是必要的：有 .indx 時可讀到不完整的記錄，.meta 只放進 header()。
資料列表只讀檔頭、檔尾的記錄 (PROBE；時間、波段、範圍)，畫圖時才完整解析 (傳統硬碟冷啟動每檔 3~5 秒 -> 約 0.3 秒)。
完整檔頭：reader.header(scene) 回傳 {'mode', 'records' (各記錄的文件欄位), 'offsets' / 'ir_pos' / 'vis_pos' (區塊位置)、
  'index_file'、'meta' (.meta 的內容)、'oa' (O&A)、'nav' (導航參數)、'calibration'}。
"""
import datetime
import json
import os
import re

import numpy as np

import stretched_format as sf
from goes_imager import sector_label
from mcidas_goes_nav import McGoesNav
from satima.core.i18n import tr
from satima.geo.reader import Channel, GeoReader, Scene

_NAME = re.compile(r'^(goes|sms)(\d\d)\.(\d{4})\.(\d{3})\.(\d{4})(?:-(\d+))?\.(AAA|AA|A)(?:-([A-Za-z]+))?\.data$', re.I)

IR_LINE_OFFSET = 3.5            # IR 線中心 = 該步第一條 VIS 線 + 3.5 (8 條 VIS 線的中央，同 GOES-PREGVAR)
IR_ELE_OFFSET = 0.
MIN_SPINS = 4                   # 一個波段至少要有幾圈資料才列入選單 (排除校正、換模式時的零星幾圈)
# 資料列表只讀檔頭、檔尾的記錄 (前 n 筆、後 m 筆)：VAS 的波段在檔頭就會全部出現 (樣本 MSI 前 14 筆、DS 前 57 筆)；
# Mode A 的波段固定，只要時間與 O&A (20 圈一組)
PROBE = {'A': (40, 16), 'AA': (160, 16), 'AAA': (160, 16)}
FILL_STEPS = 2                  # VAS 波段：沒有資料的步以最近的線補 (最多相隔幾步；DS 的小偵測器波段只涵蓋一半的步，
                                # 同 McIDAS AREA 以最近一圈的資料重複填線)


def _fill_rows(img, has, limit):
    """沒有資料的列 (has = False) 以最近 (相隔 <= limit 列) 的有資料列補上；距離相同時取上方 (北) 的列。"""
    rows = np.flatnonzero(has)
    if rows.size == 0 or rows.size == len(has):
        return
    for i in np.flatnonzero(~has):
        j = np.searchsorted(rows, i)
        cand = [rows[k] for k in (j - 1, j) if 0 <= k < rows.size]
        best = min(cand, key=lambda r: (abs(r - i), r > i))
        if abs(best - i) <= limit:
            img[i] = img[best]


class _Info:
    """一個檔的影像配置：VIS 各線來自哪一筆記錄、哪一段；IR 各波段來自哪些記錄、槽位、掃描步。"""

    def __init__(self, f):
        m = f.meta
        lay = f.layout
        frame = np.flatnonzero(m['frame'] & (m['scan'] >= 1) & (m['scan'] <= 1821))
        if frame.size >= 3:                 # 掃描數的位元錯誤：與前後 9 筆的中位數相差超過 20 的記錄不用
            s = m['scan'][frame].astype(float)
            med = np.median(np.lib.stride_tricks.sliding_window_view(np.pad(s, 4, mode='edge'), 9), axis=1)
            frame = frame[np.abs(s - med) <= 20]
        if frame.size == 0:
            raise ValueError(tr('沒有影像資料 (frame code)', 'no image data (frame code)'))
        self.frame = frame
        ref = frame[0]
        self.ref = ref
        self.yyddd, self.seconds, self.scan0 = int(m['yyddd'][ref]), float(m['seconds'][ref]), int(m['scan'][ref])
        blocks = sf.assemble(m['blk'], m['mfi'], m['oa'])
        self.oa = sf.oa_fields(blocks)
        self.params = sf.nav_params(self.oa, self.yyddd, self.seconds, self.scan0)
        # VIS 線號
        base = {'A': -3, 'AA': 5, 'AAA': -7}[lay.mode]              # 掃描 s 段 0 的線號 = 8 s + base
        vf = frame
        if lay.mode == 'AA':
            # AA DS 駐留 (submode 2) 的記錄掃描數已加 1，可見光卻仍是前一圈的線 (樣本 GOES-5 1982、GOES-6 1985 每段駐留
            # 除最後一筆 (步進) 外皆如此)：不用。AAA DS 駐留的可見光在本圈，照常使用。
            vf = frame[m['submode'][frame] != 2]
        have = f.vis_pos[vf] >= 0                                    # 不完整的記錄可能缺某幾段
        if not have.any():
            raise ValueError(tr('沒有可見光資料', 'no visible data'))
        lines = ((8 * m['scan'][vf] + base)[:, None] + np.arange(8))[have]
        recs = np.repeat(vf[:, None], 8, axis=1)[have]
        secs = np.repeat(np.arange(8)[None, :], len(vf), axis=0)[have]
        first, last = int(lines.min()), int(lines.max())
        self.vis_first = first
        self.vis_rec = np.full(last - first + 1, -1, np.int64)
        self.vis_sec = np.zeros(last - first + 1, np.int64)
        # 同一線出現多次：Mode A 取第一圈、AA / AAA 取最後一圈 (畫格末端重複的記錄；兩者皆與 SSEC AREA 逐像素相同)
        if lay.mode == 'A':
            uniq, idx = np.unique(lines, return_index=True)
        else:
            uniq, idx = np.unique(lines[::-1], return_index=True)
            idx = len(lines) - 1 - idx
        self.vis_rec[uniq - first] = recs[idx]
        self.vis_sec[uniq - first] = secs[idx]
        # IR：各波段的 (記錄, 槽位, 第一步, 步數)
        self.bands = {}
        if lay.mode == 'A':
            self.bands[0] = [(int(r), 0, int(m['scan'][r]), 1) for r in frame]
        else:
            # VAS 的 VISSR 模式 (第 24 字 bit 9 = 0；GOES-6 1987/183 的 .AAA)：IR1、IR2 都是 band 8 小偵測器，IR2 的
            # adjusted earth count 比 IR1 多 1，影像卻與南 4 圈的 IR1 相符 -> 只用 IR1 (每圈一條，同 Mode A 與 SSEC AREA)
            vissr = np.mean((m['vasmode'][frame] & 2) == 0) > .5
            for slot in range(1 if vissr else len(lay.ir_slots)):
                for r in frame:
                    if not m['present'][r, slot] or not m['slot_frame'][r, slot]:
                        continue
                    self.bands.setdefault(int(m['band'][r, slot]), []).append(
                        (int(r), slot, int(m['step'][r, slot]), int(m['nsteps'][r, slot])))
            self.bands = {b: v for b, v in sorted(self.bands.items()) if len(v) >= MIN_SPINS}
        steps = [t + i for v in self.bands.values() for _, _, t, n in v for i in range(n)]
        if not steps:
            raise ValueError(tr('沒有紅外線資料', 'no infrared data'))
        self.ir_first, self.ir_last = min(steps), max(steps)
        self.cal = None
        if lay.mode == 'AAA':
            self.cal = sf.aaa_calibration(m['blk'], m['mfi'], m['ircal'])
        # Mode AA 可見光：反算地面處理器的非線性換算 (計數值 -> 原始計數值；係數為恆等時 None)
        self.vis_lut = sf.vis_inverse(m) if lay.mode == 'AA' else None

    @property
    def time(self):
        y, d = self.yyddd // 1000, self.yyddd % 1000
        year = 1900 + y if y >= 50 else 2000 + y
        return datetime.datetime(year, 1, 1) + datetime.timedelta(days=d - 1, seconds=self.seconds)


class GoesStretchedReader(GeoReader):
    NAME = 'GOES-STRETCHED'
    SENSOR = 'VISSR'
    DESCRIPTION = 'SMS-1/2、GOES-1~7 Mode A / AA / AAA 原始資料 (goesNN.yyyy.ddd.hhmm.A.data、.AA-MSI.data、.AAA-DS.data…)'
    DESCRIPTION_EN = 'SMS-1/2, GOES-1~7 Mode A / AA / AAA raw stretched data (goesNN.yyyy.ddd.hhmm.A.data, .AA-MSI.data, ...)'
    PATTERNS = ('goes0?.????.???.????*.A*.data', 'sms0?.????.???.????*.A*.data')
    CHANNELS = {}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._files = {}
        self._infos = {}
        self._navs = {}
        self._images = {}

    # ---- 檔案 --------------------------------------------------------------
    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            path = os.path.join(folder, name)
            mode = m.group(7).upper()
            try:
                info = _Info(sf.RawFile(path, mode, probe=PROBE[mode]))       # 列表只讀檔頭、檔尾；畫圖時才完整解析
            except (OSError, ValueError) as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                continue
            n = int(m.group(2))
            sat = ('GOES-%d' if m.group(1).lower() == 'goes' else 'SMS-%d') % n
            vas = m.group(1).lower() == 'goes' and n >= 4
            sub = (m.group(8) or '').upper()
            s = Scene(sat, info.time.replace(microsecond=0), 'FD', format='MODE ' + mode + ('-' + sub if sub else ''),
                      sensor='VAS' if vas else 'VISSR')
            s.files['VIS'] = path
            if mode == 'A':
                s.files['IR'] = path
            else:
                for b in info.bands:
                    s.files['B%02d' % b] = path
            s.extra.update(mode=mode, path=path)
            if len(info.vis_rec) <= 9000:                    # 涵蓋大部分地球 (> 9000 條 VIS 線) 為 FD，否則為區域
                s.area = 'SECTOR'
                try:
                    s.extra['georange'] = self.nav_range(self._make_nav(info, 'VIS', (path, 'probe')), s.area)
                    s.region, s.auto_range = sector_label(s.extra['georange'])
                except (SystemExit, ValueError):
                    pass
            out.append(s)
        return sorted(out, key=lambda s: (s.time, s.satellite))

    def _file(self, path, mode=None):
        if path not in self._files:
            if mode is None:
                mode = _NAME.match(os.path.basename(path)).group(7).upper()
            self._files[path] = sf.RawFile(path, mode)
        return self._files[path]

    def _info(self, path, mode=None):
        if path not in self._infos:
            f = self._file(path, mode)
            if len(f.offsets) == 0:
                raise ValueError(tr('找不到 Mode %s 的記錄', 'no Mode %s records found') % f.mode)
            self._infos[path] = _Info(f)
        return self._infos[path]

    def channels(self, scene):
        out = []
        for name in scene.files:
            if name == 'VIS':
                out.append(Channel('VIS', 0.65, 1., 'vis', band='VIS', alias='Visible'))
            elif name == 'IR':
                wl = 11.17 if scene.sensor == 'VAS' else 11.5
                out.append(Channel('IR', wl, 4., 'ir', band='IR', alias='LW Window'))
            else:
                b = int(name[1:])
                wl, kind, alias = sf.VAS_BANDS[b]
                out.append(Channel(name, wl, 4., kind, band='BAND-%02d' % b, alias=alias))
        return out

    def header(self, scene, channel=None):
        f = self._file(scene.extra['path'])
        info = self._info(scene.extra['path'])
        meta = None
        path = os.path.splitext(f.path)[0] + '.meta'
        if os.path.exists(path):
            try:
                with open(path, encoding='utf-8') as fp:
                    meta = json.load(fp)
            except (OSError, ValueError):
                pass
        return {'mode': f.mode, 'records': f.meta, 'offsets': f.offsets, 'ir_pos': f.ir_pos, 'vis_pos': f.vis_pos,
                'ir_valid': f.ir_valid, 'vis_valid': f.vis_valid, 'index_file': f.index_file, 'meta': meta, 'oa': info.oa,
                'nav': info.params, 'calibration': info.cal, 'vis_inverse': info.vis_lut}

    # ---- 導航 --------------------------------------------------------------
    def nav(self, scene, channel):
        path = scene.extra['path']
        kind = 'VIS' if channel == 'VIS' else 'IR'
        key = (path, kind)
        if key not in self._navs:
            self._navs[key] = self._make_nav(self._info(path), kind, key)
        return self._navs[key]

    def _make_nav(self, info, kind, key):
        engine = str(self.cfg.navigation.engine).lower()
        if kind == 'VIS':
            return McGoesNav((len(info.vis_rec), sf.NVIS), info.params, info.vis_first, 1., 1., 1., key=key, engine=engine)
        first = (info.ir_first - 1) * 8 + 1 + IR_LINE_OFFSET
        return McGoesNav((info.ir_last - info.ir_first + 1, sf.NIR), info.params, first, 1. + IR_ELE_OFFSET, 8., 4.,
                         key=key, engine=engine)

    def satellite_position(self, scene):
        lon, lat, height = self.nav(scene, 'VIS').subpoint()
        return lon, height

    # ---- 讀取與校正 --------------------------------------------------------
    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        path = scene.extra['path']
        if channel != 'VIS':
            img = self._ir_image(path, channel)
            out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
            rr0, rr1 = max(r0, 0), min(r1, img.shape[0])
            cc0, cc1 = max(c0, 0), min(c1, img.shape[1])
            if rr1 > rr0 and cc1 > cc0:
                out[rr0 - r0:rr1 - r0, cc0 - c0:cc1 - c0] = img[rr0:rr1, cc0:cc1]
            return out
        f, info = self._file(path), self._info(path)
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        rr0, rr1 = max(r0, 0), min(r1, len(info.vis_rec))
        cc0, cc1 = max(c0, 0), min(c1, sf.NVIS)
        if rr1 <= rr0 or cc1 <= cc0:
            return out
        rows = np.arange(rr0, rr1)
        rows = rows[info.vis_rec[rows] >= 0]
        for i in range(0, len(rows), 256):                       # 分批解碼 (大範圍時控制記憶體)
            part = rows[i:i + 256]
            v = f.vis_lines(info.vis_rec[part], info.vis_sec[part], cc0, cc1)
            if info.vis_lut is not None:
                v = info.vis_lut[v]
            vals = (v / 63.) ** 2
            nv = f.vis_valid[info.vis_rec[part], info.vis_sec[part]]
            vals[np.arange(cc0, cc1)[None, :] >= nv[:, None]] = np.nan     # 截短區塊的補零部分
            out[part - r0, cc0 - c0:cc1 - c0] = vals
        return out

    def _ir_image(self, path, channel):
        key = (path, channel)
        if key in self._images:
            return self._images[key]
        f, info = self._file(path), self._info(path)
        m = f.meta
        nrow = info.ir_last - info.ir_first + 1
        if channel == 'IR':
            img = np.full((nrow, sf.NIR), np.nan, dtype=np.float32)
            recs = [r for r, _, _, _ in info.bands[0]]
            counts = f.ir_video(recs, 0)
            temp = sf.mode_a_temperature(counts)
            temp[:, f.layout.ir_valid:] = np.nan
            temp[np.arange(sf.NIR)[None, :] >= f.ir_valid[recs, 0][:, None]] = np.nan     # 截短區塊的補零部分
            temp[~counts.any(axis=1)] = np.nan                      # 整行為 0
            for (r, _, t, _), row in zip(info.bands[0][::-1], temp[::-1]):
                img[t - info.ir_first] = row                        # 同一步有多圈時取第一圈
        else:
            band = int(channel[1:])
            total = np.zeros((nrow, sf.NIR), np.float64)
            count = np.zeros((nrow, sf.NIR), np.int32)                  # 逐像素的圈數 (無效像素不計)
            for slot in (0, 1):
                items = [it for it in info.bands.get(band, ()) if it[1] == slot]    # 列表 (只讀檔頭) 有、整檔不足的波段為 NaN
                if not items:
                    continue
                recs = np.array([it[0] for it in items])
                video = f.ir_video(recs, slot).astype(np.float64)
                rad = self._radiance(f, info, recs, slot, band, video)
                rad[np.arange(sf.NIR)[None, :] >= f.ir_valid[recs, slot][:, None]] = np.nan     # 截短區塊的補零部分
                for (r, _, t, n), row in zip(items, rad):
                    if not np.isfinite(row).any():
                        continue
                    ok = np.isfinite(row)
                    for i in range(n):
                        j = t + i - info.ir_first
                        total[j] += np.where(ok, row, 0.)
                        count[j] += ok
            with np.errstate(invalid='ignore'):
                rad = np.where(count > 0, total / np.maximum(count, 1), np.nan)
            img = sf.brightness_temperature(rad, band).astype(np.float32)
            _fill_rows(img, count.any(axis=1), FILL_STEPS)
        self._images[key] = img
        while len(self._images) > 6:
            self._images.pop(next(iter(self._images)))
        return img

    @staticmethod
    def _radiance(f, info, recs, slot, band, video):
        """各記錄 (recs) 槽位 slot 的影像 -> 輻射 (n, 3822)；無法校正的行為 NaN。"""
        m = f.meta
        if f.mode == 'AA':
            yzq = m['yzq'][recs, slot].astype(np.float64)[:, None]
            dsf = m['dsf'][recs, slot] & 15
            scale = 2. ** (sf.F_AA[band - 1] + np.array([sf.DELTAF[d] for d in dsf]) - 15.)[:, None]
            rad = (32. * video - yzq) * scale
            rad[yzq[:, 0] <= 0] = np.nan
        else:
            if info.cal is None:
                raise ValueError(tr('Mode AAA 的 IR 校正參數不完整', 'Mode AAA IR calibration parameters incomplete'))
            ab, fab = info.cal
            rad = np.full(video.shape, np.nan)
            for i, r in enumerate(recs):
                k = sf.detector_index(int(m['det'][r, slot]), int(m['sector'][r, slot]) == sf.UPPER, band)
                if k == 0 or ab[k - 1, 1] == 0:
                    continue
                rad[i] = (ab[k - 1, 1] * video[i] - ab[k - 1, 0]) / 2. ** (15 - int(fab[k - 1]))
        rad[video <= 0] = np.nan
        return rad


READER = GoesStretchedReader
