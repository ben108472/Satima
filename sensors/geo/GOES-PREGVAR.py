"""GVAR 以前的 GOES 與 SMS (SMS-1/2、GOES-1~7) — VISSR / VAS 的 McIDAS AREA
檔名 (SSEC / NOAA CLASS)：goesNN.YYYY.DDD.HHMM[-n].AREA_IR / .AREA_VIS、smsNN.YYYY.DDD.HHMM.AREA_IR / .AREA_VIS
  同一衛星、時間 (含 -n 序號) 的可見光與紅外檔合成一筆資料。
格式 (結構與欄位名稱同舊版 AREA_format.py，見 sensors/pregvar_area_format.py)：
  AREA 目錄 64 words (依第 2 個 word = 4 判斷位元組順序)；影像座標以可見光線、像素為單位 (由 1 起算)，
    可見光解析度 1 / 1，紅外 8 (一圈 8 條可見光線) / 4。W34 = 影像位置、W35 = 導航區塊、W63 = 校正區塊。
  影像來源 (W52)：
    VISR  VISSR 模式 (SMS、GOES-1~7)，每像素 1 byte 的亮度值 (校正類型 BRIT)。紅外每行前置資料 = 有效碼 4 + 文件 128 bytes。
    VAS   VAS mode AA (GOES-4~7，1987/03 以前)，每像素 2 bytes、最多 3 個波段交錯；每行前置資料 =
          有效碼 4 + IR 文件 512 + 校正資訊 116 (日期、時間、掃描序號、13 組 [頻道, 自旋數, RAWDELTAF, YSUBZ]) + 波段表。
    AAA   VAS mode AAA (GOES-7)，同上，但校正資訊的 13 組為 [頻道, 自旋數, 0, 0]，係數在校正區塊 (IAB、IFAB)。
  導航區塊 'GOES' (128 words)：Kepler 軌道、自旋軸、掃描幾何、姿態 (sensors/mcidas_goes_nav.py)。
校正：
  可見光 (VISR)：6 位元計數值 x 4 (0~252)。VISSR 可見光為平方根編碼 (同型的 GMS VISSR 反照率表 = (計數值/63)^2，
    rms 0.002)，SMS / GOES 沒有絕對校正，反射率取 (B/252)^2 (舊版直接用 B/255 顯示，等同已做 gamma)。
  紅外 (VISR、BRIT)：McIDAS 亮度值 -> 溫度：B >= 176 時 418 - B，否則 330 - B/2 (K)。
  VAS mode AA：輻射 = (P - YSUBZ) x 2^(F + DELTAF - 15) (F：各波段的放大位數；DELTAF 由 RAWDELTAF 查表)。
  AAA：輻射 = (IAB2 x P/32 - IAB1) / 2^(15 - IFAB)，依該行的頻道 (偵測器) 取係數。
  輻射 -> 亮溫：T = (FK2 / ln(FK1/R + 1) - TC1) / TC2 (McIDAS 的形式；舊版寫成 ((FK2 - TC1)/TC2) / ln(...)，偏暖約 0.3 K)。
    FK1、FK2、TC1、TC2 沿用舊版的一組 VAS 係數 (未註明衛星；GOES-4~7 的光譜響應略有不同，目前無各衛星的係數)。
  計數值 0、整行為 0、波段表中沒有該波段的行 為 NaN。
導航：McIDAS 'GOES' 導航 (mcidas_goes_nav.McGoesNav)；導航區塊沒有軌道參數的檔案 (例如 1979 年的 GOES-1、GOES-3) 無法定位，略過。
完整檔頭：reader.header(scene, 頻道) 回傳 {'Control_Block', 'Navigation', 'Calibration', 'Line_prefix', 'Level_map', 'nav'}。
"""
import datetime
import os
import re

import numpy as np

import pregvar_area_format as fmt
import stretched_format as sf
from goes_imager import sector_label
from mcidas_goes_nav import McGoesNav, params as nav_params
from satima.core.i18n import tr
from satima.geo.reader import Channel, GeoReader, Scene

_NAME = re.compile(r'^(goes|sms)(\d\d)\.(\d{4})\.(\d{3})\.(\d{4})(?:-(\d+))?\.AREA_(IR|VIS)$', re.I)

# 紅外線、像素中心相對於影像座標的偏移 (可見光單位)：紅外每圈取樣一次，AREA 的紅外行號 = 該圈第一條可見光線
# (例如 817 = 8 x 102 + 1)，中心在同一圈 8 條可見光線的中央 (+3.5)，同 GVAR 紅外的 +1.5。
# (GOES-7 同時次可見光 / 紅外的雲邊比對：27 個區域皆為 +3.25~+8，中位數 +5.5；東西向 1988 年 +7、1994 年 -3，
#  不一致，屬當時的配準誤差，不修正。)
IR_LINE_OFFSET = 3.5
IR_ELE_OFFSET = 0.

# VAS 波段 (中心波長、種類、功能代稱)、輻射 -> 亮溫係數 (舊版 goes_pregvar_tools.py)、mode AA 的放大位數 F 與
# RAWDELTAF -> DELTAF 查表：與 GOES-STRETCHED 共用，定義在 sensors/stretched_format.py
_VAS, _F_AA, _DELTAF = sf.VAS_BANDS, sf.F_AA, sf.DELTAF
_brightness_temperature = sf.brightness_temperature


def _dt(dt, order):
    return dt.newbyteorder(order) if order != '=' else dt


class AreaFile:
    def __init__(self, path):
        self.path = path
        with open(path, 'rb') as f:
            raw = f.read(256)
            self.order = '<' if np.frombuffer(raw[4:8], '<i4')[0] == 4 else '>'
            d = np.frombuffer(raw, _dt(fmt._Control_Block, self.order), count=1)[0]
            self.dir = d
            f.seek(int(d['W35']))
            nraw = f.read(fmt._Navigation.itemsize)
            self.nav = np.frombuffer(nraw, _dt(fmt._Navigation, self.order), count=1)[0]
            self.navblock = np.frombuffer(nraw, self.order + 'i4').astype(np.int64)
            cal_off = int(d['W63'])
            f.seek(cal_off)
            craw = f.read(fmt._Calibration_AAA.itemsize) if cal_off > 0 else b''
        self.source = bytes(d['Image source type']).decode('ascii', 'replace').strip()
        self.caltype = bytes(d['Calibration Type']).decode('ascii', 'replace').strip()
        self.cal = (np.frombuffer(craw, _dt(fmt._Calibration_AAA, self.order), count=1)[0]
                    if len(craw) == fmt._Calibration_AAA.itemsize else None)
        self.nlines, self.neles = int(d['NLINES']), int(d['NELES'])
        self.elesiz, self.nchans, self.presiz = int(d['ELESIZ']), int(d['NCHANS']), int(d['PRESIZ'])
        self.doc, self.callen, self.levlen = int(d['W49']), int(d['W50']), int(d['W51'])
        self.rowlen = self.presiz + self.neles * self.elesiz * self.nchans
        fm = int(d['FILTER MAP']) & 0xFFFFFFFF
        self.bands = [i + 1 for i in range(32) if fm >> i & 1]
        self._prefix = None

    def rows(self, r0, r1):
        """第 r0~r1-1 行的原始資料 (含前置資料)，uint8 陣列 (n, rowlen)。"""
        r0, r1 = max(0, r0), min(self.nlines, r1)
        if r1 <= r0:
            return np.zeros((0, self.rowlen), np.uint8)
        with open(self.path, 'rb') as f:
            f.seek(int(self.dir['W34']) + r0 * self.rowlen)
            buf = f.read((r1 - r0) * self.rowlen)
        n = len(buf) // self.rowlen
        return np.frombuffer(buf[:n * self.rowlen], np.uint8).reshape(n, self.rowlen)

    def parse_prefix(self, rows):
        """前置資料 -> (有效碼, 13 組 [頻道, 自旋數, RAWDELTAF/0, YSUBZ/0], 波段表)；沒有 VAS 校正資訊時後兩者為 None。"""
        p = np.ascontiguousarray(rows[:, :self.presiz])
        valid = np.frombuffer(p[:, :4].tobytes(), self.order + 'i4') if self.presiz >= 4 else None
        if self.callen < 116:
            return valid, None, None
        c0 = 4 + self.doc
        groups = np.frombuffer(np.ascontiguousarray(p[:, c0 + 12:c0 + 116]).tobytes(),
                               self.order + 'i2').reshape(len(p), 13, 4)
        lev = p[:, c0 + self.callen:c0 + self.callen + self.levlen]
        return valid, groups, lev

    @property
    def prefix(self):
        """每行的前置資料 (VAS / AAA：fmt._Data_Block 與波段表；VISR：有效碼與文件)，第一次用到時才讀。"""
        if self._prefix is None:
            rows = self.rows(0, self.nlines)[:, :self.presiz]
            if self.presiz >= fmt._Data_Block.itemsize:
                blk = np.frombuffer(np.ascontiguousarray(rows[:, :fmt._Data_Block.itemsize]).tobytes(),
                                    _dt(fmt._Data_Block, self.order))
                lev = rows[:, fmt._Data_Block.itemsize:fmt._Data_Block.itemsize + self.levlen].copy()
            else:
                dt = np.dtype([('Validity code', 'i4'), ('Documentation', 'u1', (max(0, self.presiz - 4),))])
                blk = np.frombuffer(np.ascontiguousarray(rows).tobytes(), _dt(dt, self.order)) if self.presiz else None
                lev = None
            self._prefix = (blk, lev)
        return self._prefix


class GoesPreGvarReader(GeoReader):
    NAME = 'GOES-PREGVAR'
    SENSOR = 'VISSR'
    DESCRIPTION = 'SMS-1/2、GOES-1~7 VISSR / VAS McIDAS AREA (goesNN.yyyy.ddd.hhmm.AREA_IR / AREA_VIS)'
    DESCRIPTION_EN = 'SMS-1/2, GOES-1~7 VISSR / VAS McIDAS AREA (goesNN.yyyy.ddd.hhmm.AREA_IR / AREA_VIS)'
    PATTERNS = ('goes0?.????.???.????*.AREA_*', 'sms0?.????.???.????*.AREA_*')
    CHANNELS = {}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._areas = {}
        self._navs = {}

    # ---- 檔案 --------------------------------------------------------------
    def scenes(self, folder='.'):
        groups = {}
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            path = os.path.join(folder, name)
            try:
                a = self._open(path)
                nav_params(a.navblock)
            except (OSError, ValueError) as err:
                print(tr('  [略過] %s：%s，無法定位', '  [skip] %s: %s; cannot geolocate') % (name, err))
                continue
            n = int(m.group(2))
            sat = ('GOES-%d' if m.group(1).lower() == 'goes' else 'SMS-%d') % n
            hms = int(a.dir['HHMMSS'])
            t = datetime.datetime.strptime(m.group(3) + m.group(4), '%Y%j') + datetime.timedelta(
                hours=hms // 10000, minutes=hms // 100 % 100, seconds=hms % 100)
            key = (sat, m.group(3), m.group(4), m.group(5), m.group(6))
            if key not in groups:
                groups[key] = Scene(sat, t, 'FD', format='AREA',
                                    sensor='VAS' if m.group(1).lower() == 'goes' and n >= 4 else 'VISSR')
            s = groups[key]
            if m.group(7).upper() == 'VIS':
                s.files['VIS'] = path
            elif a.source == 'VISR':
                s.files['IR'] = path
            else:
                for b in a.bands:
                    s.files.setdefault('B%02d' % b, path)
            s.extra.setdefault('source', {})[m.group(7).upper()] = a.source
        out = sorted(groups.values(), key=lambda s: (s.time, s.satellite))
        for s in out:
            # 涵蓋大部分地球 (> 9000 條可見光線) 為 FD，否則為區域 (小範圍時地圖範圍自動帶入)
            ch = next((c for c in s.files if c != 'VIS'), 'VIS')
            a = self._area(s, ch)
            if a.nlines * int(a.dir['LINERES']) > 9000:
                continue
            s.area = 'SECTOR'
            try:
                s.extra['georange'] = self._image_range(s)
                s.region, s.auto_range = sector_label(s.extra['georange'])
            except (SystemExit, ValueError):
                pass
        return out

    def _open(self, path):
        if path not in self._areas:
            self._areas[path] = AreaFile(path)
        return self._areas[path]

    def _area(self, scene, channel):
        return self._open(scene.files[channel])

    def channels(self, scene):
        out = []
        n = int(scene.satellite.split('-')[1])
        for name in sorted(scene.files, key=lambda c: (c != 'VIS', c != 'IR', c)):
            a = self._area(scene, name)
            res = float(a.dir['ELERES'])
            if name == 'VIS':
                out.append(Channel('VIS', 0.65, res, 'vis', band='VIS', alias='Visible'))
            elif name == 'IR':
                wl = 11.17 if scene.satellite.startswith('GOES') and n >= 4 else 11.5
                out.append(Channel('IR', wl, res, 'ir', band='IR', alias='LW Window'))
            else:
                b = int(name[1:])
                wl, kind, alias = _VAS[b]
                out.append(Channel(name, wl, res, kind, band='BAND-%02d' % b, alias=alias))
        return out

    def header(self, scene, channel):
        a = self._area(scene, channel)
        blk, lev = a.prefix
        return {'Control_Block': a.dir, 'Navigation': a.nav, 'Calibration': a.cal, 'Line_prefix': blk,
                'Level_map': lev, 'nav': nav_params(a.navblock)}

    # ---- 導航 --------------------------------------------------------------
    def nav(self, scene, channel):
        path = scene.files[channel]
        if path not in self._navs:
            a = self._open(path)
            d = a.dir
            lres, eres = int(d['LINERES']), int(d['ELERES'])
            first_line = int(d['UPPERLEFTLINE']) + (IR_LINE_OFFSET if lres > 1 else 0.)
            first_ele = int(d['UPPERLEFTELE']) + (IR_ELE_OFFSET if eres > 1 else 0.)
            self._navs[path] = McGoesNav((a.nlines, a.neles), a.navblock, first_line, first_ele, lres, eres,
                                         key=(path,), engine=str(self.cfg.navigation.engine).lower())
        return self._navs[path]

    def satellite_position(self, scene):
        lon, lat, height = self.nav(scene, next(iter(scene.files))).subpoint()
        return lon, height

    # ---- 讀取與校正 --------------------------------------------------------
    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        a = self._area(scene, channel)
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        rows = a.rows(r0, r1)
        n = len(rows)
        if n == 0:
            return out
        dir_valid = int(a.dir['Validity Code'])
        valid, groups, lev = a.parse_prefix(rows)
        bad_line = np.zeros(n, bool)
        if dir_valid and valid is not None:
            bad_line |= valid != dir_valid
        if a.elesiz == 1:
            b = rows[:, a.presiz + c0:a.presiz + c1].astype(np.float32)
            bad_line |= ~rows[:, a.presiz:].any(axis=1)
            if channel == 'VIS':
                vals = (b / 252.) ** 2
            else:
                vals = np.where(b >= 176, 418. - b, 330. - b / 2.) - 273.15
            vals[bad_line] = np.nan
            out[:n, :vals.shape[1]] = vals
            return out
        band = int(channel[1:])
        data = np.frombuffer(np.ascontiguousarray(rows[:, a.presiz:]).tobytes(),
                             a.order + 'i2').reshape(n, a.neles, a.nchans)
        ch2band = fmt.channel_to_band
        for i in range(n):
            if bad_line[i] or groups is None:
                continue
            slots = np.flatnonzero(lev[i, :a.nchans] == band)
            if slots.size == 0:
                continue
            slot = int(slots[0])
            g = next((g for g in groups[i] if ch2band.get(int(g[0])) == band), None)
            if g is None:
                continue
            P = data[i, c0:c1, slot].astype(np.float64)
            if a.source == 'AAA':
                ch = int(g[0])
                if a.cal is None or not 1 <= ch <= 38:
                    continue
                iab, ifab = a.cal['IAB'][ch - 1], int(a.cal['IFAB'][ch - 1])
                rad = (iab[1] * P / 32. - iab[0]) / 2. ** (15 - ifab)
            else:
                deltaf = _DELTAF[int(g[2]) % 16]
                rad = (P - int(g[3])) * 2. ** (_F_AA[band - 1] + deltaf - 15)
            t = _brightness_temperature(rad, band)
            t[P <= 0] = np.nan
            out[i, :t.size] = t
        return out


READER = GoesPreGvarReader
