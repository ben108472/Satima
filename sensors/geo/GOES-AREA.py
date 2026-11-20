"""GOES-8~15 Imager — McIDAS AREA (NOAA CLASS 的 GVAR_IMG，例如 goes08.1998.299.181515.BAND_01)

檔案：一個波段一個檔案 (BAND_01~06 = GVAR 頻道 1~6)，同一時間的各波段合成一筆資料。
格式 (結構與欄位名稱同舊版 GVAR_AREA_format.py，見 sensors/gvar_area_format.py)：
  AREA 目錄 64 words (依第 2 個 word = 4 判斷位元組順序)：UPPERLEFTLINE / ELE、LINERES / ELERES 為以可見光
    為單位的影像座標 (紅外 4)，W34 = 影像位置、W35 = 導航區塊、W63 = 校正區塊
  導航區塊 'GVAR' (640 words)：GVAR 的 O&A 參數 (放大 10^7 的整數)
  校正區塊 (128 words)：GVAR Block 0 的可見光係數與紅外 scaling factor (Gould 浮點數原樣保存)
  每行：前置資料 PRESIZ (232) bytes (有效碼、掃描狀態、BCD 時間、GVAR 區塊頭、Line Documentation：side、偵測器、頻道…)
        + 每像素 2 bytes (10 位元計數值 x 32)
校正：同 GOES-GVAR (紅外依 conversion_table/GOES 的表；可見光依校正區塊各偵測器的係數)，計數值 0 (填充) 為 NaN。
導航：GOES I-M O&A 模型 (sensors/goes_imager.py)；IMC 關閉時 (例如 2003~2005 在 155°E 的 GOES-9) 依每行前置資料的
      掃描時間逐線計算。
完整檔頭：reader.header(scene, 'CH4') 回傳 {'Control_Block': AREA 目錄, 'Navigation': 導航區塊, 'Calibration': 校正區塊,
  'Line_prefix': 每行前置資料 (_Data_Block), 'nav': 導航參數}。
"""
import datetime
import os
import re

import numpy as np

import gvar_area_format as fmt
from goes_imager import (IR_LINE_CENTER, IR_PIXEL_CENTER, GoesNav, bcd_time, fit_line_time, gould, ir_lut,
                         make_channel, params_from_area, sector_label, vis_coefficients, vis_scale)
from satima.core.i18n import tr
from satima.geo.reader import GeoReader, Scene

_NAME =re.compile(r'^goes(\d\d)\.(\d{4})\.(\d{3})\.(\d{6})\.BAND_(\d\d)$', re.I)


def _dtype(dt, order):
    """把舊版 (原生位元組順序) 的結構改成指定位元組順序。"""
    return dt.newbyteorder(order) if order != '=' else dt


class AreaFile:
    def __init__(self, path):
        self.path = path
        with open(path, 'rb') as f:
            raw = f.read(256)
            self.order = '<' if np.frombuffer(raw[4:8], '<i4')[0] == 4 else '>'
            d = np.frombuffer(raw, _dtype(fmt._Control_Block, self.order), count=1)[0]
            self.dir = d
            f.seek(int(d['W35']))
            self.nav = np.frombuffer(f.read(fmt._Navigation.itemsize), _dtype(fmt._Navigation, self.order), count=1)[0]
            cal_off = int(d['W63'])
            f.seek(cal_off)
            craw = f.read(fmt._Calibration.itemsize) if cal_off > 0 else b''
            self.cal = np.frombuffer(craw, _dtype(fmt._Calibration, self.order), count=1)[0] \
                if len(craw) == fmt._Calibration.itemsize else None
            self.nlines, self.neles = int(d['NLINES']), int(d['NELES'])
            self.presiz, self.elesiz = int(d['PRESIZ']), int(d['ELESIZ'])
            self.rowlen = self.presiz + self.neles * self.elesiz * int(d['NCHANS'])
        fm = int(d['FILTER MAP'])
        self.band = fm.bit_length() if fm and fm & (fm - 1) == 0 else None
        self._prefix = None

    def _read_prefix(self):
        """每行的前置資料 (第一次用到時才讀)。"""
        raw = []
        with open(self.path, 'rb') as f:
            for i in range(self.nlines):
                f.seek(int(self.dir['W34']) + i * self.rowlen)
                raw.append(f.read(fmt._Data_Block.itemsize).ljust(fmt._Data_Block.itemsize, b'\0'))
        self._prefix = (np.frombuffer(b''.join(raw), _dtype(fmt._Data_Block, self.order)), raw)

    @property
    def prefix(self):
        if self._prefix is None:
            self._read_prefix()
        return self._prefix[0]

    @property
    def prefix_raw(self):
        if self._prefix is None:
            self._read_prefix()
        return self._prefix[1]

    def line_time(self, first_line, line_step):
        """每行前置資料的 BCD 時間 (bytes 10~17) -> 逐線時間模型。"""
        times = [bcd_time(p[10:18]) for p in self.prefix_raw]
        lines = first_line + line_step * np.arange(self.nlines)
        return fit_line_time(lines, times)


class GoesAreaReader(GeoReader):
    NAME = 'GOES-AREA'
    SENSOR = 'IMAGER'
    DESCRIPTION = 'GOES-8~15 Imager McIDAS AREA (NOAA CLASS GVAR_IMG .BAND_0?)'
    PATTERNS = ('goes??.????.???.??????.BAND_0?',)
    CHANNELS = {'CH%d' % n: make_channel(n, 'GOES-8' if n != 6 else 'GOES-13') for n in range(1, 7)}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._areas = {}
        self._cal = {}

    def scenes(self, folder='.'):
        groups = {}
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            sat = 'GOES-%d' % int(m.group(1))
            t = datetime.datetime.strptime(m.group(2) + m.group(3) + m.group(4), '%Y%j%H%M%S')
            key = (sat, t)
            if key not in groups:
                groups[key] = Scene(sat, t, 'AREA')
            groups[key].files['CH%d' % int(m.group(5))] = os.path.join(folder, name)
        out = sorted(groups.values(), key=lambda s: (s.time, s.satellite))
        for s in out:
            # 涵蓋全圓盤 (> 9000 條可見光線) 為 FD，否則為區域 (小範圍時地圖範圍自動帶入)
            a = self._area(s, next(iter(s.files)))
            if a.nlines * int(a.dir['LINERES']) > 9000:
                s.area = 'FD'
                continue
            s.area = 'SECTOR'
            try:
                s.extra['georange'] = self._image_range(s)
                s.region, s.auto_range = sector_label(s.extra['georange'])
            except (SystemExit, ValueError):
                pass
        return out

    def _area(self, scene, channel):
        path = scene.files[channel]
        if path not in self._areas:
            self._areas[path] = AreaFile(path)
        return self._areas[path]

    def channels(self, scene):
        out = []
        for name in sorted(scene.files, key=lambda c: int(c[2:])):
            a = self._area(scene, name)
            out.append(make_channel(int(name[2:]), scene.satellite, int(a.dir['LINERES']) if name != 'CH1' else 1))
        return out

    def header(self, scene, channel):
        a = self._area(scene, channel)
        return {'Control_Block': a.dir, 'Navigation': a.nav, 'Calibration': a.cal, 'Line_prefix': a.prefix,
                'nav': params_from_area(a.nav)}

    # ---- 導航 --------------------------------------------------------------
    def nav(self, scene, channel):
        a = self._area(scene, channel)
        d = a.dir
        lres, eres = int(d['LINERES']), int(d['ELERES'])
        first_line = int(d['UPPERLEFTLINE']) + (IR_LINE_CENTER * lres / 4. if lres > 1 else 0.)
        first_pixel = int(d['UPPERLEFTELE']) + (IR_PIXEL_CENTER if eres > 1 else 0.)
        p = params_from_area(a.nav)
        lt = None if p['imc'] else a.line_time(first_line, lres)
        return GoesNav((a.nlines, a.neles), p, first_line, first_pixel, lres, eres, key=(scene.files[channel],),
                       engine=str(self.cfg.navigation.engine).lower(), line_time=lt)

    def satellite_position(self, scene):
        st = self.nav(scene, next(iter(scene.files))).state
        return st['sublon'], st['height']

    # ---- 讀取 --------------------------------------------------------------
    def _vis_coef(self, scene, a):
        if a.cal is None:
            bias, g1, g2, alb = vis_coefficients(scene.satellite, [0] * 8, [0] * 8, [0] * 8, 0)
        else:
            c = a.cal
            get = lambda fmt_, i: gould(int(c[fmt_ % i]))
            bias, g1, g2, alb = vis_coefficients(
                scene.satellite, [get('Visible bias coefficients %d', i) for i in range(1, 9)],
                [get('Visible first order gain coefficients %d', i) for i in range(1, 9)],
                [get('Visible second order gain coefficients %d', i) for i in range(1, 9)],
                gould(int(c['Visible radiance to albedo conversion factor'])))
        return bias, g1, g2, alb * vis_scale(self.cfg, scene.satellite, scene.time)

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        a = self._area(scene, channel)
        if a.elesiz != 2 or int(a.dir['NCHANS']) != 1:
            raise ValueError(tr('只支援單一波段、每像素 2 bytes 的 GVAR AREA (ELESIZ=%d)',
                                'only single-band GVAR AREA files with 2 bytes per pixel are supported (ELESIZ=%d)')
                             % a.elesiz)
        number = int(channel[2:])
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        with open(scene.files[channel], 'rb') as f:
            f.seek(int(a.dir['W34']) + r0 * a.rowlen)
            buf = f.read((r1 - r0) * a.rowlen)
        n = len(buf) // a.rowlen
        if n == 0:
            return out
        rows = np.frombuffer(buf[:n * a.rowlen], np.uint8).reshape(n, a.rowlen)[:, a.presiz:]
        counts = np.frombuffer(np.ascontiguousarray(rows).tobytes(), a.order + 'i2').reshape(n, a.neles)
        counts = (counts[:, c0:c1].astype(np.int64) >> 5)
        pre = a.prefix[r0:r0 + n]
        if number == 1:
            bias, g1, g2, alb = self._vis_coef(scene, a)
            det = np.clip(pre['LIDET'].astype(np.int64) - 1, 0, 7)[:, None]
            vals = (bias[det] + g1[det] * counts + g2[det] * counts * counts) * alb
        else:
            vals = np.full(counts.shape, np.nan)
            table_dir = self.cfg.general.conversion_table_dir
            side = np.where(pre['LSIDE'] == 1023, 2, 1)          # 同 GOES-GVAR：1023 = side 2，其餘 (含損壞值) side 1
            det = np.where(pre['LIDET'] % 2 == 0, 2, 1)
            for s_, d_ in set(zip(side.tolist(), det.tolist())):
                key = (scene.satellite, s_, number, d_)
                if key not in self._cal:
                    self._cal[key] = ir_lut(table_dir, scene.satellite, s_, number, d_)
                lut = self._cal[key]
                if lut is None:
                    continue
                sel = (side == s_) & (det == d_)
                vals[sel] = lut[np.clip(counts[sel], 0, 1023)]
        vals = np.where(counts > 0, vals, np.nan)
        out[:n] = vals
        return out


READER = GoesAreaReader
