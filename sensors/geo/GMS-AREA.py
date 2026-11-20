"""GMS (VISSR) McIDAS AREA 檔 — 例如 GMS_IR_1_1979_007_1533.area、GMS_VIS_1_1979_007_0533.area

檔案：McIDAS AREA 格式，一個檔案一個波段；衛星與時間取自 AREA 目錄 (Memo、YYDDD、HHMMSS)。
格式：
  AREA 目錄 64 個 word (依第 2 個 word = 4 判斷位元組順序)，欄位名稱同舊版 GVAR_AREA_format._Control_Block：
    NLINES、NELES 影像大小，UPPERLEFTLINE / ELE、LINERES / ELERES 為以可見光為單位的影像座標 (紅外 1 行 = 4 個
    可見光行、1 像素 = 2 個可見光像素)，W34 = 影像資料位置、W35 = 導航區塊位置，Image source type / Calibration Type
  導航區塊 'GMSX' (McIDAS nvxgmsx / VisAD GMSXnav)：JMA 座標轉換、姿態預報 (10 筆)、軌道預報 (8 筆) 以 GMS 自訂的
    BCD 編碼存放，解碼同舊版 read_gmsx_nav.py。舊版姿態與軌道預報的位置多減了一筆 (第一筆為 0、最後一筆遺失)，這裡已更正。
  影像：每行 PRESIZ bytes 前置 + NELES 個 1 byte 亮度 (BRIT)
校正 (同舊版)：紅外 McIDAS 亮度 B -> 溫度：B >= 176 為 418 - B，否則 330 - B/2 (K)；
              可見光 B / 256 (AREA 沒有校正表，只是相對亮度)。
導航：JMA MGIVSR (sensors/vissr_nav.py)，參數全部取自 GMSX 區塊。
完整檔頭：reader.header(scene, 'IR') 回傳 {'Control_Block': AREA 目錄, 'GMSX': 導航區塊解碼結果 (同舊版 read_gmsx_nav)}。
"""
import datetime
import os
import re

import numpy as np

from satima.core.i18n import tr
from satima.geo.reader import Channel, GeoReader, Scene
from vissr_nav import EA, VissrNav

_Control_Block = [('W1', 'i4'), ('Area Format', 'i4'), ('SSS', 'i4'), ('YYDDD', 'i4'), ('HHMMSS', 'i4'),
                  ('UPPERLEFTLINE', 'i4'), ('UPPERLEFTELE', 'i4'), ('W8', 'i4'), ('NLINES', 'i4'), ('NELES', 'i4'),
                  ('ELESIZ', 'i4'), ('LINERES', 'i4'), ('ELERES', 'i4'), ('NCHANS', 'i4'), ('PRESIZ', 'i4'),
                  ('PROJ', 'i4'), ('CREATION DATE', 'i4'), ('CREATION TIME', 'i4'), ('FILTER MAP', 'i4')] + \
                 [('W%d' % i, 'i4') for i in range(20, 25)] + \
                 [('Memo', 'S32'), ('Area file number', 'i4'), ('W34', 'i4'), ('W35', 'i4'), ('Validity Code', 'i4')] + \
                 [('W%d' % i, 'i4') for i in range(37, 52)] + \
                 [('Image source type', 'S4'), ('Calibration Type', 'S4')] + [('W%d' % i, 'i4') for i in range(54, 65)]

# (名稱, 中心波長 µm, 解析度 km, 種類, 功能代稱)；GMS-1~4 VISSR
_BANDS = (('VIS', 0.63, 1.25, 'vis', 'Visible'), ('IR', 11.5, 5.0, 'ir', 'LW Window'))
_SENSOR_INDEX = {'VIS': 0, 'IR': 1}       # GMSX 參數陣列的感測器索引 (VIS、IR1、IR2、WV)


# ---- GMSX 導航區塊 -----------------------------------------------------------
def _sv0100(iword, ipos, b, off):
    """GMS 自訂的數值編碼：iword (4 或 6) bytes 的整數，最高位元為正負號，再除以 10^ipos。"""
    value = int.from_bytes(bytes([b[off] & 0x7F]) + bytes(b[off + 1:off + iword]), 'big') / 10.0 ** ipos
    return -value if b[off] & 0x80 else value


def gmsx_bytes(block):
    """導航區塊 (raw bytes) -> 參數 bytes：去掉每 128 words 開頭的文字 word ('GMSX' 與續行標記)。"""
    if block[:4] != b'GMSX':
        raise ValueError(tr('不是 GMSX 導航區塊 (%r)', 'not a GMSX navigation block (%r)') % block[:4])
    return block[4:4 * 127] + b''.join(block[4 * 128 * k:4 * (128 * k + 127)] for k in range(1, 5))


def decode_gmsx(b):
    """參數 bytes -> 導航參數 (鍵名同舊版 read_gmsx_nav.decode_oa_block)。"""
    sv = _sv0100
    nav = dict(dtims=sv(6, 8, b, 0), dspin=sv(6, 8, b, 240),
               resLin=[sv(4, 8, b, 6), sv(4, 8, b, 10), sv(4, 8, b, 10), sv(4, 8, b, 10)],
               resEle=[sv(4, 10, b, 14), sv(4, 10, b, 18), sv(4, 10, b, 18), sv(4, 10, b, 18)],
               rlic=[sv(4, 4, b, 22), sv(4, 4, b, 26), sv(4, 4, b, 110), sv(4, 4, b, 114)],
               relmfc=[sv(4, 4, b, 30), sv(4, 4, b, 34), sv(4, 4, b, 118), sv(4, 4, b, 122)],
               senssu=[sv(4, 0, b, 38), sv(4, 0, b, 42), sv(4, 0, b, 42), sv(4, 0, b, 42)],
               rline=[sv(4, 0, b, 46), sv(4, 0, b, 50), sv(4, 0, b, 50), sv(4, 0, b, 50)],
               relem=[sv(4, 0, b, 54), sv(4, 0, b, 58), sv(4, 0, b, 58), sv(4, 0, b, 58)],
               vmis=[sv(4, 10, b, 62), sv(4, 10, b, 66), sv(4, 10, b, 70)],
               subLon=sv(6, 6, b, 198), subLat=sv(6, 6, b, 204))
    nav['elmis'] = np.array([sv(4, 7, b, 74), sv(4, 10, b, 78), sv(4, 10, b, 82),
                             sv(4, 10, b, 86), sv(4, 7, b, 90), sv(4, 10, b, 94),
                             sv(4, 10, b, 98), sv(4, 10, b, 102), sv(4, 7, b, 106)]).reshape(3, 3).T
    # 姿態預報 10 筆 (每筆 48 bytes，自 256 起)：0 MJD、2 自旋軸赤經、3 赤緯、4 太陽-地球夾角 β、5 (保留)
    atit = {k: [0.0] * 10 for k in (0, 2, 3, 4, 5)}
    for i in range(10):
        j = 256 + 48 * i
        atit[0][i], atit[2][i] = sv(6, 8, b, j), sv(6, 8, b, j + 12)
        atit[3][i], atit[4][i], atit[5][i] = sv(6, 11, b, j + 18), sv(6, 8, b, j + 24), sv(6, 8, b, j + 30)
    nav['atit'] = atit
    # 軌道預報 8 筆 (每筆 200 bytes，自 736 起)：0 MJD、8~10 地固座標位置 (m)、14 格林威治恆星時、
    # 17 / 18 太陽方位 / 仰角 (度)、19~27 歲差章動矩陣 A1~A9
    spec = {0: (0, 6, 8), 8: (48, 6, 6), 9: (54, 6, 6), 10: (60, 6, 6), 14: (84, 6, 8), 17: (102, 6, 8),
            18: (108, 6, 8), 19: (128, 6, 12), 20: (134, 6, 14), 21: (140, 6, 14), 22: (146, 6, 14),
            23: (152, 6, 12), 24: (158, 6, 16), 25: (164, 6, 12), 26: (170, 6, 16), 27: (176, 6, 12)}
    orbt1 = {k: [0.0] * 8 for k in spec}
    for i in range(8):
        j = 736 + 200 * i
        for key, (rel, iword, ipos) in spec.items():
            orbt1[key][i] = sv(iword, ipos, b, j + rel)
    nav['orbt1'] = orbt1
    return nav


def read_directory(path):
    """AREA 目錄與導航區塊 -> (目錄 numpy 結構, GMSX 解碼結果)。"""
    with open(path, 'rb') as f:
        raw = f.read(256)
        order = '<' if np.frombuffer(raw[4:8], '<i4')[0] == 4 else '>'
        dtype = np.dtype([(n, order + t if t == 'i4' else t) for n, t in _Control_Block])
        d = np.frombuffer(raw, dtype, count=1)[0]
        f.seek(int(d['W35']))
        block = f.read(4 * 640)
    return d, decode_gmsx(gmsx_bytes(block))


class GmsAreaReader(GeoReader):
    NAME = 'GMS-AREA'
    SENSOR = 'VISSR'
    DESCRIPTION = 'GMS VISSR McIDAS AREA (GMSX navigation)'
    PATTERNS = ('GMS*.area',)
    CHANNELS = {b[0]: Channel(*b[:4], alias=b[4]) for b in _BANDS}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._dirs = {}

    def _directory(self, path):
        if path not in self._dirs:
            self._dirs[path] = read_directory(path)
        return self._dirs[path]

    def scenes(self, folder='.'):
        groups = {}
        for name in self.find(folder):
            path = os.path.join(folder, name)
            try:
                d, _ = self._directory(path)
            except (OSError, ValueError) as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                continue
            memo = d['Memo'].decode('ascii', 'replace').upper()
            band = 'VIS' if 'VIS' in memo or '_VIS_' in name.upper() else 'IR'
            m = re.search(r'GMS\s*-?\s*(\d)', memo) or re.search(r'GMS_[A-Z]+_(\d)', name.upper())
            sat = 'GMS-%s' % m.group(1) if m else 'GMS'
            yyddd, hhmmss = int(d['YYDDD']), int(d['HHMMSS'])
            t = (datetime.datetime(1900 + yyddd // 1000, 1, 1) + datetime.timedelta(days=yyddd % 1000 - 1) +
                 datetime.timedelta(hours=hhmmss // 10000, minutes=hhmmss // 100 % 100))
            key = (sat, t)
            if key not in groups:
                groups[key] = Scene(sat, t, 'FD')
            groups[key].files[band] = path
        return sorted(groups.values(), key=lambda s: (s.time, s.satellite))

    def header(self, scene, channel):
        d, gmsx = self._directory(scene.files[channel])
        return {'Control_Block': d, 'GMSX': gmsx}

    def nav(self, scene, channel):
        d, g = self._directory(scene.files[channel])
        n = _SENSOR_INDEX[channel]
        # 目錄的影像座標以可見光為單位 -> 該波段的 VISSR 行號、像素號
        kl, ke = g['rline'][0] / g['rline'][n], g['relem'][0] / g['relem'][n]
        scan = dict(step=g['resLin'][n], samp=g['resEle'][n], line0=g['rlic'][n], pixel0=g['relmfc'][n],
                    sensors=g['senssu'][n] or 1., elmis=g['elmis'], t0=g['dtims'], spin=g['dspin'])
        a, o = g['atit'], g['orbt1']
        npa = np.stack([o[k] for k in range(19, 28)], axis=-1).reshape(-1, 3, 3).transpose(0, 2, 1)
        return VissrNav(
            (int(d['NLINES']), int(d['NELES'])), scan,
            attitude=dict(mjd=a[0], alpha=a[2], delta=a[3], beta=a[4]),
            orbit=dict(mjd=o[0], sat=np.stack([o[8], o[9], o[10]], axis=-1), gst=o[14], sun_az=o[17], sun_el=o[18],
                       npa=npa),
            first_line=(int(d['UPPERLEFTLINE']) - 1) / kl + 1, line_step=int(d['LINERES']) / kl,
            first_pixel=(int(d['UPPERLEFTELE']) - 1) / ke + 1, pixel_step=int(d['ELERES']) / ke,
            key=(scene.files[channel], channel), engine=str(self.cfg.navigation.engine).lower())

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        d, _ = self._directory(scene.files[channel])
        nele, size, nch = int(d['NELES']), int(d['ELESIZ']), int(d['NCHANS'])
        if size != 1 or nch != 1:
            raise ValueError(tr('只支援單一波段、每像素 1 byte 的 AREA (ELESIZ=%d, NCHANS=%d)',
                                'only single-band AREA files with 1 byte per pixel are supported (ELESIZ=%d, NCHANS=%d)')
                             % (size, nch))
        rowlen = int(d['PRESIZ']) + nele
        with open(scene.files[channel], 'rb') as f:
            f.seek(int(d['W34']) + r0 * rowlen)
            buf = f.read((r1 - r0) * rowlen)
        n = len(buf) // rowlen
        b = np.frombuffer(buf[:n * rowlen], np.uint8).reshape(n, rowlen)[:, int(d['PRESIZ']) + c0:int(d['PRESIZ']) + c1]
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        b = b.astype(np.float32)
        if channel == 'VIS':
            out[:n] = b / 256.
        else:
            out[:n] = np.where(b >= 176, 418. - b, 330. - b / 2.) - 273.15
        return out

    def satellite_position(self, scene):
        _, g = self._directory(next(iter(scene.files.values())))
        sat = np.array([g['orbt1'][k][0] for k in (8, 9, 10)])
        return g['subLon'], float(np.linalg.norm(sat)) / 1000. - EA / 1000.


READER = GmsAreaReader
