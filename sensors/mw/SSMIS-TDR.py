"""DMSP SSMIS TDR 天線溫度 (NOAA CLASS 的 NPR.TDRN.*，大端二進位)

!! 這支讀取器沒有樣本檔可以驗證：格式照舊版 Satima-SSMIS 的 SSMIS_TDR.py 移植 (該程式只讀過 91 GHz)。
!! 檔案大小與預期的記錄長度不符時會直接報錯而不會畫出錯的圖；其餘頻道的欄位位置請以實際檔案確認。

檔名  NPR.TDRN.S<A-D>.D<yyddd>.S<hhmm>.E<hhmm>.B<軌道>.<站>   SA / SB / SC / SD = F16 / F17 / F18 / F19
結構  檔頭 40 bytes (第 18 byte 起的 int16 為掃描數)，之後每條掃描 9592 bytes：
        掃描檔頭 36 bytes (年 u4、日序 u2、時、分、掃描時間 u4 (當日毫秒)…)
        星曆 20 bytes x 3
        影像場景 24 bytes x 180：ch8~11 的緯度、經度、場景號、地表、降雨、ch8、ch9、ch10、ch11、
                                 ch17/18 的緯度、經度、ch17、ch18
        環境場景 20 bytes x 90 ：ch12~14 的緯度、經度、場景號、地表、ch12、ch13、ch14、
                                 ch15/16 的緯度、經度、ch15、ch16
        LAS 場景 24 bytes x 60 ：緯度、經度、場景號、地表、ch1~ch7、ch24
        UAS 場景 16 bytes x 30、輔助資料 1456 bytes (此處不讀)
      緯度、經度為 0.01 度；溫度為 0.01 °C (+273.15 -> K)
TDR 是天線溫度 (TA)，頻道附註 TA。頻道與 SSMIS-SDR 相同。
"""
import datetime
import os
import re

import numpy as np

from satima.core.i18n import tr
from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, to_seconds

_NAME = re.compile(r'^NPR\.TDRN\.S([A-D])\.D(\d{5})\.S(\d{4})\.E(\d{4})\.B(\d+)', re.I)
_SATELLITES = {'A': 'F16', 'B': 'F17', 'C': 'F18', 'D': 'F19'}
_HEAD = np.dtype([('year', '>u4'), ('jday', '>u2'), ('hour', 'u1'), ('minute', 'u1'), ('ms', '>u4'), ('rest', 'V24')])
_IMG = np.dtype([('lat8', '>i2'), ('lon8', '>i2'), ('scene', '>u2'), ('sfc', 'u1'), ('rain', 'u1'),
                 ('c8', '>i2'), ('c9', '>i2'), ('c10', '>i2'), ('c11', '>i2'),
                 ('lat', '>i2'), ('lon', '>i2'), ('c17', '>i2'), ('c18', '>i2')])
_ENV = np.dtype([('lat12', '>i2'), ('lon12', '>i2'), ('scene', 'u1'), ('sfc', 'u1'),
                 ('c12', '>i2'), ('c13', '>i2'), ('c14', '>i2'), ('lat', '>i2'), ('lon', '>i2'),
                 ('c15', '>i2'), ('c16', '>i2')])
_LAS = np.dtype([('lat', '>i2'), ('lon', '>i2'), ('scene', '>u2'), ('sfc', '>u2')] +
                [('c%d' % k, '>i2') for k in (1, 2, 3, 4, 5, 6, 7, 24)])
_SCAN = np.dtype([('head', _HEAD), ('eph', 'V60'), ('img', _IMG, 180), ('env', _ENV, 90), ('las', _LAS, 60),
                  ('uas', 'V480'), ('aux', 'V1456')])
assert _SCAN.itemsize == 9592
# 組 -> (結構中的欄位, 緯度欄位, 經度欄位, [(欄位, 頻率, 偏移, 極化), ...])
GROUPS = {
    'ENV19': ('env', 'lat12', 'lon12', [('c13', 19.35, None, 'V'), ('c12', 19.35, None, 'H'),
                                        ('c14', 22.235, None, 'V')]),
    'ENV37': ('env', 'lat', 'lon', [('c16', 37.0, None, 'V'), ('c15', 37.0, None, 'H')]),
    'IMG91': ('img', 'lat', 'lon', [('c17', 91.655, None, 'V'), ('c18', 91.655, None, 'H')]),
    'IMG150': ('img', 'lat8', 'lon8', [('c8', 150.0, None, 'H'), ('c11', 183.31, 1., 'H'), ('c10', 183.31, 3., 'H'),
                                       ('c9', 183.31, 6.6, 'H')]),
    'LAS': ('las', 'lat', 'lon', [('c1', 50.3, None, ''), ('c2', 52.8, None, ''), ('c3', 53.596, None, ''),
                                  ('c4', 54.4, None, ''), ('c5', 55.5, None, ''), ('c6', 57.29, None, ''),
                                  ('c7', 59.4, None, ''), ('c24', 60.79, None, '')]),
}
_TB_RANGE = (40., 400.)


def _times(yyddd, start, end):
    yy = int(yyddd[:2])
    day = datetime.datetime(1900 + yy if yy >= 80 else 2000 + yy, 1, 1) + datetime.timedelta(days=int(yyddd[2:]) - 1)
    a = day + datetime.timedelta(hours=int(start[:2]), minutes=int(start[2:]))
    b = day + datetime.timedelta(hours=int(end[:2]), minutes=int(end[2:]))
    if b < a:
        b += datetime.timedelta(days=1)
    return a, b


def parse(path):
    """整個檔讀成 (掃描的結構化陣列, 各掃描時間)。檔案大小不符時 ValueError。"""
    with open(path, 'rb') as f:
        buf = f.read()
    body = len(buf) - 40
    if body <= 0 or body % _SCAN.itemsize:
        raise ValueError(tr('%s 的大小不符合 SSMIS TDR 的記錄長度 (%d bytes)；這支讀取器未經實際檔案驗證，'
                            '格式可能與預期不同',
                            'the size of %s does not match the SSMIS TDR record length (%d bytes); this reader has '
                            'not been verified against real files, so the format may differ from what it expects')
                         % (os.path.basename(path), _SCAN.itemsize))
    arr = np.frombuffer(buf, _SCAN, body // _SCAN.itemsize, 40)
    head = arr['head']
    time = np.full(arr.size, np.nan)
    for k in range(arr.size):
        try:
            day = to_seconds(datetime.datetime(int(head['year'][k]), 1, 1)) + (int(head['jday'][k]) - 1) * 86400.
        except (ValueError, OverflowError):
            continue
        time[k] = day + float(head['ms'][k]) / 1000.
    return arr, time


class SSMISTDRReader(SwathReader):
    NAME = 'SSMIS-TDR'
    SENSOR = 'SSMIS'
    SOURCE = 'CLASS'
    DESCRIPTION = 'DMSP F16~F19 SSMIS TDR 天線溫度 (NOAA CLASS NPR.TDRN；未經樣本驗證)'
    DESCRIPTION_EN = ('DMSP F16~F19 SSMIS TDR antenna temperature (NOAA CLASS NPR.TDRN; not verified against '
                      'sample files)')
    PATTERNS = ('NPR.TDRN.S[A-Da-d].*',)

    def __init__(self, cfg):
        super().__init__(cfg)
        self._parsed = {}

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(os.path.basename(name))
            if not m or name.lower().endswith(('.gz', '.zip', '.tar')):
                continue
            start, end = _times(m.group(2), m.group(3), m.group(4))
            out.append(Scene(_SATELLITES[m.group(1).upper()], 'SSMIS', start, end, file=os.path.join(folder, name),
                             orbit=m.group(5), format='TDR'))
        return out

    def channels(self, scene):
        return [Channel(band_name(freq, pol, offset), freq, pol, g, key=field, note='TA')
                for g in ('ENV19', 'ENV37', 'IMG91', 'IMG150', 'LAS') for field, freq, offset, pol in GROUPS[g][3]]

    def _data(self, scene):
        if scene.file not in self._parsed:
            self._parsed.clear()
            self._parsed[scene.file] = parse(scene.file)
        return self._parsed[scene.file]

    def geo(self, scene, group):
        arr, time = self._data(scene)
        part, la, lo, _ = GROUPS[group]
        lat, lon = clean_latlon(arr[part][la] / 100., arr[part][lo] / 100.)
        return [Geo(lat, lon, time, (scene.file, group))]

    def tb(self, scene, channel):
        arr, _ = self._data(scene)
        data = arr[GROUPS[channel.group][0]][channel.key].astype(np.float32) / np.float32(100.) + np.float32(273.15)
        data[~((data > _TB_RANGE[0]) & (data < _TB_RANGE[1]))] = np.nan
        return [data]


READER = SSMISTDRReader
