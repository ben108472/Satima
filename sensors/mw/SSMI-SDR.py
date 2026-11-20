"""DMSP SSM/I SDR 亮溫與 TDR 天線溫度 (NOAA CLASS 的 NPR.SDRR.*、NPR.TDRR.*，大端二進位)

檔名  NPR.<SDRR|TDRR>.S<4-9>.D<yyddd>.S<hhmm>.E<hhmm>.B<軌道>.<站>   S4~S9 = F10~F15
結構  (2026/11 以 6 個樣本檔確認：檔案大小 = 檔頭 + 整數筆記錄)
  SDR：檔頭 3348 bytes，每筆記錄 3348 bytes，場景資料自記錄內第 16 byte 起；掃描時間為第 6 byte 的 int32 (當日秒數)
  TDR：檔頭 3604 bytes，每筆記錄 3604 bytes，場景資料自記錄內第 274 byte 起；掃描時間為第 6 byte 的 int32 (1e-4 秒)
  一筆記錄 = 一對 A / B 掃描 = 64 個低頻場景，每個低頻場景後面跟 3 個 85 GHz 場景 (共 52 bytes)：
    低頻場景 22 bytes：序號、緯度、經度、19V、19H、22V、37V、37H、85V、85H、地表、位置
    85 GHz 場景 10 bytes：緯度、經度、85V、85H、地表、位置
    緯度 = 值 / 100 - 90，經度 = 值 / 100 (0~360)，溫度 = 值 / 100 K
  每個低頻場景 (第 k 個) 的 4 個 85 GHz 取樣點：
    低頻場景本身 = A 掃描第 2k 點，第 2 個 85 GHz 場景 = A 掃描第 2k+1 點，
    第 1 個 = B 掃描第 2k 點，第 3 個 = B 掃描第 2k+1 點   (由各點經緯度的相對位置確認)
  所以 19 / 22 / 37 GHz 是 (記錄數, 64)、間距 25 km；85 GHz 把 A、B 掃描依序交錯成 (2 x 記錄數, 128)、間距 12.5 km。
TDR 是天線溫度 (TA，未做天線場型修正，與亮溫差數 K)：頻道附註 TA，仍可做 PCT 與 RGB，但數值略有偏差。
"""
import datetime
import os
import re

import numpy as np

from satima.core.i18n import tr
from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, to_seconds

_NAME = re.compile(r'^NPR\.(SDR|TDR)R\.S([4-9])\.D(\d{5})\.S(\d{4})\.E(\d{4})\.B(\d+)', re.I)
_LOW = np.dtype([('cnt', '>i2'), ('lat', '>u2'), ('lon', '>u2'), ('v19', '>i2'), ('h19', '>i2'), ('v22', '>i2'),
                 ('v37', '>i2'), ('h37', '>i2'), ('v85', '>i2'), ('h85', '>i2'), ('sfc', 'i1'), ('pos', 'u1')])
_HI = np.dtype([('lat', '>u2'), ('lon', '>u2'), ('v85', '>i2'), ('h85', '>i2'), ('sfc', 'i1'), ('pos', 'u1')])
_SCENE = np.dtype([('low', _LOW), ('hi', _HI, 3)])
# 種類 -> (檔頭長度, 記錄長度, 場景資料在記錄內的位置, 掃描時間的單位 (秒))
_LAYOUT = {'SDR': (3348, 3348, 16, 1.), 'TDR': (3604, 3604, 274, 1e-4)}
_CHANNELS = (('v19', 19.35, 'V', 'LOW'), ('h19', 19.35, 'H', 'LOW'), ('v22', 22.235, 'V', 'LOW'),
             ('v37', 37.0, 'V', 'LOW'), ('h37', 37.0, 'H', 'LOW'), ('v85', 85.5, 'V', 'HI'), ('h85', 85.5, 'H', 'HI'))
_TB_RANGE = (40., 400.)


def _times(yyddd, start, end):
    """檔名的 D<yyddd>、S<hhmm>、E<hhmm> -> (起, 訖)；跨日時結束時間加一天。"""
    yy = int(yyddd[:2])
    day = datetime.datetime(1900 + yy if yy >= 80 else 2000 + yy, 1, 1) + datetime.timedelta(days=int(yyddd[2:]) - 1)
    a = day + datetime.timedelta(hours=int(start[:2]), minutes=int(start[2:]))
    b = day + datetime.timedelta(hours=int(end[:2]), minutes=int(end[2:]))
    if b < a:
        b += datetime.timedelta(days=1)
    return a, b


def parse(path, kind, day_seconds):
    """整個檔讀成 (場景陣列 (記錄數, 64)，各記錄 A 掃描的時間 (1970 起算的秒數))。"""
    header, size, offset, unit = _LAYOUT[kind]
    with open(path, 'rb') as f:
        buf = f.read()
    nrec = (len(buf) - header) // size
    if nrec <= 0:
        raise ValueError(tr('%s 不是 SSM/I %s 檔 (太短)', '%s is not an SSM/I %s file (too short)')
                         % (os.path.basename(path), kind))
    rec = np.dtype({'names': ['time', 'scene'], 'formats': ['>i4', (_SCENE, 64)], 'offsets': [6, offset],
                    'itemsize': size})
    arr = np.frombuffer(buf, rec, nrec, header)
    sec = arr['time'].astype(np.float64) * unit
    # 掃描時間是當日秒數：檔案跨日時 (比第一筆小很多) 加一天
    sec = np.where(sec < sec[0] - 43200., sec + 86400., sec)
    ok = (sec >= 0) & (sec < 2 * 86400.)
    return arr['scene'], np.where(ok, day_seconds + sec, np.nan)


class SSMIReader(SwathReader):
    NAME = 'SSMI-SDR'
    SENSOR = 'SSM/I'
    SOURCE = 'CLASS'
    DESCRIPTION = 'DMSP F10~F15 SSM/I SDR 亮溫 / TDR 天線溫度 (NOAA CLASS NPR.SDRR、NPR.TDRR)'
    DESCRIPTION_EN = ('DMSP F10~F15 SSM/I SDR brightness temperature / TDR antenna temperature '
                      '(NOAA CLASS NPR.SDRR, NPR.TDRR)')
    PATTERNS = ('NPR.SDRR.S[4-9].*', 'NPR.TDRR.S[4-9].*')
    PARSE_CACHE = 2

    def __init__(self, cfg):
        super().__init__(cfg)
        self._parsed = {}

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(os.path.basename(name))
            if not m or name.lower().endswith(('.gz', '.zip', '.tar')):
                continue
            start, end = _times(m.group(3), m.group(4), m.group(5))
            kind = m.group(1).upper()
            out.append(Scene('F%d' % (int(m.group(2)) + 6), 'SSM/I', start, end, file=os.path.join(folder, name),
                             orbit=m.group(6), format=kind, extra={'kind': kind}))
        return out

    def channels(self, scene):
        note = 'TA' if scene.extra['kind'] == 'TDR' else ''
        return [Channel(band_name(freq, pol), freq, pol, group, key=field, note=note)
                for field, freq, pol, group in _CHANNELS]

    def _data(self, scene):
        if scene.file not in self._parsed:
            while len(self._parsed) >= self.PARSE_CACHE:
                self._parsed.pop(next(iter(self._parsed)))
            day = scene.time.replace(hour=0, minute=0, second=0, microsecond=0)
            self._parsed[scene.file] = parse(scene.file, scene.extra['kind'], to_seconds(day))
        return self._parsed[scene.file]

    @staticmethod
    def _hi(scene_arr, field):
        """85 GHz：A、B 掃描交錯成 (2 x 記錄數, 128)。"""
        n = scene_arr.shape[0]
        out = np.empty((2 * n, 128), dtype=scene_arr['hi'][field].dtype)
        out[0::2, 0::2] = scene_arr['low'][field]           # A 掃描第 2k 點
        out[0::2, 1::2] = scene_arr['hi'][field][:, :, 1]   # A 掃描第 2k+1 點
        out[1::2, 0::2] = scene_arr['hi'][field][:, :, 0]   # B 掃描第 2k 點
        out[1::2, 1::2] = scene_arr['hi'][field][:, :, 2]   # B 掃描第 2k+1 點
        return out

    def geo(self, scene, group):
        arr, time = self._data(scene)
        if group == 'LOW':
            lat, lon = arr['low']['lat'], arr['low']['lon']
        else:
            lat, lon = self._hi(arr, 'lat'), self._hi(arr, 'lon')
            step = np.nanmedian(np.diff(time)) if time.size > 1 else 3.798
            time = np.repeat(time, 2)
            time[1::2] += 0.5 * step                            # B 掃描在 A 掃描之後半個週期
        raw_bad = (lat == 0) & (lon == 0)
        lat, lon = clean_latlon(lat / 100. - 90., lon / 100.)
        lat[raw_bad] = np.nan
        lon[raw_bad] = np.nan
        return [Geo(lat, lon, time, (scene.file, group))]

    def tb(self, scene, channel):
        arr, _ = self._data(scene)
        raw = arr['low'][channel.key] if channel.group == 'LOW' else self._hi(arr, channel.key)
        data = raw.astype(np.float32) / np.float32(100.)
        data[~((data > _TB_RANGE[0]) & (data < _TB_RANGE[1]))] = np.nan
        return [data]


READER = SSMIReader
