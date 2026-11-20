"""DMSP SSMIS SDR 亮溫 (NOAA CLASS 的 NPR.SDRN.*，大端二進位)

檔名  NPR.SDRN.S<A-D>.D<yyddd>.S<hhmm>.E<hhmm>.B<軌道>.<站>   SA / SB / SC / SD = F16 / F17 / F18 / F19
結構  (2026/11 逐筆解析 15 個樣本檔，讀到的位置都與檔案大小完全吻合)
  檔頭 40 bytes (第 18 byte 起的 int16 為記錄數)，資料自 512 byte 起，每筆記錄補齊到 512 的倍數：
    掃描檔頭 360 bytes：年、日序、時、分、…；之後四種場景各有「掃描時間 (當日毫秒) x N」與「場景數 x N」：
      影像 (imager)       offset 20 起 28 個時間、132 起 28 個場景數   (每條掃描 180 點，間距 12.5 km)
      環境 (environmental) offset 160 起 24 個時間、256 起 24 個場景數  (90 點，25 km；每條掃描都有)
      LAS (低層探空)       offset 280 起 8 個時間、312 起 8 個場景數    (60 點，37.5 km；每 3 條掃描一條)
      UAS (高層探空)       offset 320 起 4 個時間、336 起 4 個場景數    (30 點；此處不讀)
    影像場景 20 bytes：緯度、經度、場景號、地表、降雨旗標、ch8、ch9、ch10、ch11、ch17、ch18
    環境場景：記錄內第 0、2、4… 條每場景 36 bytes，第 1、3、5… 條 18 bytes；前 18 bytes 相同：
              緯度、經度、場景號、海冰、地表、ch12、ch13、ch14、ch15、ch16
    LAS 場景 40 bytes：緯度、經度、ch1~ch7、ch8~ch11、ch18、ch24、…
  緯度、經度為 0.01 度的 int16；亮溫為 0.01 °C 的 int16 (+273.15 -> K)
頻道  影像：ch8 150H、ch9 183.31±6.6H、ch10 183.31±3H、ch11 183.31±1H、ch17 91.655V、ch18 91.655H
      環境：ch12 19.35H、ch13 19.35V、ch14 22.235V、ch15 37H、ch16 37V
      LAS ：ch1~ch7 (50.3、52.8、53.596、54.4、55.5、57.29、59.4 GHz) 與 ch24 (60.79 GHz)，溫度探空用
"""
import datetime
import os
import re

import numpy as np

from satima.core.i18n import tr
from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, to_seconds

_NAME = re.compile(r'^NPR\.SDRN\.S([A-D])\.D(\d{5})\.S(\d{4})\.E(\d{4})\.B(\d+)', re.I)
SATELLITES = {'A': 'F16', 'B': 'F17', 'C': 'F18', 'D': 'F19'}

_IMG = np.dtype([('lat', '>i2'), ('lon', '>i2'), ('scene', '>u2'), ('sfc', 'u1'), ('rain', 'u1'),
                 ('c8', '>i2'), ('c9', '>i2'), ('c10', '>i2'), ('c11', '>i2'), ('c17', '>i2'), ('c18', '>i2')])
_ENV = np.dtype([('lat', '>i2'), ('lon', '>i2'), ('scene', '>u2'), ('ice', 'u1'), ('sfc', 'u1'),
                 ('c12', '>i2'), ('c13', '>i2'), ('c14', '>i2'), ('c15', '>i2'), ('c16', '>i2')])
_ENV_LONG = np.dtype([('head', _ENV), ('extra', 'V18')])          # 記錄內第 0、2、4… 條環境掃描 (36 bytes)
_LAS = np.dtype([('lat', '>i2'), ('lon', '>i2')] + [('c%d' % k, '>i2') for k in (1, 2, 3, 4, 5, 6, 7)] +
                [('c8_11', '>i2', 4), ('c18', '>i2'), ('c24', '>i2'), ('rest', 'V10')])
_UAS_SIZE = 28
# 組 -> (每條掃描的點數, [(欄位, 頻率, 偏移, 極化), ...])
GROUPS = {
    'ENV': (90, [('c13', 19.35, None, 'V'), ('c12', 19.35, None, 'H'), ('c14', 22.235, None, 'V'),
                 ('c16', 37.0, None, 'V'), ('c15', 37.0, None, 'H')]),
    'IMG': (180, [('c17', 91.655, None, 'V'), ('c18', 91.655, None, 'H'), ('c8', 150.0, None, 'H'),
                  ('c11', 183.31, 1., 'H'), ('c10', 183.31, 3., 'H'), ('c9', 183.31, 6.6, 'H')]),
    'LAS': (60, [('c1', 50.3, None, ''), ('c2', 52.8, None, ''), ('c3', 53.596, None, ''), ('c4', 54.4, None, ''),
                 ('c5', 55.5, None, ''), ('c6', 57.29, None, ''), ('c7', 59.4, None, ''), ('c24', 60.79, None, '')]),
}
_TB_RANGE = (40., 400.)         # 這個範圍以外的亮溫視為無效


def file_times(yyddd, start, end):
    """檔名的 D<yyddd>、S<hhmm>、E<hhmm> -> (起, 訖)；跨日時結束時間加一天。"""
    yy = int(yyddd[:2])
    day = datetime.datetime(1900 + yy if yy >= 80 else 2000 + yy, 1, 1) + datetime.timedelta(days=int(yyddd[2:]) - 1)
    a = day + datetime.timedelta(hours=int(start[:2]), minutes=int(start[2:]))
    b = day + datetime.timedelta(hours=int(end[:2]), minutes=int(end[2:]))
    if b < a:
        b += datetime.timedelta(days=1)
    return a, b


def parse(path):
    """整個檔讀成 {組: (結構化陣列 (掃描數, 點數), 各掃描時間 (1970 起算的秒數))}。場景數不足的掃描略過。"""
    with open(path, 'rb') as f:
        buf = f.read()
    nrec = int(np.frombuffer(buf, '>i2', 1, 18)[0])
    scans = {g: [] for g in GROUPS}
    times = {g: [] for g in GROUPS}
    pos = 512
    for _ in range(nrec):
        if pos + 360 > len(buf):
            break
        year, jday = int(np.frombuffer(buf, '>u4', 1, pos + 4)[0]), int(np.frombuffer(buf, '>u2', 1, pos + 8)[0])
        hour = int(buf[pos + 10])
        try:
            day = to_seconds(datetime.datetime(year, 1, 1)) + (jday - 1) * 86400.
        except ValueError:
            day = np.nan
        t_img, n_img = np.frombuffer(buf, '>u4', 28, pos + 20), np.frombuffer(buf, 'u1', 28, pos + 132)
        t_env, n_env = np.frombuffer(buf, '>u4', 24, pos + 160), np.frombuffer(buf, 'u1', 24, pos + 256)
        t_las, n_las = np.frombuffer(buf, '>u4', 8, pos + 280), np.frombuffer(buf, 'u1', 8, pos + 312)
        n_uas = np.frombuffer(buf, 'u1', 4, pos + 336)

        def stamp(ms):
            # 掃描時間是當日毫秒；記錄跨日時 (檔頭在午夜前、掃描在午夜後) 加一天
            sec = ms / 1000.
            return day + sec + (86400. if hour >= 18 and sec < 6 * 3600. else 0.)

        q = pos + 360
        for k in range(28):
            n = int(n_img[k])
            if n == 180:
                scans['IMG'].append(np.frombuffer(buf, _IMG, n, q))
                times['IMG'].append(stamp(float(t_img[k])))
            q += _IMG.itemsize * n
        for k in range(24):
            n = int(n_env[k])
            dt = _ENV_LONG if k % 2 == 0 else _ENV
            if n == 90:
                arr = np.frombuffer(buf, dt, n, q)
                scans['ENV'].append(arr['head'] if k % 2 == 0 else arr)
                times['ENV'].append(stamp(float(t_env[k])))
            q += dt.itemsize * n
        for k in range(8):
            n = int(n_las[k])
            if n == 60:
                scans['LAS'].append(np.frombuffer(buf, _LAS, n, q))
                times['LAS'].append(stamp(float(t_las[k])))
            q += _LAS.itemsize * n
        q += _UAS_SIZE * int(n_uas.sum())
        total = q - pos
        pos += total + (512 - total % 512)
    out = {}
    for g, (npix, _) in GROUPS.items():
        if scans[g]:
            out[g] = (np.stack(scans[g]), np.array(times[g], dtype=np.float64))
    return out


class SSMISReader(SwathReader):
    NAME = 'SSMIS-SDR'
    SENSOR = 'SSMIS'
    SOURCE = 'CLASS'
    DESCRIPTION = 'DMSP F16~F19 SSMIS SDR 亮溫 (NOAA CLASS NPR.SDRN)'
    DESCRIPTION_EN = 'DMSP F16~F19 SSMIS SDR brightness temperature (NOAA CLASS NPR.SDRN)'
    PATTERNS = ('NPR.SDRN.S[A-Da-d].*',)
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
            start, end = file_times(m.group(2), m.group(3), m.group(4))
            out.append(Scene(SATELLITES[m.group(1).upper()], 'SSMIS', start, end, file=os.path.join(folder, name),
                             orbit=m.group(5), format='SDR'))
        return out

    def channels(self, scene):
        return [Channel(band_name(freq, pol, offset), freq, pol, g, key=field)
                for g in ('ENV', 'IMG', 'LAS') for field, freq, offset, pol in GROUPS[g][1]]

    def _data(self, scene):
        if scene.file not in self._parsed:
            while len(self._parsed) >= self.PARSE_CACHE:
                self._parsed.pop(next(iter(self._parsed)))
            self._parsed[scene.file] = parse(scene.file)
        return self._parsed[scene.file]

    def geo(self, scene, group):
        data = self._data(scene)
        if group not in data:
            raise ValueError(tr('%s 沒有 %s 的資料', '%s has no %s data') % (os.path.basename(scene.file), group))
        arr, time = data[group]
        lat, lon = clean_latlon(arr['lat'] / 100., arr['lon'] / 100.)
        return [Geo(lat, lon, time, (scene.file, group))]

    def tb(self, scene, channel):
        arr, _ = self._data(scene)[channel.group]
        data = arr[channel.key].astype(np.float32) / np.float32(100.) + np.float32(273.15)
        data[~((data > _TB_RANGE[0]) & (data < _TB_RANGE[1]))] = np.nan
        return [data]


READER = SSMISReader
