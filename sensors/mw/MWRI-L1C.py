"""風雲三號 MWRI 軌道 L1c 產品 (NSMC 數據同化支撐產品 DAS，二進位)

檔名  FY3D_MWRI{A|D}_ORBT_L2_MIP_MLT_NUL_<YYYYMMDD>_<hhmm>_010KM_MS.L1c   (MWRIA = 升軌、MWRID = 降軌，各約半圈)
格式  NSMC「FY-3D MWRI 轨道 L1C (升轨 / 降轨) 格式说明」(FY-3D_Product_Instruction_MWRI_DAS_Asc/Desc.pdf)：
      每個觀測點一筆 152 bytes (little-endian)，依掃描線、點序排列，每條掃描 266 點：
        char[12] 衛星名 'FY-3D'、uint32 衛星代碼 (4)、儀器代碼 (36)、掃描線序號、點序 (1~266)、年月日時分秒 (各 uint32)、
        int32 緯度、經度 (x100 度)、海陸標記 (1 陸、2 陸地水體、3 海、5 交界)、高度 (x100 m)、
        局地天頂角、局地方位角、太陽天頂角、太陽方位角 (x100 度)、衛星高度 (x100 km)、品質標記、
        亮溫 x 10 (x100 K：10.65V、10.65H、18.7V、18.7H、23.8V、23.8H、36.5V、36.5H、89V、89H)、
        雲水、降水、海面風速、海面溫度、地表類型 (反演產品，這裡不用)
      填補值 999999。V / H 順序以海面中位數確認 (V 較高)。
時間  每條掃描線的年月日時分秒 (到秒)；同一秒的掃描線依序內插出小數秒。
"""
import datetime
import os
import re

import numpy as np

from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, to_seconds

_NAME = re.compile(r'^(FY3[A-Z])_MWRI([AD])_ORBT_L2_MIP_\w+?_(\d{8})_(\d{4})_', re.I)
FOVS = 266
FILL = 999999
_FREQS = (10.65, 18.7, 23.8, 36.5, 89.0)
RECORD = np.dtype([('platform', 'S12'), ('sat_id', '<u4'), ('inst_id', '<u4'), ('scan_line', '<u4'), ('fov', '<u4'),
                   ('year', '<u4'), ('mon', '<u4'), ('day', '<u4'), ('hour', '<u4'), ('min', '<u4'), ('sec', '<u4'),
                   ('lat', '<i4'), ('lon', '<i4'), ('surface', '<i4'), ('height', '<i4'), ('lza', '<i4'),
                   ('laz', '<i4'), ('sza', '<i4'), ('saz', '<i4'), ('sat_alt', '<i4'), ('quality', '<i4'),
                   ('bt', '<i4', (10,)), ('clw', '<i4'), ('mrr', '<i4'), ('sws', '<i4'), ('sst', '<i4'),
                   ('landcover', '<i4')])


def _read(path, fields=None):
    """整個檔案 -> (掃描數, 266) 的結構陣列。"""
    size = os.path.getsize(path)
    if size % (RECORD.itemsize * FOVS):
        raise ValueError('%s: size %d is not a multiple of %d x %d' % (path, size, RECORD.itemsize, FOVS))
    a = np.fromfile(path, RECORD).reshape(-1, FOVS)
    if not (a['fov'][0] == np.arange(1, FOVS + 1)).all():
        raise ValueError('%s: unexpected FOV order' % path)
    return a


def _scan_seconds(a):
    """每條掃描線的時間 (1970 起的秒)；記錄只到秒，同一秒內的掃描線依序平均分配。"""
    first = a[:, 0]
    out = np.full(first.size, np.nan)
    for k in range(first.size):
        try:
            out[k] = to_seconds(datetime.datetime(int(first['year'][k]), int(first['mon'][k]), int(first['day'][k]),
                                                  int(first['hour'][k]), int(first['min'][k]), int(first['sec'][k])))
        except ValueError:
            pass
    ok = np.isfinite(out)
    if ok.sum() > 2:                # 秒的整數部分重複：以首尾的平均掃描週期內插 (MWRI 約 1.8 秒)
        k = np.arange(out.size)
        period = (out[ok][-1] - out[ok][0]) / max(1, k[ok][-1] - k[ok][0])
        if period > 0:
            fit = out[ok][0] + (k - k[ok][0]) * period
            out = np.where(ok & (np.abs(out - fit) < 2.), fit, out)
    return out


class MWRIL1cReader(SwathReader):
    NAME = 'MWRI-L1C'
    SENSOR = 'MWRI'
    SOURCE = 'NSMC'
    DESCRIPTION = '風雲三號 (FY-3D) MWRI 軌道 L1c 亮溫 (NSMC 數據同化支撐產品，二進位 .L1c)'
    DESCRIPTION_EN = 'Fengyun-3 (FY-3D) MWRI orbital L1c brightness temperature (NSMC DAS product, binary .L1c)'
    PATTERNS = ('FY3?_MWRI?_ORBT_L2_MIP_*.L1c',)

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(os.path.basename(name))
            if not m:
                continue
            path = os.path.join(folder, name)
            try:
                size = os.path.getsize(path)
                if size % (RECORD.itemsize * FOVS):
                    continue
                n = size // RECORD.itemsize
                head = np.fromfile(path, RECORD, count=1)[0]
                with open(path, 'rb') as fh:
                    fh.seek((n - 1) * RECORD.itemsize)
                    tail = np.frombuffer(fh.read(RECORD.itemsize), RECORD)[0]
                t = [datetime.datetime(int(r['year']), int(r['mon']), int(r['day']), int(r['hour']), int(r['min']),
                                       int(r['sec'])) for r in (head, tail)]
            except (OSError, ValueError):
                continue
            sat = 'FY-3' + m.group(1)[-1].upper()
            out.append(Scene(sat, 'MWRI', t[0], t[1], file=path, orbit=m.group(2).upper(), format='L1C'))
        return out

    def channels(self, scene):
        return [Channel(band_name(freq, pol), freq, pol, 'all', key=2 * k + i)
                for k, freq in enumerate(_FREQS) for i, pol in enumerate('VH')]

    def _data(self, scene):
        key = '_l1c'
        if scene.extra.get(key) is None:
            scene.extra[key] = _read(scene.file)
        return scene.extra[key]

    def geo(self, scene, group):
        a = self._data(scene)
        lat = np.where(a['lat'] == FILL, np.nan, a['lat'] / 100.)
        lon = np.where(a['lon'] == FILL, np.nan, a['lon'] / 100.)
        lat, lon = clean_latlon(lat, lon)
        return [Geo(lat, lon, _scan_seconds(a), (scene.file, 'all'))]

    def tb(self, scene, channel):
        raw = self._data(scene)['bt'][:, :, channel.key]
        data = raw.astype(np.float32) / np.float32(100.)
        data[(raw == FILL) | (raw <= 0) | (raw > 40000)] = np.nan
        return [data]


READER = MWRIL1cReader
