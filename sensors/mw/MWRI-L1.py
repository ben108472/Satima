"""風雲三號 MWRI Level 1 亮溫 (NSMC，HDF5)

檔名  FY3D_MWRIA_GBAL_L1_<YYYYMMDD>_<hhmm>_010KM_MS.HDF   (MWRIA = 升軌、MWRID = 降軌)
亮溫  Calibration/EARTH_OBSERVE_BT_10_to_89GHz (10, 掃描, 像素)，int16：亮溫 K = 值 x Slope (0.01) + Intercept (327.68)
      頻道順序：10.65V、10.65H、18.7V、18.7H、23.8V、23.8H、36.5V、36.5H、89V、89H
      (2026/11 以熱帶海面的中位數確認 V 在前：10.65V 167 K、10.65H 90 K…)
經緯度  Geolocation/Latitude、Longitude (所有頻道共用)
時間  Calibration/Scan_Daycnt (2000-01-01 12:00 UTC 起算的日數)、Scan_Mscnt[:, 0] (該日 12:00 起算的毫秒)
FY-3B / FY-3C 的檔案沒有 Calibration / Geolocation 群組 (資料集在最上層)，這裡兩種位置都會找；
FY-3F MWRI-II、FY-3G MWRI-RM 是不同的格式，尚未支援。
"""
import datetime
import os
import re

import h5py
import numpy as np

from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, to_seconds

_NAME = re.compile(r'^(FY3[A-E])_MWRI([AD])_GBAL_L1_(\d{8})_(\d{4})_', re.I)
_FREQS = (10.65, 18.7, 23.8, 36.5, 89.0)
_BT = 'EARTH_OBSERVE_BT_10_to_89GHz'
_EPOCH = to_seconds(datetime.datetime(2000, 1, 1, 12))


def _find(f, name, groups):
    for g in groups:
        path = '%s/%s' % (g, name) if g else name
        if path in f:
            return f[path]
    raise KeyError(name)


def _attr(f, key, default=''):
    v = f.attrs.get(key, default)
    v = v[0] if isinstance(v, np.ndarray) and v.size else v
    return v.decode('latin1').strip() if isinstance(v, bytes) else str(v).strip()


class MWRIReader(SwathReader):
    NAME = 'MWRI-L1'
    SENSOR = 'MWRI'
    SOURCE = 'NSMC'
    DESCRIPTION = '風雲三號 (FY-3B/C/D) MWRI Level 1 亮溫 (NSMC)'
    DESCRIPTION_EN = 'Fengyun-3 (FY-3B/C/D) MWRI Level 1 brightness temperature (NSMC)'
    PATTERNS = ('FY3?_MWRI?_GBAL_L1_*_010KM_MS.HDF',)

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(os.path.basename(name))
            if not m:
                continue
            path = os.path.join(folder, name)
            start = datetime.datetime.strptime(m.group(3) + m.group(4), '%Y%m%d%H%M')
            end = start + datetime.timedelta(minutes=55)
            try:
                with h5py.File(path, 'r') as f:
                    a = _attr(f, 'Observing Beginning Date') + ' ' + _attr(f, 'Observing Beginning Time')[:8]
                    b = _attr(f, 'Observing Ending Date') + ' ' + _attr(f, 'Observing Ending Time')[:8]
                    start = datetime.datetime.strptime(a, '%Y-%m-%d %H:%M:%S')
                    end = datetime.datetime.strptime(b, '%Y-%m-%d %H:%M:%S')
            except (KeyError, ValueError):
                pass                                # 檔頭沒有或格式不同：用檔名的時間
            sat = m.group(1).upper()
            out.append(Scene('FY-3' + sat[-1], 'MWRI', start, end, file=path, orbit=m.group(2).upper(), format='L1'))
        return out

    def channels(self, scene):
        return [Channel(band_name(freq, pol), freq, pol, 'all', key=2 * k + i)
                for k, freq in enumerate(_FREQS) for i, pol in enumerate('VH')]

    def geo(self, scene, group):
        with h5py.File(scene.file, 'r') as f:
            lat = _find(f, 'Latitude', ('Geolocation', ''))[:]
            lon = _find(f, 'Longitude', ('Geolocation', ''))[:]
            day = _find(f, 'Scan_Daycnt', ('Calibration', ''))[:].astype(np.float64).reshape(-1)
            ms = _find(f, 'Scan_Mscnt', ('Calibration', ''))[:].astype(np.float64)
        ms = ms[:, 0] if ms.ndim == 2 else ms
        time = np.where((day > 0) & (ms >= 0), _EPOCH + day * 86400. + ms / 1000., np.nan)
        lat, lon = clean_latlon(lat, lon)
        return [Geo(lat, lon, time, (scene.file, 'all'))]

    def tb(self, scene, channel):
        with h5py.File(scene.file, 'r') as f:
            d = _find(f, _BT, ('Calibration', ''))
            raw = d[channel.key]
            slope = float(np.ravel(d.attrs.get('Slope', 0.01))[0])
            inter = float(np.ravel(d.attrs.get('Intercept', 327.68))[0])
        data = raw.astype(np.float32) * np.float32(slope) + np.float32(inter)
        data[(raw <= -32766) | ~(data > 0.)] = np.nan
        return [data]


READER = MWRIReader
