"""GOSAT-GW AMSR3 Level 1B 亮溫 (JAXA G-Portal，netCDF-4)

檔名  GGWAM3_<YYYYMMDDhhmm><A/D><路徑>_S1BTBB…nc   (A = 升軌、D = 降軌，一個檔約半圈)
亮溫  Tb_Ch<頻段><極化>，uint16 x 0.01 K (填補值 65535)
      06 (6.925)、07 (7.3)、10u (10.25)、10 (10.65)、18 (18.7)、23 (23.8)、36 (36.42) GHz：V、H，每條掃描 243 點
      89A、89B (89.0 GHz 的 A、B 兩組掃描)：V、H，各 486 點
      165 (165.5)、183r3 (183.31±3)、183r7 (183.31±7) GHz：只有 V，243 點
經緯度  每個頻段各有 Latitude_P<頻段> / Longitude_P<頻段> (不需像 AMSR2 由 89A 推算)
時間  ScanTimeTAI93 (1993-01-01 起算、含閏秒)。每個檔前後各多 30 條與相鄰檔重疊的掃描，
      多檔合併時依掃描時間去掉重複。
89 GHz：[swath] merge_ab = true 時 A、B 合併成一個頻道 (89V、89H)；false 時分列為 89V-A、89V-B…
"""
import datetime
import os
import re

import h5py
import numpy as np

from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, tai93_to_unix

_NAME = re.compile(r'^GGWAM3_(\d{12})([AD])(\d{3})_', re.I)
# (頻段代號, 中心頻率 GHz, 偏移, 極化)
_BANDS = (('06', 6.925, None, 'VH'), ('07', 7.3, None, 'VH'), ('10u', 10.25, None, 'VH'), ('10', 10.65, None, 'VH'),
          ('18', 18.7, None, 'VH'), ('23', 23.8, None, 'VH'), ('36', 36.42, None, 'VH'),
          ('165', 165.5, None, 'V'), ('183r3', 183.31, 3., 'V'), ('183r7', 183.31, 7., 'V'))


def _attr(f, key):
    v = f.attrs[key]
    v = v[0] if isinstance(v, np.ndarray) else v
    return v.decode('latin1') if isinstance(v, bytes) else str(v)


def _time(text):
    return datetime.datetime.strptime(text[:19], '%Y-%m-%dT%H:%M:%S')


class AMSR3Reader(SwathReader):
    NAME = 'AMSR3-L1B'
    SENSOR = 'AMSR3'
    SOURCE = 'JAXA'
    DESCRIPTION = 'GOSAT-GW AMSR3 Level 1B 亮溫 (JAXA G-Portal)'
    DESCRIPTION_EN = 'GOSAT-GW AMSR3 Level 1B brightness temperature (JAXA G-Portal)'
    PATTERNS = ('GGWAM3_*1BTBB*.nc', 'GGWAM3_*1BTBB*.h5')

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(os.path.basename(name))
            if not m:
                continue
            path = os.path.join(folder, name)
            with h5py.File(path, 'r') as f:
                start, end = _time(_attr(f, 'ObservationStartDateTime')), _time(_attr(f, 'ObservationEndDateTime'))
            out.append(Scene('GOSAT-GW', 'AMSR3', start, end, file=path, orbit=m.group(3) + m.group(2).upper(),
                             format='L1B'))
        return out

    def channels(self, scene):
        key = 'channels_%s' % self.cfg.swath.merge_ab
        if key in scene.extra:
            return scene.extra[key]
        with h5py.File(scene.file, 'r') as f:
            have = set(f.keys())
        out = []

        def add(tag, freq, offset, pols, group=None, suffix=''):
            for pol in pols:
                name = 'Tb_Ch%s%s' % (tag, pol)
                if name in have:
                    out.append(Channel(band_name(freq, pol, offset, suffix), freq, pol, group or tag, key=(name,)))

        for tag, freq, offset, pols in _BANDS[:7]:
            add(tag, freq, offset, pols)
        if self.cfg.swath.merge_ab:
            for pol in 'VH':
                a, b = 'Tb_Ch89A' + pol, 'Tb_Ch89B' + pol
                if a in have and b in have:
                    out.append(Channel(band_name(89.0, pol), 89.0, pol, '89A+89B', key=(a, b)))
        else:
            add('89A', 89.0, None, 'VH', suffix='-A')
            add('89B', 89.0, None, 'VH', suffix='-B')
        for tag, freq, offset, pols in _BANDS[7:]:
            add(tag, freq, offset, pols)
        scene.extra[key] = out
        return out

    def geo(self, scene, group):
        out = []
        with h5py.File(scene.file, 'r') as f:
            time = tai93_to_unix(f['ScanTimeTAI93'][:])
            for g in group.split('+'):
                lat, lon = clean_latlon(f['Latitude_P' + g][:], f['Longitude_P' + g][:])
                out.append(Geo(lat, lon, time, (scene.file, g)))
        return out

    def tb(self, scene, channel):
        out = []
        with h5py.File(scene.file, 'r') as f:
            for name in channel.key:
                d = f[name]
                raw = d[:]
                scale = float(np.ravel(d.attrs.get('scale_factor', 0.01))[0])
                offset = float(np.ravel(d.attrs.get('add_offset', 0.))[0])
                fill = int(np.ravel(d.attrs.get('_FillValue', 65535))[0])
                data = raw.astype(np.float32) * np.float32(scale) + np.float32(offset)
                data[(raw == fill) | (raw == 0)] = np.nan
                out.append(data)
        return out


READER = AMSR3Reader
