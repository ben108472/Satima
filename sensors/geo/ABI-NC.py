"""GOES-16~19 ABI — Level 1b Radiances (NetCDF)

檔名：OR_ABI-L1b-RadC-M6C01_G19_s20263121801183_e20263121803557_c20263121803582.nc
  觀測範圍 F (全圓盤)、C (CONUS / PACUS)、M1 / M2 (mesoscale，地圖範圍自動帶入)
  同一次掃描的各波段起始時間 (s) 相同，結束 (e) 與產製 (c) 時間各波段不同。
導航：檔案中的 x / y 掃描角 (弧度) 與 goes_imager_projection (sweep = x)。
"""
import datetime
import os
import re

import numpy as np
from netCDF4 import Dataset

from satima.geo.nav import GeosNav
from satima.geo.reader import Channel, FileCache, GeoReader, Scene, to_float

_NAME = re.compile(r'OR_ABI-L1b-Rad(F|C|M1|M2)-M\dC(\d\d)_G(\d\d)_s(\d{13})\d_e\d{14}_c\d{14}\.nc$')
_AREA = {'F': 'FD', 'C': 'CONUS', 'M1': 'MESO1', 'M2': 'MESO2'}
_WAVE = (0.47, 0.64, 0.86, 1.37, 1.6, 2.2, 3.9, 6.2, 6.9, 7.3, 8.4, 9.6, 10.3, 11.2, 12.3, 13.3)   # NOAA 標稱值
_RES = {1: 1.0, 2: 0.5, 3: 1.0, 5: 1.0}
_KINDS = {1: 'vis', 2: 'vis', 3: 'nir', 4: 'nir', 5: 'nir', 6: 'nir', 7: 'swir', 8: 'wv', 9: 'wv', 10: 'wv'}
# NOAA 官方代稱 (簡寫)
_ALIAS = ('Blue', 'Red', 'Veggie', 'Cirrus', 'Snow/Ice', 'Cloud Particle Size', 'Shortwave Window',
          'Upper-level WV', 'Mid-level WV', 'Lower-level WV', 'Cloud-top Phase', 'Ozone',
          'Clean LW Window', 'LW Window', 'Dirty LW Window', 'CO2 LW')


class AbiNcReader(GeoReader):
    NAME = 'ABI-NC'
    SENSOR = 'ABI'
    DESCRIPTION = 'GOES-16~19 ABI (L1b Radiances NetCDF)'
    PATTERNS = ('OR_ABI-L1b-Rad*-M*C*_G*_s*_e*_c*.nc',)
    CHANNELS = {'C%02d' % b: Channel('C%02d' % b, _WAVE[b - 1], _RES.get(b, 2.0), _KINDS.get(b, 'ir'),
                                     'BAND-%02d' % b, _ALIAS[b - 1])
                for b in range(1, 17)}
    # ABI 沒有綠光：以藍、紅合成 0.51 µm 綠光再混入近紅外 (CIRA 作法，見 satima/truecolor/composite.py)
    TRUECOLOR = {'TrueColor': dict(red='C02', green=None, blue='C01', nir='C03', ir='C13', swir='C07',
                                   lo_blur=(0.4, 0.4))}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._ds = FileCache(Dataset)

    def scenes(self, folder='.'):
        groups = {}
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            area, band, sat, start = m.groups()
            key = (sat, area, start)
            if key not in groups:
                t = datetime.datetime.strptime(start, '%Y%j%H%M%S')
                groups[key] = Scene('GOES-%d' % int(sat), t, _AREA[area], auto_range=area.startswith('M'))
            groups[key].files['C' + band] = os.path.join(folder, name)
        return sorted(groups.values(), key=lambda s: (s.time, s.area))

    def _dataset(self, path):
        return self._ds.get(path)

    def nav(self, scene, channel):
        ds = self._dataset(scene.files[channel])
        x, y = ds.variables['x'], ds.variables['y']
        p = ds.variables['goes_imager_projection']
        return GeosNav((y.size, x.size), *self._axis(x), *self._axis(y),
                       p.longitude_of_projection_origin, p.perspective_point_height,
                       p.semi_major_axis, p.semi_minor_axis, p.sweep_angle_axis)

    @staticmethod
    def _axis(v):
        """掃描角軸 (起點, 間距)：由打包整數與 scale_factor / add_offset 以倍精度計算。
        不可用解包後的 float32 相減求間距，誤差會隨行號累積 (CONUS 東緣可達半個像素)。"""
        v.set_auto_maskandscale(False)
        i0, i1 = (float(k) for k in v[:2])
        v.set_auto_maskandscale(True)
        scale, offset = float(v.scale_factor), float(v.add_offset)
        return offset + scale * i0, scale * (i1 - i0)

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        ds = self._dataset(scene.files[channel])
        rad = to_float(ds.variables['Rad'][r0:r1, c0:c1], np.float64)
        v = ds.variables
        if self.CHANNELS[channel].reflective:
            out = rad * float(v['kappa0'][:])
        else:
            fk1, fk2, bc1, bc2 = (float(v[k][:]) for k in ('planck_fk1', 'planck_fk2', 'planck_bc1', 'planck_bc2'))
            with np.errstate(invalid='ignore', divide='ignore'):
                out = (fk2 / np.log(fk1 / rad + 1) - bc1) / bc2 - 273.15
        return out.astype(np.float32)

    def satellite_position(self, scene):
        p = self._dataset(next(iter(scene.files.values()))).variables['goes_imager_projection']
        return float(p.longitude_of_projection_origin), float(p.perspective_point_height) / 1000.


READER = AbiNcReader
