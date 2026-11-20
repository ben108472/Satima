"""MFG (Meteosat-2~7) MVIRI — EUMETSAT NetCDF (舊版 MFG.py 的格式)
檔名：W_XX-EUMETSAT-Darmstadt,VIS+IR+IMAGERY,MET7+MVIRI_C_EUMS_20050101120000.nc
變數 (同舊版)：ch1 = VIS、ch2 = IR、ch3 = WV (計數值，屬性 scale_factor / add_offset)，
  lon_vis / lat_vis (VIS 格點)、lon / lat (IR / WV 格點)。
定位：直接以檔案的經緯度 (LonLatNav)。
定標 (同舊版)：輻射 L = scale_factor x (計數值 - add_offset)；
  IR / WV 亮溫 T = B / (ln L - A)，A、B 為 EUMETSAT 公布的各衛星係數 (舊版 _sate_param.py)；
  (舊版 Meteosat-5~7 另以 conversion_table/MFG 的表內插，但表的欄位意義未註明且舊版 Meteosat-6 水氣取錯欄，
   此處一律用 A、B 公式。)
  VIS 沒有絕對校正，同舊版以 L / 255 顯示。
注意：沒有實際樣本，依舊版程式移植，尚未驗證 (有檔案後請告知以便確認變數名稱、影像方向與定標)。
"""
import datetime
import os
import re

import numpy as np

from satima.core.i18n import tr
from satima.geo.nav import LonLatNav
from satima.geo.reader import Channel, GeoReader, Scene

_NAME = re.compile(r'^W_XX-EUMETSAT-Darmstadt,VIS\+IR\+IMAGERY,MET(\d+)\+MVIRI_C_EUMS_(\d{14})\.nc$')
# EUMETSAT 亮溫係數 T = B / (ln L - A) (舊版 _sate_param.py)
_AB = {2: ((6.7631624, -1264.2322), (9.4199790, -2177.3002)), 3: ((6.9028900, -1259.9131), (9.3665346, -2164.6560)),
       4: ((6.7160127, -1269.4460), (9.0881146, -2253.0677)), 5: ((6.7210227, -1269.3131), (9.2208773, -2263.6804)),
       6: ((6.7473998, -1264.4411), (9.1090648, -2261.2885)), 7: ((6.9590805, -1252.6850), (9.2347942, -2221.2987))}
_VAR = {'VIS': 'ch1', 'IR': 'ch2', 'WV': 'ch3'}
_BANDS = (('VIS', 0.7, 2.25, 'vis', 'Visible'), ('IR', 11.5, 4.5, 'ir', 'LW Window'), ('WV', 6.4, 4.5, 'wv', 'Water Vapor'))


class MviriNcReader(GeoReader):
    NAME = 'MVIRI-NC'
    SENSOR = 'MVIRI'
    DESCRIPTION = 'MFG (Meteosat-2~7) MVIRI EUMETSAT NetCDF (W_XX-EUMETSAT-Darmstadt,VIS+IR+IMAGERY,MET*+MVIRI_C_EUMS_*.nc)'
    PATTERNS = ('W_XX-EUMETSAT-Darmstadt,VIS+IR+IMAGERY,MET*+MVIRI_C_EUMS_*.nc',)
    CHANNELS = {b[0]: Channel(*b[:4], band=b[0], alias=b[4]) for b in _BANDS}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._ds = {}
        self._navs = {}

    def _open(self, path):
        if path not in self._ds:
            import netCDF4
            ds = netCDF4.Dataset(path)
            ds.set_auto_maskandscale(False)
            self._ds[path] = ds
        return self._ds[path]

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            path = os.path.join(folder, name)
            sat = int(m.group(1))
            try:
                ds = self._open(path)
            except OSError as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                continue
            files = {c: path for c, v in _VAR.items() if v in ds.variables}
            out.append(Scene('Meteosat-%d' % sat, datetime.datetime.strptime(m.group(2), '%Y%m%d%H%M%S'), 'FD',
                             files=files, format='NC', extra={'sat': sat}))
        return sorted(out, key=lambda s: (s.time, s.satellite))

    def header(self, scene, channel=None):
        ds = self._open(next(iter(scene.files.values())))
        return {'attributes': {k: ds.getncattr(k) for k in ds.ncattrs()},
                'variables': {n: {k: v.getncattr(k) for k in v.ncattrs()} for n, v in ds.variables.items()}}

    def nav(self, scene, channel):
        key = (scene.files[channel], channel == 'VIS')
        if key not in self._navs:
            ds = self._open(scene.files[channel])
            lo, la = ('lon_vis', 'lat_vis') if channel == 'VIS' else ('lon', 'lat')
            lon = np.asarray(ds.variables[lo][:], dtype=np.float64)
            lat = np.asarray(ds.variables[la][:], dtype=np.float64)
            bad = ~np.isfinite(lon) | ~np.isfinite(lat) | (np.abs(lat) > 90.) | (np.abs(lon) > 360.)
            lon[bad], lat[bad] = np.nan, np.nan
            self._navs[key] = LonLatNav(lon, lat)
        return self._navs[key]

    def satellite_position(self, scene):
        nav = self.nav(scene, 'IR' if 'IR' in scene.files else next(iter(scene.files)))
        r, c = nav.shape[0] // 2, nav.shape[1] // 2
        lon, _ = nav.to_lonlat(np.array([r]), np.array([c]))
        return float(lon[0]), 35786.

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        ds = self._open(scene.files[channel])
        v = ds.variables[_VAR[channel]]
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        a0, b0 = max(r0, 0), max(c0, 0)
        cnt = np.asarray(v[a0:r1, b0:c1]).astype(np.float64)
        rad = float(v.getncattr('scale_factor')) * (cnt - float(v.getncattr('add_offset')))
        if channel == 'VIS':
            vals = np.where(rad > 0, rad / 255., 0.)
        else:
            a, b = _AB.get(scene.extra['sat'], _AB[7])[0 if channel == 'IR' else 1]
            with np.errstate(divide='ignore', invalid='ignore'):
                vals = np.where(rad > 0, b / (np.log(rad) - a) - 273.15, np.nan)
        out[a0 - r0:a0 - r0 + vals.shape[0], b0 - c0:b0 - c0 + vals.shape[1]] = vals
        return out


READER = MviriNcReader
