"""MSG SEVIRI (Meteosat-8~11) — EUMETSAT Level 1.5 NetCDF
檔名：W_XX-EUMETSAT-Darmstadt,VIS+IR+HRV+IMAGERY,MSG4+SEVIRI_C_EUMG_20190301120000.nc
變數：ch1~ch11 (VIS006 ~ IR_134 的 10 位元計數值；屬性 scale_factor / add_offset 即定標斜率 / 截距)，
  planned_chan_processing (輻射種類)、channel_data_visir_data_line_validity / _radiometric_quality /
  _geometric_quality (逐行品質)；全域屬性：longitude_of_SSP、equatorial_radius、north/south_polar_radius、
  vis_ir_column/line_dir_grid_step、vis_ir_grid_origin、type_of_earth_model、north/south_most_line、
  east/west_most_pixel、satellite_id、true_repeat_cycle_start_day / _mi_sec。HRV 不在此格式中。
影像方向：第 0 列最南；行的方向依 satpy 為由西向東 (與 HRIT / native 相反)。
  檔案若附有 lon / lat 變數 (舊版 MSG_NC.py 即以此定位)，會以它們檢查導航，方向不符時自動改用由東向西；
  兩種都不符時直接以檔案的經緯度定位 (LonLatNav)。
定標、亮溫、反射率同 SEVIRI-HRIT / SEVIRI-NAT (seviri_common.py)。
完整檔頭：reader.header(scene, 'IR_108') 回傳 {'attributes': 全域屬性, 'variables': {名稱: 屬性}}。
"""
import datetime
import os
import re

import numpy as np

import seviri_common as sc
from satima.core.i18n import tr
from satima.geo.nav import GeosNav, LonLatNav
from satima.geo.reader import GeoReader, Scene

_NAME = re.compile(r'^W_XX-EUMETSAT-Darmstadt,VIS\+IR\+HRV\+IMAGERY,(MSG\d)\+SEVIRI_C_EUMG_(\d{14})\.nc$')


def _attr(ds, *names, default=None):
    for n in names:
        if n in ds.ncattrs():
            return ds.getncattr(n)
    return default


def _num(v):
    if isinstance(v, (bytes, str)):
        s = v.decode() if isinstance(v, bytes) else v
        return int(s, 16) if s.lower().startswith('0x') else float(s)
    return float(np.asarray(v).ravel()[0])


class SeviriNcReader(GeoReader):
    NAME = 'SEVIRI-NC'
    SENSOR = 'SEVIRI'
    DESCRIPTION = 'MSG (Meteosat-8~11) SEVIRI Level 1.5 NetCDF (W_XX-EUMETSAT-Darmstadt,...SEVIRI_C_EUMG_*.nc)'
    PATTERNS = ('W_XX-EUMETSAT-Darmstadt,VIS+IR+HRV+IMAGERY,MSG*+SEVIRI_C_EUMG_*.nc',)
    CHANNELS = {n: sc.make_channel(i + 1) for i, n in enumerate(sc.NAMES[:11])}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._ds = {}
        self._cache = {}

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
            try:
                ds = self._open(path)
            except OSError as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                continue
            satid = int(_num(_attr(ds, 'satellite_id', default=sc.SATID.get(m.group(1), 324))))
            day, ms = _attr(ds, 'true_repeat_cycle_start_day'), _attr(ds, 'true_repeat_cycle_start_mi_sec')
            start = sc.cds(_num(day), _num(ms)) if day is not None and ms is not None else None
            start = start or datetime.datetime.strptime(m.group(2), '%Y%m%d%H%M%S')
            reduced = str(_attr(ds, 'reduced_scanning', default='F')).strip().upper().startswith('T')
            ssp = _num(_attr(ds, 'longitude_of_SSP', default=0.))
            chans = {sc.NAMES[i - 1]: path for i in range(1, 12) if 'ch%d' % i in ds.variables}
            s = Scene(sc.SATELLITES.get(satid, 'Meteosat-?'), sc.nominal_time(start, reduced), 'RSS' if reduced else 'FD',
                      files=chans, region='IODC' if ssp > 30. else '', format='NC', extra={'satid': satid, 'start': start})
            out.append(s)
        return sorted(out, key=lambda s: (s.time, s.satellite))

    def channels(self, scene):
        return [self.CHANNELS[n] for n in sc.NAMES[:11] if n in scene.files]

    def header(self, scene, channel=None):
        ds = self._open(next(iter(scene.files.values())))
        return {'attributes': {k: ds.getncattr(k) for k in ds.ncattrs()},
                'variables': {n: {k: v.getncattr(k) for k in v.ncattrs()} for n, v in ds.variables.items()}}

    def _shape(self, ds):
        v = ds.variables['ch9' if 'ch9' in ds.variables else next(n for n in ds.variables if re.match(r'ch\d+$', n))]
        return v.shape[-2], v.shape[-1]

    def nav(self, scene, channel):
        path = scene.files[channel]
        if ('nav', path) in self._cache:
            return self._cache[('nav', path)]
        ds = self._open(path)
        shape = self._shape(ds)
        ssp = _num(_attr(ds, 'longitude_of_SSP', default=0.))
        a = _num(_attr(ds, 'equatorial_radius')) * 1000.
        b = (_num(_attr(ds, 'north_polar_radius')) + _num(_attr(ds, 'south_polar_radius'))) / 2. * 1000.
        cs = _num(_attr(ds, 'vis_ir_column_dir_grid_step')) * 1000. / sc.H
        ls = _num(_attr(ds, 'vis_ir_line_dir_grid_step')) * 1000. / sc.H
        model = int(_num(_attr(ds, 'type_of_earth_model', default='0x02')))
        south = int(_num(_attr(ds, 'south_most_line', default=1)))
        east = int(_num(_attr(ds, 'east_most_pixel', default=1)))
        west = int(_num(_attr(ds, 'west_most_pixel', default=3712)))
        shift = 0.5 if model == 1 else 0.
        y0 = (south - sc.CENTER_VISIR - shift) * ls
        west_first = GeosNav(shape, (sc.CENTER_VISIR - west + shift) * cs, cs, y0, ls, ssp, sc.H, a, b, 'y')
        east_first = sc.geos_nav(shape, south, east, ls, cs, ssp, a, b, model)
        nav = west_first
        if 'lon' in ds.variables and 'lat' in ds.variables:
            nav = self._check(ds, shape, ((tr('由西向東', 'west to east'), west_first),
                                          (tr('由東向西', 'east to west'), east_first)))
        self._cache[('nav', path)] = nav
        return nav

    @staticmethod
    def _check(ds, shape, candidates):
        """以檔案的 lon / lat 變數檢查導航 (取內部數點)；都不符時以檔案經緯度定位。"""
        lon_v, lat_v = ds.variables['lon'], ds.variables['lat']
        rows = np.linspace(shape[0] * 0.2, shape[0] * 0.8, 7).astype(int)
        cols = np.linspace(shape[1] * 0.2, shape[1] * 0.8, 7).astype(int)
        R, C = np.meshgrid(rows, cols, indexing='ij')
        flon = np.array([[float(lon_v[r, c]) for c in cols] for r in rows])
        flat = np.array([[float(lat_v[r, c]) for c in cols] for r in rows])
        good = np.isfinite(flon) & np.isfinite(flat) & (np.abs(flon) <= 360) & (np.abs(flat) <= 90)
        for tag, nav in candidates:
            lon, lat = nav.to_lonlat(R, C)
            d = np.hypot(((lon - flon + 180.) % 360.) - 180., lat - flat)[good & np.isfinite(lon)]
            if d.size and np.nanmedian(d) < 0.05:
                if tag != candidates[0][0]:
                    print(tr('  [SEVIRI-NC] 依檔案經緯度，行的方向為%s',
                             '  [SEVIRI-NC] according to the lon/lat in the file, the columns run %s') % tag)
                return nav
        print(tr('  [SEVIRI-NC] 投影參數與檔案經緯度不符，改用檔案的經緯度定位',
                 '  [SEVIRI-NC] the projection parameters do not match the lon/lat in the file; '
                 'geolocating with the lon/lat in the file instead'))
        return LonLatNav(np.asarray(lon_v[:], dtype=np.float64), np.asarray(lat_v[:], dtype=np.float64))

    def satellite_position(self, scene):
        ds = self._open(next(iter(scene.files.values())))
        return _num(_attr(ds, 'longitude_of_SSP', default=0.)), sc.H / 1000.

    def _calibration(self, scene, channel):
        key = ('cal', scene.files[channel], channel)
        if key not in self._cache:
            ds = self._open(scene.files[channel])
            n = sc.NUMBER[channel]
            v = ds.variables['ch%d' % n]
            rad_type = 2
            if 'planned_chan_processing' in ds.variables:
                rad_type = int(np.asarray(ds.variables['planned_chan_processing'][:]).ravel()[n - 1])
            c = sc.Calibration(scene.extra['satid'], n, _num(v.getncattr('scale_factor')),
                               _num(v.getncattr('add_offset')), rad_type=rad_type,
                               mode='nominal' if str(self.cfg.seviri.calibration).lower() == 'gsics'
                               else self.cfg.seviri.calibration, time=scene.extra['start'])
            self._cache[key] = c
        return self._cache[key]

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        ds = self._open(scene.files[channel])
        n = sc.NUMBER[channel]
        v = ds.variables['ch%d' % n]
        counts = np.asarray(v[max(r0, 0):r1, max(c0, 0):c1])
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        vals = self._calibration(scene, channel)(counts)
        names = ['channel_data_visir_data_line_%s' % k for k in ('validity', 'radiometric_quality', 'geometric_quality')]
        if all(k in ds.variables for k in names):
            q = [np.asarray(ds.variables[k][max(r0, 0):r1, n - 1]) for k in names]
            vals[sc.bad_lines(*q)] = np.nan
        out[max(r0, 0) - r0:max(r0, 0) - r0 + vals.shape[0], max(c0, 0) - c0:max(c0, 0) - c0 + vals.shape[1]] = vals
        return out


READER = SeviriNcReader
