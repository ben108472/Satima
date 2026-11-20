"""MFG (Meteosat-2~7) MVIRI — EUMETSAT / FIDUCEO 基礎氣候資料紀錄 (FCDR) Level 1.5 NetCDF
檔名：MVIRI_FCDR-EASY_L15_MET7-E5700_201604172330_201604180000_0200.nc (或 Data Store 的同名 .zip，不必解開)、
      FIDUCEO_FCDR_L15_MVIRI_MET7-57.0_200607060600_200607060630_EASY_v2.6_fv3.1.nc；FULL 版亦可。
  檔名的 E5700 = 影像的投影經度 57.00 度；時間為掃描起訖 (例如 23:30~00:00 = 舊制的 00:00 時段)。
內容：VIS 5000x5000 (toa_bidirectional_reflectance_vis；FULL 版為 count_vis 與定標係數)，
  IR、WV 2500x2500 (count_ir / count_wv，填充 255；定標 a_*、b_*、亮溫係數 bt_a_*、bt_b_*)，
  逐像素時間 time (IR 格點)、太陽 / 衛星角度 (500x500 tie points)、品質旗標。
導航 (FCDR 的 PUG 第 2 節、舊版 MFG_FCDR.py 的 regfeo 與 satpy 相同)：標準靜止衛星投影 (sweep y)，
  視場 18 度等角分成 N 格 (VIS 5000、IR 2500)，中心 N/2 + 0.5；第 0 列最南、第 0 行最東；
  地球 a = 6378.140、b = 6356.755 km，衛星距地心 42164 km；投影經度取自檔名。
定標：IR / WV 輻射 L = a + b x 計數值，亮溫 T = bt_b / (ln L - bt_a)；
  VIS：FCDR 提供的是雙向反射因子 (已除以 cos(太陽天頂角))，乘回 cos(太陽天頂角) 成為與其他衛星一致的 TOA 反射率
  (太陽天頂角以每個像素的觀測時間計算；[visible] sun_zenith_correction 開啟時再由繪圖程式統一除以 cos)。
  FULL 版的 VIS：L = (計數值 - 太空計數值) x (a0 + a1 y + a2 y^2)，反射率 = π L d^2 / E (y = 發射後年數)。
注意：部分檔案的時間變數與屬性比檔名晚一天 (例如 MET7 2016/04/17 23:30 檔的時間變數為 04/18)，
  以檔名為準 (來源 L1.0 檔名 20160418000000 即 04/17 23:30 開始的時段)；計算太陽角時把時間變數平移整數天對齊。
完整檔頭：reader.header(scene) 回傳 {'attributes': 全域屬性, 'variables': {名稱: 屬性}}。
"""
import datetime
import os
import re
import zipfile

import numpy as np

from satima.core.i18n import tr
from satima.geo.nav import GeosNav
from satima.geo.reader import Channel, GeoReader, Scene

_NAMES = (re.compile(r'^MVIRI_FCDR-(?P<kind>EASY|FULL)_L15_MET(?P<sat>\d)-E(?P<lon>\d{4})_(?P<t0>\d{12})_(?P<t1>\d{12})_\d{4}\.(nc|zip)$', re.I),
          re.compile(r'^FIDUCEO_FCDR_L15_MVIRI_MET(?P<sat>\d)-(?P<lon>[\d.]+)_(?P<t0>\d{12})_(?P<t1>\d{12})_(?P<kind>EASY|FULL)_.*\.nc$', re.I))
_REQ, _RPOL, _DIST = 6378140.0, 6356755.0, 42164000.0
_FOV = 18.0
_BANDS = (('VIS', 0.7, 2.25, 'vis', 'Visible'), ('IR', 11.5, 4.5, 'ir', 'LW Window'), ('WV', 6.4, 4.5, 'wv', 'Water Vapor'))


def _match(name):
    for p in _NAMES:
        m = p.match(name)
        if m:
            return m
    return None


class MviriFcdrReader(GeoReader):
    NAME = 'MVIRI-FCDR'
    SENSOR = 'MVIRI'
    DESCRIPTION = 'MFG (Meteosat-2~7) MVIRI FCDR Level 1.5 NetCDF (MVIRI_FCDR-EASY/FULL_*.nc / .zip)'
    PATTERNS = ('MVIRI_FCDR-*_L15_MET*.nc', 'MVIRI_FCDR-*_L15_MET*.zip', 'FIDUCEO_FCDR_L15_MVIRI_MET*.nc')
    CHANNELS = {b[0]: Channel(*b[:4], band=b[0], alias=b[4]) for b in _BANDS}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._ds = {}
        self._cache = {}

    def _open(self, path):
        if path not in self._ds:
            import netCDF4
            if path.lower().endswith('.zip'):
                with zipfile.ZipFile(path) as z:
                    member = next(n for n in z.namelist() if n.lower().endswith('.nc'))
                    ds = netCDF4.Dataset(member, memory=z.read(member))
            else:
                ds = netCDF4.Dataset(path)
            ds.set_auto_maskandscale(False)
            self._ds[path] = ds
        return self._ds[path]

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _match(name)
            if not m:
                continue
            lon = m.group('lon')
            proj_lon = float(lon) if '.' in lon else float(lon) / 100.
            t0 = datetime.datetime.strptime(m.group('t0'), '%Y%m%d%H%M')
            path = os.path.join(folder, name)
            s = Scene('Meteosat-%s' % m.group('sat'), t0, 'FD', files={b[0]: path for b in _BANDS},
                      region='IODC' if proj_lon > 30. else '', format='FCDR ' + m.group('kind').upper(),
                      extra={'proj_lon': proj_lon, 't1': datetime.datetime.strptime(m.group('t1'), '%Y%m%d%H%M')})
            out.append(s)
        return sorted(out, key=lambda s: (s.time, s.satellite))

    def header(self, scene, channel=None):
        ds = self._open(next(iter(scene.files.values())))
        return {'attributes': {k: ds.getncattr(k) for k in ds.ncattrs()},
                'variables': {n: {k: v.getncattr(k) for k in v.ncattrs()} for n, v in ds.variables.items()}}

    # ---- 導航 --------------------------------------------------------------
    def nav(self, scene, channel):
        n = 5000 if channel == 'VIS' else 2500
        key = ('nav', scene.extra['proj_lon'], n)
        if key not in self._cache:
            step = np.deg2rad(_FOV / n)
            half = n / 2. - 0.5
            self._cache[key] = GeosNav((n, n), half * step, -step, -half * step, step, scene.extra['proj_lon'],
                                       _DIST - _REQ, _REQ, _RPOL, sweep='y')
        return self._cache[key]

    def satellite_position(self, scene):
        return scene.extra['proj_lon'], (_DIST - _REQ) / 1000.

    # ---- 讀取與定標 --------------------------------------------------------
    @staticmethod
    def _scalar(ds, name):
        return float(np.asarray(ds.variables[name][:]).ravel()[0])

    def _pixel_time(self, scene, ds, r0, r1, c0, c1, factor):
        """視窗內每個像素的觀測時間 (datetime64)；時間變數與檔名差整數天時平移對齊。"""
        v = ds.variables['time']
        off = float(v.getncattr('add_offset')) if 'add_offset' in v.ncattrs() else 0.
        fill = int(v.getncattr('_FillValue')) if '_FillValue' in v.ncattrs() else None
        rows = np.arange(r0, r1) // factor
        cols = np.arange(c0, c1) // factor
        sub = np.asarray(v[rows.min():rows.max() + 1, cols.min():cols.max() + 1]).astype(np.float64)
        sub = sub[np.ix_(rows - rows.min(), cols - cols.min())]
        if fill is not None:
            sub[sub == fill] = np.nan
        sec = sub + off
        if 'day_shift' not in scene.extra:
            good = np.asarray(v[::50, ::50]).astype(np.float64)
            good = good[good != fill] + off if fill is not None else good.ravel() + off
            mid = (scene.time + (scene.extra['t1'] - scene.time) / 2 - datetime.datetime(1970, 1, 1)).total_seconds()
            scene.extra['day_shift'] = round((np.median(good) - mid) / 86400.) if good.size else 0
            if scene.extra['day_shift']:
                print(tr('  [MVIRI-FCDR] 時間變數與檔名相差 %d 天，以檔名為準',
                         '  [MVIRI-FCDR] the time variable differs from the file name by %d day(s); using the file name')
                      % scene.extra['day_shift'])
        sec = sec - scene.extra['day_shift'] * 86400.
        fallback = (scene.time - datetime.datetime(1970, 1, 1)).total_seconds()
        sec = np.where(np.isfinite(sec), sec, fallback)
        return (np.datetime64('1970-01-01T00:00:00') + (sec * 1000).astype('timedelta64[ms]'))

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        n = 5000 if channel == 'VIS' else 2500
        a0, a1, b0, b1 = max(r0, 0), min(r1, n), max(c0, 0), min(c1, n)
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        if a1 <= a0 or b1 <= b0:
            return out
        ds = self._open(scene.files[channel])
        if channel in ('IR', 'WV'):
            k = channel.lower()
            cnt = np.asarray(ds.variables['count_' + k][a0:a1, b0:b1]).astype(np.float64)
            fill = int(ds.variables['count_' + k].getncattr('_FillValue'))
            rad = self._scalar(ds, 'a_' + k) + self._scalar(ds, 'b_' + k) * cnt
            with np.errstate(divide='ignore', invalid='ignore'):
                t = self._scalar(ds, 'bt_b_' + k) / (np.log(rad) - self._scalar(ds, 'bt_a_' + k)) - 273.15
            vals = np.where((cnt != fill) & (rad > 0), t, np.nan)
        else:
            if 'toa_bidirectional_reflectance_vis' in ds.variables:
                v = ds.variables['toa_bidirectional_reflectance_vis']
                raw = np.asarray(v[a0:a1, b0:b1]).astype(np.float64)
                brf = np.where(raw != int(v.getncattr('_FillValue')),
                               raw * float(v.getncattr('scale_factor')) + float(v.getncattr('add_offset')), np.nan)
            else:                                                     # FULL 版：由計數值計算
                cnt = np.asarray(ds.variables['count_vis'][a0:a1, b0:b1]).astype(np.float64)
                y = self._scalar(ds, 'years_since_launch')
                gain = self._scalar(ds, 'a0_vis') + self._scalar(ds, 'a1_vis') * y + self._scalar(ds, 'a2_vis') * y * y
                rad = (cnt - self._scalar(ds, 'mean_count_space_vis')) * gain
                refl = np.pi * self._scalar(ds, 'distance_sun_earth') ** 2 * rad / self._scalar(ds, 'solar_irradiance_vis')
                out[a0 - r0:a1 - r0, b0 - c0:b1 - c0] = np.where(rad > 0, refl, np.nan)
                return out
            from satima.geo.truecolor.astronomy import cos_sun_zenith
            R, C = np.meshgrid(np.arange(a0, a1), np.arange(b0, b1), indexing='ij')
            lon, lat = self.nav(scene, 'VIS').to_lonlat(R, C)
            mu = cos_sun_zenith(self._pixel_time(scene, ds, a0, a1, b0, b1, 2), lon, lat)
            vals = np.where(np.isfinite(mu), brf * np.clip(mu, 0., 1.), np.nan)
        out[a0 - r0:a1 - r0, b0 - c0:b1 - c0] = vals
        return out


READER = MviriFcdrReader
