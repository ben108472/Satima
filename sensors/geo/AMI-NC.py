"""GK-2A AMI — Level 1B (NetCDF)

檔名：gk2a_ami_le1b_vi006_fd005ge_202611080400.nc
  波段 (vi004 ... ir133)、觀測範圍 (fd 全圓盤、la 區域、ko 韓國域...)、解析度 (005 = 0.5 km)、
  投影 (ge = 靜止衛星投影、lc = 蘭伯特正形圓錐投影)
  全圓盤以外的範圍，地圖範圍自動帶入。
導航：ge 為 CGMS 標準靜止衛星投影 (行列號由 1 起算；檔案屬性 image_upperleft_x/y 與此一致)；
      KMA 檔案的 lfac 為負值，換成 JMA 慣例 (正值) 使用，同 satpy。
校正：計數值取有效位元後查校正表 (conversion_table/GK-2A)。可見光的表與檔案屬性的係數完全相同；
      紅外的表使用精確的有效波數，檔案只有標稱中心波長 (誤差約 1.6 K)，所以沿用校正表。
"""
import datetime
import os
import re

import numpy as np
from netCDF4 import Dataset

from satima.geo.nav import GeosNav, ProjGridNav
from satima.geo.reader import Channel, FileCache, GeoReader, Scene

_NAME = re.compile(r'gk2a_ami_le1b_([a-z]{2}\d{3})_([a-z]+?)(\d{3})(ge|lc)_(\d{12})\.nc$')
# (代號, 中心波長 µm, 解析度 km, 種類, 功能代稱)；波長取 KMA GK2A Fact Sheet 表 3-1 (取到小數 2 位)，
# 順序即 KMA 的頻道編號 1~16 (顯示為 BAND-01~BAND-16)
_BANDS = (('VI004', 0.47, 1.0, 'vis', 'Blue'), ('VI005', 0.51, 1.0, 'vis', 'Green'),
          ('VI006', 0.64, 0.5, 'vis', 'Red'), ('VI008', 0.86, 1.0, 'nir', 'Veggie'),
          ('NR013', 1.38, 2.0, 'nir', 'Cirrus'), ('NR016', 1.61, 2.0, 'nir', 'Snow/Ice'),
          ('SW038', 3.83, 2.0, 'swir', 'Shortwave Window'), ('WV063', 6.24, 2.0, 'wv', 'Upper-level WV'),
          ('WV069', 6.95, 2.0, 'wv', 'Mid-level WV'), ('WV073', 7.34, 2.0, 'wv', 'Lower-level WV'),
          ('IR087', 8.59, 2.0, 'ir', 'Cloud-top Phase'), ('IR096', 9.63, 2.0, 'ir', 'Ozone'),
          ('IR105', 10.4, 2.0, 'ir', 'Clean LW Window'), ('IR112', 11.21, 2.0, 'ir', 'LW Window'),
          ('IR123', 12.36, 2.0, 'ir', 'Dirty LW Window'), ('IR133', 13.31, 2.0, 'ir', 'CO2 LW'))
# 計數值的有效位元數
_BITS = {'VI004': 11, 'VI005': 11, 'NR016': 11, 'VI006': 12, 'NR013': 12, 'WV063': 12, 'SW038': 14}


def _open(path):
    ds = Dataset(path)
    ds.variables['image_pixel_values'].set_auto_maskandscale(False)     # 計數值自行取有效位元
    return ds


class AmiNcReader(GeoReader):
    NAME = 'AMI-NC'
    SENSOR = 'AMI'
    DESCRIPTION = 'GK-2A AMI (Level 1B NetCDF)'
    PATTERNS = ('gk2a_ami_le1b_*_*_*.nc',)
    CHANNELS = {b[0]: Channel(b[0], b[1], b[2], b[3], 'BAND-%02d' % i, b[4]) for i, b in enumerate(_BANDS, 1)}
    TRUECOLOR = {'TrueColor': dict(red='VI006', green='VI005', blue='VI004', nir='VI008', ir='IR105', swir='SW038',
                                   lo_blur=(0.6, 0.4))}   # AMI 1 km 波段南北向較模糊

    def __init__(self, cfg):
        super().__init__(cfg)
        self._ds = FileCache(_open)
        self._tables = {}

    def scenes(self, folder='.'):
        groups = {}
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            band, area, _, proj, stamp = m.groups()
            key = (area, proj, stamp)
            if key not in groups:
                t = datetime.datetime.strptime(stamp, '%Y%m%d%H%M')
                groups[key] = Scene('GK-2A', t, area.upper(), auto_range=(area != 'fd'), extra={'proj': proj})
            groups[key].files[band.upper()] = os.path.join(folder, name)
        return sorted(groups.values(), key=lambda s: (s.time, s.area))

    def platform(self, scene):
        return 'GEO-KOMPSAT-2A'

    def _dataset(self, path):
        return self._ds.get(path)

    def nav(self, scene, channel):
        ds = self._dataset(scene.files[channel])
        a = {k: ds.getncattr(k) for k in ds.ncattrs()}
        shape = (int(a['number_of_lines']), int(a['number_of_columns']))
        if scene.extra['proj'] == 'ge':
            req, rpol = float(a['earth_equatorial_radius']), float(a['earth_polar_radius'])
            return GeosNav.from_cgms(shape, a['cfac'], -float(a['lfac']), a['coff'], a['loff'],
                                     np.rad2deg(float(a['sub_longitude'])),
                                     float(a['nominal_satellite_height']) - req, req, rpol)
        # lc：蘭伯特正形圓錐投影的規則格點 (第 0 列為北邊)；與舊版 GK2A_NC 相同的定義 (尚無樣本資料驗證)
        params = dict(proj='lcc', lon_0=float(a['central_meridian']), lat_0=float(a['origin_latitude']),
                      lat_1=float(a['standard_parallel1']), lat_2=float(a['standard_parallel2']))
        x0, x1 = float(a['lower_left_easting']), float(a['lower_right_easting'])
        y_top, y_bottom = float(a['upper_left_northing']), float(a['lower_left_northing'])
        return ProjGridNav(shape, params, x0, (x1 - x0) / (shape[1] - 1), y_top, (y_bottom - y_top) / (shape[0] - 1))

    def _table(self, channel):
        if channel not in self._tables:
            kind = 'alb' if self.CHANNELS[channel].reflective else 'bt'
            path = os.path.join(self.cfg.general.conversion_table_dir, 'GK-2A', '%s_con_%s.txt' % (channel, kind))
            table = np.loadtxt(path, dtype=np.float64)
            self._tables[channel] = (table if kind == 'alb' else table - 273.15).astype(np.float32)
        return self._tables[channel]

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        counts = self._dataset(scene.files[channel]).variables['image_pixel_values'][r0:r1, c0:c1]
        counts = np.bitwise_and(counts.astype(np.int64), (1 << _BITS.get(channel, 13)) - 1)
        table = self._table(channel)
        out = np.full(counts.shape, np.nan, dtype=np.float32)
        ok = counts < table.size
        out[ok] = table[counts[ok]]
        return out

    def satellite_position(self, scene):
        a = self._dataset(next(iter(scene.files.values())))
        try:
            return (float(np.rad2deg(a.getncattr('sub_longitude'))),
                    (float(a.getncattr('nominal_satellite_height')) - float(a.getncattr('earth_equatorial_radius'))) / 1000.)
        except AttributeError:
            return 128.2, 35786.0


READER = AmiNcReader
