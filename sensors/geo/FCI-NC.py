"""MTG-I FCI — Level 1c (NetCDF chunks in zip)

檔名：W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI1+FCI-1C-RRAD-FDHSI-FD--x-x---x_C_EUMT_..._IDPFI_OPE_<開始>_<結束>_N__O_0081_0000.zip
  FDHSI：16 波段 (可見光/近紅外 1 km、紅外 2 km)
  HRFI ：vis_06_hr、nir_22_hr (0.5 km)、ir_38_hr、ir_105_hr (1 km)
  FDHSI 與 HRFI 分成兩個場景，各有一個 TrueColor：
    FDHSI 的 TrueColor 為 1 km；HRFI 的 TrueColor 以 HRFI 0.5 km 紅光銳化 (藍、綠、近紅外取同時次的 FDHSI)。
資料直接從 zip 讀進記憶體 (不解壓縮)；每個 zip 有 40 個 BODY chunk，由南到北排列。
顯示名稱：頻道依波長順序編為 BAND-01~BAND-16 (EUMETSAT 頻道表順序；注意 BAND-13 為 9.7 µm 臭氧、
  BAND-14 才是 10.5 µm 窗區)，HRFI 波段沿用同波段的編號；標註與檔名再接上資料集標記，
  例如 ir_105 -> BAND-14 FDHSI、ir_105_hr -> BAND-14 HRFI (檔名 BAND-14_FDHSI、BAND-14_HRFI)。
導航：x / y 的封裝值為由 1 起算的行列號，以 scale_factor / add_offset 換算成掃描角 (弧度)；
      x 向西為正 (與 proj 慣例相反)，y 向北為正。
"""
import datetime
import os
import re
import zipfile

import numpy as np
from netCDF4 import Dataset

from satima.geo.nav import GeosNav
from satima.geo.reader import Channel, GeoReader, Scene, to_float

_NAME = re.compile(r'W_XX-EUMETSAT-Darmstadt,IMG\+SAT,MTI(\d)\+FCI-1C-RRAD-(FDHSI|HRFI)-(\w+?)--.*?'
                   r'_IDPFI_[A-Z]+_(\d{14})_(\d{14})_.*\.zip$')
# (代號, 中心波長 µm, 解析度 km, 種類, 功能代稱)；波長取 EUMETSAT FCI 頻道表
_FDHSI = (('vis_04', 0.444, 1.0, 'vis', 'Blue'), ('vis_05', 0.51, 1.0, 'vis', 'Green'),
          ('vis_06', 0.64, 1.0, 'vis', 'Red'), ('vis_08', 0.865, 1.0, 'nir', 'Veggie'),
          ('vis_09', 0.914, 1.0, 'nir', 'Near-IR WV'), ('nir_13', 1.38, 1.0, 'nir', 'Cirrus'),
          ('nir_16', 1.61, 1.0, 'nir', 'Snow/Ice'), ('nir_22', 2.25, 1.0, 'nir', 'Cloud Particle Size'),
          ('ir_38', 3.8, 2.0, 'swir', 'Shortwave Window'), ('wv_63', 6.3, 2.0, 'wv', 'Upper-level WV'),
          ('wv_73', 7.35, 2.0, 'wv', 'Lower-level WV'), ('ir_87', 8.7, 2.0, 'ir', 'Cloud-top Phase'),
          ('ir_97', 9.66, 2.0, 'ir', 'Ozone'), ('ir_105', 10.5, 2.0, 'ir', 'Clean LW Window'),
          ('ir_123', 12.3, 2.0, 'ir', 'Dirty LW Window'), ('ir_133', 13.3, 2.0, 'ir', 'CO2 LW'))
_HRFI = (('vis_06_hr', 0.64, 0.5, 'vis', 'Red'), ('nir_22_hr', 2.25, 0.5, 'nir', 'Cloud Particle Size'),
         ('ir_38_hr', 3.8, 1.0, 'swir', 'Shortwave Window'), ('ir_105_hr', 10.5, 1.0, 'ir', 'Clean LW Window'))
# 頻道編號：FDHSI 依波長順序 1~16，HRFI 沿用同波段的編號
_NUMBER = {b[0]: i for i, b in enumerate(_FDHSI, 1)}
_NUMBER.update({h[0]: _NUMBER[h[0][:-3]] for h in _HRFI})
_AU_KM = 149597870.7
_DS_CACHE = 12          # 同時保留在記憶體中的 chunk 數


class FciNcReader(GeoReader):
    NAME = 'FCI-NC'
    SENSOR = 'FCI'
    DESCRIPTION = 'MTG-I FCI (Level 1c FDHSI / HRFI, zip)'
    PATTERNS = ('W_XX-EUMETSAT-Darmstadt,IMG+SAT,MTI*+FCI-1C-RRAD-*.zip',)
    CHANNELS = {b[0]: Channel(b[0], b[1], b[2], b[3], 'BAND-%02d' % _NUMBER[b[0]], b[4]) for b in _FDHSI + _HRFI}
    _TC = dict(green='vis_05', blue='vis_04', nir='vis_08', ir='ir_105', swir='ir_38')
    # 兩個配置在選單上都叫 TrueColor (FDHSI 場景用前者、HRFI 場景用後者，標註/檔名靠資料集標記區分)
    TRUECOLOR = {'TrueColor': dict(_TC, red='vis_06'),
                 # HRFI 的 0.5 km 紅光銳化；FDHSI 1 km 波段比 HRFI 紅光 2x2 平均模糊 (同波段 vis_06 擬合)
                 'TrueColorHR': dict(_TC, red='vis_06_hr', lo_blur=(0.6, 0.6), rsr={'vis_06_hr': 'vis_06'})}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._zips = {}
        self._ds = {}
        self._rows = {}

    def scenes(self, folder='.'):
        """FDHSI 與 HRFI 各成一個場景；HRFI 場景另外帶有同時次 FDHSI 的檔案 (真實色用，選單不列出)。"""
        groups = {}
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            sat, product, area, start, _ = m.groups()
            key = (sat, area, start, product)
            if key not in groups:
                t = datetime.datetime.strptime(start, '%Y%m%d%H%M%S')
                groups[key] = Scene('MTG-I%s' % sat, t, '%s %s' % (area, product), tag=product)
            for b in (_FDHSI if product == 'FDHSI' else _HRFI):
                groups[key].files[b[0]] = os.path.join(folder, name)
        for (sat, area, start, product), scene in groups.items():
            fdhsi = groups.get((sat, area, start, 'FDHSI'))
            if product == 'HRFI' and fdhsi:
                for band, path in fdhsi.files.items():
                    scene.files.setdefault(band, path)
        return sorted(groups.values(), key=lambda s: (s.time, s.area))

    def channels(self, scene):
        """選單只列出該資料集自己的波段 (HRFI 場景附帶的 FDHSI 波段不列出)。"""
        hr = scene.tag == 'HRFI'
        return [ch for name, ch in self.CHANNELS.items() if name in scene.files and name.endswith('_hr') == hr]

    def truecolor_products(self, scene):
        """FDHSI：1 km TrueColor；HRFI：HRFI 紅光銳化的 TrueColor (需要同時次的 FDHSI)。"""
        recipe = self.TRUECOLOR['TrueColorHR' if scene.tag == 'HRFI' else 'TrueColor']
        need = [recipe[k] for k in ('red', 'green', 'blue', 'nir')]
        return self.with_geocolor({'TrueColor': recipe} if all(b in scene.files for b in need) else {}, scene)

    # ---- zip / chunk -------------------------------------------------------
    def _chunks(self, path):
        if path not in self._zips:
            z = zipfile.ZipFile(path)
            names = sorted((n for n in z.namelist() if n.endswith('.nc') and '-CHK-BODY-' in n),
                           key=lambda n: n.split('_')[-1])
            self._zips[path] = (z, names)
        return self._zips[path]

    def _open(self, path, name):
        key = (path, name)
        if key not in self._ds:
            z, _ = self._chunks(path)
            self._ds[key] = Dataset(name, memory=z.read(name))
            while len(self._ds) > _DS_CACHE:
                self._ds.pop(next(iter(self._ds))).close()
        else:
            self._ds[key] = self._ds.pop(key)
        return self._ds[key]

    def _row_range(self, path, name, channel):
        """chunk 在該波段格點中的列範圍 [start, end) (陣列索引)。"""
        key = (path, name, channel)
        if key not in self._rows:
            m = self._open(path, name).groups['data'][channel]['measured']
            self._rows[key] = (int(m.variables['start_position_row'][:]) - 1, int(m.variables['end_position_row'][:]))
        return self._rows[key]

    def _chunks_for(self, path, channel, r0, r1):
        """與列範圍 [r0, r1) 重疊的 chunk；chunk 由南到北排列，以二分搜尋只需開啟少數 chunk。"""
        _, names = self._chunks(path)
        lo, hi = 0, len(names)
        while lo < hi:
            mid = (lo + hi) // 2
            if self._row_range(path, names[mid], channel)[1] > r0:
                hi = mid
            else:
                lo = mid + 1
        out = []
        for name in names[lo:]:
            start, end = self._row_range(path, name, channel)
            if start >= r1:
                break
            out.append((name, start, end))
        return out

    # ---- 讀取 --------------------------------------------------------------
    def nav(self, scene, channel):
        path = scene.files[channel]
        _, names = self._chunks(path)
        ds = self._open(path, names[len(names) // 2])
        m = ds.groups['data'][channel]['measured']
        x, y = m.variables['x'], m.variables['y']
        sx, ox, sy, oy = float(x.scale_factor), float(x.add_offset), float(y.scale_factor), float(y.add_offset)
        shape = (int(y.valid_range[1]), int(x.valid_range[1]))
        p = ds.groups['data'].variables['mtg_geos_projection']
        # 陣列索引 i 對應封裝值 i+1；x 向西為正 -> proj 慣例取負號
        return GeosNav(shape, -(sx + ox), -sx, sy + oy, sy, p.longitude_of_projection_origin,
                       p.perspective_point_height, p.semi_major_axis, p.semi_minor_axis, p.sweep_angle_axis)

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        path = scene.files[channel]
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        for name, start, end in self._chunks_for(path, channel, r0, r1):
            a, b = max(r0, start), min(r1, end)
            ds = self._open(path, name)
            m = ds.groups['data'][channel]['measured']
            rad = to_float(m.variables['effective_radiance'][a - start:b - start, c0:c1], np.float64)
            out[a - r0:b - r0] = self._calibrate(ds, m, channel, rad)
        return out

    def _calibrate(self, ds, m, channel, rad):
        v = m.variables
        if self.CHANNELS[channel].reflective:
            d = float(np.mean(ds.groups['state'].groups['celestial'].variables['earth_sun_distance'][:])) / _AU_KM
            return (np.pi * rad * d ** 2 / float(v['channel_effective_solar_irradiance'][:])).astype(np.float32)
        a, b, wnc, c1, c2 = (float(v[k][:]) for k in (
            'radiance_to_bt_conversion_coefficient_a', 'radiance_to_bt_conversion_coefficient_b',
            'radiance_to_bt_conversion_coefficient_wavenumber', 'radiance_to_bt_conversion_constant_c1',
            'radiance_to_bt_conversion_constant_c2'))
        with np.errstate(invalid='ignore', divide='ignore'):
            return ((c2 * wnc) / (a * np.log(1 + c1 * wnc ** 3 / rad)) - b / a - 273.15).astype(np.float32)

    def satellite_position(self, scene):
        path = next(iter(scene.files.values()))
        _, names = self._chunks(path)
        p = self._open(path, names[0]).groups['data'].variables['mtg_geos_projection']
        return float(p.longitude_of_projection_origin), float(p.perspective_point_height) / 1000.


READER = FciNcReader
