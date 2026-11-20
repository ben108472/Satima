"""COMS-1 MI (Meteorological Imager) — Level 1B HDF-EOS5 (.he5)

檔名：coms_mi_le1b_zzz_fd000ge_201309290145.he5 (一個檔案含全部 5 個頻道)
  觀測範圍 (檔名第 18~19 字)：fd 全圓盤 (FD)、cn (ENH)、en (EN)，同舊版 COMS_HDF.py；全圓盤以外地圖範圍自動帶入
  頻道：VIS 0.675 µm (1 km)、SWIR 3.75、WV 6.75、IR1 10.8、IR2 12.0 µm (4 km)，名稱沿用 KMA
導航：標準靜止衛星投影 (CGMS)。COFF / CFAC / LOFF / LFAC 為 HDFEOS/POINTS/Navigation Result/Data/NAVI MI
      的屬性 (Visible、IR 各一組)，LFAC 為負值 (KMA 慣例，同 GK-2A)，取絕對值使用。
校正：conversion_table/COMS/COMS-1.txt (計數值 -> IR 亮溫 K、VIS 反照率 %)，同舊版；
      表中為 0 的計數值 (IR 約 950 以後、VIS 901 以後) 為無效值。

注意：本讀取器依舊版 COMS_HDF.py 的資料結構改寫，撰寫時沒有樣本資料，以下各點待有實際資料時確認：
  * 行列號起算：依 CGMS 規範由 1 起算 (同其他衛星)。舊版以 0 起算的角點格網繪圖，與本讀取器差半個像素
    (IR 約 2 km、VIS 約 0.5 km)。
  * 衛星經度屬性 (Satellite Definition / Longitude) 的單位：舊版遇到小於 120 的值一律改用 128.2，
    推測為弧度；本讀取器把 2π 以內的值當弧度換算，換算後不在 100~160 度時用 128.2。
  * 'Altitude at Scene Center Time' 的意義與單位：舊版當作地心距離 (km)；本讀取器依數值大小判斷
    (公尺 / 離地高度 / 地心距離)，不合理時用 42164 km。
"""
import datetime
import math
import os
import re

import h5py
import numpy as np

from satima.geo.nav import GeosNav
from satima.geo.reader import Channel, FileCache, GeoReader, Scene

_NAME = re.compile(r'coms_mi_le1b_([a-z0-9]+)_([a-z]{2})[a-z0-9]*_(\d{12})\.he5$')
_AREA = {'fd': 'FD', 'cn': 'ENH', 'en': 'EN'}                 # 同舊版
# (名稱, 中心波長 µm, 解析度 km, 種類, 功能代稱)
_BANDS = (('VIS', 0.675, 1.0, 'vis', 'Visible'), ('IR1', 10.8, 4.0, 'ir', 'LW Window'),
          ('IR2', 12.0, 4.0, 'ir', 'Dirty LW Window'), ('WV', 6.75, 4.0, 'wv', 'Water Vapor'),
          ('SWIR', 3.75, 4.0, 'swir', 'Shortwave Window'))
_DATASET = {'VIS': 'HDFEOS/GRIDS/Visible Image Pixel Value/Data Fields/Visible Image Pixel Values',
            'IR1': 'HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/IR1 band Image Pixel Values',
            'IR2': 'HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/IR2 band Image Pixel Values',
            'WV': 'HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/WV band Image Pixel Values',
            'SWIR': 'HDFEOS/GRIDS/IR Image Pixel Value/Data Fields/SWIR band Image Pixel Values'}
_NAV = 'HDFEOS/POINTS/Navigation Result/Data/NAVI MI'
_SAT_DEF = 'HDFEOS/POINTS/Satellite Status/Data/Satellite Definition'
_ORBIT = 'HDFEOS/POINTS/Satellite Status/Data/Satellite_Orbit'
_EARTH = 'HDFEOS/POINTS/Geometric Processing/Data/Earth Model'
_TABLE_COLUMN = {'IR1': 1, 'IR2': 2, 'WV': 3, 'SWIR': 4, 'VIS': 5}   # COMS-1.txt 的欄位
# CGMS 標準值 (檔案中沒有或數值不合理時使用)
_REQ, _RPOL, _DISTANCE, _SUB_LON = 6378.169, 6356.5838, 42164.0, 128.2


def _attr(f, path, name):
    """HDF5 屬性的第一個數值；沒有時回傳 None。"""
    try:
        return float(np.ravel(f[path].attrs[name])[0])
    except (KeyError, IndexError, TypeError, ValueError):
        return None


class MiHdfReader(GeoReader):
    NAME = 'MI-HDF'
    SENSOR = 'MI'
    DESCRIPTION = 'COMS-1 MI (Level 1B HDF-EOS5 .he5)'
    PATTERNS = ('coms_mi_le1b_*_*.he5',)
    CHANNELS = {b[0]: Channel(*b[:4], alias=b[4]) for b in _BANDS}
    TRUECOLOR = {}                  # 只有一個寬頻可見光，無法合成真實色

    def __init__(self, cfg):
        super().__init__(cfg)
        self._files = FileCache(lambda path: h5py.File(path, 'r'), maxsize=16)
        self._tables = {}

    def scenes(self, folder='.'):
        scenes = []
        for name in self.find(folder):
            m = _NAME.search(name)
            if not m:
                continue
            _, area, stamp = m.groups()
            t = datetime.datetime.strptime(stamp, '%Y%m%d%H%M')
            path = os.path.join(folder, name)
            scenes.append(Scene('COMS-1', t, area.upper(), files={ch: path for ch in _DATASET},
                                auto_range=(area != 'fd'), region=_AREA.get(area, '')))
        return sorted(scenes, key=lambda s: (s.time, s.area))

    def channels(self, scene):
        """檔案中實際有的頻道 (一個檔案含全部頻道，個別頻道可能缺)。"""
        try:
            f = self._files.get(next(iter(scene.files.values())))
        except OSError:
            return []
        return [ch for name, ch in self.CHANNELS.items() if _DATASET[name] in f]

    # ---- 導航 --------------------------------------------------------------
    def _geometry(self, f):
        """(衛星經度 度, 赤道半徑 km, 極半徑 km, 地心距離 km)。"""
        req = _attr(f, _EARTH, 'Equitorial Radius of Earth Ellipsoid')      # 屬性名稱原文拼字 (Equitorial)
        rpol = _attr(f, _EARTH, 'Polar Radius of Earth Ellipsoid')
        req = req / 1000. if req and req > 1e5 else (req or _REQ)
        rpol = rpol / 1000. if rpol and rpol > 1e5 else (rpol or _RPOL)
        dist = _attr(f, _ORBIT, 'Altitude at Scene Center Time')
        if dist and dist > 1e6:
            dist /= 1000.                       # 公尺 -> km
        if dist and 30000. < dist < 40000.:
            dist += req                         # 離地高度 -> 地心距離
        if not dist or not 41000. < dist < 43500.:
            dist = _DISTANCE
        lon = _attr(f, _SAT_DEF, 'Longitude')
        if lon is not None and abs(lon) <= 2 * math.pi + 0.01:
            lon = math.degrees(lon)             # 弧度 -> 度
        if lon is None or not 100. <= lon <= 160.:
            lon = _SUB_LON                      # 同舊版：不合理時用 128.2
        return lon, req, rpol, dist

    def nav(self, scene, channel):
        f = self._files.get(scene.files[channel])
        coff, cfac, loff, lfac = (float(v) for v in np.ravel(f[_NAV].attrs['Visible' if channel == 'VIS' else 'IR'])[:4])
        lon, req, rpol, dist = self._geometry(f)
        # KMA 的 LFAC 為負值 (舊版乘 -1)，取絕對值即 JMA / CGMS 慣例
        return GeosNav.from_cgms(f[_DATASET[channel]].shape, abs(cfac), abs(lfac), coff, loff, lon,
                                 (dist - req) * 1000., req * 1000., rpol * 1000.)

    # ---- 讀取 --------------------------------------------------------------
    def _table(self, channel):
        if channel not in self._tables:
            path = os.path.join(self.cfg.general.conversion_table_dir, 'COMS', 'COMS-1.txt')
            col = np.loadtxt(path)[:, _TABLE_COLUMN[channel]]
            lut = np.full(65536, np.nan, dtype=np.float32)
            valid = np.flatnonzero(col)
            last = int(valid[-1])
            if channel == 'VIS':
                lut[:last + 1] = col[:last + 1] / 100.          # 反照率 % -> 0~1 (同舊版 /100)
            else:
                lut[valid] = col[valid] - 273.15                # 亮溫 K -> °C；表中為 0 的計數值為無效
            self._tables[channel] = lut
        return self._tables[channel]

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        counts = np.asarray(self._files.get(scene.files[channel])[_DATASET[channel]][r0:r1, c0:c1]).astype(np.int64)
        counts[(counts < 0) | (counts > 65535)] = 65535                    # 範圍外 -> 無效值
        return self._table(channel)[counts]

    def satellite_position(self, scene):
        lon, req, _, dist = self._geometry(self._files.get(next(iter(scene.files.values()))))
        return lon, dist - req


READER = MiHdfReader
