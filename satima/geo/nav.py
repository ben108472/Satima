"""導航 (經緯度 <-> 影像行列號)

核心 (裁切、取樣、真實色對齊) 只透過 Nav 介面使用導航，各衛星依自己的資料選一種實作：
  GeosNav      CGMS 標準靜止衛星投影 (AHI、AMI、ABI、FCI、SEVIRI、AGRI...)，只需填參數
  ProjGridNav  規則投影格點 (例如 GK2A 的 LCC 區域產品)
  LonLatNav    只有每個像素的經緯度陣列時使用 (萬用備援，例如以軌道/姿態計算經緯度的老衛星)
若衛星有特殊之處 (例如 SEVIRI HRV 的上下視窗)，在該衛星的檔案中繼承 GeosNav 覆寫即可。

行列號一律為「陣列索引」：由 0 起算、像素中心位於整數 (像素 i 涵蓋 i-0.5 ~ i+0.5)。
看不到的位置 (太空、地球背面) 回傳 NaN。
"""
import numpy as np


class Nav:
    shape = (0, 0)          # 影像大小 (列, 行)

    def to_rowcol(self, lon, lat):
        """經緯度 (度) -> (列, 行) 浮點陣列。"""
        raise NotImplementedError

    def to_lonlat(self, row, col):
        """(列, 行) -> 經緯度 (度)。"""
        raise NotImplementedError

    def key(self):
        """相同 key 的導航可共用取樣結果 (例如同解析度的波段)。"""
        return id(self)


def _finite(*arrays):
    out = []
    for a in arrays:
        a = np.asarray(a, dtype=np.float64)
        out.append(np.where(np.isfinite(a) & (np.abs(a) < 1e20), a, np.nan))
    return out


class GeosNav(Nav):
    """CGMS 標準靜止衛星投影 (normalized geostationary projection)。

    影像行列號與掃描角 (弧度，採 proj 慣例：x 向東為正、y 向北為正) 為線性關係：
        x = x0 + col * dx        y = y0 + row * dy
    sweep：掃描軸，'y' = Meteosat / Himawari / GK2A / MTG 慣例，'x' = GOES-R 慣例。
    height：衛星離地 (橢球面) 高度 (m)；a、b：地球長、短半軸 (m)。
    """

    def __init__(self, shape, x0, dx, y0, dy, sub_lon, height, a=6378137.0, b=6356752.31414, sweep='y'):
        from pyproj import Proj
        self.shape = (int(shape[0]), int(shape[1]))
        self.x0, self.dx, self.y0, self.dy = float(x0), float(dx), float(y0), float(dy)
        self.sub_lon, self.height = float(sub_lon), float(height)
        self.a, self.b, self.sweep = float(a), float(b), str(sweep)
        self._proj = Proj(proj='geos', h=self.height, lon_0=self.sub_lon, a=self.a, b=self.b, sweep=self.sweep)

    @classmethod
    def from_cgms(cls, shape, cfac, lfac, coff, loff, sub_lon, height, a, b, first_line=1, first_col=1, sweep='y'):
        """CGMS LRIT/HRIT Global Specification 4.4.4：
            c = COFF + nint(x · 2^-16 · CFAC)     l = LOFF + nint(y · 2^-16 · LFAC)
        c、l 由 1 起算 (4.2.3.1：左上角像素為 (1,1))，x、y 為度，y 向南為正 (LFAC > 0)。
        陣列第 0 列、第 0 行分別對應行號 first_line、first_col。"""
        kx = np.deg2rad(2.0 ** 16) / float(cfac)
        ky = np.deg2rad(2.0 ** 16) / float(lfac)
        x0 = (first_col - float(coff)) * kx
        y0 = -(first_line - float(loff)) * ky       # y 向南為正 -> proj 慣例取負號
        return cls(shape, x0, kx, y0, -ky, sub_lon, height, a, b, sweep)

    def key(self):
        return ('geos', self.shape, self.x0, self.dx, self.y0, self.dy, self.sub_lon, self.height, self.sweep)

    def scan_angles(self, row, col):
        return self.x0 + np.asarray(col, dtype=np.float64) * self.dx, self.y0 + np.asarray(row, dtype=np.float64) * self.dy

    def to_lonlat(self, row, col):
        x, y = self.scan_angles(row, col)
        lon, lat = self._proj(x * self.height, y * self.height, inverse=True)
        return _finite(lon, lat)

    def to_rowcol(self, lon, lat):
        with np.errstate(invalid='ignore'):
            x, y = self._proj(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        x, y = _finite(x, y)
        return (y / self.height - self.y0) / self.dy, (x / self.height - self.x0) / self.dx


class ProjGridNav(Nav):
    """規則投影格點：proj 平面座標與行列號為線性關係 x = x0 + col*dx、y = y0 + row*dy (公尺)。"""

    def __init__(self, shape, proj_params, x0, dx, y0, dy):
        from pyproj import Proj
        self.shape = (int(shape[0]), int(shape[1]))
        self.params = dict(proj_params)
        self.x0, self.dx, self.y0, self.dy = float(x0), float(dx), float(y0), float(dy)
        self._proj = Proj(**self.params)

    def key(self):
        return ('proj', self.shape, tuple(sorted(self.params.items())), self.x0, self.dx, self.y0, self.dy)

    def to_lonlat(self, row, col):
        x = self.x0 + np.asarray(col, dtype=np.float64) * self.dx
        y = self.y0 + np.asarray(row, dtype=np.float64) * self.dy
        return _finite(*self._proj(x, y, inverse=True))

    def to_rowcol(self, lon, lat):
        x, y = _finite(*self._proj(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)))
        return (y - self.y0) / self.dy, (x - self.x0) / self.dx


class LonLatNav(Nav):
    """只有經緯度陣列時的導航：以 KD-tree 找最近的像素 (回傳整數行列號)。
    max_km：最近像素超過此距離即視為看不到 (預設為像素間距的 2 倍)。"""

    def __init__(self, lon, lat, max_km=None):
        from scipy.spatial import cKDTree
        lon, lat = _finite(lon, lat)
        self.shape = lon.shape
        self._lon, self._lat = lon, lat
        ok = np.isfinite(lon) & np.isfinite(lat)
        self._index = np.flatnonzero(ok)
        self._tree = cKDTree(self._xyz(lon.ravel()[self._index], lat.ravel()[self._index]))
        if max_km is None:
            step = np.nanmedian(np.hypot(np.diff(lat, axis=1), np.diff(lon, axis=1) * np.cos(np.deg2rad(lat[:, 1:]))))
            max_km = 2.0 * 111.2 * float(step)
        self.max_km = max_km

    @staticmethod
    def _xyz(lon, lat):
        lo, la = np.deg2rad(lon), np.deg2rad(lat)
        return np.column_stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])

    def to_rowcol(self, lon, lat):
        lon, lat = _finite(lon, lat)
        shape = lon.shape
        lon, lat = lon.ravel(), lat.ravel()
        rows = np.full(lon.size, np.nan)
        cols = np.full(lon.size, np.nan)
        ok = np.isfinite(lon) & np.isfinite(lat)
        dist, k = self._tree.query(self._xyz(lon[ok], lat[ok]), distance_upper_bound=self.max_km / 6371.0)
        hit = np.isfinite(dist)
        r, c = np.unravel_index(self._index[k[hit]], self.shape)
        sel = np.flatnonzero(ok)[hit]
        rows[sel], cols[sel] = r, c
        return rows.reshape(shape), cols.reshape(shape)

    def to_lonlat(self, row, col):
        r = np.clip(np.rint(row).astype(int), 0, self.shape[0] - 1)
        c = np.clip(np.rint(col).astype(int), 0, self.shape[1] - 1)
        return self._lon[r, c], self._lat[r, c]
