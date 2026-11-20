"""規則 tie-point 網格內插 (共用)

適用：tie point 第 k 點位於第 k×interval 像元的規則網格 (例如 SGLI 每 10 像元一點)。
  interp_geo   經緯度：轉成地心直角座標 (x, y, z) 各自做樣條內插再轉回，換日線與極區不會出錯
               (同 satpy sgli_l1b / python-geotiepoints 的 GeoSplineInterpolator，kx = ky = 2)
  interp_field 純量場 (例如太陽天頂角)

SGLI 實測 (tie point 間距加倍成 20 像元，內插回原 tie point 與真值比較)：
  griddata 線性 (原作法)      平均 4.2 m、99% 21.5 m；全解析度 7820×5000 約 34 秒
  直角座標二次樣條 (此模組)   平均 0.3 m、99%  0.8 m；全解析度約 3 秒
"""
import numpy as np
from scipy.interpolate import RectBivariateSpline

from . import swath


def _tie_axes(tie_shape, interval):
    ri, ci = (interval, interval) if np.isscalar(interval) else interval
    return np.arange(tie_shape[0]) * float(ri), np.arange(tie_shape[1]) * float(ci)


def _lonlat_to_xyz(lons, lats):
    lat_r = np.deg2rad(lats)
    lon_r = np.deg2rad(lons)
    return np.cos(lat_r) * np.cos(lon_r), np.cos(lat_r) * np.sin(lon_r), np.sin(lat_r)


def pixel_coords(n_out, n_ref):
    """輸出網格 (n_out 個像元) 的像元中心在參考網格 (n_ref 個像元，tie point 以此為單位) 上的座標。
    同解析度為 0, 1, 2, ...；例如 1 km 對 250 m (n_ref = 4 n_out) 為 1.5, 5.5, 9.5, ... (4 個 250 m 像元的中心)。"""
    s = float(n_ref) / n_out
    return np.arange(n_out) * s + (s - 1.) / 2.


def _eval_axes(shape, rows, cols):
    rows = np.arange(shape[0], dtype=np.float64) if rows is None else np.asarray(rows, dtype=np.float64)
    cols = np.arange(shape[1], dtype=np.float64) if cols is None else np.asarray(cols, dtype=np.float64)
    return rows, cols


def interp_geo(lats, lons, interval, shape, rows=None, cols=None, k=2, rows_per_chunk=1024, dtype=np.float32):
    """tie-point 經緯度 -> 全解析度 (lats, lons)，shape = (列, 行)；經度範圍 -180~180。
    rows / cols：各輸出列/行在 tie 參考網格上的座標 (預設 0, 1, 2, ...，見 pixel_coords)。"""
    lats = swath.fill_invalid(np.asarray(lats, dtype=np.float64))
    lons = swath.fill_invalid(np.asarray(lons, dtype=np.float64))
    tie_r, tie_c = _tie_axes(lats.shape, interval)
    rows, cols = _eval_axes(shape, rows, cols)
    splines = [RectBivariateSpline(tie_r, tie_c, v, kx=k, ky=k) for v in _lonlat_to_xyz(lons, lats)]
    out_lat = np.empty(shape, dtype=dtype)
    out_lon = np.empty(shape, dtype=dtype)
    for r0 in range(0, shape[0], rows_per_chunk):
        r = rows[r0:r0 + rows_per_chunk]
        x, y, z = (s(r, cols) for s in splines)
        out_lat[r0:r0 + r.size] = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
        out_lon[r0:r0 + r.size] = np.rad2deg(np.arctan2(y, x))
    return out_lat, out_lon


def interp_field(field, interval, shape, rows=None, cols=None, k=1, rows_per_chunk=1024, dtype=np.float32):
    """tie-point 純量場 -> 全解析度 (預設線性)；無效值先以鄰近值補上。rows / cols 同 interp_geo。"""
    field = swath.fill_invalid(np.asarray(field, dtype=np.float64))
    tie_r, tie_c = _tie_axes(field.shape, interval)
    rows, cols = _eval_axes(shape, rows, cols)
    spline = RectBivariateSpline(tie_r, tie_c, field, kx=k, ky=k)
    out = np.empty(shape, dtype=dtype)
    for r0 in range(0, shape[0], rows_per_chunk):
        r = rows[r0:r0 + rows_per_chunk]
        out[r0:r0 + r.size] = spline(r, cols)
    return out
