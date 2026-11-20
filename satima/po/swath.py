"""軌道 (swath) 資料的共用處理：經度展開、範圍裁切、座標補值、颱風眼溫度

經度約定：所有經度先依中心經度 clon 展開到 [clon-180, clon+180)，
因此換日線兩側 (例如中心 180°、範圍 170°~190°) 不需要另外處理 0~360 / -180~180 的轉換，
繪圖時地圖的中央經線也設在 clon，投影後座標是連續的。
"""
import numpy as np
import pandas as pd

from ..core.i18n import tr


def as_float(a):
    """轉成浮點 ndarray，masked 值改為 NaN。"""
    if np.ma.isMaskedArray(a):
        if a.dtype.kind != 'f':
            a = a.astype(np.float64)
        return a.filled(np.nan)
    a = np.asarray(a)
    if a.dtype.kind != 'f':
        a = a.astype(np.float64)
    return a


def wrap_lon(lons, clon):
    """經度展開到 [clon-180, clon+180)。"""
    lons = as_float(lons)
    return (lons - clon + 180.) % 360. - 180. + clon


def swath_box(georange, lats, lons):
    """包住 georange 的最小矩形 (列、行皆含端點)；lons 需已用 wrap_lon 展開。
    回傳 (yi, yj, xi, xj)，切片用 [yi:yj, xi:xj]；範圍內沒有資料時回傳 None。"""
    latmin, latmax, lonmin, lonmax = georange
    lats = as_float(lats)
    with np.errstate(invalid='ignore'):
        inside = (lats > latmin) & (lats < latmax) & (lons > lonmin) & (lons < lonmax)
    rows = np.flatnonzero(inside.any(axis=1))
    if rows.size == 0:
        return None
    cols = np.flatnonzero(inside.any(axis=0))
    return int(rows[0]), int(rows[-1]) + 1, int(cols[0]), int(cols[-1]) + 1


def crop(box, arr):
    yi, yj, xi, xj = box
    return arr[yi:yj, xi:xj]


def fill_invalid(a):
    """NaN 以沿軌 (axis 0) 前後值補上，再以沿掃描 (axis 1) 補上。
    matplotlib 3.5 的 pcolormesh 不接受含 NaN 的座標，繪圖前必須先補 (VIIRS fill_gaps 的通用版)。"""
    a = as_float(a)
    if np.isfinite(a).all():
        return a
    df = pd.DataFrame(np.where(np.isfinite(a), a, np.nan)).ffill().bfill()
    if df.isna().values.any():
        df = df.ffill(axis=1).bfill(axis=1)
    return df.values


def angular_distance(lats, lons, clat, clon):
    """各點到 (clat, clon) 的大圓角距離 (度，haversine)；經緯度無效處為 masked。"""
    lats = np.ma.masked_invalid(np.ma.asarray(lats))
    lons = np.ma.masked_invalid(np.ma.asarray(lons))
    dlat = np.deg2rad(lats - clat)
    dlon = np.deg2rad((lons - clon + 180.0) % 360.0 - 180.0)
    hav = (np.sin(dlat / 2.0) ** 2
           + np.cos(np.deg2rad(clat)) * np.cos(np.deg2rad(lats))
           * np.sin(dlon / 2.0) ** 2)
    return np.rad2deg(2.0 * np.arcsin(np.sqrt(np.clip(hav, 0.0, 1.0))))


def eyetemp(data, lats, lons, clat, clon, radius_deg=0.5):
    """中心點半徑 radius_deg (大圓角距離) 內的最高亮溫 (MODIS / METimage 版本)。"""
    data = np.ma.masked_invalid(np.ma.asarray(data))
    distance = angular_distance(lats, lons, clat, clon)
    inside = np.ma.filled(distance <= radius_deg, False)
    eye_data = np.ma.masked_where(~inside, np.ma.masked_outside(data, -150., 100.))

    print('====EYETEMP====')
    if eye_data.count():
        print(tr('最高亮溫 (中心 %.2f 度內)：%.5f C', 'MAX (within %.2f deg): %.5f C') % (radius_deg, eye_data.max()))
    else:
        print(tr('中心 %.2f 度內沒有有效的資料。', 'No valid data within %.2f deg of center.') % radius_deg)
    print('====END====')


def eyebase(data, t):
    return data[data > t].shape[0]


def eyetemp_counts(data):
    """舊版 eyetemp：印出暖於 9~24°C 的像元數 (扣除 >100°C 的無效值)。"""
    data = np.ma.filled(np.ma.masked_invalid(np.ma.asarray(data, dtype=float)), np.nan)
    with np.errstate(invalid='ignore'):
        print('====EYETEMP (counts)====')
        minus = eyebase(data, 100)
        for i in [9, 12, 15, 18, 20, 21, 22, 23, 24]:
            print(i, eyebase(data, i) - minus)
    print('MAX', np.nanmax(data))
    print('MIN', np.nanmin(data))
    print('====END====')
