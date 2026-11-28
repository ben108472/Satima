"""微光通道 (VIIRS DNB、FY-3E MERSI-LL) 共用的縮放

  erf_dnb            ERF Dynamic DNB (CIRA/SSEC，移植自 satpy ERFDNB / Polar2Grid dynamic_dnb)：
                     依太陽天頂角、月亮天頂角與月相動態決定縮放範圍；輻亮度單位為 W cm-2 sr-1 (VIIRS SDR 原生單位)
  lunar_phase_factor 月相 -> ERF 的月光項 (Allen 月亮星等公式)
  moon_phase_angle   月亮相位角 (度，0 = 滿月)：Meeus《Astronomical Algorithms》第 48 章的低精度公式 (誤差約 0.1 度)，
                     不需要 ephem；給檔案裡沒有月相的格式 (MERSI-LL) 用
VIIRS-HDF5 的 erf_dnb 與這裡的演算法相同 (VIIRS 用 GDNBO 檔裡的月相)。
"""
import datetime

import numpy as np

DNB_SAT_PCT = 0.005    # 飽和修正：允許超過上限的像元比例
DNB_SAT_STEP = 1.1     # 飽和修正：上限每次放大倍率


def lunar_phase_factor(phase):
    """相位角 phase (度) 時月光比滿月暗 dm 星等，換算成 ERF 月光項；上限 0.7 為 ERF 原本的新月 (氣輝) 設定。"""
    p = abs(phase)
    dm = 0.026 * p + 4e-9 * p ** 4
    return min(0.2 * dm, 0.7)


def erf_dnb(rad, sza, lza, moon_frac, saturation=False, moon_phase=None):
    """rad：W cm-2 sr-1；sza / lza：太陽 / 月亮天頂角 (度)；moon_frac：月亮照明比例 (0~1)；
    moon_phase：相位角 (度)，有的話月光項依相位角計算。回傳 0~1。"""
    from scipy.special import erf
    rad = np.where(rad < -900., np.nan, rad)
    if moon_phase is not None:
        moon_factor1 = lunar_phase_factor(moon_phase)
    else:
        moon_factor1 = 0.7 * (1.0 - moon_frac)
    moon_factor2 = 0.0022 * lza
    erf_portion = 1 + erf((sza - 95.0) / (5.0 * np.sqrt(2.0)))
    max_val = np.power(10, -1.7 - (2.65 + moon_factor1 + moon_factor2) * erf_portion)
    min_val = np.power(10, -4.0 - (2.95 + moon_factor2) * erf_portion)
    if saturation:
        valid = np.isfinite(rad) & np.isfinite(max_val)
        nvalid = max(valid.sum(), 1)
        while np.sum(rad[valid] > max_val[valid]) / nvalid > DNB_SAT_PCT:
            max_val = max_val * DNB_SAT_STEP
    with np.errstate(invalid='ignore'):
        inner = (rad - min_val) / (max_val - min_val)
        inner = np.where(inner < 0, 0., inner)
    return np.clip(np.sqrt(inner), 0., 1.)


def moon_phase_angle(dt):
    """月亮相位角 (度，0 = 滿月、180 = 新月) 與照明比例 (0~1)。dt：UTC datetime。"""
    jd = (dt - datetime.datetime(2000, 1, 1, 12)).total_seconds() / 86400. + 2451545.0
    t = (jd - 2451545.0) / 36525.
    d = np.deg2rad(297.8501921 + 445267.1114034 * t)       # 月亮平距角
    m = np.deg2rad(357.5291092 + 35999.0502909 * t)        # 太陽平近點角
    mp = np.deg2rad(134.9633964 + 477198.8675055 * t)      # 月亮平近點角
    i = (180. - np.rad2deg(d) - 6.289 * np.sin(mp) + 2.100 * np.sin(m) - 1.274 * np.sin(2 * d - mp)
         - 0.658 * np.sin(2 * d) - 0.214 * np.sin(2 * mp) - 0.110 * np.sin(d))
    i = abs((i + 180.) % 360. - 180.)
    return float(i), float((1. + np.cos(np.deg2rad(i))) / 2.)
