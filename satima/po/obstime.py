"""觀測時間：各格式的時間欄位換成 UTC (numpy datetime64[us]，缺值為 NaT)

各衛星檔案使用的時間尺度不同 (以實際資料對照檔名、全域屬性與產品內的太陽天頂角確認)：
  UTC：AVHRR (EPS GRH、NetCDF record_start_time)、METimage scan_start_time_utc
  TAI：VIIRS IET (1958-01-01 起的微秒)、MODIS / SGLI TAI93 (1993-01-01 起的秒，含閏秒)
  GPS：Sentinel-3 OLCI / SLSTR time_stamp (2000-01-01 起的微秒；GPS = UTC + 18 s，自 2017 年起)
"""
import datetime

import numpy as np

from ..core.i18n import tr

US ='datetime64[us]'
NAT = np.datetime64('NaT', 'us')

# TAI − UTC (秒) 與生效日 (UTC)；之後若再有閏秒，於此表最後加一列
TAI_UTC = [
    ('1972-01-01', 10), ('1972-07-01', 11), ('1973-01-01', 12), ('1974-01-01', 13), ('1975-01-01', 14),
    ('1976-01-01', 15), ('1977-01-01', 16), ('1978-01-01', 17), ('1979-01-01', 18), ('1980-01-01', 19),
    ('1981-07-01', 20), ('1982-07-01', 21), ('1983-07-01', 22), ('1985-07-01', 23), ('1988-01-01', 24),
    ('1990-01-01', 25), ('1991-01-01', 26), ('1992-07-01', 27), ('1993-07-01', 28), ('1994-07-01', 29),
    ('1996-01-01', 30), ('1997-07-01', 31), ('1999-01-01', 32), ('2006-01-01', 33), ('2009-01-01', 34),
    ('2012-07-01', 35), ('2015-07-01', 36), ('2017-01-01', 37),
]
_LEAP_DATES = np.array([d for d, _ in TAI_UTC], dtype=US)
_LEAP_VALUES = np.array([v for _, v in TAI_UTC], dtype=np.float64)
TAI_GPS = 19.   # TAI − GPS (秒，固定)


def tai_minus_utc(t):
    """t (datetime64) 當時的 TAI − UTC 秒數。"""
    idx = np.searchsorted(_LEAP_DATES, np.asarray(t, dtype=US), side='right') - 1
    return np.where(idx >= 0, _LEAP_VALUES[np.clip(idx, 0, None)], 0.)


def _delta(values, unit):
    """數值 (整數 / 浮點，可為 masked、NaN) -> timedelta64[us]；無效值為 NaT。"""
    a = np.ma.asarray(values)
    scale = {'s': 1_000_000, 'ms': 1_000, 'us': 1}[unit]
    if a.dtype.kind in 'iu':
        out = (a.filled(0).astype(np.int64) * scale).astype('timedelta64[us]')
        bad = np.ma.getmaskarray(a)
    else:
        v = a.astype(np.float64).filled(np.nan)
        bad = ~np.isfinite(v)
        out = np.round(np.where(bad, 0., v) * scale).astype(np.int64).astype('timedelta64[us]')
    out = np.asarray(out)
    out[bad] = np.timedelta64('NaT', 'us')
    return out


def from_utc(epoch, values, unit='s'):
    """以 UTC 計數的時間 (epoch 起的 values 個單位，不含閏秒)。"""
    return np.datetime64(epoch, 'us') + _delta(values, unit)


def from_tai(epoch, values, unit='s', epoch_tai_utc=0.):
    """以 TAI 計數的時間 (含閏秒)；epoch_tai_utc 為 epoch 當時的 TAI − UTC。"""
    naive = np.datetime64(epoch, 'us') + _delta(values, unit)
    return naive - _delta(tai_minus_utc(naive) - epoch_tai_utc, 's')


def tai93(seconds):
    """MODIS / SGLI：1993-01-01 00:00:00 UTC 起的 TAI 秒數 (當時 TAI − UTC = 27)。"""
    return from_tai('1993-01-01', seconds, 's', epoch_tai_utc=27.)


def iet(microseconds):
    """VIIRS IET：1958-01-01 起的 TAI 微秒數。"""
    return from_tai('1958-01-01', microseconds, 'us', epoch_tai_utc=0.)


def gps2000(microseconds):
    """Sentinel-3 time_stamp：2000-01-01 起的 GPS 時間微秒數 (UTC = GPS − (TAI−UTC − 19 s))。"""
    naive = np.datetime64('2000-01-01', 'us') + _delta(microseconds, 'us')
    return naive - _delta(tai_minus_utc(naive) - TAI_GPS, 's')


def per_scan(times, nrows):
    """掃描式感測器：每個掃描帶一個時間 -> 每一列 (同一掃描帶的列同一時間；缺值保留 NaT)。"""
    times = np.asarray(times, dtype=US)
    if times.size == 0 or nrows % times.size:
        raise ValueError(tr('列數 %d 不是掃描數 %d 的整數倍', '%d rows is not a multiple of %d scans')
                         % (nrows, times.size))
    return np.repeat(times, nrows // times.size)


def resample_rows(times, nrows):
    """逐列時間換到不同解析度的 nrows 列 (例如 SGLI 250 m 的列時間 -> 1 km 頻道)，依列中心線性內插。"""
    times = np.asarray(times, dtype=US)
    if nrows == times.size:
        return times
    ok = ~np.isnat(times)
    if not ok.any():
        return np.full(nrows, NAT)
    base = times[ok][0]
    us = (times[ok] - base) / np.timedelta64(1, 'us')
    s = times.size / nrows
    pos = np.arange(nrows) * s + (s - 1.) / 2.        # 輸出列中心在原列號上的座標
    out = np.interp(pos, np.flatnonzero(ok), us)
    return base + np.round(out).astype(np.int64).astype('timedelta64[us]')


def linear(start, end, nrows):
    """沒有逐列時間時的備援：起訖時間之間平均分配 (原腳本的作法)。"""
    start, end = np.datetime64(start, 'us'), np.datetime64(end, 'us')
    frac = (np.arange(nrows) + 0.5) / nrows
    span = (end - start) / np.timedelta64(1, 'us')
    return start + np.round(frac * span).astype(np.int64).astype('timedelta64[us]')


def to_datetime(t):
    """datetime64 -> datetime.datetime (UTC，naive)。"""
    return np.datetime64(t, 'us').astype(datetime.datetime)
