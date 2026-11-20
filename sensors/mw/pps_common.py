"""NASA PPS (GPM) 的 HDF5 / netCDF-4 swath 檔共用讀取 (PPS-L1C.py、GMI-L1B.py 使用)

結構：S1、S2… 每個 swath 一組：Latitude / Longitude (掃描, 像素)、亮溫 (掃描, 像素, 頻道)、ScanTime/*。
L1C 的亮溫為 Tc，各 swath 的頻道寫在 Tc 的 long_name，例如
  "Intercalibrated Tb for channels 1) 19.35 GHz V-Pol 2) 19.35 GHz H-Pol and 3) 22.235 GHz V-Pol"
頻道表由這段文字解析；1B (GMI) 的亮溫為 Tb，沒有這段文字，用固定的頻道表 (TABLE)。
AMSR2 / AMSR3 的 89 GHz 有 A-Scan、B-Scan 兩個 swath：[swath] merge_ab = true 時合併成一個頻道
(內插時取兩組掃描線的聯集，見 satima/swath.py)，false 時分列為 89V-A、89V-B。
"""
import datetime
import os
import re

import h5py
import numpy as np

from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon

SATELLITES = {'GCOMW1': 'GCOM-W1', 'GOSATGW': 'GOSAT-GW', 'METOPA': 'MetOp-A', 'METOPB': 'MetOp-B',
              'METOPC': 'MetOp-C', 'NPP': 'S-NPP', 'WSFM1': 'WSF-M1', 'MT1': 'Megha-Tropiques', 'AQUA': 'Aqua'}
SENSORS = {'SSMI': 'SSM/I', 'AMSRE': 'AMSR-E', 'AMSUB': 'AMSU-B'}       # 檔名中的儀器代號 -> 顯示名稱
# 頻率後面的 +/- 只在吸收線附近的頻道才是頻道的一部分 (183.31±7、118.75±…、50~60 GHz)；
# AMSU-B 的「89.0 +/- 0.9」「150.0 +/- 0.9」是窗區頻道的雙邊帶，名稱只取中心頻率
_LINES = ((45., 70.), (113., 124.), (175., 192.))
_ITEM = re.compile(r'(?:\band\s+)?\d+\)\s*')
_POL = re.compile(r'\b(Q?[VH])-?Pol', re.I)
_SCAN = re.compile(r'\b([AB])-Scan', re.I)
_NUM = r'(\d+(?:\.\d+)?)'
_FREQ = re.compile(_NUM + r'\s*(?:GHz)?\s*(?:\+\s*/?\s*-|±)\s*' + _NUM + r'|' + _NUM)


def text(v):
    return v.decode('latin1') if isinstance(v, bytes) else str(v)


def satellite_name(code):
    code = code.upper()
    if code in SATELLITES:
        return SATELLITES[code]
    m = re.fullmatch(r'NOAA(\d+)', code)
    return 'NOAA-%s' % m.group(1) if m else code


def parse_channels(long_name):
    """Tc 的 long_name -> [(頻率, 偏移, 極化, A/B 掃描)]；解析不出來時回傳 None。"""
    s = re.sub(r'\s+', ' ', long_name)
    k = s.find('1)')
    out = []
    for item in _ITEM.split(s[k:] if k >= 0 else s):
        item = item.strip()
        if not item:
            continue
        f = _FREQ.search(item)
        if not f:
            return None
        freq = float(f.group(1) or f.group(3))
        offset = float(f.group(2)) if f.group(2) else None
        if offset and not any(lo <= freq <= hi for lo, hi in _LINES):
            offset = None
        pol = _POL.search(item)
        scan = _SCAN.search(item)
        out.append((freq, offset, pol.group(1).upper() if pol else '', scan.group(1).upper() if scan else ''))
    return out


def granule_times(date, start, end):
    """檔名的 YYYYMMDD、S / E 的 HHMMSS -> (起, 訖)；跨日時結束時間加一天。"""
    a = datetime.datetime.strptime(date + start, '%Y%m%d%H%M%S')
    b = datetime.datetime.strptime(date + end, '%Y%m%d%H%M%S')
    if b < a:
        b += datetime.timedelta(days=1)
    return a, b


class PPSReader(SwathReader):
    VAR = 'Tc'          # 亮溫的資料集名稱
    TABLE = {}          # 沒有 long_name 時的頻道表 {swath: [(頻率, 偏移, 極化, 掃描), ...]}
    FILENAME = None     # 檔名的正規表示式 (具名群組 sat、sensor、date、s、e、orbit；可有 variant：接在格式名稱後面)
    FORMAT = ''
    SOURCE = 'PPS'

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = self.FILENAME.match(os.path.basename(name))
            if not m:
                continue
            start, end = granule_times(m['date'], m['s'], m['e'])
            sensor = m['sensor'].upper()
            out.append(Scene(satellite_name(m['sat']), SENSORS.get(sensor, sensor), start, end,
                             file=os.path.join(folder, name), orbit=m['orbit'] or '',
                             format=self.FORMAT + (m.groupdict().get('variant') or '').upper()))
        return out

    # ---- 頻道 --------------------------------------------------------------
    def channels(self, scene):
        key = 'channels_%s' % self.cfg.swath.merge_ab
        if key not in scene.extra:
            scene.extra[key] = self._channels(scene)
        return scene.extra[key]

    def _channels(self, scene):
        raw = []                                    # (swath, 索引, 頻率, 偏移, 極化, 掃描)
        with h5py.File(scene.file, 'r') as f:
            groups = sorted((g for g in f if re.fullmatch(r'S\d+', g) and self.VAR in f[g]), key=lambda g: int(g[1:]))
            for g in groups:
                tc = f[g][self.VAR]
                n = tc.shape[2] if tc.ndim == 3 else 1
                items = self.TABLE.get(g) or parse_channels(text(tc.attrs.get('long_name', '')))
                if not items or len(items) != n:
                    items = [(0., None, '', '')] * n          # 解析不出來：以 swath 與序號命名
                for k, (freq, offset, pol, scan) in enumerate(items):
                    raw.append((g, k, freq, offset, pol, scan))
        merge = bool(self.cfg.swath.merge_ab)
        out, done = [], set()
        for g, k, freq, offset, pol, scan in raw:
            if (g, k) in done:
                continue
            if not freq:
                out.append(Channel('%s-%d' % (g, k + 1), 0., pol, g, key=((g, k),)))
                continue
            name = band_name(freq, pol, offset)
            if scan and merge:
                # 同頻率、同極化的另一組掃描 (A <-> B)
                mate = [r for r in raw if r[2:5] == (freq, offset, pol) and r[5] and r[5] != scan]
                if mate:
                    a, b = ((g, k), mate[0][:2]) if scan == 'A' else (mate[0][:2], (g, k))
                    done.update((a, b))
                    out.append(Channel(name, freq, pol, a[0] + '+' + b[0], key=(a, b)))
                    continue
            out.append(Channel(name + ('-' + scan if scan else ''), freq, pol, g, key=((g, k),)))
        return out

    # ---- 資料 --------------------------------------------------------------
    def geo(self, scene, group):
        out = []
        with h5py.File(scene.file, 'r') as f:
            for g in group.split('+'):
                lat, lon = clean_latlon(f[g]['Latitude'][:], f[g]['Longitude'][:])
                out.append(Geo(lat, lon, scan_time(f[g]['ScanTime']), (scene.file, g)))
        return out

    def tb(self, scene, channel):
        out = []
        mask = bool(self.cfg.swath.mask_bad_quality)
        with h5py.File(scene.file, 'r') as f:
            for g, k in channel.key:
                tc = f[g][self.VAR]
                data = (tc[:, :, k] if tc.ndim == 3 else tc[:]).astype(np.float32)
                bad = ~(data > 0.)
                if mask and 'Quality' in f[g]:
                    bad |= f[g]['Quality'][:] < 0
                data[bad] = np.nan
                out.append(data)
        return out


def scan_time(st):
    """ScanTime 群組 (Year、Month、DayOfMonth、SecondOfDay…) -> 各掃描 1970 起算的秒數 (無效為 NaN)。"""
    y, mo, d = (st[k][:].astype(np.int64) for k in ('Year', 'Month', 'DayOfMonth'))
    sec = st['SecondOfDay'][:].astype(np.float64) if 'SecondOfDay' in st else (
        st['Hour'][:] * 3600. + st['Minute'][:] * 60. + st['Second'][:] + st['MilliSecond'][:] / 1000.)
    ok = (y > 1900) & (mo >= 1) & (mo <= 12) & (d >= 1) & (d <= 31) & (sec >= 0)
    y, mo, d = np.where(ok, y, 1970), np.where(ok, mo, 1), np.where(ok, d, 1)
    days = ((y - 1970).astype('datetime64[Y]').astype('datetime64[M]') + (mo - 1).astype('timedelta64[M]')
            ).astype('datetime64[D]') + (d - 1).astype('timedelta64[D]')
    return np.where(ok, days.astype(np.int64) * 86400. + sec, np.nan)
