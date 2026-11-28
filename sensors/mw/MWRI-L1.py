"""風雲三號 MWRI Level 1 亮溫 (NSMC，HDF5)：FY-3A~D MWRI、FY-3F/H MWRI-II、FY-3G MWRI-RM

檔名  FY3A~D：FY3X_MWRI{A|D}_GBAL_L1_<YYYYMMDD>_<hhmm>_010KM_MS.HDF       (MWRIA = 升軌、MWRID = 降軌)
      FY3F~H：FY3X_MWRI-_ORB{A|D}_L1_<YYYYMMDD>_<hhmm>_{010KM|7000M}_Vn.HDF
亮溫  int16：亮溫 K = 值 x Slope (0.01) + Intercept (327.68)
  MWRI     EARTH_OBSERVE_BT_10_to_89GHz (10, 掃描, 像素)：10.65、18.7、23.8、36.5、89 GHz 的 V、H (V 在前)
           FY-3D 在 Calibration/、FY-3C 在 Data/、FY-3A/B 在最上層；經緯度在 Geolocation/ 或最上層
  MWRI-II  Window Channel/Calibration/EARTH_OBSERVE_BT (掃描, 520, 10)：同上 10 個窗區頻道
           Sounding Channel/Calibration/EARTH_OBSERVE_BT (掃描, 520, 12)：50.3、52.61、53.24、53.75 GHz 的 V、H，
           118.7503±3.2、±2.1、±1.4、±1.2 GHz (V)；兩組各有自己的 Geolocation/ (窗區以 89V、探測以 50.3V 定位)
  MWRI-RM  S1/Data/EARTH_OBSERVE_BT_10_to_89GHz (掃描, 492, 10)：同上 10 個窗區頻道
           S2/Data/EARTH_OBSERVE_BT_50_to_183GHz (掃描, 492, 16)：MWRI-II 的 12 個探測頻道，加 165.5±0.75、
           183.31±2、±3.4、±7 GHz (V)；S1、S2 各有 Geolocation/
  (V / H 的順序 2026/11 以各衛星樣本的海面中位數確認：每一對都是 V 較高)
時間  2000-01-01 12:00 UTC 起的日數 + 當日 12:00 起的毫秒 (x Slope；MWRI-II / RM 為 0.1)：
      FY-3C/D、MWRI-II、MWRI-RM：Scan_Daycnt + Scan_Mscnt (FY-3C 為小寫 Scan_daycnt / Scan_mscnt)
      FY-3A/B：Scan_Time_and_Period 的第 1、2 欄 (日數、毫秒)
格式說明：NSMC 各衛星的 MWRI「一級產品使用說明」、FY-3A/B L1 Product Format (satellite.nsmc.org.cn 文件頁)。
FY-3A 沒有樣本，照格式說明寫 (每條掃描 240 點)，未驗證。
"""
import datetime
import os
import re

import h5py
import numpy as np

from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, to_seconds

_NAME = re.compile(r'^(FY3[A-Z])_MWRI([AD\-X]?)_(GBAL|ORB[ADT]?)_L1_(\d{8})_(\d{4})_', re.I)
_EPOCH = to_seconds(datetime.datetime(2000, 1, 1, 12))
_SENSOR = {'F': 'MWRI-II', 'G': 'MWRI-RM', 'H': 'MWRI-II'}

# (頻率 GHz, 偏移, 極化)，依資料集內的順序
_WINDOW = tuple((f, None, p) for f in (10.65, 18.7, 23.8, 36.5, 89.0) for p in 'VH')
_SOUNDING = (tuple((f, None, p) for f in (50.3, 52.61, 53.24, 53.75) for p in 'VH')
             + tuple((118.7503, d, 'V') for d in (3.2, 2.1, 1.4, 1.2))
             + ((165.5, 0.75, 'V'),) + tuple((183.31, d, 'V') for d in (2.0, 3.4, 7.0)))

# 各種檔案的分組：(群組名, 亮溫資料集, 經緯度所在群組, 時間所在群組, 頻道在最後一維, 頻道表)
_LAYOUTS = {
    'MWRI': (('all', ('Calibration/EARTH_OBSERVE_BT_10_to_89GHz', 'Data/EARTH_OBSERVE_BT_10_to_89GHz',
                      'EARTH_OBSERVE_BT_10_to_89GHz'), ('Geolocation', ''), ('Calibration', 'Data', ''), False,
              _WINDOW),),
    'MWRI-II': (('win', ('Window Channel/Calibration/EARTH_OBSERVE_BT',), ('Window Channel/Geolocation',),
                 ('Window Channel/Calibration',), True, _WINDOW),
                ('snd', ('Sounding Channel/Calibration/EARTH_OBSERVE_BT',), ('Sounding Channel/Geolocation',),
                 ('Sounding Channel/Calibration',), True, _SOUNDING)),
    'MWRI-RM': (('S1', ('S1/Data/EARTH_OBSERVE_BT_10_to_89GHz',), ('S1/Geolocation',), ('S1/Data',), True, _WINDOW),
                ('S2', ('S2/Data/EARTH_OBSERVE_BT_50_to_183GHz',), ('S2/Geolocation',), ('S2/Data',), True,
                 _SOUNDING)),
}


def _find(f, name, groups):
    """在各群組找資料集 (名稱不分大小寫：FY-3C 為 Scan_daycnt)。找不到回傳 None。"""
    for g in groups:
        if g and g not in f:
            continue
        grp = f[g] if g else f
        if not isinstance(grp, h5py.Group):
            continue
        if name in grp:
            return grp[name]
        low = name.lower()
        for k in grp.keys():
            if k.lower() == low and isinstance(grp[k], h5py.Dataset):
                return grp[k]
    return None


def _first(f, paths):
    for p in paths:
        if p in f:
            return f[p]
    return None


def _attr(f, key, default=''):
    v = f.attrs.get(key, default)
    v = v[0] if isinstance(v, np.ndarray) and v.size else v
    return v.decode('latin1').strip() if isinstance(v, bytes) else str(v).strip()


def _slope(ds, default=1.):
    v = ds.attrs.get('Slope')
    return float(np.ravel(v)[0]) if v is not None and np.size(v) else default


def _layout(f):
    if 'Window Channel' in f:
        return 'MWRI-II'
    if 'S1' in f:
        return 'MWRI-RM'
    return 'MWRI'


class MWRIReader(SwathReader):
    NAME = 'MWRI-L1'
    SENSOR = 'MWRI'
    SOURCE = 'NSMC'
    DESCRIPTION = '風雲三號 (FY-3A~D MWRI、FY-3F/H MWRI-II、FY-3G MWRI-RM) Level 1 亮溫 (NSMC)'
    DESCRIPTION_EN = 'Fengyun-3 (FY-3A~D MWRI, FY-3F/H MWRI-II, FY-3G MWRI-RM) Level 1 brightness temperature (NSMC)'
    PATTERNS = ('FY3?_MWRI?_*_L1_*.HDF',)

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(os.path.basename(name))
            if not m:
                continue
            path = os.path.join(folder, name)
            start = datetime.datetime.strptime(m.group(4) + m.group(5), '%Y%m%d%H%M')
            end = start + datetime.timedelta(minutes=55)
            layout = None
            try:
                with h5py.File(path, 'r') as f:
                    layout = _layout(f)
                    a = _attr(f, 'Observing Beginning Date') + ' ' + _attr(f, 'Observing Beginning Time')[:8]
                    b = _attr(f, 'Observing Ending Date') + ' ' + _attr(f, 'Observing Ending Time')[:8]
                    start = datetime.datetime.strptime(a, '%Y-%m-%d %H:%M:%S')
                    end = datetime.datetime.strptime(b, '%Y-%m-%d %H:%M:%S')
            except (KeyError, ValueError):
                pass                                # 檔頭沒有或格式不同：用檔名的時間
            except OSError:
                continue                            # 不是 HDF5 (下載不完整等)
            sat = 'FY-3' + m.group(1)[-1].upper()
            orbit = m.group(2).upper() if m.group(2).upper() in ('A', 'D') else m.group(3)[3:].upper()
            sensor = _SENSOR.get(sat[-1], 'MWRI') if layout is None else layout
            out.append(Scene(sat, sensor, start, end, file=path, orbit=orbit, format='L1', extra={'layout': layout}))
        return out

    def channels(self, scene):
        if 'channels' in scene.extra:
            return scene.extra['channels']
        out = []
        with h5py.File(scene.file, 'r') as f:
            layout = scene.extra.get('layout') or _layout(f)
            for group, bt_paths, _, _, last, table in _LAYOUTS[layout]:
                ds = _first(f, bt_paths)
                if ds is None:
                    continue
                n = ds.shape[-1] if last else ds.shape[0]
                for k, (freq, offset, pol) in enumerate(table[:n]):
                    out.append(Channel(band_name(freq, pol, offset), freq, pol, group, key=k))
        scene.extra['layout'] = layout
        scene.extra['channels'] = out
        return out

    def _spec(self, scene, group):
        for spec in _LAYOUTS[scene.extra.get('layout') or 'MWRI']:
            if spec[0] == group:
                return spec
        raise KeyError(group)

    def geo(self, scene, group):
        _, _, geo_groups, time_groups, _, _ = self._spec(scene, group)
        with h5py.File(scene.file, 'r') as f:
            lat = _find(f, 'Latitude', geo_groups)[:]
            lon = _find(f, 'Longitude', geo_groups)[:]
            dc = _find(f, 'Scan_Daycnt', time_groups)
            ms = _find(f, 'Scan_Mscnt', time_groups)
            if dc is not None and ms is not None:
                day = dc[:].astype(np.float64).reshape(-1)
                msec = ms[:].astype(np.float64)
                msec = msec[:, 0] if msec.ndim == 2 else msec
                if scene.extra.get('layout') != 'MWRI':     # MWRI-II / RM 的毫秒以 0.1 ms 存 (Slope)；FY-3C/D 照舊不乘
                    msec = msec * _slope(ms)
            else:                                   # FY-3A/B：Scan_Time_and_Period (日數, 毫秒, …)
                tp = _find(f, 'Scan_Time_and_Period', time_groups)
                if tp is None:
                    raise KeyError('Scan_Daycnt / Scan_Time_and_Period')
                tp = tp[:].astype(np.float64)
                day, msec = tp[:, 0], tp[:, 1]
        time = np.where((day > 0) & (msec >= 0), _EPOCH + day * 86400. + msec / 1000., np.nan)
        lat, lon = clean_latlon(lat, lon)
        return [Geo(lat, lon, time, (scene.file, group))]

    def tb(self, scene, channel):
        _, bt_paths, _, _, last, _ = self._spec(scene, channel.group)
        with h5py.File(scene.file, 'r') as f:
            d = _first(f, bt_paths)
            raw = d[:, :, channel.key] if last else d[channel.key]
            slope = _slope(d, 0.01)
            inter = float(np.ravel(d.attrs.get('Intercept', 327.68))[0])
        data = raw.astype(np.float32) * np.float32(slope) + np.float32(inter)
        data[(raw <= -32766) | ~(data > 0.) | (data > 400.)] = np.nan
        return [data]


READER = MWRIReader
