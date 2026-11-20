"""JAXA AMSR 系列 Level 1B 亮溫的共用讀取 (AMSR2-L1B.py、AMSRE-L1B.py 使用)

GCOM-W1 AMSR2 的 L1B (HDF5) 與 Aqua AMSR-E 重新處理成「AMSR2 格式」的 L1B 結構相同：
亮溫  'Brightness Temperature (<頻率>GHz,<極化>)'，uint16 x 0.01 K
      低頻每條掃描 243 點；89.0 GHz 的 A、B 兩組掃描各 486 點
經緯度  檔案裡只有 89A / 89B 的經緯度。低頻各頻段的經緯度由 89A 推算 (JAXA 的 co-registration)：
        取 89A 第 2i、2i+1 點 (0 起算) 的地心單位向量 p1、p2，θ 為兩者夾角，
          ex = p1，ez = (p1 x p2) / |p1 x p2|，ey = ez x ex
          低頻第 i 點 = cos(A2 θ) [cos(A1 θ) ex + sin(A1 θ) ey] + sin(A2 θ) ez
        A1、A2 為檔頭 CoRegistrationParameterA1 / A2 中該頻段的值。
時間  'Scan Time' 為 TAI93 (1993-01-01 起算、含閏秒)。每個檔前後各多 OverlapScans (30) 條與相鄰檔重疊的掃描，
      多檔合併時依掃描時間去掉重複。
89 GHz：[swath] merge_ab = true 時 A、B 合併成一個頻道 (89V、89H)；false 時分列為 89V-A、89V-B…
頻道  依檔頭 AMSRChannel 列出的頻率 (AMSR-E 沒有 7.3 GHz：檔案中的 7.3 GHz 資料集只是 6.9 GHz 的複本，不列)。
"""
import datetime
import os
import re

import h5py
import numpy as np

from satima.core.i18n import tr
from satima.mw.reader import Channel, Geo, Scene, SwathReader, band_name, clean_latlon, tai93_to_unix
from satima.mw.swath import unit

# (資料集中的頻率字樣, 中心頻率 GHz, 地理定位組 = co-registration 參數的名稱)
_LOW = (('6.9', 6.925, '6G'), ('7.3', 7.3, '7G'), ('10.7', 10.65, '10G'), ('18.7', 18.7, '18G'),
        ('23.8', 23.8, '23G'), ('36.5', 36.5, '36G'))
_FILL = 65534           # 65534、65535 為無效值


def _attr(f, key):
    v = f.attrs[key]
    v = v[0] if isinstance(v, np.ndarray) else v
    return v.decode('latin1') if isinstance(v, bytes) else str(v)


def _time(text):
    return datetime.datetime.strptime(text[:19], '%Y-%m-%dT%H:%M:%S')


def coregister(lat89, lon89, a1, a2):
    """由 89A 的經緯度推算低頻頻段的經緯度 (見模組說明)。回傳 (lat, lon)，形狀 (掃描, 89A 點數 // 2)。"""
    p = unit(lat89, lon89)
    p1, p2 = p[:, 0::2], p[:, 1::2]
    n = min(p1.shape[1], p2.shape[1])
    p1, p2 = p1[:, :n], p2[:, :n]
    with np.errstate(invalid='ignore', divide='ignore'):
        ez = np.cross(p1, p2)
        ez /= np.linalg.norm(ez, axis=-1, keepdims=True)
        ey = np.cross(ez, p1)
        th = np.arccos(np.clip((p1 * p2).sum(axis=-1), -1., 1.))[..., None]
        q = np.cos(a2 * th) * (np.cos(a1 * th) * p1 + np.sin(a1 * th) * ey) + np.sin(a2 * th) * ez
        lat = np.rad2deg(np.arcsin(np.clip(q[..., 2], -1., 1.)))
        lon = np.rad2deg(np.arctan2(q[..., 1], q[..., 0]))
    return lat.astype(np.float32), lon.astype(np.float32)


class AMSRReader(SwathReader):
    """子類別設定 NAME、SENSOR、SATELLITE、DESCRIPTION、PATTERNS、FILENAME (群組 1~3 = 時間、路徑、A/D)。"""
    SATELLITE = ''
    SOURCE = 'JAXA'
    FILENAME = None

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = self.FILENAME.match(os.path.basename(name))
            if not m:
                continue
            path = os.path.join(folder, name)
            with h5py.File(path, 'r') as f:
                start, end = _time(_attr(f, 'ObservationStartDateTime')), _time(_attr(f, 'ObservationEndDateTime'))
                listed = [c.strip() for c in _attr(f, 'AMSRChannel').split(',')] if 'AMSRChannel' in f.attrs else None
            out.append(Scene(self.SATELLITE, self.SENSOR, start, end, file=path, orbit=m.group(2) + m.group(3).upper(),
                             format='L1B', extra={'listed': listed}))
        return out

    def channels(self, scene):
        out = []
        ds = 'Brightness Temperature (%sGHz,%s)'
        listed = scene.extra.get('listed')
        for tag, freq, grp in _LOW:
            if listed is not None and not any(c.startswith(('%g' % freq, tag)) for c in listed):
                continue                        # 檔頭沒有列的頻率 (AMSR-E 的 7.3 GHz)
            for pol in 'VH':
                out.append(Channel(band_name(freq, pol), freq, pol, grp, key=(ds % (tag, pol),)))
        for pol in 'VH':
            a = 'Brightness Temperature (89.0GHz-A,%s)' % pol
            b = 'Brightness Temperature (89.0GHz-B,%s)' % pol
            if self.cfg.swath.merge_ab:
                out.append(Channel(band_name(89.0, pol), 89.0, pol, '89A+89B', key=(a, b)))
            else:
                out.append(Channel(band_name(89.0, pol, suffix='-A'), 89.0, pol, '89A', key=(a,)))
                out.append(Channel(band_name(89.0, pol, suffix='-B'), 89.0, pol, '89B', key=(b,)))
        return out

    def geo(self, scene, group):
        out = []
        with h5py.File(scene.file, 'r') as f:
            time = tai93_to_unix(f['Scan Time'][:])
            for g in group.split('+'):
                src = g if g in ('89A', '89B') else '89A'
                lat = f['Latitude of Observation Point for %s' % src][:]
                lon = f['Longitude of Observation Point for %s' % src][:]
                lat, lon = clean_latlon(lat, lon)
                if g != src:
                    a1 = self._param(f, 'CoRegistrationParameterA1', g)
                    a2 = self._param(f, 'CoRegistrationParameterA2', g)
                    lat, lon = coregister(lat, lon, a1, a2)
                out.append(Geo(lat, lon, time, (scene.file, g)))
        return out

    @staticmethod
    def _param(f, key, band):
        """'6G-1.51820,7G-1.40210,...' -> 該頻段的值 (值本身可能為負：'6G--0.14760')。"""
        for item in _attr(f, key).split(','):
            name, _, value = item.strip().partition('-')
            if name == band:
                return float(value)
        raise ValueError(tr('%s 沒有 %s 的參數', '%s has no parameter for %s') % (key, band))

    def tb(self, scene, channel):
        out = []
        with h5py.File(scene.file, 'r') as f:
            for name in channel.key:
                d = f[name]
                raw = d[:]
                scale = float(np.ravel(d.attrs.get('SCALE FACTOR', 0.01))[0])
                data = raw.astype(np.float32) * np.float32(scale)
                data[(raw >= _FILL) | (raw == 0)] = np.nan
                out.append(data)
        return out

