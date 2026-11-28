"""
VIRR-HDF5 — 風雲三號 A/B/C 可見光紅外掃描輻射計 VIRR L1 (NSMC，HDF5)
* 沒有樣本，照 NSMC 的格式說明 (FY-3A/B L1 Product Format、FY-3C 扫描辐射计 L1 数据格式说明) 與 satpy virr_l1b 寫，未驗證。

所需檔案 (同一目錄；一個片段 5 分鐘 = 1800 條掃描，每條 2048 點，1.1 km)
  FY3X_VIRRX_GBAL_L1_YYYYMMDD_HHmm_1000M_MS.HDF    資料 (FY-3A/B 檔內也有經緯度、角度、時間)
  FY3X_VIRRX_GBAL_L1_YYYYMMDD_HHmm_GEOXX_MS.HDF    FY-3C 的定位檔 (Geolocation/、Timedata/)

頻道 (10 個)
  反射  EV_RefSB (7)：1 (0.63)、2 (0.865)、6 (1.6)、7 (0.455)、8 (0.505)、9 (0.555)、10 (1.36 µm)
        反射率 % = DN x k + b (全域屬性 RefSB_Cal_Coefficients，每頻道 k、b 一組；FY-3B 檔內的值有誤，用 satpy 的係數)
  紅外  EV_Emissive (3)：3 (3.74)、4 (10.8)、5 (12.0 µm)
        輻亮度 = DN x Emissive_Radiance_Scales + Emissive_Radiance_Offsets (逐掃描線，mW m-2 sr-1 cm)；
        亮溫 = Planck^-1(輻亮度, Emissive_Centroid_Wave_Number)；FY-3A/B 的屬性名稱拼成 Emmisive_…
        Emissive_BT_Coefficients (3 x 2) 若能判斷出 A、B 的排列 (A 約 1、|B| < 5) 就以 Tbb = A Te + B 修正
時間  Day_Count + Millisecond_Count (Msec_Count)：2000-01-01 12:00 UTC 起的日數與毫秒

產品  單頻道 (反射率灰階 0~1、亮溫色表) / 真實色 (1、9、7 頻道；樣式同 MERSI)
"""
import datetime
import os
import re

import h5py
import numpy as np

import fy3_common as fc
from satima.core import segments as seg
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr
from satima.po import SatReader, Segment, ui, swath

_NAME = re.compile(r'^(FY3[A-Z])_VIRRX_GBAL_L1_(\d{8})_(\d{4})_([0-9A-Z]{5})_MS', re.I)
GEO_GROUPS = ('Geolocation', '')
TIME_GROUPS = ('Timedata', 'Geolocation', '')
# 頻道：(編號, 波長, 種類, 資料集內的索引)
CHANNELS = ((1, 0.63, 'ref', 0), (2, 0.865, 'ref', 1), (3, 3.74, 'bt', 0), (4, 10.8, 'bt', 1), (5, 12.0, 'bt', 2),
            (6, 1.6, 'ref', 2), (7, 0.455, 'ref', 3), (8, 0.505, 'ref', 4), (9, 0.555, 'ref', 5), (10, 1.36, 'ref', 6))
TC = {'R': 1, 'G': 9, 'B': 7}
# satpy virr_l1b：FY-3B 檔內的反射通道定標係數有誤，改用這組 (k, b)
FY3B_REF_COEFFS = (0.1264, -1.4320, 0.1353, -1.6236, 0.0919, -2.4821, 0.0938, -1.1494, 0.0857, -1.0280,
                   0.0803, -0.9636, 0.0630, -0.7628)
NOMINAL_WN = (2673.8, 925.9, 833.3)          # 備援：標稱中心波數 (3.74、10.8、12.0 µm)


def _text(c):
    n, wl, kind, idx = c
    return 'CH%02d(%s)' % (n, ('%.3f' % wl).rstrip('0').rstrip('.'))


class VirrPart:
    def __init__(self, name):
        m = _NAME.match(os.path.basename(name))
        if not m:
            raise ValueError(tr('不認得的 VIRR 檔名：%s', 'Unexpected VIRR file name: %s') % name)
        self.fullname = name
        self.satellite = 'FY-3' + m.group(1)[-1].upper()
        self.kind = m.group(4).upper()
        self.time = datetime.datetime.strptime(m.group(2) + m.group(3), '%Y%m%d%H%M')


class VirrGranule:
    def __init__(self, satellite, time):
        self.satellite = satellite
        self.time = time
        self.end = time + datetime.timedelta(minutes=5)
        self.files = {}

    def __lt__(self, other):
        return self.time < other.time

    def geo_file(self):
        """有經緯度的檔：GEOXX，沒有時資料檔本身 (FY-3A/B)。"""
        for k in ('GEOXX', '1000M'):
            if k in self.files:
                with h5py.File(self.files[k], 'r') as f:
                    if fc.find_ds(f, 'Latitude', GEO_GROUPS) is not None:
                        return self.files[k]
        return None


class VirrReader(SatReader):
    NAME = 'VIRR-HDF5'
    SENSOR = 'VIRR'
    SOURCE = 'NSMC'
    DESCRIPTION = '風雲三號 A/B/C VIRR L1 (NSMC，HDF5；未以真實檔案驗證)'
    DESCRIPTION_EN = 'Fengyun-3A/B/C VIRR L1 (NSMC, HDF5; not yet verified with real files)'
    PATTERNS = ('FY3?_VIRRX_GBAL_L1_*.HDF',)
    GRANULE_SECONDS = 300
    MULTI_FILE = True

    def make_file(self, name):
        return VirrPart(name)

    def segments(self, parts):
        grans = {}
        for p in sorted(parts, key=lambda p: p.fullname):
            g = grans.setdefault((p.satellite, p.time), VirrGranule(p.satellite, p.time))
            g.files[p.kind] = p.fullname
        by_sat = {}
        for g in grans.values():
            if '1000M' in g.files:
                by_sat.setdefault(g.satellite, []).append(g)
        out = []
        for sat, items in by_sat.items():
            for members, missing in seg.chain(items, lambda g: g.time, lambda g: g.end):
                out.append(Segment(sat, members[0].time, max(m.end for m in members), files=members,
                                   slots=len(members), missing=missing, product='1km', short=True))
        return out

    def merge_key(self, segment):
        return (segment.satellite,)

    def begin(self, segment):
        self.granules = sorted(segment.files)
        self.files = self.granules

    # ---- 選單 --------------------------------------------------------------
    def select_band(self):
        refs = [c for c in CHANNELS if c[2] == 'ref']
        bts = [c for c in CHANNELS if c[2] == 'bt']
        options = refs + bts + [None]
        groups = [('Reflectance', list(range(1, len(refs) + 1))),
                  ('Brightness temperature', list(range(len(refs) + 1, len(refs) + len(bts) + 1))),
                  ('RGB', [len(options)])]
        print('Products:')
        ui.print_columns([(t, ['%2d.%s' % (k, _text(options[k - 1]) if options[k - 1] else 'TrueColor (1km)')
                               for k in ks]) for t, ks in groups])
        k = ui.ask_int(tr('輸入產品編號：', 'Input the product: '), range(1, len(options) + 1))
        self.channel = options[k - 1]
        self.res = '1km'

    def select_variants(self):
        c = self.channel
        if c is None:
            return ui.ask_menu('Style:', fc.TC_STYLE_DESC, tr('選擇樣式 (可複選，0 = 全部)：', 'Select style: '))
        if c[2] == 'bt':
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    # ---- 讀檔 --------------------------------------------------------------
    def _rows(self, g):
        """這個片段在地圖內的列範圍 (r0, r1) 與行範圍 (c0, c1)；不在地圖內 None。"""
        path = g.geo_file()
        if path is None:
            raise IOError(tr('%s %s 找不到 VIRR 定位檔 (GEOXX)', '%s %s: no VIRR geolocation file (GEOXX)')
                          % (g.satellite, g.time.strftime('%Y/%m/%d %H:%M')))
        with h5py.File(path, 'r') as f:
            lat_ds, lon_ds = fc.find_ds(f, 'Latitude', GEO_GROUPS), fc.find_ds(f, 'Longitude', GEO_GROUPS)
            shape = lat_ds.shape
            lat, lon = fc.clean_latlon(lat_ds[::10, ::10], lon_ds[::10, ::10])
        latmin, latmax, lonmin, lonmax = self.georange
        m = 0.3
        lonw = swath.wrap_lon(lon, self.center[1])
        with np.errstate(invalid='ignore'):
            inside = (lat > latmin - m) & (lat < latmax + m) & (lonw > lonmin - m) & (lonw < lonmax + m)
        rows = np.flatnonzero(inside.any(axis=1))
        if rows.size == 0:
            return None
        cols = np.flatnonzero(inside.any(axis=0))
        return (max(0, (rows[0] - 1) * 10), min(shape[0], (rows[-1] + 2) * 10),
                max(0, (cols[0] - 1) * 10), min(shape[1], (cols[-1] + 2) * 10), shape)

    def _read_geo(self, g, r0, r1, c0, c1):
        with h5py.File(g.geo_file(), 'r') as f:
            lat = fc.find_ds(f, 'Latitude', GEO_GROUPS)[r0:r1, c0:c1]
            lon = fc.find_ds(f, 'Longitude', GEO_GROUPS)[r0:r1, c0:c1]
        return fc.clean_latlon(lat, lon)

    def _read_times(self, g, r0, r1):
        for path in (g.files.get('GEOXX'), g.files.get('1000M')):
            if not path:
                continue
            with h5py.File(path, 'r') as f:
                dc = fc.find_ds(f, 'Day_Count', TIME_GROUPS)
                ms = fc.find_ds(f, 'Millisecond_Count', TIME_GROUPS)
                if ms is None:
                    ms = fc.find_ds(f, 'Msec_Count', TIME_GROUPS)
                if dc is not None and ms is not None:
                    t = fc.day_ms_times(dc[:], ms[:].astype(np.float64) * fc.scaled(ms)[0])
                    return t[r0:r1]
        start = np.datetime64(g.time, 'us')
        return start + (np.arange(r0, r1) * (300e6 / 1800)).astype(np.int64).astype('timedelta64[us]')

    def _read_angles(self, g, r0, r1, c0, c1, lat, lon, times):
        out = {}
        with h5py.File(g.geo_file(), 'r') as f:
            for key, name in (('sza', 'SolarZenith'), ('saz', 'SolarAzimuth'), ('vza', 'SensorZenith'),
                              ('vaz', 'SensorAzimuth')):
                ds = fc.find_ds(f, name, GEO_GROUPS)
                if ds is not None:
                    raw = ds[r0:r1, c0:c1]
                    slope, inter = fc.scaled(ds)
                    v = raw.astype(np.float64) * slope + inter
                    fv = fc.attr_array(ds, 'FillValue')
                    if fv is None:
                        fv = fc.attr_array(ds, '_FillValue')
                    if fv is not None and fv.size:
                        v[raw == fv[0]] = np.nan
                    out[key] = v
            h = fc.find_ds(f, 'Height', GEO_GROUPS)
            if h is None:
                h = fc.find_ds(f, 'DEM', GEO_GROUPS)
            if h is not None:
                out['height'] = np.nan_to_num(np.clip(h[r0:r1, c0:c1].astype(np.float64), 0., 9000.), nan=0.)
        if 'sza' not in out or 'saz' not in out:
            out['sza'], out['saz'] = fc.solar_angles(times, lat, lon)
        if 'vza' not in out or 'vaz' not in out:
            out['vza'] = np.zeros_like(out['sza'])
            out['vaz'] = np.zeros_like(out['sza'])
        return out

    @staticmethod
    def _bt_coeffs(f, n):
        v = None
        for key in ('Emissive_BT_Coefficients', 'Emmisive_BT_Coefficients', 'Emisive_BT_Coefficients'):
            v = fc.attr_array(f, key)
            if v is not None:
                break
        if v is None or v.size < 2 * n:
            return None
        ok = lambda a, b: np.all((a > 0.95) & (a < 1.05)) and np.all(np.abs(b) < 5.)
        if ok(v[0:2 * n:2], v[1:2 * n:2]):
            return v[0:2 * n:2], v[1:2 * n:2]
        if ok(v[:n], v[n:2 * n]):
            return v[:n], v[n:2 * n]
        return None

    def _read_channel(self, g, c, r0, r1, c0, c1):
        n, wl, kind, idx = c
        with h5py.File(g.files['1000M'], 'r') as f:
            name = 'EV_RefSB' if kind == 'ref' else 'EV_Emissive'
            ds = fc.find_ds(f, name)
            raw = ds[idx, r0:r1, c0:c1]
            vr = fc.attr_array(ds, 'valid_range')
            bad = (raw > (vr[1] if vr is not None and vr.size == 2 else 32767)) | (raw == 0)
            if kind == 'ref':
                coef = FY3B_REF_COEFFS if g.satellite == 'FY-3B' else fc.attr_array(f, 'RefSB_Cal_Coefficients')
                if coef is None or len(coef) < 14:
                    raise ValueError(tr('%s 沒有 RefSB_Cal_Coefficients', '%s has no RefSB_Cal_Coefficients')
                                     % g.files['1000M'])
                k, b = coef[2 * idx], coef[2 * idx + 1]
                out = (raw.astype(np.float64) * (k if k else 1.) + b) / 100.
            else:
                sc = fc.find_ds(f, 'Emissive_Radiance_Scales')[r0:r1, idx].astype(np.float64)[:, None]
                of = fc.find_ds(f, 'Emissive_Radiance_Offsets')[r0:r1, idx].astype(np.float64)[:, None]
                sc = np.where(sc == 0, 1., sc)
                rad = raw.astype(np.float64) * sc + of
                wn = None
                for key in ('Emissive_Centroid_Wave_Number', 'Emmisive_Centroid_Wave_Number',
                            'Emisive_Centroid_Wave_Number'):
                    v = fc.attr_array(f, key)
                    if v is not None and v.size >= 3:
                        wn = float(v[idx])
                        break
                if wn is None or not (0.8 < wn * wl / 1e4 < 1.2):
                    wn = NOMINAL_WN[idx]
                out = fc.planck_bt(rad, wn)
                ab = self._bt_coeffs(f, 3)
                if ab is not None:
                    out = ab[0][idx] * out + ab[1][idx]
                out = out - 273.15
        out[bad] = np.nan
        return out

    # ---- 資料 --------------------------------------------------------------
    def prepare(self):
        picks = []
        for g in self.granules:
            box = self._rows(g)
            if box is not None:
                picks.append((g, box))
        if not picks:
            return False
        c0 = min(b[2] for _, b in picks)
        c1 = max(b[3] for _, b in picks)
        lats, lons, times = [], [], []
        for g, (r0, r1, _, _, _) in picks:
            la, lo = self._read_geo(g, r0, r1, c0, c1)
            lats.append(la)
            lons.append(lo)
            times.append(self._read_times(g, r0, r1))
        lats, lons, times = np.concatenate(lats), np.concatenate(lons), np.concatenate(times)
        if self.slice_swath(lats, lons) is None:
            return False
        self.time = self.center_time(times)
        if self.channel is None:
            refl = {}
            byn = {c[0]: c for c in CHANNELS}
            for col, num in TC.items():
                refl[col] = self.crop(self.concat([self._read_channel(g, byn[num], b[0], b[1], c0, c1)
                                                   for g, b in picks]))
            angs, k = [], 0
            for g, b in picks:
                nr = b[1] - b[0]
                angs.append(self._read_angles(g, b[0], b[1], c0, c1, lats[k:k + nr], lons[k:k + nr], times[k:k + nr]))
                k += nr
            ang = {key: self.crop(np.concatenate([a[key] for a in angs])) for key in ('sza', 'saz', 'vza', 'vaz')}
            height = self.crop(np.concatenate([a['height'] for a in angs])) if all('height' in a for a in angs) else None
            with h5py.File(picks[0][0].files['1000M'], 'r') as f:
                des = fc.attr_array(f, 'EarthSun Distance Ratio')
            des = float(des[0]) if des is not None and 0.9 < des[0] < 1.1 else 1.
            wls = {col: byn[num][1] for col, num in TC.items()}
            self._tc = fc.TrueColor(refl, wls, des, (ang['sza'], ang['saz'], ang['vza'], ang['vaz']), height)
        else:
            self.data = self.crop(self.concat([self._read_channel(g, self.channel, b[0], b[1], c0, c1)
                                               for g, b in picks]))
            if self.channel[2] == 'bt':
                self.report_eyetemp(self.data)
        return True

    def plot(self, canvas, variant):
        if self.channel is None:
            canvas.rgb(self.lons, self.lats, self._tc.rgb(variant))
            return dict(band='TrueColor', colorbar=fc.TC_STYLES[variant], label='TrueColor')
        name = 'CH%02d' % self.channel[0]
        if self.channel[2] == 'ref':
            canvas.mesh(self.lons, self.lats, np.ma.masked_greater(np.ma.masked_invalid(self.data), 1.5),
                        cmap='gray', vmin=0., vmax=1.)
            return dict(band=name, colorbar='', label=name)
        cm = canvas.bt_mesh(self.lons, self.lats, self.data, variant)
        return dict(band=name, colorbar=cm, label=name)


READER = VirrReader
