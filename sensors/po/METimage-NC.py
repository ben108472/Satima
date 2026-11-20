"""
METimage-NC — Metop-SG METimage (VII) L1B (NetCDF，zip 壓縮)
原始腳本：SatimaV3-METOPSG-loop.py

所需檔案
  W_XX-EUMETSAT-Darmstadt,SAT,SGA1-VII-1B-RAD_C_EUMT_<處理>_G_O_<開始>_<結束>_C_N_T__.zip (或解壓後的 .nc)

產品
  1–20 單頻道 (1–11 反射率、12–20 亮溫；15/16 水氣色表)
  21 真實色 (VII-12/VII-8/VII-4)：1.General (只做 sunz，無 Rayleigh) 2.CIRA 3.NASA Worldview

與原腳本的差異
  - 頻道選單由頻道表產生 (原選單最後一個頻道誤印為 18)
  - 頻道 11 (VII-25, 2.25 µm 反射率) 原本被當成紅外線而出現色表選單 (畫出 13 張相同的灰階圖)，改為 12 起才有色表
  - 地圖改由共用繪圖 (預設 Basemap，設定檔可改回 Cartopy)
  - 檔名時間改用正規表示式解析；可直接讀解壓後的 .nc
  - CREFL 共用部分移到 satima/truecolor.py
  - 中心時間：取最接近輸入中心點的像元，用該掃描的 scan_start_time_utc
    (原本以檔名起訖時間依裁切範圍的中間列線性推算)
"""
import re
from datetime import datetime
from zipfile import ZipFile

import netCDF4 as nc
import numpy as np

from satima.po import SatReader, ui, obstime, truecolor as tc
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr

WAVELENGTHS_NM = (
    443, 555, 668, 752, 763, 865, 914, 1240, 1375, 1630, 2250,
    3740, 3959, 4050, 6725, 7325, 8540, 10690, 12020, 13345,
)

# 選單編號: (頻道名稱, 中心波長 µm, 類別)
CHANNELS = {
    1: ('VII-4', '0.443', 'VIS'), 2: ('VII-8', '0.555', 'VIS'), 3: ('VII-12', '0.668', 'VIS'),
    4: ('VII-15', '0.7515', 'NIR'), 5: ('VII-16', '0.7627', 'NIR'), 6: ('VII-17', '0.865', 'NIR'),
    7: ('VII-20', '0.914', 'NIR'), 8: ('VII-22', '1.240', 'NIR'),
    9: ('VII-23', '1.375', 'SWIR'), 10: ('VII-24', '1.63', 'SWIR'), 11: ('VII-25', '2.25', 'SWIR'),
    12: ('VII-26', '3.74', 'MWIR'), 13: ('VII-28', '3.959', 'MWIR'), 14: ('VII-30', '4.05', 'MWIR'),
    15: ('VII-33', '6.725', 'WV'), 16: ('VII-34', '7.325', 'WV'),
    17: ('VII-35', '8.54', 'TIR'), 18: ('VII-37', '10.69', 'TIR'), 19: ('VII-39', '12.02', 'TIR'),
    20: ('VII-40', '13.345', 'TIR'),
}
GROUP_LABEL = (('VIS', 'VIS Channel  '), ('NIR', 'NIR Channel  '), ('SWIR', 'SWIR Channel '),
               ('MWIR', 'MWIR Channel '), ('WV', 'WV Channel   '), ('TIR', 'TIR Channel  '))
REFL_MAX_CHANNEL = 11            # 1–11 反射率，12–20 亮溫
WV_CHANNELS = (15, 16)
_NAME_TIMES = re.compile(r'_G_O_(\d{14})_(\d{14})')
_NAME_SAT = re.compile(r'SG([AB]\d)')

# ---- VII L1B geolocation (EPS-SG VII Level 1B Product Format Specification) ----
# latitude/longitude/solar_zenith are stored on tie points (num_tie_points_alt, num_tie_points_act).
# Each scan has SCAN_ALT_TIE_POINTS tie rows; tie points are every TIE_POINTS_FACTOR pixels/lines
# (zone_size_act / zone_size_alt = 8). Interpolation is done per scan (no interpolation across the
# scan boundary, because of the bow-tie effect), same as the reference geotiepoints/satpy reader.
SCAN_ALT_TIE_POINTS = 4
TIE_POINTS_FACTOR = 8
MEAN_EARTH_RADIUS = 6371008.7714  # [m], IUGG


def _interp_axis(arr, n_out, factor, axis):
    """Linear interpolation from tie grid (0, factor, 2*factor, ...) to pixel grid 0..n_out-1."""
    n_tie = arr.shape[axis]
    pos = np.arange(n_out) / factor
    i0 = np.clip(np.floor(pos).astype(int), 0, n_tie - 2)
    w = pos - i0
    shape = [1] * arr.ndim
    shape[axis] = n_out
    w = w.reshape(shape)
    return np.take(arr, i0, axis=axis) * (1 - w) + np.take(arr, i0 + 1, axis=axis) * w


def tie_points_interpolation(tie, scan_alt_tie_points=SCAN_ALT_TIE_POINTS, factor=TIE_POINTS_FACTOR):
    """tie: (num_tie_points_alt, num_tie_points_act) -> (num_lines, num_pixels)."""
    n_tie_alt, n_tie_act = tie.shape
    if n_tie_alt % scan_alt_tie_points != 0:
        raise ValueError(tr('num_tie_points_alt 不是每個掃描 tie point 數的整數倍',
                            'num_tie_points_alt is not a multiple of tie points per scan'))
    n_scans = n_tie_alt // scan_alt_tie_points
    n_pixels = (n_tie_act - 1) * factor
    n_lines_scan = (scan_alt_tie_points - 1) * factor
    # across track
    act = _interp_axis(tie, n_pixels, factor, axis=1)
    # along track, scan by scan
    act = act.reshape(n_scans, scan_alt_tie_points, n_pixels)
    alt = _interp_axis(act, n_lines_scan, factor, axis=1)
    return alt.reshape(n_scans * n_lines_scan, n_pixels)


def tie_points_angle_interpolation(angle_tie):
    """Interpolate azimuth angles (deg) via sin/cos to avoid the 0/360 wrap."""
    r = np.deg2rad(angle_tie)
    s = tie_points_interpolation(np.sin(r))
    c = tie_points_interpolation(np.cos(r))
    return np.rad2deg(np.arctan2(s, c)) % 360.0


def tie_points_geo_interpolation(lon_tie, lat_tie):
    """Interpolate lon/lat in Cartesian space (safe near poles and across the dateline)."""
    lon_r, lat_r = np.deg2rad(lon_tie), np.deg2rad(lat_tie)
    x = np.cos(lat_r) * np.cos(lon_r)
    y = np.cos(lat_r) * np.sin(lon_r)
    z = np.sin(lat_r)
    x, y, z = (tie_points_interpolation(v) for v in (x, y, z))
    lon = np.rad2deg(np.arctan2(y, x)) % 360.0
    lat = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
    return lon, lat


# ---- True color ----
# Bands: R = vii_668 (index 2), G = vii_555 (index 1), B = vii_443 (index 0)
# Corrections (same chain as satpy's VII 'true_color' composite):
#   1. radiance -> TOA reflectance (Earth-Sun distance, band-averaged solar irradiance)
#   2. sun-zenith correction (satpy sunz_corrected: 1/cos(SZA), limited at 88 deg, fades out to 95 deg)
#   3. Rayleigh + gaseous absorption correction (NASA CREFL algorithm, 6S-based, as in satpy _crefl_utils)
# Styles: General = steps 1-2 only, linear (Himawari-style simple composite, keeps the Rayleigh-blue sky/ocean)
#         CIRA / Worldview = steps 1-3 + non-linear stretch
# CREFL has no official VII coefficients:
#   - Rayleigh optical depth: Hansen & Travis (1974) formula weighted by the VII ISRF stored in the L1B file
#   - O3 absorption: interpolated from MODIS/VIIRS CREFL values for neighbouring wavelengths
#   - H2O absorption: VIIRS CREFL coefficients of the matching band (M04 for 555, M05 for 668)
#   - sea level (no terrain-height correction of the Rayleigh optical depth)
TRUE_COLOR_CHANNEL = 21
TRUE_COLOR_BANDS = (2, 1, 0)
TC_STYLES = {1: 'General', 2: 'CIRA', 3: 'Worldview'}
TC_STYLE_DESC = {1: 'General(no Rayleigh)', 2: 'CIRA', 3: 'NASA Worldview'}
VII_CREFL_GAS = {   # band index: (aH2O, bH2O, aO3)
    0: (0.0, 0.0, 0.0030),               # vii_443
    1: (1.78644e-05, 0.867785, 0.0870),  # vii_555  (VIIRS M04 H2O)
    2: (0.000406601, 0.812659, 0.0460),  # vii_668  (VIIRS M05 H2O)
}
UO3_CLIM = 0.285   # climatological ozone [atm-cm] (CREFL VIIRS)
UH2O_CLIM = 2.93   # climatological water vapour [g/cm2] (CREFL)
MAXAIRMASS = tc.MAXAIRMASS
REFLMIN, REFLMAX = tc.REFLMIN, tc.REFLMAX

# 共用版本 (satima/truecolor.py)，保留原名稱
rayleigh_optical_depth = tc.rayleigh_tau
_csalbr = tc.csalbr
_chand = tc.chand


def crefl_correct(refl, mus, muv, phi, taur, ah2o, bh2o, ao3):
    """CREFL Rayleigh + O3/H2O correction of sun-zenith-corrected TOA reflectance (sea level)."""
    airmass = np.minimum(1.0 / mus + 1.0 / muv, MAXAIRMASS)
    t_o3 = tc.o3_transmittance(airmass, UO3_CLIM, ao3)
    t_h2o = tc.h2o_transmittance_viirs(airmass, UH2O_CLIM, ah2o, bh2o)
    corr = tc.crefl_core(refl, mus, muv, phi, taur, t_o3, t_h2o)
    return np.clip(corr, REFLMIN, REFLMAX)


def sunz_correction_factor(sza, limit=88.0, max_sza=95.0):
    """satpy sunz_corrected: 1/cos(SZA) up to `limit`, then fades to 0 at `max_sza`."""
    return tc.sunz_corr(sza, limit, max_sza)


def enhance_true_color(rgb, style):
    """rgb: output of MetimageFile.true_color(), (..., 6) = Rayleigh-corrected R,G,B + sunz-only R,G,B.
    style 1 = General, 2 = CIRA, 3 = NASA Worldview. Returns (..., 3) floats in 0~1."""
    rgb = np.nan_to_num(np.asarray(rgb, dtype=np.float64), nan=0.0)
    if style == 1:
        # General: sun-zenith corrected reflectance WITHOUT Rayleigh correction, plain linear 0~1
        return np.clip(rgb[..., 3:6], 0, 1)
    rgb = rgb[..., 0:3]
    if style == 2:
        # CIRA geostationary (ABI/AHI) true color: logarithmic stretch adapted to human vision
        out = tc.cira_stretch(rgb)
    else:
        # MODIS Rapid Response / NASA Worldview "Corrected Reflectance (True Color)"
        out = tc.worldview_stretch(rgb)
    return np.clip(out, 0, 1)


def parse_name(name):
    """檔名 -> (衛星, 開始時間, 結束時間)。"""
    m = _NAME_TIMES.search(name)
    if not m:
        raise ValueError(tr('不認得的 METimage 檔名：%s', 'Unexpected METimage file name: %s') % name)
    s = _NAME_SAT.search(name)
    satellite = 'METOP-SG-%s' % (s.group(1) if s else '?')
    return satellite, datetime.strptime(m.group(1), '%Y%m%d%H%M%S'), datetime.strptime(m.group(2), '%Y%m%d%H%M%S')


# ===== 單一檔案 =====
class MetimageFile:

    def __init__(self, fullname, channel=None, time=None):
        self.fullname = fullname
        self.channel = channel
        self.satellite, self.start, self.end = parse_name(fullname)
        self.time = datetime.strptime(time, '%Y%m%d%H%M%S') if time else self.start

    def __lt__(self, other):
        return self.time < other.time

    def extract(self):
        if self.fullname.lower().endswith('.zip'):
            archive = ZipFile(self.fullname)
            nc_names = [name for name in archive.namelist() if name.lower().endswith('.nc')]
            raw = archive.read(nc_names[0])
            archive.close()
            self.ncfile = nc.Dataset('memory.nc', mode='r', memory=raw)
        else:
            self.ncfile = nc.Dataset(self.fullname, mode='r')

        self.geolocate()

        if self.channel == TRUE_COLOR_CHANNEL:
            data = self.true_color()
        elif self.channel <= REFL_MAX_CHANNEL:
            data = self.rad_to_ref(self.channel - 1)
        else:
            data = self.rad_to_bt(self.channel - 1)

        return data

    def rad_to_ref(self, band):
        calib_data = self.ncfile.groups['data']['calibration_data']
        measure_data = self.ncfile.groups['data']['measurement_data']
        sate_data = self.ncfile.groups['status']['satellite']
        rad = measure_data.variables[f'vii_{WAVELENGTHS_NM[band]}'][:]
        des = sate_data.variables['earth_sun_distance_ratio'][:]
        solar_irrad = calib_data.variables['band_averaged_solar_irradiance'][:][band]
        with np.errstate(invalid='ignore', divide='ignore'):
            ref = np.pi * rad / (np.cos(np.deg2rad(self.solar_zenith)) * solar_irrad) * np.square(des)
        return ref

    def rad_to_bt(self, band):
        calib_data = self.ncfile.groups['data']['calibration_data']
        measure_data = self.ncfile.groups['data']['measurement_data']
        rad = measure_data.variables[f'vii_{WAVELENGTHS_NM[band]}'][:] * 10 ** 6

        A = calib_data.variables['bt_conversion_a'][:][band - 11]
        B = calib_data.variables['bt_conversion_b'][:][band - 11]
        central_wavelength = calib_data.variables['channel_cw_thermal'][:][band - 11]

        central_wavelength = central_wavelength * 1e-6

        c = 299792458
        h = 6.626069 * 10 ** -34
        k = 1.38065 * 10 ** -23

        with np.errstate(invalid='ignore', divide='ignore'):
            Tb = (h * c) / ((k * central_wavelength) * np.log(1 + 2 * h * c ** 2 / (central_wavelength ** 5 * rad))) * A + B

        return Tb - 273.15

    def geolocate(self, orthorectify=True):
        """Compute pixel lats/lons (0~360) and solar zenith from VII L1B tie points."""
        md = self.ncfile.groups['data']['measurement_data']
        lat_tie = np.ma.filled(md.variables['latitude'][:].astype(np.float64), np.nan)
        lon_tie = np.ma.filled(md.variables['longitude'][:].astype(np.float64), np.nan)
        sza_tie = np.ma.filled(md.variables['solar_zenith'][:].astype(np.float64), np.nan)

        lons, lats = tie_points_geo_interpolation(lon_tie, lat_tie)

        if orthorectify:
            # DEM terrain correction: delta (m) between orthorectified and ellipsoid position
            dn = np.ma.filled(md.variables['delta_lat_N_dem'][:].astype(np.float64), 0.0)
            de = np.ma.filled(md.variables['delta_lon_E_dem'][:].astype(np.float64), 0.0)
            lats = lats + np.rad2deg(dn / MEAN_EARTH_RADIUS)
            lons = (lons + np.rad2deg(de / (MEAN_EARTH_RADIUS * np.cos(np.deg2rad(lats))))) % 360.0

        self.lats, self.lons = lats, lons
        self.solar_zenith = tie_points_interpolation(sza_tie)

    def true_color(self):
        """Shape (lines, pixels, 6): [0:3] sunz + CREFL (Rayleigh/O3/H2O) corrected R,G,B,
        [3:6] sunz-corrected only R,G,B (for the General style)."""
        md = self.ncfile.groups['data']['measurement_data']
        calib = self.ncfile.groups['data']['calibration_data']
        des = float(self.ncfile.groups['status']['satellite'].variables['earth_sun_distance_ratio'][:])
        solar_irrad = calib.variables['band_averaged_solar_irradiance'][:]

        def tie(name):
            return np.ma.filled(md.variables[name][:].astype(np.float64), np.nan)

        sza = self.solar_zenith
        vza = tie_points_interpolation(tie('observation_zenith'))
        saa = tie_points_angle_interpolation(tie('solar_azimuth'))
        vaa = tie_points_angle_interpolation(tie('observation_azimuth'))   # pixel -> satellite
        mus = np.cos(np.deg2rad(sza))
        mus = np.where(mus > 0, mus, np.nan)
        muv = np.cos(np.deg2rad(vza))
        phi = saa - vaa
        sunz = sunz_correction_factor(sza)

        rgb, rgb_sunz = [], []
        for band in TRUE_COLOR_BANDS:
            wl = WAVELENGTHS_NM[band]
            # ISRF-weighted Rayleigh optical depth
            try:
                w = np.asarray(calib.variables[f'ISRF_{wl}'][:], dtype=np.float64)
                lam = np.asarray(calib.variables[f'ISRF_{wl}_wavelength'][:], dtype=np.float64)
                taur = np.sum(w * rayleigh_optical_depth(lam)) / np.sum(w)
            except KeyError:
                taur = rayleigh_optical_depth(wl / 1000.0)
            rad = np.ma.filled(md.variables[f'vii_{wl}'][:].astype(np.float64), np.nan)
            refl = np.pi * rad * des ** 2 / solar_irrad[band] * sunz
            with np.errstate(invalid='ignore', divide='ignore'):
                rgb.append(crefl_correct(refl, mus, muv, phi, taur, *VII_CREFL_GAS[band]))
            rgb_sunz.append(refl)
        return np.dstack(rgb + rgb_sunz).astype(np.float32)

    def geocoord(self):
        return self.lats, self.lons

    def plot(self, canvas, lats, lons, data, variant):
        """variant：真實色為樣式編號，紅外線/水氣為色表。回傳樣式/色表名稱。"""
        print(np.nanmax(lats), np.nanmax(lons), np.nanmin(lats), np.nanmin(lons))
        if data.size == 0:
            return ''
        if self.channel == TRUE_COLOR_CHANNEL:
            canvas.rgb(lons, lats, enhance_true_color(data, variant))
            return TC_STYLES[variant]
        if self.channel <= REFL_MAX_CHANNEL:
            canvas.mesh(lons, lats, data, cmap='gray', vmin=0.0, vmax=1.0)
            return ''
        return canvas.bt_mesh(lons, lats, data, variant)

    def row_times(self, nrows):
        """每一列的觀測時間 (UTC)：逐掃描的 scan_start_time_utc (2020-01-01 起的秒數)，
        同一掃描的列同一時間 (須先呼叫 extract)。"""
        md = self.ncfile.groups['data']['measurement_data']
        return obstime.per_scan(obstime.from_utc('2020-01-01', md.variables['scan_start_time_utc'][:], 's'), nrows)


# ===== 讀取器 =====
class MetimageReader(SatReader):
    NAME = 'METimage-NC'
    SENSOR = 'METimage'
    DESCRIPTION = 'Metop-SG METimage (VII) L1B (NetCDF, zip)'
    SOURCE = 'EUMETSAT'
    PATTERNS = ('*VII-1B*.zip', '*VII-1B*.nc')
    MULTI_FILE = True

    @staticmethod
    def sort_key(name):
        m = _NAME_TIMES.search(name)
        return m.group(1) if m else name

    @classmethod
    def find(cls, folder=None):
        """同一 granule 同時有 .zip 與解壓後的 .nc 時只列 .nc。"""
        names = super().find(folder)
        have = set(names)
        return [n for n in names if not (n.lower().endswith('.zip') and n[:-4] + '.nc' in have)]

    def make_file(self, name):
        return MetimageFile(name, channel=None)

    def listing(self, f):
        return (f.satellite, 'METimage', f.time)

    def select_band(self):
        for key, label in GROUP_LABEL:
            items = ['%d.%s(%sµm)' % (i, n, wl) for i, (n, wl, g) in CHANNELS.items() if g == key]
            print('%s: %s' % (label, ' '.join(items)))
        print('RGB          : %d.True Color (VII-12/VII-8/VII-4, sunz + Rayleigh corrected)' % TRUE_COLOR_CHANNEL)
        self.channel = ui.ask_int(tr('輸入頻道編號：', 'Channel : '), list(CHANNELS) + [TRUE_COLOR_CHANNEL])
        self.channel_name = 'TrueColor' if self.channel == TRUE_COLOR_CHANNEL else CHANNELS[self.channel][0]
        for f in self.files:
            f.channel = self.channel

    def select_variants(self):
        if self.channel == TRUE_COLOR_CHANNEL:
            return ui.ask_menu('', TC_STYLE_DESC, tr('選擇樣式 (可複選，0 = 全部)：', 'Select style: '))
        if self.channel in WV_CHANNELS:
            return ask_colormaps(self.cfg, 'wv')
        if self.channel > REFL_MAX_CHANNEL:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare(self):
        datas, lats, lons, rowtime = [], [], [], []
        for f in self.files:
            datas.append(f.extract())
            la, lo = f.geocoord()
            lats.append(la)
            lons.append(lo)
            rowtime.append(f.row_times(la.shape[0]))
        data = self.concat(datas)
        lats = np.concatenate(lats, axis=0)
        lons = np.concatenate(lons, axis=0)
        rowtime = np.concatenate(rowtime)
        if len(self.files) > 1:
            print(tr('已合併 %d 個片段：', 'Merged %d granules:') % len(self.files), lats.shape)
        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        self.data = self.crop(data)
        self.time = self.center_time(rowtime)
        if REFL_MAX_CHANNEL < self.channel < TRUE_COLOR_CHANNEL:
            self.report_eyetemp(self.data)
        return True

    def slice(self, georange, lats, lons):
        return self.slice_swath(lats, lons)

    def plot(self, canvas, variant):
        cm = self.files[0].plot(canvas, self.lats, self.lons, self.data, variant)
        return dict(band=self.channel_name, colorbar=cm, label=self.channel_name)


READER = MetimageReader
