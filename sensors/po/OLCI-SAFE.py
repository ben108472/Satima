"""
OLCI-SAFE — Sentinel-3A/B OLCI L1B (EFR / ERR)，SAFE 格式
原始腳本：SatimaV3-S3-OLCI.py

所需檔案
  SAFE 產品 S3A_OL_1_EFR____<開始>_<結束>_..._.SEN3 (資料夾內為 xfdumanifest.xml + 各頻道 NetCDF)
  可為 .SEN3.tar、.SEN3.zip 或已解壓的 .SEN3 資料夾

產品
  1–21 單頻道 (TOA 反射率 πL/(F0·cosθs))
  0 真實色 (Oa08/Oa06/Oa03)：1.RGB 原始合成 2.CORR (sunz + CREFL + 銳化 + CIRA) 3.WV (NASA Worldview)

與原腳本的差異
  - 不再解壓縮：直接從 zip / tar 讀取 NetCDF (原本每次整包解壓到目前目錄)，也接受已解壓的 .SEN3 資料夾
  - CREFL / sunz / 增強曲線的共用部分移到 satima/truecolor.py (OLCI 特有：ECMWF 臭氧、
    海平面氣壓與地形高度的 Rayleigh 縮放、unsharp mask 銳化、CIRA 白點正規化，均保留)
  - 中心時間：取最接近輸入中心點的像元，用該列的 time_stamp (GPS 時間，已換成 UTC)
    (原本以檔名起訖時間依裁切範圍的中間列線性推算)
"""
from datetime import datetime

import numpy as np
from scipy.interpolate import RectBivariateSpline
from scipy.ndimage import gaussian_filter

from satima.po import SatReader, ui, archive, obstime, truecolor as tc
from satima.core.i18n import tr

# ============================================================
# True Color 設定
#   R = Oa08 (665 nm), G = Oa06 (560 nm), B = Oa03 (442.5 nm)
# ============================================================
TC_BANDS = (8, 6, 3)
TC_WAVELENGTH = (0.665, 0.560, 0.4425)          # um
TC_O3_COEF = (0.050, 0.095, 0.0035)             # 臭氧吸收係數 (cm-atm)^-1
TC_VERSIONS = {1: 'RGB', 2: 'CORR', 3: 'WV'}
TC_VERSION_DESC = {1: 'RGB (raw composite)', 2: 'CORR (sunz + CREFL + sharpening + CIRA stretch)',
                   3: 'WV (NASA Worldview style)'}
TC_DESC = {'RGB': 'RGB', 'CORR': 'CORR', 'WV': 'WORLDVIEW'}
# True color stack 的圖層索引
L_R, L_G, L_B, L_SZA, L_VZA, L_RAA, L_TSCALE, L_O3 = range(8)
SCALE_HEIGHT = tc.SCALE_HEIGHT     # m, CREFL 所用 Rayleigh 標高
DU_KGM2 = 2.1415e-5                # 1 DU = 2.1415e-5 kg/m2


# ------------------------------------------------------------
# 校正 / 增強函式 (共用版本見 satima/truecolor.py，保留原名稱)
# ------------------------------------------------------------
rayleigh_tau = tc.rayleigh_tau     # Rayleigh 光學厚度 (Hansen & Travis 1974)，海平面標準大氣
_csalbr = tc.csalbr                # Rayleigh 大氣球面反照率
_chand = tc.chand                  # Rayleigh 路徑反射率


def sunz_corr(sza, limit=88., max_sza=95.):
    '''太陽天頂角校正因子 (同 satpy sunz_corrected)：1/cos(SZA)，
    超過 limit 後平滑遞減，到 max_sza 時歸零，避免晨昏線附近爆亮。'''
    return tc.sunz_corr(sza, limit, max_sza)


def crefl(refl, sza, vza, raa, tscale, o3, wl, o3coef):
    '''CREFL Rayleigh + 臭氧校正 (MODIS CREFL 演算法套用於 OLCI 波段)。
    refl: 已做 sunz 的 TOA 反射率；tscale: 地面氣壓 / 1013.25 (含地形高度)；o3: cm-atm'''
    mus = np.clip(np.cos(np.radians(sza)), 1e-3, 1.)
    muv = np.clip(np.cos(np.radians(vza)), 1e-3, 1.)
    taur = rayleigh_tau(wl) * tscale
    to3 = np.exp(-(1. / mus + 1. / muv) * o3 * o3coef)
    return tc.crefl_core(refl, mus, muv, raa, taur, t_o3=to3)


def sharpen(rgb, sigma=1.0, amount=0.6):
    '''亮度比值銳化 (unsharp mask 套在亮度上，再按比例回乘 RGB，避免色邊)。
    OLCI 各波段同為 300 m，無高解析波段可做 ratio sharpening，故以 unsharp mask 代替。'''
    with np.errstate(invalid='ignore'):
        lum = np.mean(rgb, axis=-1)
    valid = np.isfinite(lum)
    filled = np.where(valid, lum, 0.)
    weight = gaussian_filter(valid.astype(np.float32), sigma)
    blur = gaussian_filter(filled, sigma) / np.maximum(weight, 1e-6)
    lum_sharp = lum + amount * (lum - blur)
    with np.errstate(invalid='ignore', divide='ignore'):
        ratio = np.where(lum > 1e-4, lum_sharp / lum, 1.)
    ratio = np.clip(np.nan_to_num(ratio, nan=1.), 0.5, 2.)
    return rgb * ratio[..., None]


def cira_stretch(refl):
    '''CIRA 對數拉伸 (同 satpy cira_stretch)，輸入反射率 0-1。
    標準 CIRA 把反射率 1.0 映射到 0.83，雲頂偏灰；這裡再除以該值，使反射率 1.0 = 滿白。'''
    return tc.cira_stretch(refl, normalize_white=True)


def worldview_stretch(refl):
    '''NASA Worldview / MODIS CREFL 真實色非線性拉伸：
    反射率 0-1.1 線性映射到 0-255，再以分段線性 [0,30,60,120,190,255]->[0,110,160,210,240,255]。'''
    return tc.worldview_stretch(refl)


def make_truecolor(stack, version):
    '''由 true color stack (rows, cols, 8) 產生 0-1 RGB 影像。
    RGB : 純合成 (pi*L/F0，未做任何校正，線性 0-1)
    CORR: sunz + CREFL + sharpening + CIRA stretch
    WV  : sunz + CREFL + sharpening + NASA Worldview 拉伸'''
    toa = stack[..., L_R:L_B + 1].astype(np.float64)
    sza = stack[..., L_SZA]
    night = ~(sza < 90.)

    if version == 'RGB':
        rgb = np.clip(toa, 0., 1.)
    else:
        refl = toa * sunz_corr(sza)[..., None]
        corr = np.empty_like(refl)
        with np.errstate(invalid='ignore', divide='ignore'):
            for k in range(3):
                corr[..., k] = crefl(refl[..., k], sza, stack[..., L_VZA], stack[..., L_RAA],
                                     stack[..., L_TSCALE], stack[..., L_O3], TC_WAVELENGTH[k], TC_O3_COEF[k])
        corr = np.clip(corr, 0., 1.2)
        corr = sharpen(corr)
        rgb = cira_stretch(corr) if version == 'CORR' else worldview_stretch(corr)

    rgb = np.array(rgb, dtype=np.float64)
    rgb[night] = np.nan
    return rgb


# ===== 單一檔案 =====
class OlciFile:

    def __init__(self, fullname, channel=None, time=None):
        self.fullname = fullname
        self.channel = channel
        self.satellite = 'Sentinel-%s' % (self.fullname[1:3])
        if time:
            self.time = datetime.strptime(time, '%Y%m%dT%H%M%S')
        else:
            self.time = datetime.strptime(self.fullname[16:31], '%Y%m%dT%H%M%S')
        try:                    # 檔名第 2 個時間是結束時間：S3A_xx_1_xxx____<起>_<訖>_…
            self.end = datetime.strptime(self.fullname[32:47], '%Y%m%dT%H%M%S')
        except ValueError:
            self.end = None

    def __lt__(self, other):
        return self.time < other.time

    def interpolate_matrix(self, matrix, new_shape, order=3):
        # 計算原始矩陣的網格座標
        x = np.arange(matrix.shape[1])
        y = np.arange(matrix.shape[0])

        # 建立插值函數
        interpolator = RectBivariateSpline(y, x, matrix, kx=order, ky=order)

        # 產生新的網格座標
        new_y = np.linspace(0, matrix.shape[0] - 1, new_shape[0])
        new_x = np.linspace(0, matrix.shape[1] - 1, new_shape[1])

        # 執行插值
        interpolated_matrix = interpolator(new_y, new_x, grid=True)

        return interpolated_matrix

    def _dataset(self, key):
        '''開啟 SAFE 產品中名稱包含 key 的 NetCDF (zip / tar 直接讀進記憶體，不解壓縮)。'''
        if getattr(self, 'safe', None) is None:
            self.safe = archive.SafeProduct(self.fullname)
        return self.safe.dataset(key)

    def _tie(self, ncobj, var, shape, fill):
        '''讀取 tie-point 變數並線性內插到全解析度 (float32)。'''
        arr = np.ma.filled(ncobj.variables[var][:].astype(np.float64), np.nan)
        arr = np.where(np.isfinite(arr), arr, fill)
        return self.interpolate_matrix(arr, shape, order=1).astype(np.float32)

    def _toa(self, channel, inst):
        '''pi * L / F0 (未除 cos(SZA))，F0 依 detector_index 逐像素查表。'''
        nc_obj = self._dataset('Oa%02d_radiance' % channel)
        rad = np.ma.filled(nc_obj.variables['Oa%02d_radiance' % channel][:].astype(np.float32), np.nan)
        nc_obj.close()
        det = np.ma.filled(inst.variables['detector_index'][:], -1).astype(np.int32)
        flux = np.ma.filled(inst.variables['solar_flux'][channel - 1, :].astype(np.float32), np.nan)
        f0 = np.where(det >= 0, flux[np.clip(det, 0, None)], np.nan)
        return (rad * np.pi / f0).astype(np.float32)

    def extract(self):
        print(self.fullname)
        if self.channel == 0:
            return self.extract_truecolor()

        print('Oa%02d_radiance' % (self.channel))
        inst = self._dataset('instrument_data')
        toa = self._toa(int(self.channel), inst)
        inst.close()

        # Correct (sun zenith)
        geom = self._dataset('tie_geometries')
        sun_zenith = np.ma.filled(geom.variables['SZA'][:].astype(np.float64), np.nan)
        geom.close()
        sun_zenith = self.interpolate_matrix(sun_zenith, toa.shape)

        with np.errstate(invalid='ignore', divide='ignore'):
            data = toa / np.cos(np.radians(sun_zenith))

        print(tr('最大反射率', 'max reflectance'), np.nanmax(data))
        return data

    def extract_truecolor(self):
        '''回傳 (rows, cols, 8) float32：R/G/B 的 pi*L/F0、SZA、VZA、相對方位角、
        Rayleigh 氣壓縮放因子、臭氧 (cm-atm)。校正延後到裁切後才做以節省運算。'''
        print(tr('真實色頻道：Oa%02d / Oa%02d / Oa%02d', 'True color bands: Oa%02d / Oa%02d / Oa%02d') % TC_BANDS)
        inst = self._dataset('instrument_data')
        bands = [self._toa(ch, inst) for ch in TC_BANDS]
        inst.close()
        shape = bands[0].shape

        geom = self._dataset('tie_geometries')
        sza = self._tie(geom, 'SZA', shape, 90.)
        vza = self._tie(geom, 'OZA', shape, 0.)
        saa = np.radians(np.ma.filled(geom.variables['SAA'][:].astype(np.float64), np.nan))
        oaa = np.radians(np.ma.filled(geom.variables['OAA'][:].astype(np.float64), np.nan))
        geom.close()
        # 方位角差以 cos/sin 分量內插，避免 0/360 跳躍
        dphi = np.nan_to_num(saa - oaa, nan=0.)
        raa = np.degrees(np.arctan2(self.interpolate_matrix(np.sin(dphi), shape, order=1),
                                    self.interpolate_matrix(np.cos(dphi), shape, order=1))).astype(np.float32)

        meteo = self._dataset('tie_meteo')
        slp = self._tie(meteo, 'sea_level_pressure', shape, 1013.25)
        o3 = self._tie(meteo, 'total_ozone', shape, 0.319 * 1000. * DU_KGM2) / (1000. * DU_KGM2)
        meteo.close()

        geo = self._dataset('geo_coordinates.nc')
        alt = np.ma.filled(geo.variables['altitude'][:].astype(np.float32), 0.)
        geo.close()
        tscale = slp / 1013.25 * np.exp(-np.clip(alt, 0, None) / SCALE_HEIGHT)

        return np.dstack(bands + [sza, vza, raa, tscale.astype(np.float32), o3.astype(np.float32)])

    def geocoord(self):
        nc_obj = self._dataset('geo_coordinates.nc')
        lats = nc_obj.variables['latitude'][:]
        lons = nc_obj.variables['longitude'][:]
        nc_obj.close()
        print(np.amax(lats), np.amin(lats), np.amax(lons), np.amin(lons))
        return lats, lons

    def plot_rgb(self, canvas, lats, lons, rgb):
        '''畫 RGB (NaN 為透明)，實作見 satima/canvas.py。'''
        canvas.rgb(lons, lats, rgb)

    def plot(self, canvas, lats, lons, data, version=None):
        print(data.shape)
        if data.size == 0:
            return ''
        if self.channel == 0:
            print(tr('製作真實色：%s', 'Making true color: %s') % version)
            self.plot_rgb(canvas, lats, lons, make_truecolor(data, version))
            return version
        data = np.where(data > 100, np.nanmax(data), data)
        canvas.mesh(lons, lats, data, cmap='gray', vmin=0.0, vmax=1.0)
        return ''

    def row_times(self):
        """每一列的觀測時間 (UTC)：time_coordinates.nc 的 time_stamp。
        time_stamp 為 GPS 時間 (比 UTC 快 18 s；與全域屬性 start_time 及產品內太陽天頂角對照確認)。"""
        nc_obj = self._dataset('time_coordinates.nc')
        t = obstime.gps2000(nc_obj.variables['time_stamp'][:])
        nc_obj.close()
        return t


# ===== 讀取器 =====
class OlciReader(SatReader):
    NAME = 'OLCI-SAFE'
    SENSOR = 'OLCI'
    DESCRIPTION = 'Sentinel-3 OLCI L1B (SAFE .SEN3 tar/zip/資料夾)'
    DESCRIPTION_EN = 'Sentinel-3 OLCI L1B (SAFE .SEN3 tar/zip/folder)'
    SOURCE = 'Copernicus'
    PATTERNS = ('S*_OL*tar', 'S*_OL*zip', 'S*_OL*.SEN3')
    MULTI_FILE = True

    @classmethod
    def find(cls, folder=None):
        return archive.prefer_archives(super().find(folder))

    def make_file(self, name):
        # S3A_OL_1_EFR____20181001T211249_20181001T211549_20181003T063018_0179_036_214_2700_MAR_O_NT_003.SEN3.tar
        return OlciFile(name, channel=None, time=name[16:31])

    def listing(self, f):
        return (f.satellite, 'OLCI', f.time)

    def select_band(self):
        self.channel = ui.ask_int(tr('輸入頻道編號 (0 TRUE COLOR | 1-21 VIS (0.4-1.02))：',
                                     '0 TRUE COLOR | 1-21 VIS (0.4-1.02) : '), range(0, 22))
        for f in self.files:
            f.channel = self.channel

    def select_variants(self):
        if self.channel == 0:
            return [TC_VERSIONS[k] for k in ui.ask_menu('', TC_VERSION_DESC, tr('選擇真實色版本 (可複選，0 = 全部)：',
                                                                                'Select true color version: '))]
        return [None]

    def prepare(self):
        # 讀檔 + 裁切只做一次，各版本共用
        datas, lats, lons = [], [], []
        for f in self.files:
            datas.append(f.extract())
            la, lo = f.geocoord()
            lats.append(la)
            lons.append(lo)
        data = self.concat(datas)
        lats = self.concat(lats)
        lons = self.concat(lons)
        rowtime = np.concatenate([f.row_times() for f in self.files])
        if len(self.files) > 1:
            print(tr('已合併 %d 個片段：', 'Merged %d granules:') % len(self.files), np.shape(lats))
        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        self.data = self.crop(data)
        del data
        self.time = self.center_time(rowtime)
        return True

    def slice(self, georange, lats, lons):
        return self.slice_swath(lats, lons)

    def plot(self, canvas, version):
        cm = self.files[0].plot(canvas, self.lats, self.lons, self.data, version)
        if self.channel == 0:
            return dict(band='TrueColor', colorbar=cm, label='TRUE COLOR %s' % TC_DESC[version])
        return dict(band='B%02d' % self.channel, colorbar=cm, label='BAND%02d' % self.channel)


READER = OlciReader
