"""
SGLI-HDF5 — GCOM-C SGLI L1B (HDF5)
原始腳本：SatimaV3-GCOMC-loop.py

所需檔案
  GC1SG1_<YYYYMMDDhhmm><D/N><路徑>_1BSG_<VNR|IRS><D/N><解析度碼>_<版本>.h5
  紅外線亮溫另需轉換表 conversion_table/sgli_r2tfl_pfm_400.dat
  (依序在資料目錄、Satima-PO.py 所在目錄尋找)

頻道
  VNR：VN01–VN11 (反射率，除以 cos(太陽天頂角))
  IRS：SW01–SW04 (短波紅外反射率，同上)、TI01 (10.8)、TI02 (12.0) (亮溫)
  {res} 欄位為該頻道的實際解析度 (例如 250M；250 m 產品中的 SW01/SW02/SW04 為 1000M)
  注意：VN07、VN10 為高增益 (海色) 頻道，雲上約在反射率 0.16 / 0.13 飽和；VN08 約 0.49 飽和；
        畫雲建議用 VN11 (868.5 nm，飽和約 1.1)

真實色 (VNR 檔選 0；比照 MODIS / OLCI 提供 1.RGB 2.CORR 3.WV)
  頻道  R = VN08 (673.5 nm)、G = VN05 (530 nm)、B = VN03 (443 nm)
        (JAXA JASMES SGLI NRT 規格與 SGLI 雲演算法論文使用的 L1B RGB 組合)
  RGB   TOA 反射率 rt 直接線性 0–1，不做任何校正 (同 MODIS / OLCI 的 RGB)
  CORR  太陽天頂角校正 (satima.truecolor.sunz_corr) + Rayleigh / 臭氧校正 + CIRA 增強
  WV    校正同 CORR，增強改為 NASA Worldview
  Rayleigh 光學厚度與臭氧係數取自 JAXA SGLI 大氣校正 ATBD (Table 3.1、Table 3-1，Bodhaine 1999，
  已考慮光譜響應)；Rayleigh 路徑反射率以 CREFL (6S) 公式計算 (官方使用查表)，與其他衛星相同。
  L1B 沒有臭氧與地面氣壓，使用臭氧氣候值 300 DU、海平面氣壓 (無地形修正)。
  感測器天頂角 / 方位角使用各頻道的 Sensor_zenith_VNxx / Sensor_azimuth_VNxx。
  飽和：VN08 在反射率約 0.49、VN05 約 0.67 飽和，亮雲上會先飽和而使雲偏青色；
        飽和像元以 VN11 (868.5 nm 低增益，雲在可見光–近紅外約為灰色) 補值，且不低於飽和值
        (此為 Satima 的處理，非 JAXA 官方作法)

與原腳本的差異
  - 經緯度內插改為地心直角座標二次樣條 (satima/tiepoint.py)：誤差約為原 griddata 的 1/14，速度快約 10 倍；
    太陽天頂角改為規則網格雙線性內插 (與原 griddata 線性相近，但快很多)，並排除 Error_DN
  - 缺值 (原始值 = Error_DN 65535，或 DN = 16383) 設為 NaN：原本 65535 & Mask 後變成 16383 被當成最大值，
    掃描範圍外約 14% 的像元 (VNR) 會畫成一片亮灰/白色
  - SW01–SW04 原本套用熱紅外的亮溫轉換表 (沒有物理意義)，且未套用 Mask (86% 像元帶旗標位元)；
    改為與 VN 相同的反射率
  - 與 tie 網格解析度不同的頻道 (250 m 產品中的 1 km SW 頻道、500 m TIR 產品) 依像元中心換算 tie 座標，
    原本直接用像元編號，定位會錯
  - 轉換表不限定在執行目錄
  - 中心時間：取最接近輸入中心點的像元，用該列的 Line_tai93 (已扣閏秒) 加上該頻道的焦平面時間差
    (原本以 Scene_start/end_time 依裁切範圍的中間列線性推算)
  - 注意：TI01 與 TI02 目前共用同一個轉換表 (同原程式)
"""
import os
from datetime import datetime
from glob import glob

import h5py
import numpy as np

from satima.po import SatReader, ui, obstime, tiepoint, truecolor as tc
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr

CONVERSION_TABLE = os.path.join('conversion_table', 'sgli_r2tfl_pfm_400.dat')
MISSING_DN = 16383        # Lt 的 DN (bit 0–13)：16383 = 缺值、16382 = 飽和 (見各 Lt 的屬性說明)
SATURATION_DN = 16382
VNR_CHANNELS = {i: 'VN%02d' % i for i in range(1, 12)}
IRS_CHANNELS = {1: 'SW01', 2: 'SW02', 3: 'SW03', 4: 'SW04', 5: 'TI01', 6: 'TI02'}
IRS_MENU = '1.SW01 2.SW02 3.SW03 4.SW04 5.TI01(10.85) 6.TI02(12.0)'

# ===== True Color =====
TRUE_COLOR = 0                                   # VNR 頻道選單的真實色編號
TC_BANDS = ('VN08', 'VN05', 'VN03')              # R, G, B
TC_FILL_BAND = 'VN11'                            # 飽和像元補值用 (868.5 nm 低增益，飽和約 1.1)
# JAXA SGLI 大氣校正 ATBD Table 3.1：標準氣壓 (1013.25 hPa) 下的 Rayleigh 光學厚度
TC_TAUR = {'VN08': 0.04265, 'VN05': 0.1132, 'VN03': 0.2361}
# 同 ATBD Table 3-1：臭氧光學厚度係數 Koz [DU-1]，τoz = Koz × DU
TC_KOZ = {'VN08': 4.42e-05, 'VN05': 6.79e-05, 'VN03': 3.74e-06}
TC_OZONE_DU = 300.                               # 臭氧氣候值 (L1B 無臭氧資料)
TC_MODES = {1: 'RGB', 2: 'CORR', 3: 'WV'}
TC_MODE_DESC = {1: 'RGB (raw composite)', 2: 'CORR (sunz + Rayleigh/O3 + CIRA stretch)',
                3: 'WV (NASA Worldview style)'}
TC_DESC = {'RGB': 'RGB', 'CORR': 'CORR', 'WV': 'WORLDVIEW'}
# true color stack 的圖層：R/G/B 反射率、太陽天頂角、各頻道的感測器天頂角與相對方位角
L_R, L_G, L_B, L_SZA, L_VZA, L_RAA = 0, 1, 2, 3, 4, 7


def make_truecolor(stack, mode):
    """由 true color stack (rows, cols, 10) 產生 0–1 RGB。
    RGB : TOA 反射率 rt (未除 cosθs) 線性 0–1
    CORR: sunz + Rayleigh / 臭氧校正 (JAXA 係數，CREFL 公式) + CIRA 增強
    WV  : 校正同 CORR，NASA Worldview 增強"""
    toa = stack[..., L_R:L_B + 1].astype(np.float64)
    sza = stack[..., L_SZA].astype(np.float64)
    night = ~(sza < 90.)

    if mode == 'RGB':
        rgb = np.clip(toa, 0., 1.)
    else:
        refl = toa * tc.sunz_corr(sza)[..., None]
        with np.errstate(invalid='ignore', divide='ignore'):
            mus = np.cos(np.deg2rad(sza))
            mus = np.where(mus > 0, mus, np.nan)
            corr = np.empty_like(refl)
            for k, band in enumerate(TC_BANDS):
                muv = np.cos(np.deg2rad(stack[..., L_VZA + k].astype(np.float64)))
                air_mass = np.minimum(1. / mus + 1. / muv, tc.MAXAIRMASS)
                t_o3 = tc.o3_transmittance(air_mass, TC_OZONE_DU, TC_KOZ[band])
                corr[..., k] = tc.crefl_core(refl[..., k], mus, muv, stack[..., L_RAA + k], TC_TAUR[band], t_o3)
        corr = np.clip(corr, tc.REFLMIN, tc.REFLMAX)
        rgb = tc.cira_stretch(corr) if mode == 'CORR' else tc.worldview_stretch(corr)

    rgb = np.array(rgb, dtype=np.float64)
    rgb[night] = np.nan
    return rgb


# ===== 單一檔案 =====
class SgliFile:

    table_dirs = ['.']      # 轉換表搜尋目錄 (由讀取器設定)

    def __init__(self, fullname, rescode=None, channel=None, time=None):
        self.fullname = fullname
        self.channel = channel
        self.satellite = 'GCOM-C'
        self.kind = self.fullname[31:34]          # VNR / IRS
        rescode = rescode or self.fullname[35:36]

        # 'GC1SG1_202310120107P04111_1BSG_VNRDQ_2008.h5'
        self.res = ('unknown', ['', ''])
        if self.kind == 'VNR':
            if rescode == 'K':
                self.res = 'VIS/NIR:1000m', ['1000M', '1000M']
            elif rescode == 'Q':
                self.res = 'VIS/NIR:250m', ['250M', '250M']
        elif self.kind == 'IRS':
            if rescode == 'K':
                self.res = 'SWI:1000m/TIR:1000m', ['1000M', '1000M']
            elif rescode == 'H':
                self.res = 'SWI:1000m/TIR:500m', ['1000M', '500M']
            elif rescode == 'Y':
                self.res = 'SWI:1000m/TIR:250m', ['1000M', '250M']
            elif rescode == 'X':
                self.res = 'SWI:250m/TIR:1000m', ['250M', '1000M']
            elif rescode == 'M':
                self.res = 'SWI:250m/TIR:500m', ['250M', '500M']
            elif rescode == 'Q':
                self.res = 'SWI:250m/TIR:250m', ['250M', '250M']
            elif rescode == 'L':
                self.res = 'SWI:1000m/TIR:1000m (re-sampling)', ['1000M', '1000M']
        if time:
            self.time = datetime.strptime(time, '%Y%m%d%H%M')
        else:
            self.time = datetime.strptime(self.fullname[7:19], '%Y%m%d%H%M')

    def __lt__(self, other):
        return self.time < other.time

    def planckInvert(self, wavel, Llambda):
        """input wavelength in microns and Llambda in W/m^2/micron/sr, output
        output brightness temperature in K  (note that we've remove the factor
        of pi because we are working with radiances, not fluxes)
        """
        c = 2.99792458e+08  # m/s -- speed of light in vacumn
        h = 6.62606876e-34  # J s  -- Planck's constant
        kb = 1.3806503e-23  # J/K  -- Boltzman's constant
        c1 = 2. * h * c ** 2.
        c2 = h * c / kb

        Llambda = Llambda * 1.e6  # convert to W/m^2/m/sr
        wavel = wavel * 1.e-6  # convert wavelength to m
        Tbright = c2 / (wavel * np.log(c1 / (wavel ** 5. * Llambda) + 1.))

        return Tbright

    def conversion_table(self):
        for folder in self.table_dirs:
            path = os.path.join(folder, CONVERSION_TABLE)
            if os.path.isfile(path):
                return path
        raise FileNotFoundError(tr('找不到 SGLI 轉換表 %s (已找過：%s)', 'SGLI conversion table %s not found in: %s')
                                % (CONVERSION_TABLE, ', '.join(self.table_dirs)))

    @staticmethod
    def _attr(attrs, key, default=None):
        return np.atleast_1d(attrs[key])[0] if key in attrs else default

    def extract(self, channel):
        """VN / SW：反射率 (除以 cos 太陽天頂角)；TI：亮溫 (°C)。
        DN = 原始值 & Mask；DN 16383 (或原始值 = Error_DN) 為缺值 -> NaN；DN 16382 為飽和 (保留飽和值)。"""
        self._h5 = h5py.File(self.fullname, 'r')
        lt = self._h5[f'Image_data/Lt_{channel}']
        attrs = lt.attrs
        raw = lt[:]
        dn = np.bitwise_and(raw, int(self._attr(attrs, 'Mask', 16383)))
        missing = (raw == self._attr(attrs, 'Error_DN', 65535)) | (dn == MISSING_DN)
        saturated = (dn == SATURATION_DN) & ~missing
        self.shape = raw.shape
        self.spatial_res = self._attr(attrs, 'Spatial_resolution')
        print(tr('%s %s：缺值 %.1f%%，飽和 %.1f%%', '%s %s: missing %.1f%%, saturated %.1f%%')
              % (channel, raw.shape, 100. * missing.mean(), 100. * saturated.mean()))

        if not channel.startswith('TI'):
            # VN01–VN11、SW01–SW04：反射率 rt = Lt·π/(F0/d²) = DN·Slope_reflectance + Offset_reflectance
            slope = float(self._attr(attrs, 'Slope_reflectance'))
            offset = float(self._attr(attrs, 'Offset_reflectance'))
            data = dn.astype(np.float32) * slope + offset
            data[missing] = np.nan

            # 太陽天頂角 (tie point 每 Resampling_interval 像元一點) 線性內插到全解析度
            print(tr('內插太陽天頂角 .........', 'Interpolating SZ .........'))
            sz = self.geometry_field('Solar_zenith', data.shape)

            cos_theta_0 = np.cos(np.deg2rad(sz))
            data = data / cos_theta_0

            print(np.nanmax(data), np.nanmin(data), np.nanmedian(data))
            return data

        else:
            # TI01–TI02：輻射 Lt = DN·Slope + Offset，再以轉換表換成亮溫
            slope = float(self._attr(attrs, 'Slope'))
            offset = float(self._attr(attrs, 'Offset'))
            data = dn.astype(np.float64) * slope + offset
            data[missing] = np.nan

            from scipy.interpolate import CubicSpline
            import pandas as pd
            table = pd.read_csv(self.conversion_table(), header=None)
            R = table.values[:, 0]
            TB = table.values[:, 1]
            cs = CubicSpline(R, TB)

            bt = cs(data) - 273.15
            print(bt.shape)

            return bt

    def bilin_2d(self, data, interval, lon_mode=False):
        # (目前未使用，保留)
        data = data.copy()

        if interval == 1:
            return data

        if lon_mode is True:
            max_diff = np.nanmax(np.abs(data[:, :-1] - data[:, 1:]))
            if max_diff > 180.:
                data[data < 0] = 360. + data[data < 0]

        data = np.concatenate((data, data[-1].reshape(1, -1)), axis=0)
        data = np.concatenate((data, data[:, -1].reshape(-1, 1)), axis=1)

        ratio_horizontal = np.tile(np.linspace(0, (interval - 1) / interval, interval, dtype=np.float32),
                                   (data.shape[0] * interval, data.shape[1] - 1))
        ratio_vertical = np.tile(np.linspace(0, (interval - 1) / interval, interval, dtype=np.float32).reshape(-1, 1),
                                 (data.shape[0] - 1, (data.shape[1] - 1) * interval))
        repeat_data = np.repeat(data, interval, axis=0)
        repeat_data = np.repeat(repeat_data, interval, axis=1)

        horizontal_interp = (1. - ratio_horizontal) * repeat_data[:, :-interval] + ratio_horizontal * repeat_data[:, interval:]
        ret = (1. - ratio_vertical) * horizontal_interp[:-interval, :] + ratio_vertical * horizontal_interp[interval:, :]

        if lon_mode is True:
            ret[ret > 180.] = ret[ret > 180.] - 360.

        return ret

    def _resampling_interval(self, ds):
        return int(np.atleast_1d(ds.attrs.get('Resampling_interval', [10]))[0])

    def _pixel_coords(self, h, shape):
        """影像像元中心在 tie 參考網格 (Image_data 的 Number_of_lines × Number_of_pixels) 上的座標。
        同一檔案內有不同解析度的頻道 (例如 250 m 產品中的 SW01/SW02/SW04 為 1 km)，不能直接用像元編號。"""
        im = h['Image_data'].attrs
        n_lines = int(np.atleast_1d(im.get('Number_of_lines', [shape[0]]))[0])
        n_pixels = int(np.atleast_1d(im.get('Number_of_pixels', [shape[1]]))[0])
        return tiepoint.pixel_coords(shape[0], n_lines), tiepoint.pixel_coords(shape[1], n_pixels)

    def _geometry(self, h, name, rows, cols, angle=False):
        """Geometry_data 的欄位內插到 (rows, cols) 座標 (tie 參考網格)；Error_DN 視為無效。
        angle=True：方位角以 sin / cos 分量內插，避免 0/360 跳躍。"""
        ds = h['Geometry_data/%s' % name]
        raw = ds[:]
        slope = float(np.atleast_1d(ds.attrs.get('Slope', [1.]))[0])
        offset = float(np.atleast_1d(ds.attrs.get('Offset', [0.]))[0])
        err = ds.attrs.get('Error_DN')
        interval = self._resampling_interval(ds)
        val = raw.astype(np.float64) * slope + offset
        if err is not None:
            val[raw == np.atleast_1d(err)[0]] = np.nan
        shape = (len(rows), len(cols))
        if not angle:
            return tiepoint.interp_field(val, interval, shape, rows, cols)
        rad = np.deg2rad(val)
        s = tiepoint.interp_field(np.sin(rad), interval, shape, rows, cols)
        c = tiepoint.interp_field(np.cos(rad), interval, shape, rows, cols)
        return np.rad2deg(np.arctan2(s, c)).astype(np.float32)

    def geometry_field(self, name, shape):
        """Geometry_data 中的純量場 (例如 Solar_zenith) 線性內插到 shape；Error_DN 視為無效。"""
        with h5py.File(self.fullname, 'r') as h:
            rows, cols = self._pixel_coords(h, shape)
            return self._geometry(h, name, rows, cols)

    def band_shape(self, band):
        with h5py.File(self.fullname, 'r') as h:
            return h['Image_data/Lt_%s' % band].shape

    def _band_reflectance(self, h, band, rows=slice(None), cols=slice(None)):
        """Lt_<band> 的 TOA 反射率 rt (未除 cosθs) 與飽和遮罩；可只讀取 [rows, cols] 範圍。"""
        ds = h['Image_data/Lt_%s' % band]
        attrs = ds.attrs
        raw = ds[rows, cols]
        dn = np.bitwise_and(raw, int(self._attr(attrs, 'Mask', 16383)))
        missing = (raw == self._attr(attrs, 'Error_DN', 65535)) | (dn == MISSING_DN)
        refl = dn.astype(np.float32) * float(self._attr(attrs, 'Slope_reflectance')) \
            + float(self._attr(attrs, 'Offset_reflectance'))
        refl[missing] = np.nan
        return refl, (dn == SATURATION_DN) & ~missing

    def truecolor_stack(self, rows=slice(None), cols=slice(None)):
        """真實色所需資料 (只讀取 [rows, cols] 範圍)，shape (ny, nx, 10) float32：
        R/G/B 的 TOA 反射率 rt (飽和像元以 VN11 補值)、太陽天頂角、
        R/G/B 各自的感測器天頂角、R/G/B 各自的相對方位角 (太陽方位角 − 感測器方位角)。"""
        with h5py.File(self.fullname, 'r') as h:
            fill = self._band_reflectance(h, TC_FILL_BAND, rows, cols)[0]
            bands = []
            for band in TC_BANDS:
                refl, sat = self._band_reflectance(h, band, rows, cols)
                if sat.any():
                    # 飽和值是下限：以 VN11 補值，但不低於飽和值
                    refl = np.where(sat, np.fmax(fill, refl), refl)
                print(tr('%s 飽和 %.1f%% (以 %s 補值)', '%s saturated %.1f%% (filled with %s)')
                      % (band, 100. * sat.mean(), TC_FILL_BAND))
                bands.append(refl)

            full = h['Image_data/Lt_%s' % TC_BANDS[0]].shape
            r_all, c_all = self._pixel_coords(h, full)
            r, c = r_all[rows], c_all[cols]
            geom = h['Geometry_data']
            sza = self._geometry(h, 'Solar_zenith', r, c)
            saa = self._geometry(h, 'Solar_azimuth', r, c, angle=True)
            vza, raa = [], []
            for band in TC_BANDS:
                zen = 'Sensor_zenith_%s' % band if 'Sensor_zenith_%s' % band in geom else 'Sensor_zenith'
                azi = 'Sensor_azimuth_%s' % band if 'Sensor_azimuth_%s' % band in geom else 'Sensor_azimuth'
                vza.append(self._geometry(h, zen, r, c))
                raa.append(saa - self._geometry(h, azi, r, c, angle=True))
        return np.dstack(bands + [sza] + vza + raa).astype(np.float32)

    def geocoord(self):
        """經緯度 tie point (每 Resampling_interval 像元一點) 內插到全解析度；須先呼叫 extract (取得影像大小)。
        於地心直角座標做二次樣條內插 (同 satpy / geotiepoints)，比原本的 griddata 準確且快約 10 倍，
        換日線與極區也不會出錯 (見 satima/tiepoint.py 的實測)。原作法保留在 geocoord_griddata。"""
        shape = self.shape
        with h5py.File(self.fullname, 'r') as geo:
            lat_ds, lon_ds = geo['Geometry_data/Latitude'], geo['Geometry_data/Longitude']
            lats, lons = lat_ds[:], lon_ds[:]
            interval = self._resampling_interval(lat_ds)
            rows, cols = self._pixel_coords(geo, shape)
        print(tr('內插經緯度 %s -> %s (間隔 %d) .........', 'Interpolating lat/lon %s -> %s (interval %d) .........')
              % (lats.shape, shape, interval))
        return tiepoint.interp_geo(lats, lons, interval, shape, rows, cols)

    def geocoord_griddata(self):
        """(原作法，保留對照) 經緯度格點 (每 10 像元) 以 griddata 線性內插到全解析度。"""
        with h5py.File(self.fullname, 'r') as geo:
            lats, lons = geo['Geometry_data/Latitude'][:], geo['Geometry_data/Longitude'][:]

        line, pixel = lats.shape
        print(line, pixel)

        pixel = np.arange(0, (pixel) * 10, 10)
        line = np.arange(0, (line) * 10, 10)
        pixel_new, line_new = np.meshgrid(pixel, line)
        points = np.array((pixel_new.flatten(), line_new.flatten())).T

        line, pixel = self.shape

        xnew = np.arange(0, pixel)
        ynew = np.arange(0, line)
        xnew, ynew = np.meshgrid(xnew, ynew)
        from scipy.interpolate import griddata
        print(tr('內插中 .........', 'Interpolating .........'))
        method = 'linear'
        values = lats.flatten()
        print(lats.shape)
        lats = griddata(points, values, (xnew, ynew), method=method)
        values = lons.flatten()
        lons = griddata(points, values, (xnew, ynew), method=method)

        return lats, lons

    def plot(self, canvas, lats, lons, data, variant):
        if data.size == 0:
            return ''
        if self.kind == 'VNR' and self.channel == TRUE_COLOR:
            print(tr('製作真實色：%s', 'Making true color: %s') % variant)
            canvas.rgb(lons, lats, make_truecolor(data, variant))
            return variant
        if variant is None:          # VN01–VN11、SW01–SW04：反射率
            print(np.nanmax(data), np.nanmin(data), np.nanmedian(data))
            canvas.mesh(lons, lats, data, cmap='gray', vmin=0.0, vmax=1.0)
            return ''
        return canvas.bt_mesh(lons, lats, data, variant)

    @staticmethod
    def _obs_hours(ds):
        """Geometry_data/Obs_time* (當日 UTC 小時，tie 網格)；缺值為 NaN。"""
        raw = ds[:]
        val = raw.astype(np.float64) * float(np.atleast_1d(ds.attrs.get('Slope', [0.001]))[0]) \
            + float(np.atleast_1d(ds.attrs.get('Offset', [0.]))[0])
        val[raw == -32768] = np.nan
        return val

    def row_times(self, nrows, bands=()):
        """每一列的觀測時間 (UTC)：Image_data/Line_tai93 (逐參考列，TAI93 秒數)，依列中心換到 nrows 列
        (1 km 頻道的一列對應 4 條 250 m 參考列)。
        VN 各頻道在焦平面上的位置不同，同一地點的觀測時間相差約 −5～+5 s：加上 bands 各頻道
        Obs_time_<band> 與 Obs_time 的平均差 (多個頻道取平均，例如真實色的 R/G/B)。"""
        with h5py.File(self.fullname, 'r') as h:
            if 'Line_tai93' not in h['Image_data']:
                ga = h['Global_attributes'].attrs
                start, end = (datetime.strptime(ga[k][0].decode('utf-8')[:-4], '%Y%m%d %H:%M:%S')
                              for k in ('Scene_start_time', 'Scene_end_time'))
                return obstime.linear(start, end, nrows)
            line = obstime.tai93(np.ma.masked_less_equal(h['Image_data/Line_tai93'][:], 0))
            geo = h['Geometry_data']
            offsets = []
            if 'Obs_time' in geo:
                ref = self._obs_hours(geo['Obs_time'])
                for band in bands:
                    name = 'Obs_time_%s' % band
                    if name in geo and geo[name].shape == geo['Obs_time'].shape:
                        diff = (self._obs_hours(geo[name]) - ref + 12.) % 24. - 12.     # 跨午夜
                        offsets.append(np.nanmean(diff) * 3600.)
        t = obstime.resample_rows(line, nrows)
        if offsets and np.isfinite(np.mean(offsets)):
            t = t + np.timedelta64(int(round(np.mean(offsets) * 1e6)), 'us')
        return t


# ===== 讀取器 =====
class SgliReader(SatReader):
    NAME = 'SGLI-HDF5'
    SENSOR = 'SGLI'
    DESCRIPTION = 'GCOM-C SGLI L1B (HDF5)'
    SOURCE = 'JAXA'
    PATTERNS = ('GC1SG1*1B*h5',)
    MULTI_FILE = True

    def make_file(self, name):
        # GC1SG1_202309070242D08403_1BSG_IRSNQ_2008.h5
        return SgliFile(name, channel=None, rescode=name[35:36], time=name[7:19])

    def listing(self, f):
        return (f.satellite, 'SGLI', f.time, f.res[0])

    def select_band(self):
        kinds = set(f.kind for f in self.files)
        if len(kinds) > 1:
            print(tr('注意：選到的檔案混有 VNR 與 IRS，只使用 %s 的檔案。',
                     'WARNING: VNR and IRS files mixed; only %s files are used.') % self.files[0].kind)
            self.files = [f for f in self.files if f.kind == self.files[0].kind]
        self.multifilesflag = len(self.files) > 1
        SgliFile.table_dirs = [os.getcwd(), self.cfg._home]
        f0 = self.files[0]
        print(f0.fullname)
        if f0.kind == 'VNR':
            print('%d.TRUE COLOR (VN08/VN05/VN03) ' % TRUE_COLOR + ' '.join('%d.%s' % (k, v) for k, v in VNR_CHANNELS.items()))
            channel = ui.ask_int(tr('輸入頻道編號：', 'Channel: '), [TRUE_COLOR] + list(VNR_CHANNELS))
            self.channel_name = 'TrueColor' if channel == TRUE_COLOR else VNR_CHANNELS[channel]
        else:
            print(IRS_MENU)
            channel = ui.ask_int(tr('輸入頻道編號：', 'Channel: '), IRS_CHANNELS)
            self.channel_name = IRS_CHANNELS[channel]
        self.channel = channel
        for f in self.files:
            f.channel = channel
        self.is_tir = f0.kind == 'IRS' and channel >= 5
        self.is_tc = f0.kind == 'VNR' and channel == TRUE_COLOR
        self.res = f0.res[1][1] if self.is_tir else f0.res[1][0]     # 讀檔後依頻道實際解析度更新

    def select_variants(self):
        if self.is_tc:
            return [TC_MODES[k] for k in ui.ask_menu('', TC_MODE_DESC, tr('選擇真實色模式 (可複選，0 = 全部)：',
                                                                          'Select true color mode: '))]
        if self.is_tir:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare_tc(self):
        """真實色：先由經緯度決定裁切範圍，再只讀取範圍內的 VN08/VN05/VN03/VN11 與幾何角度。"""
        lats, lons, nrows, rowtime = [], [], [], []
        for f in self.files:
            f.shape = f.band_shape(TC_BANDS[0])
            la, lo = f.geocoord()
            lats.append(la)
            lons.append(lo)
            nrows.append(f.shape[0])
            rowtime.append(f.row_times(f.shape[0], TC_BANDS))
        lats = np.concatenate(lats, axis=0)
        lons = np.concatenate(lons, axis=0)
        rowtime = np.concatenate(rowtime)
        if self.multifilesflag:
            print(tr('已合併 %d 個片段：', 'Merged %d granules:') % len(self.files), lats.shape)
        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        yi, yj, xi, xj = self.box
        stacks, r0 = [], 0
        for f, n in zip(self.files, nrows):
            a, b = max(yi, r0), min(yj, r0 + n)          # 此檔案在裁切範圍內的列
            if a < b:
                stacks.append(f.truecolor_stack(slice(a - r0, b - r0), slice(xi, xj)))
            r0 += n
        self.data = np.concatenate(stacks, axis=0)
        self.time = self.center_time(rowtime)
        return True

    def prepare(self):
        if self.is_tc:
            return self.prepare_tc()
        datas, lats, lons, rowtime = [], [], [], []
        for f in self.files:
            data = f.extract(self.channel_name)
            la, lo = f.geocoord()
            datas.append(data)
            lats.append(la)
            lons.append(lo)
            rowtime.append(f.row_times(data.shape[0], [self.channel_name]))
        data = self.concat(datas)
        lats = np.concatenate(lats, axis=0)
        lons = np.concatenate(lons, axis=0)
        rowtime = np.concatenate(rowtime)
        if self.multifilesflag:
            print(tr('已合併 %d 個片段：', 'Merged %d granules:') % len(self.files), lats.shape)
        res = getattr(self.files[0], 'spatial_res', None)
        if res:
            self.res = '%dM' % int(float(res))      # 該頻道的實際解析度 (例如 250 m 產品中的 SW01 為 1000M)
        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        self.data = self.crop(data)
        self.time = self.center_time(rowtime)
        if self.is_tir:
            self.report_eyetemp(self.data)
        return True

    def slice(self, georange, lats, lons):
        return self.slice_swath(lats, lons)

    def plot(self, canvas, variant):
        cm = self.files[0].plot(canvas, self.lats, self.lons, self.data, variant)
        if self.is_tc:
            return dict(band='TrueColor', colorbar=cm, label='TRUE COLOR %s' % TC_DESC[variant])
        return dict(band=self.channel_name, colorbar=cm, label=self.channel_name)


READER = SgliReader
