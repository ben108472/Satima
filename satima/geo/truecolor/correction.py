'''可見光反射率校正：太陽天頂角校正、Rayleigh 散射校正'''
import logging

import numpy as np

from ...core.i18n import tr
from . import astronomy


def sunzen_corr_cos(data, cos_zen, limit=88., max_sza=95.):
    '''太陽天頂角校正 (反射率除以 cos(SZA))，與 satpy 的 sunz_corrected 相同。

    SZA 超過 limit 後校正量不再增加，並在 limit~max_sza 之間逐漸降為 0，
    避免晨昏線附近過度放大。max_sza 為 None 時則固定使用 limit 的校正量。
    '''
    limit_rad = np.deg2rad(limit)
    limit_cos = np.cos(limit_rad)
    with np.errstate(invalid='ignore', divide='ignore'):
        corr = 1. / cos_zen
        if max_sza is not None:
            grad = (np.arccos(cos_zen) - limit_rad) / (np.deg2rad(max_sza) - limit_rad)
            grad = np.clip(1. - np.log(grad + 1.) / np.log(2.), 0., None)
        else:
            grad = 1.
        corr = np.where(cos_zen > limit_cos, corr, grad / limit_cos)
    corr = np.where(np.isnan(cos_zen), 0., corr)
    return data * corr


def sunzen_reduction(data, sunz, limit=80., max_sza=90., strength=1.3):
    '''晨昏線淡出 (與 satpy 的 sunz_reduced 相同)：SZA 由 limit 到 max_sza 讓訊號平滑降到 0。

    接近晨昏線時 1/cos(SZA) 會把偏藍的曙暮光放大數十倍，而 Rayleigh 查表在 SZA 87.7° 以上已飽和，
    無法完全扣除，會留下一條藍紫色帶；淡出後晨昏線會自然地轉暗。
    '''
    factor = np.clip((sunz - limit) / (max_sza - limit), 0., 1.)
    with np.errstate(invalid='ignore', divide='ignore'):
        factor = 1. - np.log2(factor + 1.)
        factor = factor ** strength / (factor ** strength + (1. - factor) ** strength)
    factor = np.where(sunz < limit, 1., factor)
    factor = np.where(np.isnan(sunz), 0., factor)
    return data * factor


def sun_zenith_correction(refl, utc_time, lon, lat, limit=88., max_sza=95.):
    '''依觀測時間與經緯度計算太陽天頂角後做校正'''
    cos_zen = astronomy.cos_sun_zenith(utc_time, lon, lat)
    return sunzen_corr_cos(refl, cos_zen, limit, max_sza)


def cloud_height_factor(bt_c, warm=283., cold=233., max_reduction=0.7):
    '''CIRA 的雲頂 Rayleigh 削減係數 (Miller et al. 2016; 2020, JTECH)。

    以 10.3 µm 亮溫 (°C) 代表雲頂高度：高雲上方的大氣路徑較短，
    若照晴空校正會過度扣除，使雲偏紅。
    BT >= 283 K 不調整；BT <= 233 K 削減 70%；中間線性內插。
    '''
    bt_k = np.asarray(bt_c, dtype=np.float32) + 273.15
    reduction = max_reduction * np.clip((warm - bt_k) / (warm - cold), 0., 1.)
    return np.where(np.isfinite(reduction), 1. - reduction, 1.).astype(np.float32)


def rayleigh_optical_depth(wavelength):
    '''Rayleigh 光學厚度 (Hansen & Travis 1974，海平面標準大氣)，波長單位 µm'''
    w = float(wavelength)
    return 0.008569 * w ** -4 * (1. + 0.0113 * w ** -2 + 0.00013 * w ** -4)


def rayleigh_transmittance(mu, tau):
    '''Rayleigh 大氣的總穿透率 (直射 + 漫射，二流近似；與 CREFL / 6S 相同)。mu 為天頂角餘弦'''
    return ((2. / 3. + mu) + (2. / 3. - mu) * np.exp(-tau / mu)) / (4. / 3. + tau)


def rayleigh_spherical_albedo(tau):
    '''Rayleigh 大氣的球面反照率 (與 CREFL 相同)；tau 可為純量 (回傳 float) 或陣列'''
    from scipy.special import exp1
    tau = np.asarray(tau, dtype=np.float64)
    e3 = (np.exp(-tau) * (1. - tau) + tau * tau * exp1(tau)) / 2.     # 三階指數積分 E3(tau)
    out = (3. * tau - e3 * (4. + 2. * tau) + 2. * np.exp(-tau)) / (4. + 3. * tau)
    return float(out) if out.ndim == 0 else out


# 臭氧 Chappuis 帶的吸收係數 (每 atm-cm；1 atm-cm = 1000 DU)：[波長 µm, 係數]
# Leckner (1978) 的表列值 (Bird & Riordan 1986 SPCTRAL2 的 Table 2-1)；已與 pvlib-python 的
# pvlib/spectrum/spectrl2.py (_SPECTRL2_COEFFS['ozone_absorption']) 逐點核對，0.40~1.0 µm 完全相同
# (原表中係數為 0 的波長只保留兩端)。以 MODIS / VIIRS 各波段的通帶平均後，
# 與 CREFL 1.7.1 的臭氧係數 (satpy 的 _crefl_utils) 比對：0.46~0.68 µm 的 12 個波段相差 0.01 /atm-cm 以內。
OZONE_TABLE = np.array([
    [0.40, 0.000], [0.44, 0.000], [0.45, 0.003], [0.46, 0.006], [0.47, 0.009], [0.48, 0.014], [0.49, 0.021],
    [0.50, 0.030], [0.51, 0.040], [0.52, 0.048], [0.53, 0.063], [0.54, 0.075], [0.55, 0.085], [0.57, 0.120],
    [0.593, 0.119], [0.61, 0.120], [0.63, 0.090], [0.656, 0.065], [0.6676, 0.051], [0.69, 0.028], [0.71, 0.018],
    [0.718, 0.015], [0.7244, 0.012], [0.74, 0.010], [0.7525, 0.008], [0.7575, 0.007], [0.7625, 0.006],
    [0.7675, 0.005], [0.78, 0.000], [1.00, 0.000]])
OZONE_HEIGHT = 22. / 6370.   # 臭氧層高度 22 km / 地球半徑 (SPCTRAL2 Eq. 2-10)


def ozone_absorption(wavelength):
    '''臭氧吸收係數 (每 atm-cm)，波長單位 µm'''
    return np.interp(wavelength, OZONE_TABLE[:, 0], OZONE_TABLE[:, 1])


def ozone_airmass(zenith):
    '''臭氧層的斜程倍數 (球面幾何，SPCTRAL2 Eq. 2-10)：天頂角 0° 為 1，90° 約 12 (平面平行的 1/cos 會趨近無限大)'''
    mu = np.cos(np.deg2rad(np.clip(zenith, 0., 90.)))
    return (1. + OZONE_HEIGHT) / np.sqrt(mu * mu + 2. * OZONE_HEIGHT)


LUT_MAX_WAVELENGTH = 0.79    # pyspectral 的 Rayleigh 查表只涵蓋 0.4~0.8 µm


class RayleighCorrector:
    '''以 pyspectral 查表 (rayleigh_only) 計算 Rayleigh 路徑反射率並校正。

    mode
      'full'      地表反射率 = y / (1 + S·y)，y = (TOA - 路徑反射) / (T↓·T↑)
                  T 為 Rayleigh 上下行穿透率、S 為大氣球面反照率 (SeaDAS / CREFL 的作法，CIRA 改自 SeaDAS)：
                  暗的地表與海面除以穿透率後較亮、較藍；亮雲的淨校正量幾乎為 0，不會偏黃
      'subtract'  只扣除路徑反射 (舊版作法)
    '''

    def __init__(self, platform, sensor, atmosphere='us-standard', mode='full'):
        from pyspectral.rayleigh import Rayleigh
        # 以波長呼叫時 pyspectral 每次都會印出警告，這裡將其關閉
        logging.getLogger('pyspectral').setLevel(logging.ERROR)
        if mode not in ('full', 'subtract'):
            raise ValueError(tr('Rayleigh 校正方式須為 full 或 subtract：%r',
                                'the Rayleigh correction mode must be full or subtract: %r') % mode)
        self.platform = platform
        self.sensor = sensor
        self.mode = mode
        self._rayleigh = Rayleigh(platform, sensor, atmosphere=atmosphere,
                                  aerosol_type='rayleigh_only')
        self._wavelengths = {}
        self._rsr = None

    def _response(self, band_name, nominal):
        '''找出波段的光譜響應：先以名稱 (pyspectral 的名稱多為 ch1、ch2…，與讀取器不同)，
        找不到時取中心波長最接近標稱波長 (相差 0.05 µm 以內) 的波段。'''
        if self._rsr is None:
            from pyspectral.rsr_reader import RelativeSpectralResponse
            self._rsr = RelativeSpectralResponse(self.platform, self.sensor).rsr
        if band_name in self._rsr:
            return self._rsr[band_name]['det-1']
        best = min(self._rsr, key=lambda b: abs(self._rsr[b]['det-1'].get('central_wavelength', 99.) - nominal))
        det = self._rsr[best]['det-1']
        if abs(det.get('central_wavelength', 99.) - nominal) > 0.05:
            raise KeyError(band_name)
        return det

    def effective_wavelength(self, band_name, nominal):
        '''以光譜響應函數 (RSR) 計算 Rayleigh 有效波長 (λ^-4 加權，與 pyspectral 相同)。
        找不到 RSR 時改用標稱中心波長。'''
        if band_name not in self._wavelengths:
            try:
                from pyspectral.utils import get_central_wave
                rsr = self._response(band_name, float(nominal))
                wvl = rsr['wavelength']
                self._wavelengths[band_name] = float(get_central_wave(wvl, rsr['response'], weight=1. / wvl ** 4))
            except Exception:
                print(tr('  注意：找不到 {} {} {} 的光譜響應，Rayleigh 校正改用標稱波長 {} µm',
                         '  NOTE: spectral response of {} {} {} not found; the Rayleigh correction uses the nominal '
                         'wavelength {} µm instead').format(
                    self.platform, self.sensor, band_name, nominal))
                self._wavelengths[band_name] = float(nominal)
        return self._wavelengths[band_name]

    def ozone_coefficient(self, band_name, nominal):
        '''波段的臭氧吸收係數 (每 atm-cm)：OZONE_TABLE 以光譜響應加權平均；找不到 RSR 時取標稱波長處的值。'''
        try:
            rsr = self._response(band_name, float(nominal))
            wvl, resp = np.asarray(rsr['wavelength'], float), np.asarray(rsr['response'], float)
            trapz = getattr(np, 'trapezoid', None) or np.trapz        # numpy 2 改名為 trapezoid
            return float(trapz(ozone_absorption(wvl) * resp, wvl) / trapz(resp, wvl))
        except Exception:
            return float(ozone_absorption(float(nominal)))

    def path_reflectance(self, sunz, satz, ssadiff, wvl, redband=None):
        '''Rayleigh 路徑反射率 (0~1)。查表範圍外 (近紅外) 光學很薄，依光學厚度比例由 0.78 µm 換算。'''
        red_pct = None if redband is None else redband * 100.
        if wvl > LUT_MAX_WAVELENGTH:
            ref = self._rayleigh.get_reflectance(sunz, satz, ssadiff, 0.78, red_pct) / 100.
            return ref * (rayleigh_optical_depth(wvl) / rayleigh_optical_depth(0.78))
        return self._rayleigh.get_reflectance(sunz, satz, ssadiff, wvl, red_pct) / 100.

    def terms(self, sunz, satz, ssadiff, band_name, nominal_wavelength, redband=None, factor=None, scale_tau=False):
        '''校正用的各項 (路徑反射, 上下行穿透率乘積, 球面反照率)；後兩項在 subtract 模式為 None。
        同一波長、同一幾何的像素 (例如銳化用的高解析紅光與低解析紅光) 可共用。
        factor：路徑反射的削減係數 (cloud_height_factor 的結果)
        scale_tau：穿透率與球面反照率的光學厚度也乘上 factor (雲頂上方只剩這個比例的大氣)。只削減路徑反射、
          穿透率仍用整層大氣時，高雲的藍光被除以過小的穿透率而偏藍，太陽越低越明顯。'''
        wvl = self.effective_wavelength(band_name, nominal_wavelength)
        ray = self.path_reflectance(sunz, satz, ssadiff, wvl, redband)
        if factor is not None:
            ray = ray * factor
        if self.mode != 'full':
            return ray, None, None
        tau = rayleigh_optical_depth(wvl)
        if scale_tau and factor is not None:
            tau = tau * np.asarray(factor, dtype=np.float64)
        mus = np.clip(np.cos(np.deg2rad(sunz)), 0.02, 1.)
        muv = np.clip(np.cos(np.deg2rad(satz)), 0.02, 1.)
        return ray, rayleigh_transmittance(mus, tau) * rayleigh_transmittance(muv, tau), rayleigh_spherical_albedo(tau)

    @staticmethod
    def apply(refl, terms):
        '''以 terms() 的結果校正反射率。'''
        ray, trans, sph = terms
        out = refl - ray
        if trans is not None:
            out = out / trans
            out = out / (1. + sph * out)
        return out.astype(np.float32)

    def __call__(self, refl, sunz, satz, ssadiff, band_name, nominal_wavelength, redband=None, bt=None):
        '''回傳 Rayleigh 校正後的反射率。

        refl、redband：0~1 的反射率 (已做太陽天頂角校正)
        sunz、satz、ssadiff：太陽天頂角、衛星天頂角、相對方位角 (度)
        redband：給定時依 pyspectral 預設，紅光反射率 > 20% 處逐漸減少路徑反射
        bt：10.3 µm 亮溫 (°C)，給定時改用 CIRA 的雲頂高度削減路徑反射
        '''
        factor = None if bt is None else cloud_height_factor(bt)
        return self.apply(refl, self.terms(sunz, satz, ssadiff, band_name, nominal_wavelength, redband, factor))
