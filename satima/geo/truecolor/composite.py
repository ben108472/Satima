"""CIRA 風格真實色 (SHAC True Color：Synthetic/Hybrid Atmospherically Corrected)

依 Miller et al. (2016, BAMS, doi:10.1175/BAMS-D-15-00154.1) 與
Miller et al. (2020, JTECH, GeoColor, doi:10.1175/JTECH-D-19-0134.1) 的白天真實色作法：

  1. 太陽天頂角校正；臭氧 (Chappuis 帶) 吸收校正 (ozone：反射率除以臭氧穿透率，臭氧量固定)
  2. 藍、綠、紅、近紅外的 Rayleigh 散射校正 (pyspectral 路徑反射查表；CIRA 改自 SeaDAS，
     扣除路徑反射後再除以上下行穿透率並考慮大氣球面反照率，見 correction.RayleighCorrector)，
     並依 10.3 µm 亮溫削減雲頂上方的路徑反射 (283 K 不減，233 K 減 70%)；
     cloud_rayleigh_transmittance：穿透率與球面反照率的光學厚度也乘上同一個削減係數
     (這兩項是本套件加的，CIRA 論文沒有：未校正時深對流雲的 藍/紅 由太陽天頂角 < 50° 的 1.16
     升到 80° 的 1.42，日出日落時雲頂偏藍；校正後為 1.06~1.10，再由 band_gain 修正到 1.00~1.03)
  3. 綠光
       0.51 µm 綠光 (AHI / AMI / FCI)：未落在植被 0.55 µm 反射峰，混入部分 0.86 µm (hybrid green)。
         CIRA 論文為約 7%；與 CIRA 實際圖片比對 (GOES-19、Himawari-9、GK-2A、MTG-I1)，9% 時植被色調最接近 (預設)
       沒有綠光 (ABI)：CIRA 先以 AHI 建立的查表合成 0.51 µm 綠光，再同樣混入 0.86 µm。
         查表未公開；這裡照同樣作法，以 Himawari-9 (2026-10-21 03 UTC、11-09 05 UTC，北緯 0~50 度，
         Rayleigh 校正後的 B01/B02/B03) 擬合 G0.51 = 0.976 x B^0.668 x R^0.347 (AHI 自身誤差中位數 < 1 個色階)
     (0. 校正前先乘上各衛星的波段相對增益 band_gain：Himawari-9 的藍、綠光比其他衛星亮，雲會偏青綠，見下)
  4. 以高解析紅光比例銳化 (輸出細於低解析波段時)；銳化前先把低解析紅光模糊到與其他波段相同的 MTF
  5. 晨昏線淡出 (SZA 80°~90°，同 satpy sunz_reduced)；CIRA 原本以夜間圖層覆蓋這一段
  6. 色調曲線 (tone)
       "cira" (預設)：反射率截斷於 [0.025, 1.2]，取 log10 後正規化於 [-1.6, 0.176] (Miller et al. 2020 原文)。
         與 CIRA SLIDER 相同，但有兩個先天限制：反射率截在 1.2、白點卻在 1.5，輸出最亮只到 241 (厚雲約 212~230，偏灰)；
         低於 0.025 一律為 0，晴空海面的紅光 (約 0.02) 有六成像素變成 0，海面出現塊狀的藍黑色。
       "enhanced"：同一條 log 曲線，中間調幾乎不變，只改兩端 (見 enhanced_stretch)：
         白點降到 white_point (預設 1.25) 並在 shoulder 以上平滑壓縮到 1，暗部以平滑下限 (shadow_toe) 取代硬截斷。
  7. 飽和度 (saturation，預設 1 = 不調整)：拉伸後以亮度為中心放大色差，亮部 (雲) 逐漸不放大，雲不會染色 (見 vibrance)

與 CIRA SLIDER GeoColor 逐像素比對 (GOES-19、Himawari-9、GK-2A、MTG-I1，見 README)：
拉伸完全相同；雲與陸地色調一致，各波段色階差的絕對值中位數為 2~8。

波段相對增益 (band_gain，只用於真實色 / GeoColor)：
  深對流雲 (10.4 µm < -50°C、紅光反射率 > 0.8，太陽與衛星天頂角 < 40°) 在可見光幾乎是白色，各波段反射率應接近。
  Himawari-9 三個時次 (2023-07-01 03、2026-10-21 01、2026-11-09 05 UTC) 的 藍/紅 為 1.11~1.12、綠/紅 1.04~1.06，
  Himawari-8 (2020-09-04 03 UTC) 為 1.035、1.008，GOES-19 (2026-10-28 13 UTC) 藍/紅 1.006；
  三顆衛星的紅光 (0.86~0.90) 與近紅外 (0.90~0.93) 一致。檔頭的更新定標係數也無法消除這個差異。
  (以上為 TOA 反射率的比值。) 加入臭氧與高雲穿透率校正後，深對流雲校正後的 藍/紅、綠/紅 為
  Himawari-9 1.06、1.02，Himawari-8 0.98、0.98，GOES-19 0.95，GOES-16 (2019) 1.06；
  Himawari-9 的 B01 乘 0.94、B02 乘 0.98 後兩個比值都是 1.00 (兩項校正關閉時為 0.92、0.95)。

各波段以導航直接取樣到輸出格點 (satima/resample.py)，區塊平均對齊各解析度的格點，波段之間不會錯位；
之後在每個輸出像素計算太陽與衛星幾何並做校正。

recipe (讀取器的 TRUECOLOR 設定)：
  red / green / blue / nir   波段名稱 (green 可為 None，改用藍/紅/近紅外合成)
  ir                         10.3 µm 附近的紅外波段 (雲頂 Rayleigh 削減用，可省略)
  lo_blur                    (南北, 東西) 高斯 sigma，單位為低解析像素：低解析波段比高解析紅光 2x2 平均模糊的程度
  rsr                        {波段: pyspectral 波段名稱}，預設與波段名稱相同
"""
import threading
import time

import numpy as np

from .. import resample
from ...core.i18n import tr
from . import astronomy
from .correction import RayleighCorrector, cloud_height_factor, ozone_airmass, sunzen_corr_cos, sunzen_reduction

DEFAULTS = {
    'hybrid_green_fraction': 0.09,
    'abi_green': 'ahi',
    'rayleigh_correction': 'full',
    'cloud_rayleigh_reduction': True,
    'cloud_rayleigh_transmittance': True,   # 高雲的穿透率與球面反照率也依雲頂削減係數縮小光學厚度
    'ozone': 0.30,               # 臭氧吸收校正的臭氧量 (atm-cm；0 / false = 不校正)
    'sharpen': True,
    'sza_limit': 88.,
    'sza_max': 95.,
    'terminator_fade_start': 80.,
    'terminator_fade_end': 90.,
    'atmosphere': 'us-standard',
    'engine': 'numba',
    'band_gain': {'Himawari-9': {'B01': 0.94, 'B02': 0.98}},
    'tone': 'cira',              # 色調曲線：cira (論文原式) / enhanced (改良兩端)
    'white_point': 1.25,         # enhanced：對應全白的反射率 (CIRA 原式為 1.5)
    'shoulder': 0.95,            # enhanced：顯示值超過此值後平滑壓縮到 1 (1 = 不壓縮，超過白點直接截斷)
    'shadow_toe': 0.006,         # enhanced：暗部平滑下限的寬度 (反射率；0 = 同 CIRA 硬截斷於 0.025)
    'saturation': 1.0,           # 飽和度 (1 = 不調整)
}

REFL_RANGE = (0.025, 1.2)    # CIRA 反射率截斷範圍
LOG_RANGE = (-1.6, 0.176)    # CIRA log10 反射率正規化範圍
VIBRANCE_HIGH = 0.85         # 飽和度：顯示亮度高於此值不放大 (雲維持中性)
VIBRANCE_SPAN = 0.35         # 由此亮度往下 0.35 之內漸增到完整倍數
LUMA = (0.2126, 0.7152, 0.0722)
CHUNK = 2000000              # 每次處理的像素數 (控制記憶體用量)
AHI_GREEN = (0.976, 0.668, 0.347)   # 以 AHI 擬合的 0.51 µm 合成綠光：k x 藍^a x 紅^b (Rayleigh 校正後)
AHI_GREEN_FLOOR = 0.003             # 取次方前的反射率下限 (擬合時相同)
AHI_GREEN_MIN_BR = 0.3              # 藍色下限 = 0.3 x 紅色 (擬合資料的範圍，見 synthetic_green)


def synthetic_green(blue, red):
    """ABI 的 0.51 µm 合成綠光 (以 AHI 擬合，見檔頭說明)。

    藍色以 0.3 倍紅色為下限：擬合資料中 藍/紅 < 0.3 只佔約 0.4% (AHI 實際 綠/紅 約 0.40，此式在下限處為 0.42)；
    晨昏線附近藍光幾乎被 Rayleigh 校正扣光，不設下限時綠色會跟著趨近 0，只剩紅色 (CIRA 的查表不會如此)。"""
    k, a, b = AHI_GREEN
    red = np.maximum(red, AHI_GREEN_FLOOR)
    blue = np.maximum(np.maximum(blue, AHI_GREEN_FLOOR), AHI_GREEN_MIN_BR * red)
    return k * blue ** a * red ** b


def cira_stretch(refl):
    """CIRA 真實色的對數拉伸，模擬人眼對亮度的響應"""
    refl = np.clip(refl, *REFL_RANGE)
    out = (np.log10(refl) - LOG_RANGE[0]) / (LOG_RANGE[1] - LOG_RANGE[0])
    return np.clip(out, 0., 1.)


def enhanced_stretch(refl, white=1.25, knee=0.95, toe=0.006):
    """CIRA log 拉伸的改良版：黑點同為 0.025、中間調斜率相近，只改兩端。

      暗部：x2 = 0.025 + s(x - 0.025)，s(d) = (d + sqrt(d^2 + 4 toe^2)) / 2 為平滑的 max(d, 0)，
            低於 0.025 的反射率仍有層次 (不會整片變成 0)；toe = 0 時與 CIRA 相同為硬截斷
      亮部：v = log10(x2 / 0.025) / log10(white / 0.025)；v 超過 knee 的部分以指數平滑壓縮，漸近 1
    """
    black = REFL_RANGE[0]
    d = np.asarray(refl, dtype=np.float64) - black
    d = 0.5 * (d + np.sqrt(d * d + 4. * toe * toe)) if toe > 0 else np.maximum(d, 0.)
    v = np.log10((black + d) / black) / np.log10(white / black)
    if knee < 1.:
        span = 1. - knee
        v = np.where(v > knee, knee + span * (1. - np.exp(-(v - knee) / span)), v)
    return np.clip(v, 0., 1.)


def vibrance(rgb, saturation):
    """飽和度：以亮度 L 為中心放大色差，out = L + k (c - L)。
    k 在 L <= 0.50 時為 saturation，往上線性遞減，L >= 0.85 時為 1：陸地與海面變鮮明，雲 (本來就接近中性) 不被染色。"""
    lum = (LUMA[0] * rgb[..., 0] + LUMA[1] * rgb[..., 1] + LUMA[2] * rgb[..., 2])[..., None]
    k = 1. + (saturation - 1.) * np.clip((VIBRANCE_HIGH - lum) / VIBRANCE_SPAN, 0., 1.)
    return np.clip(lum + k * (rgb - lum), 0., 1.)


def sharpening_ratio(high, low):
    """比例銳化係數 (與 satpy RatioSharpenedRGB 相同)"""
    with np.errstate(divide='ignore', invalid='ignore'):
        ratio = high / low
    ratio[~np.isfinite(ratio) | (ratio < 0)] = 1.
    return np.clip(ratio, 0., 1.5)


class TrueColor:
    """真實色：建立時決定各波段的區塊平均倍數並讀取資料 (每個波段只讀一次)、準備 Rayleigh 查表；
    compute(g) 對整張格點或其中一條 (grid.strips) 計算。各像素的計算互相獨立，分條結果與整張完全相同。
    low_memory：讀取範圍以稀疏取樣估計 (見 resample.Source)，不需要整張的行列號。"""

    def __init__(self, reader, scene, recipe, grid, options=None, low_memory=False):
        self.t0 = time.time()
        opts = dict(DEFAULTS)
        opts.update(options or {})
        self.opts, self.recipe, self.scene = opts, recipe, scene
        ch = reader.CHANNELS
        ref_res = ch[recipe['blue']].res
        red_res = ch[recipe['red']].res
        k_ref = max(1, int(grid.pixel_km / ref_res))
        f_red = max(1, int(round(ref_res / red_res)))
        self.sharpen = bool(opts['sharpen']) and f_red > 1 and grid.pixel_km < 0.9 * ref_res
        print(tr('  輸出 %dx%d 像素，約 %.2f km/像素%s', '  Output %dx%d pixels, about %.2f km/pixel%s')
              % (grid.shape[1], grid.shape[0], grid.pixel_km,
                 tr('，使用紅光銳化', ', sharpened with the red band') if self.sharpen else ''))

        def factor(name):
            return max(1, int(round(k_ref * ref_res / ch[name].res)))

        self.factor = factor            # GeoColor 的紅外波段沿用同一套區塊倍數
        self.low_memory = low_memory

        # (角色, 波段, 區塊平均倍數, 高斯模糊)；紅光只讀一次：低解析 (對齊其他波段的區塊平均，並匹配 MTF) 與高解析 (銳化用)
        self.k_lo = k_ref * f_red
        self.blur = recipe.get('lo_blur') if (f_red > 1 and k_ref == 1) else None
        self.sources = {}
        for role in ('blue', 'green', 'nir'):
            if recipe.get(role):
                self.sources[role] = (resample.Source(reader, scene, recipe[role], grid, factor(recipe[role]),
                                                      coarse=low_memory), factor(recipe[role]))
        self.sources['red'] = (resample.Source(reader, scene, recipe['red'], grid, self.k_lo,
                                               resample.blur_margin(self.blur), coarse=low_memory), self.k_lo)
        # 波段相對增益 {角色: 倍數} (見檔頭說明)；紅光的增益同時用於低解析與高解析
        table = (opts.get('band_gain') or {}).get(scene.satellite) or {}
        self.gain = {role: float(table[recipe[role]]) for role in ('blue', 'green', 'red', 'nir')
                     if recipe.get(role) and recipe[role] in table and float(table[recipe[role]]) != 1.}
        if self.gain:
            print(tr('  %s 波段相對增益：%s', '  %s relative band gains: %s') % (scene.satellite, tr('、', ', ').join(
                '%s x %.3g' % (recipe[role], g) for role, g in self.gain.items())))
        ir = recipe.get('ir')
        if opts['cloud_rayleigh_reduction'] and ir and ir in scene.files:
            self.sources['bt'] = (resample.Source(reader, scene, ir, grid, factor(ir), coarse=low_memory), factor(ir))

        self.sat_lon, self.sat_alt = reader.satellite_position(scene)
        sensor = recipe.get('sensor', reader.SENSOR.lower())
        self.rayleigh = RayleighCorrector(reader.platform(scene), sensor, opts['atmosphere'],
                                          str(opts['rayleigh_correction']).lower())
        rsr = dict(recipe.get('rsr', {}))
        names = {r: recipe[r] for r in ('red', 'green', 'blue', 'nir') if recipe.get(r)}
        self.wave = {r: (rsr.get(n, n), ch[n].wavelength) for r, n in names.items()}
        for band_name, nominal in self.wave.values():          # 先算好有效波長 (平行計算時不必重複查光譜響應)
            self.rayleigh.effective_wavelength(band_name, nominal)
        # 高雲的穿透率與球面反照率：與路徑反射用同一個削減係數 (需要亮溫與 full 校正)
        self.cloud_tau = (bool(opts.get('cloud_rayleigh_transmittance', True)) and 'bt' in self.sources
                          and self.rayleigh.mode == 'full')
        # 臭氧吸收校正：{角色: 臭氧量 x 波段吸收係數}，反射率乘上 exp(該值 x 太陽與衛星的臭氧斜程倍數)
        amount = float(opts.get('ozone') or 0.)
        self.ozone = {}
        if amount > 0.:
            coef = {r: self.rayleigh.ozone_coefficient(*self.wave[r]) for r in self.wave}
            self.ozone = {r: amount * k for r, k in coef.items()}
            print(tr('  臭氧吸收校正：%.2f atm-cm (吸收係數 %s /atm-cm)',
                     '  Ozone absorption correction: %.2f atm-cm (absorption coefficients %s /atm-cm)')
                  % (amount, tr('、', ', ').join(
                      '%s %.3f' % (names[r], coef[r]) for r in ('blue', 'green', 'red', 'nir') if r in coef)))
        self.valid = 0
        self._lock = threading.Lock()
        # 逐像素合成：numba 版 (預設，快數倍且可多核心) 或 numpy 版
        self.kernel = None
        if str(opts.get('engine', 'numba')).lower() == 'numba':
            from . import kernel
            if kernel.HAVE_NUMBA:
                self.kernel = kernel.prepare(self)
            else:
                print(tr('  注意：沒有安裝 numba，真實色改用 numpy 計算 (較慢)',
                         '  NOTE: numba is not installed; true color is computed with numpy instead (slower)'))

    def day(self, values, lon, lat):
        """一批有效像素的白天真實色 (n, 3)。"""
        if self.kernel is not None:
            from . import kernel
            return kernel.compose(values, lon, lat, self.kernel)
        return _compose(values, lon, lat, self.scene.time, self.sat_lon, self.sat_alt,
                        self.wave, bool(self.recipe.get('green')), self.opts, self.rayleigh,
                        self.ozone, self.cloud_tau)

    def sample_bands(self, grid):
        """各波段取樣到 grid：{'blue', 'green', 'nir', 'red_lo', 'red_hi', 'bt'}，可見光已乘上波段相對增益。"""
        bands = {}
        for role, (src, k) in self.sources.items():
            if role == 'red':
                bands['red_lo'] = src.sample(grid, k, self.blur)
                if self.sharpen:
                    bands['red_hi'] = src.sample(grid, 1)
            else:
                bands[role] = src.sample(grid, k)
        for role, g in self.gain.items():
            for key in (('red_lo', 'red_hi') if role == 'red' else (role,)):
                if key in bands:
                    bands[key] *= np.float32(g)
        return bands

    def compute(self, grid):
        """回傳 (列, 行, 3) 的 float32 陣列，數值 0~1 (無資料為 0)，第一列為北邊。"""
        bands = self.sample_bands(grid)

        valid = np.ones(grid.shape, dtype=bool)
        for role in ('blue', 'green', 'nir', 'red_lo', 'red_hi'):
            if role in bands:
                valid &= np.isfinite(bands[role])
        index = np.flatnonzero(valid)
        del valid
        with self._lock:
            self.valid += index.size

        flat = {role: band.ravel() for role, band in bands.items()}
        rgb = np.zeros((grid.size, 3), dtype=np.float32)
        for start in range(0, index.size, CHUNK):
            sel = index[start:start + CHUNK]
            values = {role: band[sel] for role, band in flat.items()}
            lon_t, lat_t = grid.lonlat_at(sel)
            rgb[sel] = self.day(values, lon_t, lat_t)
        return rgb.reshape(grid.shape + (3,))

    def finish(self):
        if self.valid == 0:
            raise ValueError(tr('指定範圍內沒有有效的可見光資料。', 'No valid visible data in the requested range.'))
        print(tr('  真實色計算完成 (%.1f 秒)', '  True color done (%.1f s)') % (time.time() - self.t0))


def make_truecolor(reader, scene, recipe, grid, options=None):
    """產生真實色影像 (整張)：回傳 (列, 行, 3) 的 float32 陣列，數值 0~1，第一列為北邊。"""
    tc = TrueColor(reader, scene, recipe, grid, options)
    rgb = tc.compute(grid)
    tc.finish()
    return rgb


def _compose(values, lon, lat, utc_time, sat_lon, sat_alt, wave, has_green, opts, rayleigh, ozone=None,
             cloud_tau=False):
    """對一批有效像素做校正與合成，回傳 (n, 3) 的 RGB
    ozone：{角色: 臭氧量 x 吸收係數} (TrueColor.ozone；空 = 不校正)；cloud_tau：見 RayleighCorrector.terms 的 scale_tau"""
    sunz, suna = astronomy.sun_angles(utc_time, lon, lat)
    satz, sata = astronomy.satellite_angles(sat_lon, sat_alt, utc_time, lon, lat)
    ssadiff = astronomy.relative_azimuth(suna, sata).astype(np.float32)
    cos_zen = np.cos(np.deg2rad(sunz))
    if ozone:           # 臭氧在平流層，位於雲與 Rayleigh 散射之上：先除以臭氧穿透率，再做 Rayleigh 校正
        airmass = ozone_airmass(sunz) + ozone_airmass(satz)
    sunz = sunz.astype(np.float32)
    satz = satz.astype(np.float32)
    del suna, sata

    refl = {role: sunzen_corr_cos(band, cos_zen, opts['sza_limit'], opts['sza_max']).astype(np.float32)
            for role, band in values.items() if role != 'bt'}
    if ozone:
        for role in refl:
            k = ozone.get('red' if role.startswith('red') else role, 0.)
            if k:
                refl[role] = (refl[role] * np.exp(k * airmass)).astype(np.float32)
    bt = values.get('bt')
    factor = None if bt is None else cloud_height_factor(bt)
    terms = {}

    def rayleigh_corrected(role, color, red_role='red_lo'):
        # 路徑反射等各項只依波長與幾何 (及削減依據) 而定：高解析紅光與低解析紅光共用
        key = (color, red_role if bt is None else '')
        if key not in terms:
            redband = refl[red_role] if bt is None else None
            band_name, nominal = wave[color]
            terms[key] = rayleigh.terms(sunz, satz, ssadiff, band_name, nominal, redband=redband, factor=factor,
                                        scale_tau=cloud_tau)
        return rayleigh.apply(refl[role], terms[key])

    red = rayleigh_corrected('red_lo', 'red')
    blue = rayleigh_corrected('blue', 'blue')
    full = rayleigh.mode == 'full'
    # 近紅外：CIRA 也做 Rayleigh 校正 (0.86 µm 路徑反射約 0.01)；舊版 (subtract) 不校正
    nir = rayleigh_corrected('nir', 'nir') if full else refl['nir']
    frac = float(opts['hybrid_green_fraction'])
    if has_green:
        green = (1. - frac) * rayleigh_corrected('green', 'green') + frac * nir
    elif isinstance(opts['abi_green'], str):    # 'ahi'：合成 0.51 µm 綠光後同樣混入近紅外
        green = (1. - frac) * synthetic_green(blue, red) + frac * nir
    else:                                       # [藍, 紅, 近紅外] 線性權重
        fb, fr, fn = (float(f) for f in opts['abi_green'])
        green = fb * blue + fr * red + fn * nir

    if 'red_hi' in refl:
        red_hi = rayleigh_corrected('red_hi', 'red', red_role='red_hi')
        ratio = sharpening_ratio(red_hi, red)
        red, green, blue = red_hi, green * ratio, blue * ratio

    fade_start, fade_end = float(opts['terminator_fade_start']), float(opts['terminator_fade_end'])
    if fade_end > fade_start:
        red, green, blue = (sunzen_reduction(band, sunz, fade_start, fade_end) for band in (red, green, blue))

    if str(opts.get('tone', 'cira')).lower() == 'enhanced':
        white, knee, toe = float(opts['white_point']), float(opts['shoulder']), float(opts['shadow_toe'])
        rgb = np.stack([enhanced_stretch(band, white, knee, toe) for band in (red, green, blue)], axis=-1)
    else:
        rgb = np.stack([cira_stretch(red), cira_stretch(green), cira_stretch(blue)], axis=-1)
    rgb[~np.isfinite(rgb)] = 0.
    saturation = float(opts.get('saturation', 1.))
    if saturation != 1.:
        rgb = vibrance(rgb, saturation)
    return rgb.astype(np.float32)
