"""靜止衛星 (GEO) 自己的設定：預設值與檢查

共同設定 (目錄、地圖、海岸線、圖片大小、標題、標註、檔名 …) 的預設值與讀取方式見 satima/core/config.py；
這裡的 DEFAULTS 只列 GEO 才有的幾節，以及與共同預設值不同的項目。設定檔為 Satima.toml + Satima-GEO.toml。
"""
from ..core.config import (HOME, MAX_SCALE, _is_num, check_colormaps, effective, fail, image_size,  # noqa: F401
                           unwrap)

KIND = 'geo'
PROJECTIONS = ('merc', 'cyl', 'lcc', 'fd')      # fd = 全圓盤 (衛星視角)

DEFAULTS = {
    # ---- 與共同預設值不同的項目 ----
    'coastline': {'color': '#FF9900'},
    'borderline': {'linewidth': 0.5},
    'image': {
        'interpolation': 'bilinear',         # 放大時的取值："bilinear" (雙線性) / "nearest" (最近鄰)
    },
    'annotation': {'text': '{time}    {satellite} {band} {wavelength} {colorbar}'},
    # ---- GEO 才有的設定 ----
    'visible': {'gamma': False, 'sun_zenith_correction': False,
                'scale': 1.0,               # 單波段可見光反射率的顯示倍數 (例如 MTSAT 同舊版 1.4)
                'post_launch_calibration': True,    # GOES-8~15：套用 NOAA 發射後校正 (感測器衰減)
                # GMS-TAPE：可見光感測器故障而缺的掃描線 (每 4 條缺 1 條) 以上下兩條的平均補上；false = 保留為無資料
                'fill_missing_lines': True},
    'report': {'max_tbb': True, 'min_tbb': True},
    # MSG SEVIRI (SEVIRI-HRIT / SEVIRI-NAT / SEVIRI-NC)
    'seviri': {
        'calibration': 'nominal',            # nominal (檔頭係數) / gsics (紅外 GSICS 修正) / meirink (可見光 Meirink 2023)
        'xrit_decompress': '',               # 壓縮 HRIT (-C_) 的 EUMETSAT xRITDecompress 執行檔 (沒有 pyPublicDecompWT 時)
    },
    # GMS (VISSR) 的 JMA 導航
    'navigation': {
        'engine': 'numba',                   # 計算方式：numba / cython / numpy (結果相同)
        'earth_edge_correction': True,       # 以各掃描線實際的地球東西邊緣修正導航偏移 (MSC Technical Note No.26)
    },
    'truecolor': {
        'hybrid_green_fraction': 0.09,
        'abi_green': 'ahi',                  # "ahi" (CIRA 作法) 或 [藍, 紅, 近紅外] 線性權重
        'rayleigh_correction': 'full',       # "full" (SeaDAS / CIRA) 或 "subtract" (舊版，只扣路徑反射)
        'cloud_rayleigh_reduction': True,
        'cloud_rayleigh_transmittance': True,    # 高雲的穿透率與球面反照率也依雲頂削減係數縮小光學厚度
        'ozone': 0.30,                       # 臭氧吸收校正的臭氧量 (atm-cm；0 / false = 不校正)
        'sharpen': True,
        'sza_limit': 88.0, 'sza_max': 95.0,
        'terminator_fade_start': 80.0, 'terminator_fade_end': 90.0,
        'atmosphere': 'us-standard',
        'engine': 'numba',                   # 逐像素合成：numba (快數倍、可多核心) / numpy
        # 各衛星可見光波段的相對增益 (只用於真實色)：Himawari-9 的藍、綠光偏亮，雲會偏青綠；false = 不調整
        'band_gain': {'Himawari-9': {'B01': 0.94, 'B02': 0.98}},
        # 色調曲線：cira (CIRA 論文原式，最亮只到 241、反射率 < 0.025 為 0) / enhanced (白點與暗部改良，中間調不變)
        'tone': 'cira',
        'white_point': 1.25,                 # enhanced：對應全白的反射率 (CIRA 原式 1.5)
        'shoulder': 0.95,                    # enhanced：顯示值超過此值後平滑壓縮到 1
        'shadow_toe': 0.006,                 # enhanced：暗部平滑下限的寬度 (0 = 硬截斷)
        'saturation': 1.0,                   # 飽和度 (1 = 不調整；亮雲不受影響)
    },
    'geocolor': {
        'nightlights': 'nightlights',            # 城市燈光 GeoTIFF (檔案或資料夾)；"" = 不畫
        'twilight_cos': [0.0, 0.18],             # 日夜混合的 cos(太陽天頂角) 區間 (CIRA 目前的圖；論文 [0.1, 0.3])
        'twilight_power': 1.75,                  # (論文 1.5)
        'ir_range': [178.0, 301.0],              # 高雲圖層的亮溫範圍 K (論文 [200, 280])
        'low_cloud_color': [0.55, 0.75, 0.98],
        'land_color': [0.06, 0.03, 0.13],
        'water_color': [0.0, 0.0, 0.0],
    },
    # 海陸底圖 (GRAY 的第 2 層)：灰階雲疊在陸地 / 海洋底色上 (仿 CWA 官網「紅外線彩色雲圖」；紅外與可見光)
    'landsea': {
        'name': 'CWA',                           # 選單與檔名 {colorbar} 的名稱
        'mask': 'landsea/gshhs_f_30s.tif',       # 海陸遮罩；不存在時由 basemap 的 GSHHS 全解析度海岸線自動產生
        'land_color': [124, 145, 85],            # 陸地 (0~255)
        'ocean_color': [31, 40, 118],            # 海洋與湖泊
        # 紅外灰階：[亮溫 °C, 灰階 0~255] 節點之間線性內插
        'ir_gray': [[-95, 255], [-80, 235], [-60, 207], [-40, 176], [-20, 144], [0, 113], [20, 86], [80, 0]],
        'ir_alpha': [26.1, -44.7],               # 雲由全透明 -> 不透明的亮溫 (°C)
        'vis_alpha': [0.10, 0.50],               # 可見光：由全透明 -> 不透明的反射率
        'alpha_power': 0.8,                      # 透明度漸變的次方
        'lines': False,                          # 是否畫海岸線、國界、經緯線 (依 [coastline] 等設定)
    },
    # 夜光底圖 (GRAY 的第 2 層，只有紅外)：雲的透明度、顏色與城市燈光依 [geocolor]，海陸遮罩依 [landsea] mask
    'nightmap': {
        'name': 'LIGHTS',                        # 選單與檔名 {colorbar} 的名稱
        'lines': False,
    },
    # Blue Marble 底圖 NASA (GRAY 的第 2 層；紅外與可見光)：白雲依亮溫 / 反射率半透明，疊在 NASA Blue Marble 上
    # (仿 NICT Himawari-8 Real-time Web；透明度以 NICT 的 B13 / B03 雲圖對照同時次 HSD 反推)。
    # 雲的設定 (cloud_color、ir_alpha、vis_alpha) 也用於 [nictmap]
    'bluemarble': {
        'name': 'NASA',                          # 選單與檔名 {colorbar} 的名稱
        # 底圖：Blue Marble Next Generation 2004；不存在時由 NASA GIBS 自動下載 (約 1 km，67.5 MB)
        'image': 'bluemarble/BlueMarble_NextGeneration_GIBS_1km.tif',
        'cloud_color': [255, 255, 255],          # 雲的顏色 (0~255)
        'ir_alpha': [[10.5, 0], [-80.0, 230.5]],  # 紅外：兩點 [亮溫 °C, 不透明度 0~255]，與 10.4 µm 輻射量線性
        'vis_alpha': [[0.0, 0], [1.0, 255]],     # 可見光：兩點 [反射率, 不透明度 0~255]，線性
        'lines': False,                          # 是否畫海岸線、國界、經緯線
    },
    # Blue Marble 底圖 NICT：同 [bluemarble]，底圖改為 NICT 網站所用的 Blue Marble 2002 (較亮、海洋為藍紫色)
    'nictmap': {
        'name': 'NICT',
        # 底圖：Blue Marble 2002 land_shallow_topo；不存在時由 Wikimedia Commons 自動下載並轉檔 (約 1.85 km，44.9 MB)
        'image': 'bluemarble/BlueMarble_2002_land_shallow_topo_2km.tif',
        'lines': False,
    },
    'colormap': {
        'ir': [{'name': n.upper(), 'file': n} for n in
               ('bd', 'color', 'cc', 'rammb', 'nhc', 'ott', 'rbtop', 'avn', 'rb', 'mb', 'ca', 'cwbm')],
        'wv': [{'name': n.upper(), 'file': n} for n in ('nrl', 'ssd', 'noaa')],
    },
}


def _rgb(v, top):
    return isinstance(v, list) and len(v) == 3 and all(_is_num(x) and 0 <= x <= top for x in v)


def _pair(v):
    return isinstance(v, list) and len(v) == 2 and all(_is_num(x) for x in v)


def validate(cfg):
    """GEO 各節的檢查 (共同的部分已由 core.config.validate 檢查)。"""
    if 'native_scale' in cfg.image:
        fail('[image] native_scale 已併入 size，例如 size = "x%s" (原始解析度的 %s 倍)',
             '[image] native_scale was merged into size, e.g. size = "x%s" (%s times the native resolution)',
             cfg.image.native_scale, cfg.image.native_scale)
    if cfg.image.interpolation not in ('bilinear', 'nearest'):
        fail('[image] interpolation 須為 "bilinear" 或 "nearest" (目前為 %r)',
             '[image] interpolation must be "bilinear" or "nearest" (got %r)', cfg.image.interpolation)
    check_colormaps(cfg, ('ir', 'wv'))
    scale = cfg.visible.scale
    if not _is_num(scale) or scale <= 0:
        fail('[visible] scale 須為正數 (目前為 %r)', '[visible] scale must be a positive number (got %r)', scale)
    if not isinstance(cfg.visible.fill_missing_lines, bool):
        fail('[visible] fill_missing_lines 須為 true 或 false (目前為 %r)',
             '[visible] fill_missing_lines must be true or false (got %r)', cfg.visible.fill_missing_lines)

    tc = cfg.truecolor
    if 'abi_green_fractions' in tc:
        fail('[truecolor] abi_green_fractions 已改名為 abi_green。預設 "ahi" 為 CIRA 作法；'
             '要沿用線性權重請改寫成 abi_green = %s',
             '[truecolor] abi_green_fractions was renamed to abi_green. The default "ahi" is the CIRA method; '
             'to keep the linear weights write abi_green = %s', list(tc.abi_green_fractions))
    green = tc.abi_green
    if not ((isinstance(green, str) and green.strip().lower() == 'ahi') or
            (isinstance(green, list) and len(green) == 3 and all(_is_num(v) for v in green))):
        fail('[truecolor] abi_green 須為 "ahi" 或 3 個數值 [藍, 紅, 近紅外] (目前為 %r)',
             '[truecolor] abi_green must be "ahi" or 3 numbers [blue, red, near-IR] (got %r)', green)
    if str(tc.rayleigh_correction).strip().lower() not in ('full', 'subtract'):
        fail('[truecolor] rayleigh_correction 須為 "full" 或 "subtract" (目前為 %r)',
             '[truecolor] rayleigh_correction must be "full" or "subtract" (got %r)', tc.rayleigh_correction)
    if isinstance(green, str):
        tc.abi_green = green.strip().lower()
    tc.rayleigh_correction = str(tc.rayleigh_correction).strip().lower()
    if str(tc.engine).strip().lower() not in ('numba', 'numpy'):
        fail('[truecolor] engine 須為 "numba" 或 "numpy" (目前為 %r)',
             '[truecolor] engine must be "numba" or "numpy" (got %r)', tc.engine)
    tc.engine = str(tc.engine).strip().lower()
    if str(tc.tone).strip().lower() not in ('cira', 'enhanced'):
        fail('[truecolor] tone 須為 "cira" 或 "enhanced" (目前為 %r)',
             '[truecolor] tone must be "cira" or "enhanced" (got %r)', tc.tone)
    tc.tone = str(tc.tone).strip().lower()
    for key, lo, hi in (('white_point', 0.3, 3.0), ('shoulder', 0.5, 1.0), ('shadow_toe', 0.0, 0.05),
                        ('saturation', 0.0, 3.0)):
        v = tc[key]
        if not _is_num(v) or not lo <= v <= hi:
            fail('[truecolor] %s 須為 %g~%g (目前為 %r)', '[truecolor] %s must be %g~%g (got %r)', key, lo, hi, v)
    if tc.ozone is not False and not (_is_num(tc.ozone) and 0. <= tc.ozone <= 1.):
        fail('[truecolor] ozone 須為 0~1 的臭氧量 (atm-cm，例如 0.30 = 300 DU) 或 false (目前為 %r)',
             '[truecolor] ozone must be an ozone amount 0~1 (atm-cm, e.g. 0.30 = 300 DU) or false (got %r)', tc.ozone)
    gain = tc.get('band_gain')
    if gain not in (None, False):
        ok = isinstance(gain, dict) and all(
            isinstance(bands, dict) and all(_is_num(g) and g > 0 for g in bands.values()) for bands in gain.values())
        if not ok:
            fail('[truecolor] band_gain 須為 false 或 { 衛星 = { 波段 = 正數 } }，'
                 '例如 { Himawari-9 = { B01 = 0.94, B02 = 0.98 } } (目前為 %r)',
                 '[truecolor] band_gain must be false or { satellite = { band = positive number } }, '
                 'e.g. { Himawari-9 = { B01 = 0.94, B02 = 0.98 } } (got %r)', unwrap(gain))

    ls = cfg.landsea
    for key in ('land_color', 'ocean_color'):
        if not _rgb(ls[key], 255):
            fail('[landsea] %s 須為 3 個 0~255 的數值 [R, G, B] (目前為 %r)',
                 '[landsea] %s must be 3 numbers 0~255 [R, G, B] (got %r)', key, ls[key])
    nodes = ls.ir_gray
    if not (isinstance(nodes, list) and len(nodes) >= 2 and all(_pair(n) for n in nodes)
            and len({float(n[0]) for n in nodes}) == len(nodes)):
        fail('[landsea] ir_gray 須為至少 2 個 [亮溫 °C, 灰階 0~255]，亮溫不可重複 (目前為 %r)',
             '[landsea] ir_gray must be at least 2 nodes [brightness temperature °C, gray 0~255] with distinct '
             'temperatures (got %r)', unwrap(nodes))
    if not (_pair(ls.ir_alpha) and ls.ir_alpha[0] != ls.ir_alpha[1]):
        fail('[landsea] ir_alpha 須為 [透明的亮溫, 不透明的亮溫] 兩個不同的數值 (目前為 %r)',
             '[landsea] ir_alpha must be two different numbers [transparent temperature, opaque temperature] '
             '(got %r)', ls.ir_alpha)
    if not (_pair(ls.vis_alpha) and ls.vis_alpha[0] != ls.vis_alpha[1]):
        fail('[landsea] vis_alpha 須為 [透明的反射率, 不透明的反射率] 兩個不同的數值 (目前為 %r)',
             '[landsea] vis_alpha must be two different numbers [transparent reflectance, opaque reflectance] '
             '(got %r)', ls.vis_alpha)
    if not (_is_num(ls.alpha_power) and ls.alpha_power > 0):
        fail('[landsea] alpha_power 須為正數 (目前為 %r)', '[landsea] alpha_power must be a positive number (got %r)',
             ls.alpha_power)

    bm = cfg.bluemarble
    if not _rgb(bm.cloud_color, 255):
        fail('[bluemarble] cloud_color 須為 3 個 0~255 的數值 [R, G, B] (目前為 %r)',
             '[bluemarble] cloud_color must be 3 numbers 0~255 [R, G, B] (got %r)', bm.cloud_color)
    for key, zh, en in (('ir_alpha', '亮溫 °C', 'brightness temperature °C'), ('vis_alpha', '反射率', 'reflectance')):
        v = bm[key]
        if not (isinstance(v, list) and len(v) == 2 and all(_pair(n) for n in v) and float(v[0][0]) != float(v[1][0])):
            fail('[bluemarble] %s 須為兩個 [%s, 不透明度 0~255]，兩點的第一個數值不可相同 (目前為 %r)',
                 '[bluemarble] %s must be two nodes [%s, opacity 0~255] with different first values (got %r)',
                 key, zh if _zh() else en, unwrap(v))
    if any(float(t) <= -273.15 for t, _ in bm.ir_alpha):
        fail('[bluemarble] ir_alpha 的亮溫須高於絕對零度 (目前為 %r)',
             '[bluemarble] ir_alpha temperatures must be above absolute zero (got %r)', unwrap(bm.ir_alpha))
    for sec in ('bluemarble', 'nictmap'):
        if not str(cfg[sec].image or '').strip():
            fail('[%s] image 不可為空白', '[%s] image must not be empty', sec)
    names = {}
    for sec in ('landsea', 'nightmap', 'bluemarble', 'nictmap'):
        name = str(cfg[sec].name or '').strip()
        if not name:
            fail('[%s] name 不可為空白', '[%s] name must not be empty', sec)
        if name.upper() == 'NONE':
            fail('[%s] name 不可為 NONE (選單的「無底圖」)', '[%s] name must not be NONE (the "no basemap" menu item)', sec)
        if name.upper() in names:
            fail('[%s] name 不可與 [%s] name 相同 (目前為 %r)', '[%s] name must differ from [%s] name (got %r)',
                 sec, names[name.upper()], name)
        names[name.upper()] = sec

    gc = cfg.geocolor
    tw = gc.twilight_cos
    if not (_pair(tw) and 0 <= tw[0] < tw[1] <= 1):
        fail('[geocolor] twilight_cos 須為 [下限, 上限]，0 <= 下限 < 上限 <= 1 (目前為 %r)',
             '[geocolor] twilight_cos must be [low, high] with 0 <= low < high <= 1 (got %r)', tw)
    for key in ('low_cloud_color', 'land_color', 'water_color'):
        if not _rgb(gc[key], 1):
            fail('[geocolor] %s 須為 3 個 0~1 的數值 [R, G, B] (目前為 %r)',
                 '[geocolor] %s must be 3 numbers 0~1 [R, G, B] (got %r)', key, gc[key])
    if not (_is_num(gc.twilight_power) and gc.twilight_power > 0):
        fail('[geocolor] twilight_power 須為正數 (目前為 %r)',
             '[geocolor] twilight_power must be a positive number (got %r)', gc.twilight_power)
    ir = gc.ir_range
    if not (_pair(ir) and 100 < ir[0] < ir[1] < 350):
        fail('[geocolor] ir_range 須為 [冷端, 暖端] (K)，例如 [178, 301] (目前為 %r)',
             '[geocolor] ir_range must be [cold end, warm end] (K), e.g. [178, 301] (got %r)', ir)
    if str(cfg.seviri.calibration).lower() not in ('nominal', 'gsics', 'meirink'):
        fail('[seviri] calibration 須為 "nominal"、"gsics" 或 "meirink" (目前為 %r)',
             '[seviri] calibration must be "nominal", "gsics" or "meirink" (got %r)', cfg.seviri.calibration)
    if str(cfg.navigation.engine).lower() not in ('numba', 'cython', 'numpy'):
        fail('[navigation] engine 須為 "numba"、"cython" 或 "numpy" (目前為 %r)',
             '[navigation] engine must be "numba", "cython" or "numpy" (got %r)', cfg.navigation.engine)


def _zh():
    from ..core import i18n
    return i18n.language() == 'zh'
