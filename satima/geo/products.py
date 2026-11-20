"""產品：單波段 (灰階)、灰階的第 2 層底圖 (CWA 海陸、LIGHTS 夜光、NASA / NICT Blue Marble)、色表強化 (紅外/水氣)、
真實色 (TrueColor；GeoColor 為加上夜間圖層的第 2 種風格)

單波段的顯示範圍與舊版相同：反射率 0~1 (gray)，亮溫 -80~50 °C (gray_r)。
底圖的灰階與透明度見 satima/landsea.py (CWA、LIGHTS)、satima/bluemarble.py (NASA、NICT) 與設定檔
[landsea]、[nightmap]、[bluemarble]、[nictmap]。
"""
import os
from dataclasses import dataclass

import numpy as np

from . import resample
from ..core import colormaps, parallel, ui
from ..core.i18n import tr


@dataclass
class Product:
    kind: str                   # 'band' / 'map' / 'enhanced' / 'truecolor'
    name: str                   # 波段名稱或真實色產品名稱
    channel: object = None      # Channel (band / map / enhanced)
    cmap: object = None         # CmapEntry (enhanced)
    recipe: dict = None         # 真實色配置 (truecolor)
    tag: str = ''               # 沒有色表時的 {colorbar} (map：底圖名稱，各底圖設定段的 name)
    basemap: str = ''           # map 的底圖 (設定段)：'landsea' (CWA 海陸) / 'nightmap' (LIGHTS 夜光) /
                                # 'bluemarble' (NASA Blue Marble 2004) / 'nictmap' (NICT Blue Marble 2002)

    @property
    def colorbar(self):
        return self.cmap.name if self.cmap is not None else self.tag

    @property
    def display(self):
        """波段顯示名稱，例如 BAND-13、TrueColor (資料集標記 FDHSI / HRFI 由主流程接上)"""
        return self.channel.display if self.channel is not None else self.name

    @property
    def wavelength(self):
        """{wavelength}：例如 10.4µm (真實色為空)"""
        return self.channel.wl_text if self.channel is not None else ''

    @property
    def label(self):
        """{product}：例如 BAND-13-OTT、BAND-03、TrueColor"""
        return self.display + ('-' + self.colorbar if self.colorbar else '')

    def native_res(self, reader, scene=None):
        """原始解析度 (km)：決定「原始解析度」的圖寬。
        原始像素不是方的 (GMS、GOES-1~15 等，東西向取樣比南北向的掃描線間距密) 時取較細的方向：
        輸出像素是方的，以較粗的方向為準會丟掉較細方向的原始像素 (GOES 紅外丟掉 43% 的行)；
        以較細的方向為準則所有原始像素都保留，較粗方向的每個像素畫成多格。"""
        if self.kind == 'truecolor':
            return reader.CHANNELS[self.recipe['red']].res
        res = self.channel.res
        if scene is not None:
            from ..core import maps
            km = maps.source_pixel_km(reader, scene, self.channel)
            if km and min(km) < 0.9 * res:
                ew, ns = km
                print(tr('  原始像素 東西 %.2f km、南北 %.2f km：以較細的 %.2f km 為原始解析度 (%s向每個原始像素約畫 %.2f 格)',
                         '  Native pixel E-W %.2f km, N-S %.2f km: using the finer %.2f km as the native resolution '
                         '(%s: each native pixel is drawn as about %.2f output pixels)') % (
                    ew, ns, min(km), tr('南北', 'N-S') if ew < ns else tr('東西', 'E-W'), max(km) / min(km)))
                res = min(km)
        return res


_GROUPS = (('Visible / Near-IR', ('vis', 'nir')), ('Infrared / Water Vapor', ('swir', 'wv', 'ir')))


def choose(reader, scene, cfg):
    """產品選單 (英文)：依類型分組，波段與真實色可複選；之後各波段再選色表 (可複選)：
    可見光 / 近紅外  只有灰階 (不詢問)
    紅外             1 = 灰階、2 起為 [colormap] ir 的色表
    水氣             1 = 灰階、2 起為 [colormap] wv 的色表
    選了灰階再選底圖 (第 2 層，可複選)：NONE (無底圖)、[landsea] name (預設 CWA)、
    [nightmap] name (預設 LIGHTS，只有紅外)、[bluemarble] name (預設 NASA)、[nictmap] name (預設 NICT)，依此順序編號。
    水氣看不到地表，沒有底圖 (不詢問)"""
    chans = reader.channels(scene)
    tc = reader.truecolor_products(scene)
    sections = [(title, [('band', ch) for ch in chans if ch.kind in kinds]) for title, kinds in _GROUPS]
    sections.append(('RGB', [('truecolor', name) for name in tc]))
    sections = [(title, items) for title, items in sections if items]
    options = [opt for _, items in sections for opt in items]

    dw = max((len(ch.display) for ch in chans), default=0)
    ww = max((len(ch.wl_text) for ch in chans), default=0)
    number = iter(range(1, len(options) + 1))

    def text(kind, item):
        if kind == 'truecolor':
            return '%2d.%s' % (next(number), item)
        return ('%2d.%-*s %-*s %s' % (next(number), dw, item.display, ww, item.wl_text, item.alias)).rstrip()

    print('Products:')
    ui.print_columns([(title, [text(*opt) for opt in items]) for title, items in sections])
    picked = ui.ask_choices(tr('輸入產品編號 (可複選，例如 13 17)：', 'Input the product numbers (several allowed, e.g. 13 17): '),
                            range(1, len(options) + 1))

    # 底圖：(設定段, 名稱, 適用的色表類型)；None = 可見光 / 近紅外。水氣看不到地表，夜光底圖只用於紅外
    bases = [('', 'NONE', (None, 'ir', 'wv')),
             ('landsea', str(cfg.landsea.name).strip(), (None, 'ir')),
             ('nightmap', str(cfg.nightmap.name).strip(), ('ir',)),
             ('bluemarble', str(cfg.bluemarble.name).strip(), (None, 'ir')),
             ('nictmap', str(cfg.nictmap.name).strip(), (None, 'ir'))]
    products = []
    for k in picked:
        kind, item = options[k - 1]
        if kind == 'truecolor':
            products.append(Product('truecolor', item, recipe=tc[item]))
            continue
        cmap_kind = {'ir': 'ir', 'swir': 'ir', 'wv': 'wv'}.get(item.kind)
        ents = colormaps.entries(cfg, cmap_kind) if cmap_kind else []     # 可見光 / 近紅外沒有色表
        menu = {1: 'GRAY'}
        menu.update({i + 2: e.name for i, e in enumerate(ents)})
        keys = ui.ask_menu('%s colormaps:' % item.display, menu,
                           tr('  選擇色表 (可複選，0 = 全部)：', '  Select colormaps (several allowed, 0 = all): '))
        for key in keys:
            if key > 1:
                products.append(Product('enhanced', item.name, channel=item, cmap=ents[key - 2]))
                continue
            # 灰階：第 2 層選底圖 (水氣看不到地表，沒有底圖；夜光底圖只用於紅外)
            avail = [(sec, name) for sec, name, kinds in bases if cmap_kind in kinds]
            for b in ui.ask_menu('%s GRAY basemaps:' % item.display, {i: n for i, (_, n) in enumerate(avail, 1)},
                                 tr('  選擇底圖 (可複選，0 = 全部)：', '  Select basemaps (several allowed, 0 = all): ')):
                sec, name = avail[b - 1]
                products.append(Product('map', item.name, channel=item, tag=name, basemap=sec) if sec else
                                Product('band', item.name, channel=item))
    return products


def render(product, reader, scene, canvas, width, cfg):
    """回傳 8 位元影像：單波段 (列, 行, 4) RGBA (無資料為透明)，真實色 (列, 行, 3) RGB；第一列為北邊。

    上色在這裡完成 (colormaps.colorize，以 float64 依色階定義計算)，不交給 matplotlib：
    結果與整張或分條處理無關，matplotlib 也不必保留浮點影像。
    [memory] low_memory：輸出影像分條處理 (每一條取樣、合成、上色後即釋放)，結果與整張處理逐像素相同。
    [memory] workers：多核心時各條由多個執行緒同時處理 (satima/parallel.py)，結果與單核心相同。
    [image] interpolation：單波段與 GeoColor 夜間紅外，輸出比原始解析度細或相近 (factor = 1) 時的取值方式。"""
    grid = canvas.grid(width)
    low = bool(cfg.memory.low_memory)
    nwork = parallel.workers(cfg)
    ranges = grid.strip_ranges(resample.strip_pixels(grid.size, low, nwork))
    if product.kind == 'truecolor':
        if product.recipe.get('style') == 'geocolor':
            from .config import HOME
            from .truecolor.geocolor import GeoColor
            tc = GeoColor(reader, scene, product.recipe, grid, dict(cfg.truecolor), dict(cfg.geocolor),
                          low_memory=low, home=HOME, interp=cfg.image.interpolation)
        else:
            from .truecolor import TrueColor
            tc = TrueColor(reader, scene, product.recipe, grid, dict(cfg.truecolor), low_memory=low)
        out = np.empty(grid.shape + (3,), dtype=np.uint8)

        def rgb_job(i0, i1):
            out[i0:i1] = colormaps.rgb8(tc.compute(grid.strip(i0, i1)))

        parallel.run([lambda a=a, b=b: rgb_job(a, b) for a, b in ranges], nwork)
        tc.finish()
        return out

    ch = product.channel
    factor = max(1, int(grid.pixel_km / ch.res))
    interp = cfg.image.interpolation
    src = resample.Source(reader, scene, ch.name, grid, factor, coarse=low, interp=interp)
    if product.kind == 'map' and product.basemap in ('bluemarble', 'nictmap'):
        # NASA (Blue Marble 2004) / NICT (Blue Marble 2002)：只有底圖檔不同，雲的合成都依 [bluemarble]
        from . import bluemarble
        from .config import HOME
        c = cfg.bluemarble
        path = os.path.normpath(os.path.join(HOME, os.path.expanduser(str(cfg[product.basemap].image))))
        image = bluemarble.open_image(path, grid.pixel_km, 'gibs' if product.basemap == 'bluemarble' else 'bm2002')
        out = np.empty(grid.shape + (4,), dtype=np.uint8)

        def bm_job(i0, i1):
            g = grid.strip(i0, i1)
            data = src.sample(g, factor)
            if ch.reflective:
                data = _visible(data, scene, g, cfg)
            bluemarble.compose(data, g.lon, g.lat, image, ch.reflective, c, out[i0:i1])

        parallel.run([lambda a=a, b=b: bm_job(a, b) for a, b in ranges], nwork)
        return out
    if product.kind == 'map':
        from . import landsea
        from .config import HOME
        c = cfg.landsea
        night = None
        if product.basemap == 'nightmap':
            from .truecolor import background
            lights = background.night_lights(cfg.geocolor.nightlights, HOME, grid.pixel_km)
            if lights is None:
                raise ValueError(tr('%s 底圖沒有城市燈光資料 ([geocolor] nightlights)',
                                    'the %s basemap has no city lights data ([geocolor] nightlights)') % product.tag)
            night = (dict(cfg.geocolor), lights)
        mask = landsea.open_mask(os.path.normpath(os.path.join(HOME, os.path.expanduser(str(c.mask)))), grid.pixel_km)
        out = np.empty(grid.shape + (4,), dtype=np.uint8)

        def map_job(i0, i1):
            g = grid.strip(i0, i1)
            data = src.sample(g, factor)
            if ch.reflective:
                data = _visible(data, scene, g, cfg)
            landsea.compose(data, g.lon, g.lat, mask, ch.reflective, c, out[i0:i1], night)

        parallel.run([lambda a=a, b=b: map_job(a, b) for a, b in ranges], nwork)
        return out
    if product.kind == 'enhanced':
        cmap, vmin, vmax = product.cmap.load()
    elif ch.reflective:
        cmap, vmin, vmax = 'gray', 0., 1.
    else:
        cmap, vmin, vmax = 'gray_r', -80., 50.
    out = np.empty(grid.shape + (4,), dtype=np.uint8)
    colormaps.colorize(np.zeros(1, np.float32), cmap, vmin, vmax)     # 色表查表在平行處理前先初始化

    def band_job(i0, i1):
        g = grid.strip(i0, i1)
        data = src.sample(g, factor)
        rng = (np.inf, -np.inf)
        if product.kind == 'enhanced':
            # 最高/最低雲頂溫度以原始像素值為準 (內插會把單一像素的極值拉平)
            near = data if interp == 'nearest' or factor > 1 else src.sample(g, factor, interp='nearest')
            if np.isfinite(near).any():
                rng = (float(np.nanmin(near)), float(np.nanmax(near)))
            del near
        elif ch.reflective:
            data = _visible(data, scene, g, cfg)
        colormaps.colorize(data, cmap, vmin, vmax, out=out[i0:i1])
        return rng

    ranges_out = parallel.run([lambda a=a, b=b: band_job(a, b) for a, b in ranges], nwork)
    lo, hi = min(r[0] for r in ranges_out), max(r[1] for r in ranges_out)
    if product.kind == 'enhanced' and lo <= hi:
        if cfg.report.max_tbb:
            print(tr('  最高雲頂溫度：%.2f °C', '  Maximum cloud-top temperature: %.2f °C') % hi)
        if cfg.report.min_tbb:
            print(tr('  最低雲頂溫度：%.2f °C', '  Minimum cloud-top temperature: %.2f °C') % lo)
    return out


def _visible(data, scene, grid, cfg):
    """單波段可見光的太陽天頂角校正、顯示倍數、gamma (逐像素)。"""
    if cfg.visible.sun_zenith_correction:
        from .truecolor import correction
        data = correction.sun_zenith_correction(data, scene.time, grid.lon, grid.lat, limit=75.)
    if float(cfg.visible.scale) != 1.:
        data = data * float(cfg.visible.scale)
    if cfg.visible.gamma:
        data = np.sqrt(np.clip(data, 0, None))
    return data
