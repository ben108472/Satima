"""產品：單頻道亮溫、PCT (極化修正亮溫)、RGB 合成

  band  單頻道亮溫 (K)：GRAY 或 [colormap] 的色表
  pct   PCT = (1 + θ) x V - θ x H：去掉海面低放射率造成的低亮溫，凸顯降水 / 冰晶散射。
        θ 見 [pct] (預設 Cecil and Chronis 2018)；同頻率同時有 V、H 的頻道才有
  rgb   三個色槍各取 PCT / V / H，線性拉伸到 [rgb.<名稱>] 的範圍 (可反相、gamma)；配方見設定檔
        (NRL-37、NRL-89 為 NRL 的 37 / 89 GHz color composite；RAMMB-37 依 RAMMB TC Realtime 的說明。
        RAMMB 的 85~92 GHz 產品不是 RGB，是單頻道亮溫配色階：色表 RAMMB)
各頻道先各自取樣到輸出格點 (swath.py)，再逐像素計算；同一組幾何的頻道 (V、H) 共用小數座標。
"""
import re
from dataclasses import dataclass

import numpy as np

from . import colormaps, config, grid as gridmod
from ..core import parallel, ui
from ..core.i18n import tr
from .reader import freq_text, kind_of
from .swath import Resampler

# RGB 配方的頻率類別 -> 標稱頻率 (同類別有多組 V/H 時取最接近的)
_NOMINAL = {'tb10': 10.65, 'tb19': 19.0, 'tb37': 37.0, 'tb89': 89.0}


@dataclass
class Product:
    kind: str                   # 'band' / 'pct' / 'rgb'
    name: str                   # 顯示名稱：89V、PCT89、NRL-37
    channels: tuple = ()        # band：(頻道,)；pct / rgb：(V, H)
    cmap: object = None         # CmapEntry；None = 灰階 (band / pct)
    theta: float = None         # pct / rgb 的 θ
    recipe: dict = None         # rgb 配方

    @property
    def colorbar(self):
        return self.cmap.name if self.cmap is not None else ''

    @property
    def display(self):
        return self.name

    @property
    def frequency(self):
        """{frequency}：例如 89GHz"""
        return self.channels[0].freq_label if self.channels else ''

    @property
    def label(self):
        return self.display + ('-' + self.colorbar if self.colorbar else '')

    def native_res(self, reader, scene):
        """原始取樣間距 (km)：決定「原始解析度」的圖寬。取沿掃描與沿軌較細的一個；
        A / B 掃描合併時沿軌間距減半。"""
        geos = reader.geo_cached((scene.parts or [scene])[0], self.channels[0].group)
        across, along = geos[0].spacing()
        vals = [v for v in (across, along / len(geos)) if np.isfinite(v) and v > 0]
        return min(vals) if vals else None


_SPEC = re.compile(r'^\s*(\d+(?:\.\d+)?)\s*(?:(?:±|\+/?-)\s*(\d+(?:\.\d+)?))?\s*(Q?[VH])?\s*$', re.I)
_OFFSET = re.compile(r'±(\d+(?:\.\d+)?)')
TOLERANCE = 0.12            # 「任意頻道」配方：頻道的頻率與指定值相差在這個比例以內才算符合


def find_channel(chans, spec):
    """「任意頻道」RGB 配方的頻道指定 -> 最符合的頻道 (沒有為 None)。
    spec：'89'、'150'、'183±7'、'89H'、'23.8'…… 頻率 (GHz)，可加 ±偏移、極化 (V / H)。
    規則：同一個頻率類別 (reader.kind_of)、頻率相差 TOLERANCE 以內；有 ± 時取偏移最接近的
    (名稱沒有 ± 的頻道以 |頻率 - 指定頻率| 當偏移：MHS 的 190.31 GHz 視為 183±7)，否則取頻率最接近的；
    沒指定極化時 H 優先 (同 NRL 的 89H)。"""
    m = _SPEC.match(str(spec))
    if not m:
        return None
    freq, offset, pol = float(m.group(1)), (float(m.group(2)) if m.group(2) else None), (m.group(3) or '').upper()
    kind = kind_of(freq)
    best = None
    for ch in chans:
        if ch.kind != kind or not ch.freq or abs(ch.freq - freq) > TOLERANCE * freq:
            continue
        if pol and not ch.pol.upper().endswith(pol[-1]):
            continue
        named = _OFFSET.search(ch.name)
        if offset is None and named and kind != 'tb183':
            continue                            # 沒要求偏移時不取探空的旁帶頻道 (57.29±0.217…)
        mine = float(named.group(1)) if named else abs(ch.freq - freq)
        rank = {'H': 0, '': 1, 'V': 2}.get(ch.pol.upper()[-1:], 1)
        if offset is not None:
            key = (abs(mine - offset), abs(ch.freq - freq), rank)
        else:
            key = (abs(ch.freq - freq), mine if named else 0., rank)
        if best is None or key < best[0]:
            best = (key, ch)
    return best[1] if best else None


def pct_name(ch):
    return 'PCT' + freq_text(ch.freq)


def available(reader, scene, cfg):
    """這個場景可做的產品 (不含色表)：[(分組標題, [Product, ...]), ...]"""
    chans = reader.channels(scene)
    pairs = reader.pairs(scene)
    bands = [Product('band', ch.name, (ch,)) for ch in chans]
    pcts = [Product('pct', pct_name(v), (v, h), theta=float(cfg.pct[kind]))
            for kind in ('tb10', 'tb19', 'tb37', 'tb89') for v, h in pairs.get(kind, [])]
    rgbs = []
    for name, rec in cfg.rgb.items():
        if rec.get('enabled', True) is False:
            continue
        if 'band' not in rec:
            # 「任意頻道」配方：各色槍指定頻率 (可再減另一個頻道)；三個色槍都找得到頻道才列出
            used, guns = [], {}
            for gun in ('red', 'green', 'blue'):
                idx = []
                for key in ('channel', 'minus'):
                    ch = find_channel(chans, rec[gun][key]) if key in rec[gun] else None
                    if ch is not None and ch not in used:
                        used.append(ch)
                    idx.append(used.index(ch) if ch is not None else None)
                guns[gun] = tuple(idx)
            if all(guns[g][0] is not None and (guns[g][1] is not None or 'minus' not in rec[g]) for g in guns):
                rgbs.append(Product('rgb', name, tuple(used), recipe=dict(rec, _guns=guns)))
            continue
        if rec.band not in pairs:
            continue
        v, h = min(pairs[rec.band], key=lambda p: abs(p[0].freq - _NOMINAL[rec.band]))
        rgbs.append(Product('rgb', name, (v, h), theta=float(rec.get('theta', cfg.pct[rec.band])), recipe=dict(rec)))
    sections = [('Brightness temperature', bands), ('PCT (polarization-corrected)', pcts), ('RGB', rgbs)]
    return [(title, items) for title, items in sections if items]


def choose(reader, scene, cfg):
    """產品選單 (可複選)；單頻道與 PCT 再選色表 (1 = GRAY，2 起為 [colormap] 的色表，可複選)。"""
    sections = available(reader, scene, cfg)
    options = [p for _, items in sections for p in items]
    number = iter(range(1, len(options) + 1))
    width = max(len(p.name) for p in options)

    def text(p):
        note = p.channels[0].note if p.kind == 'band' else ''
        return ('%2d.%-*s %s' % (next(number), width, p.name, note)).rstrip()

    print('Products:')
    ui.print_columns([(title, [text(p) for p in items]) for title, items in sections])
    picked = ui.ask_choices(tr('輸入產品編號 (可複選，例如 13 17)：',
                               'Input the product numbers (several allowed, e.g. 13 17): '),
                            range(1, len(options) + 1))
    products = []
    for k in picked:
        prod = options[k - 1]
        if prod.kind == 'rgb':
            products.append(prod)
            continue
        ents = colormaps.entries(cfg, colormaps.group_of(prod.channels[0].kind, prod.kind == 'pct'))
        menu = {1: 'GRAY'}
        menu.update({i + 2: e.name for i, e in enumerate(ents)})
        for key in ui.ask_menu('%s colormaps:' % prod.name, menu,
                               tr('  選擇色表 (可複選，0 = 全部)：', '  Select colormap (several allowed, 0 = all): ')):
            products.append(Product(prod.kind, prod.name, prod.channels, cmap=ents[key - 2] if key > 1 else None,
                                    theta=prod.theta))
    return products


def resolve(product, reader, scene, cfg):
    """把依第一筆資料選的產品對應到另一個場景 (依名稱)；該場景做不出來時回傳 None。"""
    for _, items in available(reader, scene, cfg):
        for p in items:
            if p.kind == product.kind and p.name == product.name:
                return Product(p.kind, p.name, p.channels, cmap=product.cmap, theta=p.theta, recipe=p.recipe)
    return None


def resampler(canvas, width, cfg):
    """某張畫布、某個寬度的取樣器 (各組幾何的定位器在同一張地圖的多個產品間共用)。"""
    cache = canvas.__dict__.setdefault('_resamplers', {})
    key = (width, cfg.swath.interpolation, str(cfg.swath.gauss_sigma))
    if key not in cache:
        sigma, factor = config.gauss_sigma(cfg.swath.gauss_sigma)
        cache[key] = Resampler(canvas.grid(width), cfg.swath.interpolation, sigma, factor or 1.)
    return cache[key]


def _stretch(x, lo, hi, inverse=False, gamma=1.):
    y = np.clip((x - lo) / (hi - lo), 0., 1.)
    if inverse:
        y = 1. - y
    if gamma != 1.:
        y = y ** (1. / gamma)
    return y


def render(product, reader, scene, canvas, width, cfg):
    """回傳 8 位元 RGBA 影像 (列, 行, 4)，第一列為北邊；無資料處透明 (顯示底色)。
    輸出影像分條處理 ([memory] low_memory)，多核心時各條同時處理 ([memory] workers)，結果與整張處理相同。"""
    grid = canvas.grid(width)
    rs = resampler(canvas, width, cfg)
    nwork = parallel.workers(cfg)
    ranges = grid.strip_ranges(gridmod.strip_pixels(grid.size, bool(cfg.memory.low_memory), nwork))
    swaths = [reader.swath(scene, ch) for ch in product.channels]
    out = np.zeros(grid.shape + (4,), dtype=np.uint8)
    qworkers = 1 if nwork > 1 and len(ranges) > 1 else -1          # KD-tree 查詢：各條已平行時不再開執行緒

    def values(i0, i1, cache, method=None):
        data = [rs.sample(sw, i0, i1, cache, qworkers, method) for sw in swaths]
        if product.kind == 'band':
            return data[0]
        v, h = data
        return (1. + product.theta) * v - product.theta * h

    if product.kind == 'rgb':
        rec = product.recipe

        def rgb_job(i0, i1):
            cache = {}
            data = [rs.sample(sw, i0, i1, cache, qworkers) for sw in swaths]
            ok = np.isfinite(data[0])
            for d in data[1:]:
                ok &= np.isfinite(d)
            if '_guns' not in rec:
                v, h = data
                src = {'v': v, 'h': h, 'pct': (1. + product.theta) * v - product.theta * h}
            for k, gun in enumerate(('red', 'green', 'blue')):
                g = rec[gun]
                if '_guns' in rec:
                    i, j = rec['_guns'][gun]
                    value = data[i] - data[j] if j is not None else data[i]
                else:
                    value = src[str(g['source']).lower()]
                x = _stretch(value, float(g['range'][0]), float(g['range'][1]),
                             bool(g.get('inverse', False)), float(g.get('gamma', 1.)))
                out[i0:i1, :, k] = colormaps.rgb8(np.nan_to_num(x))
            out[i0:i1, :, 3] = np.where(ok, 255, 0)
            out[i0:i1][~ok] = 0

        parallel.run([lambda a=a, b=b: rgb_job(a, b) for a, b in ranges], nwork)
        return out

    if product.cmap is not None:
        cmap, vmin, vmax = product.cmap.load()
    else:
        lo, hi = (float(v) for v in cfg.gray.range)
        cmap, vmin, vmax = ('gray_r' if cfg.gray.reverse else 'gray'), min(lo, hi), max(lo, hi)
    colormaps.colorize(np.zeros(1, np.float32), cmap, vmin, vmax)     # 色表查表在平行處理前先初始化
    report = bool(cfg.report.min_tb or cfg.report.max_tb)

    def job(i0, i1):
        cache = {}
        data = values(i0, i1, cache)
        colormaps.colorize(data, cmap, vmin, vmax, out=out[i0:i1])
        rng = (np.inf, -np.inf)
        if report:
            # 極值以原始取樣值為準 (內插不會超出取樣值，但會把單點的極值拉平)
            near = data if rs.method == 'nearest' else values(i0, i1, cache, 'nearest')
            if np.isfinite(near).any():
                rng = (float(np.nanmin(near)), float(np.nanmax(near)))
        return rng

    res = parallel.run([lambda a=a, b=b: job(a, b) for a, b in ranges], nwork)
    lo, hi = min(r[0] for r in res), max(r[1] for r in res)
    if lo <= hi:
        if cfg.report.min_tb:
            print(tr('  範圍內最低：%.2f K', '  Minimum inside the map range: %.2f K') % lo)
        if cfg.report.max_tb:
            print(tr('  範圍內最高：%.2f K', '  Maximum inside the map range: %.2f K') % hi)
    return out
