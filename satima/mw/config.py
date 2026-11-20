"""微波 (MW) 自己的設定：預設值與檢查

共同設定 (目錄、地圖、海岸線、圖片大小、標題、標註、檔名 …) 的預設值與讀取方式見 satima/core/config.py；
這裡的 DEFAULTS 只列 MW 才有的幾節。設定檔為 Satima.toml + Satima-MW.toml。
"""
from ..core.config import (HOME, MAX_SCALE, _is_num, check_colormaps, effective, fail, image_size,  # noqa: F401
                           unwrap)

KIND = 'mw'
PROJECTIONS = ('merc', 'cyl', 'lcc')
BANDS = ('tb10', 'tb19', 'tb37', 'tb89')            # 可做 PCT / RGB 的頻率類別 (同時有 V、H)
CMAP_GROUPS = ('tb', 'tb37', 'tb89', 'tb150', 'pct', 'pct37', 'pct89')

DEFAULTS = {
    'swath': {
        'interpolation': 'gauss',            # "gauss" (高斯距離加權) / "cubic" (三次，限制斜率) / "bilinear" / "nearest" (關閉內插)
        'gauss_sigma': 'x0.7',               # gauss 的 σ："xN" = 較粗方向取樣間距的 N 倍；數字 = 固定的公里數 (NRL 用 10)
        'merge_ab': True,                    # AMSR2 / AMSR3 的 89 GHz：A、B 掃描合併成一個頻道；false = 分開列出
        'merge': 'ask',                      # 選到同衛星、時間相接的多個檔時合併成一張："ask" / true / false
        'merge_gap': 120.0,                  # 兩個檔的間隔在幾秒內視為相接
        'mask_bad_quality': False,           # true = L1C Quality < 0 的取樣點整點不用 (含同組其他正常的頻道)；
                                             # 出錯頻道的亮溫 PPS 已設成無資料，一般不必再遮
        'min_coverage': 0.002,               # 掃描帶涵蓋地圖的比例低於此值時略過
    },
    # 單頻道 / PCT 的灰階 (選單的 GRAY)
    'gray': {'range': [150.0, 300.0], 'reverse': True},      # reverse = true：冷 (冰晶散射) 為白
    # PCT = (1 + θ) x V - θ x H 的 θ，預設為 Cecil and Chronis (2018)；85 GHz 的 Spencer et al. (1989) 為 0.818
    'pct': {'tb10': 1.5, 'tb19': 1.4, 'tb37': 1.15, 'tb89': 0.7},
    # RGB 合成：band = 頻率類別 (tb37 / tb89 ...)；theta = 該配方 PCT 的 θ (不寫 = [pct] 的值)；
    # 各色槍 source = "pct" / "v" / "h"，range = [最小, 最大] (K)，inverse = 反相，gamma (預設 1)
    'rgb': {
        'NRL-37': {
            'band': 'tb37', 'theta': 1.181,
            'red': {'source': 'pct', 'range': [260.0, 280.0], 'inverse': True},
            'green': {'source': 'v', 'range': [180.0, 300.0]},
            'blue': {'source': 'h', 'range': [160.0, 300.0]},
        },
        'RAMMB-37': {
            'band': 'tb37',
            'red': {'source': 'pct', 'range': [260.0, 280.0], 'inverse': True},
            'green': {'source': 'v', 'range': [160.0, 300.0]},
            'blue': {'source': 'h', 'range': [180.0, 310.0]},
        },
        'NRL-89': {
            'band': 'tb89', 'theta': 0.818,
            'red': {'source': 'pct', 'range': [220.0, 310.0], 'inverse': True},
            'green': {'source': 'h', 'range': [240.0, 300.0]},
            'blue': {'source': 'v', 'range': [270.0, 290.0]},
        },
        # RAMMB 的 85~92 GHz 產品不是 RGB，而是單頻道亮溫配色階 (色表 rammb89)
        # ---- 試作 (非官方)：「任意頻道」配方，給沒有雙極化的探空儀 (ATMS、MHS、AMSU-B) ----
        'RGB-150': {
            'red': {'channel': '89', 'range': [220.0, 280.0], 'inverse': True},
            'green': {'channel': '150', 'range': [180.0, 380.0]},
            'blue': {'channel': '183±7', 'range': [200.0, 340.0]},
        },
    },
    'report': {'min_tb': True, 'max_tb': False},     # 色表產品：印出地圖範圍內的最低 / 最高值
    # 色表選單 (順序即選單編號；GRAY 固定為 1)。分組：tb89 (80~95 GHz)、tb37 (30~40 GHz)、tb150 (140~200 GHz)、
    # tb (其他頻率)；pct89 / pct37 / pct (其他頻率的 PCT)。file = colormap 目錄中的 .cmap 檔名，
    # 或 mpl = matplotlib 色表名稱 (需 vmin / vmax)；單位 K
    'colormap': {
        'tb89': [{'name': 'NRL', 'file': 'nrl89'}, {'name': 'NRL2', 'file': '89h'}, {'name': 'RAMMB', 'file': 'rammb89'},
                 {'name': 'MWNRL', 'file': 'mwnrl'}],
        'tb37': [{'name': 'NRL', 'file': '37h'}, {'name': 'MWNRL', 'file': 'mwnrl'}],
        'tb150': [{'name': 'NRL', 'file': '150h'}, {'name': 'MWNRL', 'file': 'mwnrl'}],
        'tb': [{'name': 'MWNRL', 'file': 'mwnrl'}, {'name': 'NRL37', 'file': '37h'}],
        'pct89': [{'name': 'NRL', 'file': 'nrlpct'}, {'name': 'NRL2', 'file': '89pct'}],
        'pct37': [{'name': 'NRL', 'file': '37pct'}],
        'pct': [{'name': 'NRL37', 'file': '37pct'}, {'name': 'NRL89', 'file': 'nrlpct'}],
    },
}


def _pair(v):
    return isinstance(v, list) and len(v) == 2 and all(_is_num(x) for x in v) and v[0] != v[1]


def gauss_sigma(value):
    """[swath] gauss_sigma -> (固定的 σ 公里數 或 None, 取樣間距的倍數 或 None)。"""
    if isinstance(value, bool):
        raise ValueError(value)
    text = str(value).strip().lower()
    if text.startswith('x'):
        k = float(text[1:])                 # 不是數字時 float() 會丟出 ValueError
        if 0.3 <= k <= 3:
            return None, k
        raise ValueError(value)
    km = float(text[:-2] if text.endswith('km') else text)     # 10、"10"、"10km" 都可以
    if 1 <= km <= 200:
        return km, None
    raise ValueError(value)


def validate(cfg):
    """MW 各節的檢查 (共同的部分已由 core.config.validate 檢查)。"""
    from .swath import ALL_METHODS
    sw = cfg.swath
    sw.interpolation = str(sw.interpolation).strip().lower()
    if sw.interpolation not in ALL_METHODS:
        fail('[swath] interpolation 須為 "gauss"、"cubic"、"bilinear" 或 "nearest" (目前為 %r)',
             '[swath] interpolation must be "gauss", "cubic", "bilinear" or "nearest" (got %r)', sw.interpolation)
    try:
        gauss_sigma(sw.gauss_sigma)
    except ValueError:
        fail('[swath] gauss_sigma 須為 "x0.3"~"x3" (較粗方向取樣間距的倍數；有 x) 或 1~200 的公里數 (沒有 x，例如 10) '
             '(目前為 %r)',
             '[swath] gauss_sigma must be "x0.3"~"x3" (multiples of the coarser sampling interval; with the x) '
             'or 1~200 km (without the x, e.g. 10) (got %r)', sw.gauss_sigma)
    for key in ('merge_ab', 'mask_bad_quality'):
        if not isinstance(sw[key], bool):
            fail('[swath] %s 須為 true 或 false (目前為 %r)', '[swath] %s must be true or false (got %r)', key, sw[key])
    if not (isinstance(sw.merge, bool) or str(sw.merge).strip().lower() == 'ask'):
        fail('[swath] merge 須為 "ask"、true 或 false (目前為 %r)', '[swath] merge must be "ask", true or false (got %r)',
             sw.merge)
    if not (_is_num(sw.merge_gap) and sw.merge_gap >= 0):
        fail('[swath] merge_gap 須為 0 以上的秒數 (目前為 %r)', '[swath] merge_gap must be seconds >= 0 (got %r)',
             sw.merge_gap)
    if not (_is_num(sw.min_coverage) and 0 <= sw.min_coverage < 1):
        fail('[swath] min_coverage 須為 0~1 (目前為 %r)', '[swath] min_coverage must be 0~1 (got %r)', sw.min_coverage)

    if not _pair(cfg.gray.range):
        fail('[gray] range 須為 [最小, 最大] (K) (目前為 %r)', '[gray] range must be [min, max] (K) (got %r)',
             cfg.gray.range)
    for band in BANDS:
        if not _is_num(cfg.pct.get(band)):
            fail('[pct] %s 須為數值 (PCT 的 θ) (目前為 %r)', '[pct] %s must be a number (the PCT θ) (got %r)',
                 band, cfg.pct.get(band))
    for name, rec in cfg.rgb.items():
        if not isinstance(rec, dict):
            fail('[rgb.%s] 須為一節設定', '[rgb.%s] must be a table', name)
        free = 'band' not in rec                # 「任意頻道」配方：沒有 band，各色槍用 channel 指定頻率
        if not free and rec.get('band') not in BANDS:
            fail('[rgb.%s] band 須為 %s 之一 (目前為 %r)', '[rgb.%s] band must be one of %s (got %r)',
                 name, ' / '.join(BANDS), rec.get('band'))
        if 'theta' in rec and not _is_num(rec['theta']):
            fail('[rgb.%s] theta 須為數值 (目前為 %r)', '[rgb.%s] theta must be a number (got %r)', name, rec['theta'])
        for gun in ('red', 'green', 'blue'):
            g = rec.get(gun)
            if free:
                from .products import _SPEC
                if not (isinstance(g, dict) and _SPEC.match(str(g.get('channel', ''))) and _pair(g.get('range'))
                        and ('minus' not in g or _SPEC.match(str(g['minus'])))):
                    fail('[rgb.%s] %s 須為 { channel = "89" / "150" / "183±7" …, range = [最小, 最大] }，'
                         '可加 minus = "150" (兩個頻道相減) (目前為 %r)',
                         '[rgb.%s] %s must be { channel = "89" / "150" / "183±7" ..., range = [min, max] }, '
                         'optionally with minus = "150" (difference of two channels) (got %r)', name, gun, unwrap(g))
            elif not (isinstance(g, dict) and str(g.get('source', '')).lower() in ('pct', 'v', 'h')
                      and _pair(g.get('range'))):
                fail('[rgb.%s] %s 須為 { source = "pct" / "v" / "h", range = [最小, 最大] } (目前為 %r)',
                     '[rgb.%s] %s must be { source = "pct" / "v" / "h", range = [min, max] } (got %r)',
                     name, gun, unwrap(g))
            if 'gamma' in g and not (_is_num(g['gamma']) and g['gamma'] > 0):
                fail('[rgb.%s] %s 的 gamma 須為正數 (目前為 %r)', '[rgb.%s] %s: gamma must be a positive number (got %r)',
                     name, gun, g['gamma'])
    check_colormaps(cfg, CMAP_GROUPS, reserved=('GRAY',))
