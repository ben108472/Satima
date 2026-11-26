"""設定檔 (TOML) 讀取：共同設定 + 各類衛星的設定

檔案
  Satima.toml        共同設定 (語言、目錄、地圖、海岸線、圖片大小、標題、標註、檔名 …)
  Satima-GEO.toml    靜止衛星的設定 (真實色、底圖、紅外色表 …)，也可以覆寫共同設定的任何一節
  Satima-MW.toml     微波的設定 (內插、PCT、RGB、色表 …)
  Satima-PO.toml     極軌成像儀的設定 (颱風眼溫度、色表 …)
優先順序：程式預設值 < Satima.toml < Satima-<類別>.toml < [readers.<格式>] < 命令列
  * 每個設定只需要寫在一處；同一個鍵兩邊都有時，類別檔蓋過共同檔。設定檔可以只寫想改的部分。
  * [readers.<格式名稱>] 底下可以覆寫任何一節，只套用在該格式 (寫在共同檔或類別檔都可以)，例如
        [readers.AHI-HSD.general]
        data_dir = "D:/data/himawari"
搜尋順序：命令列 --config 指定共同設定檔 > 執行目錄的 Satima.toml > Satima.py 所在目錄的 Satima.toml；
類別檔先找共同設定檔所在的目錄，再找執行目錄與 Satima.py 所在目錄。
目錄 (data_dir、output_dir、colormap_dir、conversion_table_dir) 可寫絕對路徑或相對路徑；
相對路徑以 Satima.py 所在目錄為準，與從哪裡執行無關。
"""
import copy
import os
import re
import tomllib

from . import i18n
from .i18n import tr

COMMON_NAME = 'Satima.toml'
FAMILY_NAMES = {'geo': 'Satima-GEO.toml', 'mw': 'Satima-MW.toml', 'po': 'Satima-PO.toml'}
HOME = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))   # Satima.py 所在目錄
DIR_KEYS = ('data_dir', 'output_dir', 'colormap_dir', 'conversion_table_dir')
MODES = ('native', 'dpi')       # [image] mode：原始解析度 / 依 dpi 輸出
MAX_SCALE = 8                   # [image] size = "xN" 的倍數上限
BACKENDS = ('basemap', 'cartopy')
RESOLUTIONS = ('c', 'l', 'i', 'h', 'f', '110m', '50m', '10m')

# 共同的預設值；各類衛星自己的預設值 (含覆寫這裡的值) 在 satima/<類別>/config.py 的 DEFAULTS
DEFAULTS = {
    'general': {
        'language': 'zh',                    # 介面語言："zh" (繁體中文) / "en" (English)
        'data_dir': 'data',                  # 衛星資料目錄
        'output_dir': 'plot',                # 圖片輸出目錄 (可含檔名欄位)
        'colormap_dir': 'colormap',          # 色表目錄
        'conversion_table_dir': 'conversion_table',   # 校正表目錄
    },
    'map': {
        'backend': 'basemap',                # basemap 或 cartopy
        'resolution': 'i',                   # 海岸線解析度 c / l / i / h / f 或 110m / 50m / 10m (自動對應)
        'projection': 'merc',                # merc / cyl / lcc (GEO 另有 fd 全圓盤)；執行時直接 Enter 的預設
        'ask_projection': True,              # 執行時詢問投影；false = 一律用 projection
        'range_input': 'bounds',             # 地圖範圍的輸入方式："bounds" = 南 北 西 東；"center" = 中心 + 半寬
        'lcc_parallels': [],                 # LCC 標準緯線 [緯度1, 緯度2]；[] = 依地圖範圍自動
        'lcc_central_longitude': 'auto',     # LCC 中央經線 (度)；"auto" = 地圖範圍中央
    },
    'coastline': {'draw': True, 'color': '#FFFFFF', 'linewidth': 0.2},
    'borderline': {'draw': False, 'color': '#FFFFFF', 'linewidth': 0.3},
    'latlon': {'draw': False, 'color': '#FFFFFF', 'linewidth': 0.1, 'step': 1.0},
    'image': {
        'mode': 'native',                    # "native" = 原始解析度；"dpi" = 依 [[output.image]] 的 dpi 輸出
        # ---- mode = "native" ----
        'size': 'native',                    # "native" = 原始解析度、"x2"~"x8" = 其倍數；或固定寬度 (像素)
        'max_size': 8000,                    # 長邊的上限 (像素)
        'min_size': 2000,                    # 寬度的下限 (像素)
        'min_scale': 4,                      # 寬度不足 min_size 時至少放大的倍數
        'suffix': '',                        # 填入檔名的 {suffix}
        # ---- 兩種模式共用 ----
        'background': '#222222',             # 無資料處 (太空、掃描帶外) 的底色
    },
    'memory': {
        'low_memory': True,                  # 節省記憶體：輸出影像分條處理 (結果與整張處理相同)
        'workers': 'auto',                   # 多核心：同時處理的條數 ("auto" = CPU 邏輯核心數，1 = 不平行)
    },
    'title': {
        'ask': False,                        # 是否詢問圖片左下角的標題
        'default': '',
        'draw': True,                        # 標題畫在圖片左下角；false = 只填入檔名的 {title}
        'fontsize': 12, 'family': 'Arial', 'weight': 'normal',
        'family_nonascii': 'Source Han Sans CN', 'weight_nonascii': 'normal',
        'box_alpha': 0.5, 'box_pad': 4.0,
    },
    'annotation': {
        'draw': True,                        # 圖片下方中央的時間 / 衛星 / 產品文字
        'text': '{time}    {satellite} {sensor} {band} {colorbar}',
        'time_format': '%Y/%m/%d %H%MZ',
        'upper': True,
        'fontsize': 6, 'family': 'Courier New', 'weight': 'bold',
        'color': 'black', 'box_color': 'white', 'box_alpha': 0.7,
        'box_pad': 4.0,
    },
    'output': {
        'filename': '{title}.{satellite}.{sensor}.{band}.{colorbar}.{projection}.{time}.{suffix}',
        'time_format': '%Y%m%d%H%MZ',
        'format': 'png',
        'space': '_',
        # mode = "dpi"：每個 [[output.image]] 輸出一張 (dpi、suffix、annotation 可省略)
        'image': [{'dpi': 1000, 'suffix': ''}, {'dpi': 500, 'suffix': '500DPI'}],
    },
    'readers': {},
}

_FAMILIES = {}                  # 類別代號 -> 該類別的 config 模組 (DEFAULTS、PROJECTIONS、validate)


class Section(dict):
    """可以用屬性存取的 dict：cfg.coastline.color 等同 cfg['coastline']['color']。"""

    def __getattr__(self, key):
        try:
            return self[key]
        except KeyError:
            raise AttributeError(key) from None

    def __setattr__(self, key, value):
        self[key] = value


def wrap(obj):
    if isinstance(obj, dict):
        return Section({k: wrap(v) for k, v in obj.items()})
    if isinstance(obj, list):
        return [wrap(v) for v in obj]
    return obj


def unwrap(obj):
    if isinstance(obj, dict):
        return {k: unwrap(v) for k, v in obj.items()}
    if isinstance(obj, list):
        return [unwrap(v) for v in obj]
    return obj


def merge(base, over):
    """遞迴合併：over 的值覆蓋 base；list 整個取代。"""
    out = copy.deepcopy(unwrap(base))
    for key, value in over.items():
        if isinstance(value, dict) and isinstance(out.get(key), dict):
            out[key] = merge(out[key], value)
        else:
            out[key] = copy.deepcopy(value)
    return out


# ---- 讀檔 ------------------------------------------------------------------
def find_common(explicit=None):
    if explicit:
        if not os.path.isfile(explicit):
            raise SystemExit(tr('找不到設定檔：%s', 'Config file not found: %s') % explicit)
        return os.path.abspath(explicit)
    for folder in (os.getcwd(), HOME):
        path = os.path.join(folder, COMMON_NAME)
        if os.path.isfile(path):
            return path
    return None


def find_family(kind, common_path=None):
    folders = [os.path.dirname(common_path)] if common_path else []
    for folder in folders + [os.getcwd(), HOME]:
        path = os.path.join(folder, FAMILY_NAMES[kind])
        if os.path.isfile(path):
            return path
    return None


def read(path):
    """讀一個 TOML 檔 (utf-8-sig：記事本等編輯器存成帶 BOM 的 UTF-8 也能讀)。"""
    if not path:
        return {}
    with open(path, 'r', encoding='utf-8-sig') as f:
        try:
            return tomllib.loads(f.read())
        except tomllib.TOMLDecodeError as err:
            raise SystemExit(tr('設定檔格式錯誤 %s：%s', 'Config error in %s: %s') % (path, err))


def peek_language(explicit=None):
    """只取 [general] language (設定檔的錯誤訊息也要用對語言，所以在正式讀取之前先看)。找不到時 None。"""
    try:
        path = find_common(explicit)
        if path:
            with open(path, 'r', encoding='utf-8-sig') as f:
                return i18n.normalize(tomllib.loads(f.read()).get('general', {}).get('language', ''))
    except (SystemExit, OSError, tomllib.TOMLDecodeError, AttributeError):
        pass
    return None


def load(kind, family, explicit=None, overrides=None):
    """某一類衛星的設定 (路徑尚未解析；各格式實際使用的設定由 effective() 產生)。
    family：該類別的 config 模組，提供 DEFAULTS (類別自己的預設值)、PROJECTIONS (可用的投影)、validate(cfg)。
    overrides：命令列覆寫的設定 (例如 {'image': {'mode': 'dpi'}})，優先於設定檔與 [readers.*]。"""
    _FAMILIES[kind] = family
    common = find_common(explicit)
    own = find_family(kind, common)
    cfg = merge(DEFAULTS, family.DEFAULTS)
    for path in (common, own):
        cfg = merge(cfg, read(path))
    cfg = wrap(cfg)
    cfg['_kind'] = kind
    cfg['_paths'] = [p for p in (common, own) if p]
    cfg['_overrides'] = overrides or {}
    validate(wrap(merge(_public(cfg), cfg['_overrides'])), kind)
    return cfg


def _public(cfg):
    return {k: v for k, v in cfg.items() if not k.startswith('_')}


def effective(cfg, name, data_dir=None, output_dir=None):
    """某個格式實際使用的設定：預設值 < 設定檔 < [readers.<name>] < 命令列。"""
    kind = cfg['_kind']
    over = cfg.get('readers', {}).get(name) or {}
    data = {k: v for k, v in _public(cfg).items() if k != 'readers'}
    new = wrap(merge(merge(data, over), cfg['_overrides']))
    if data_dir:
        new.general.data_dir = data_dir
    if output_dir:
        new.general.output_dir = output_dir
    resolve_paths(new)
    validate(new, kind)
    new['_kind'] = kind
    new['_home'] = HOME
    return new


def resolve_paths(cfg):
    """所有目錄：絕對路徑照用；相對路徑以 Satima.py 所在目錄為準，留空 = 該目錄本身。
    output_dir 可含檔名欄位，存檔時才代入。
    (命令列 -d / -o 在 main() 已轉成絕對路徑：命令列上的相對路徑依一般慣例以目前目錄為準。)"""
    g = cfg.general
    for key in DIR_KEYS:
        g[key] = os.path.normpath(os.path.join(HOME, os.path.expanduser(str(g[key] or '').strip())))


# ---- 共用的解析與檢查 --------------------------------------------------------
def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def image_size(size):
    """[image] size -> (固定寬度或 None, 原始解析度的倍數)：
    "native" -> (None, 1)；"x4" (或 "4x") -> (None, 4)；正整數 2000 -> (2000, 1)。不合法時 ValueError。"""
    if _is_int(size):
        if size > 0:
            return size, 1
        raise ValueError(size)
    text = str(size).strip().lower()
    if text == 'native':
        return None, 1
    m = re.fullmatch(r'[x×*]\s*(\d+)|(\d+)\s*[x×*]', text)
    if m:
        k = int(m.group(1) or m.group(2))
        if 1 <= k <= MAX_SCALE:
            return None, k
    raise ValueError(size)


def fail(zh, en, *args):
    """設定檔錯誤：印出訊息後結束。"""
    raise SystemExit(tr('設定檔錯誤：', 'Config error: ') + tr(zh, en) % args)


def validate(cfg, kind):
    """共同設定的檢查，再交給該類別的 validate()。"""
    family = _FAMILIES[kind]
    g, m, img = cfg.general, cfg.map, cfg.image
    if i18n.normalize(g.language) is None:
        fail('[general] language 須為 "zh" 或 "en" (目前為 %r)', '[general] language must be "zh" or "en" (got %r)',
             g.language)
    if str(m.backend).lower() not in BACKENDS:
        fail('[map] backend 須為 "basemap" 或 "cartopy" (目前為 %r)',
             '[map] backend must be "basemap" or "cartopy" (got %r)', m.backend)
    if str(m.resolution).strip().lower() not in RESOLUTIONS:
        fail('[map] resolution 須為 %s 之一 (目前為 %r)', '[map] resolution must be one of %s (got %r)',
             ' / '.join(RESOLUTIONS), m.resolution)
    projections = family.PROJECTIONS
    if str(m.projection).strip().lower() not in projections:
        fail('[map] projection 須為 %s 之一 (目前為 %r)', '[map] projection must be one of %s (got %r)',
             ' / '.join(projections), m.projection)
    if str(m.range_input).strip().lower() not in ('bounds', 'center'):
        fail('[map] range_input 須為 "bounds" (南 北 西 東) 或 "center" (中心 + 範圍) (目前為 %r)',
             '[map] range_input must be "bounds" (south north west east) or "center" (center + half-width) (got %r)',
             m.range_input)
    par = m.lcc_parallels
    if not (isinstance(par, list) and len(par) <= 2 and all(_is_num(v) and -90 < v < 90 for v in par)):
        fail('[map] lcc_parallels 須為 [] 或 1~2 個緯度，例如 [10, 40] (目前為 %r)',
             '[map] lcc_parallels must be [] or 1~2 latitudes, e.g. [10, 40] (got %r)', par)
    lon0 = m.lcc_central_longitude
    if not (_is_num(lon0) or str(lon0).strip().lower() == 'auto'):
        fail('[map] lcc_central_longitude 須為 "auto" 或經度 (目前為 %r)',
             '[map] lcc_central_longitude must be "auto" or a longitude (got %r)', lon0)

    if not isinstance(cfg.title.draw, bool):
        fail('[title] draw 須為 true 或 false (目前為 %r)', '[title] draw must be true or false (got %r)', cfg.title.draw)

    mode = str(img.mode).strip().lower()
    if mode not in MODES:
        fail('[image] mode 須為 "native" 或 "dpi" (目前為 %r)', '[image] mode must be "native" or "dpi" (got %r)',
             img.mode)
    img.mode = mode
    try:
        image_size(img.size)
    except ValueError:
        fail('[image] size 須為 "native"、"x2"~"x%d" 或正整數寬度 (目前為 %r)',
             '[image] size must be "native", "x2"~"x%d" or a positive integer width (got %r)', MAX_SCALE, img.size)
    if not (_is_int(img.min_size) and _is_int(img.max_size) and 0 < img.min_size <= img.max_size):
        fail('[image] 須滿足 0 < min_size <= max_size 且為整數 (目前 %r、%r)',
             '[image] needs integers 0 < min_size <= max_size (got %r, %r)', img.min_size, img.max_size)
    if not (_is_int(img.min_scale) and 1 <= img.min_scale <= MAX_SCALE):
        fail('[image] min_scale 須為 1~%d 的整數 (目前為 %r)', '[image] min_scale must be an integer 1~%d (got %r)',
             MAX_SCALE, img.min_scale)
    images = cfg.output.get('image')
    if not (isinstance(images, list) and images):
        fail('至少需要一個 [[output.image]]', 'at least one [[output.image]] is required')
    for item in images:
        dpi = item.get('dpi') if isinstance(item, dict) else None
        if not (_is_num(dpi) and 10 <= dpi <= 5000):
            fail('每個 [[output.image]] 都要有 dpi (10~5000) (目前為 %r)',
                 'every [[output.image]] needs a dpi (10~5000) (got %r)', unwrap(item))
        if 'annotation' in item and not isinstance(item['annotation'], bool):
            fail('[[output.image]] 的 annotation 須為 true 或 false', '[[output.image]] annotation must be true or false')

    if not isinstance(cfg.memory.low_memory, bool):
        fail('[memory] low_memory 須為 true 或 false (目前為 %r)', '[memory] low_memory must be true or false (got %r)',
             cfg.memory.low_memory)
    w = cfg.memory.workers
    if not ((isinstance(w, str) and w.strip().lower() == 'auto') or (_is_int(w) and w >= 1)):
        fail('[memory] workers 須為 "auto" 或正整數 (目前為 %r)',
             '[memory] workers must be "auto" or a positive integer (got %r)', w)
    family.validate(cfg)


def check_colormaps(cfg, groups, reserved=()):
    """[colormap] 各分組的項目都要有 file 或 mpl；名稱不可用 reserved (選單保留的名稱)。"""
    for kind in groups:
        items = cfg.colormap.get(kind)
        if not isinstance(items, list):
            fail('[colormap] 缺少 %s 列表', '[colormap] is missing the %s list', kind)
        for item in items:
            if not (item.get('file') or item.get('mpl')):
                fail('colormap.%s 的 %r 需要 "file" 或 "mpl"', 'colormap.%s entry %r needs "file" or "mpl"',
                     kind, unwrap(item))
            name = str(item.get('name') or item.get('file') or item.get('mpl')).upper()
            if name in reserved:
                fail('colormap.%s 的名稱不可為 %s (選單保留的名稱)', 'colormap.%s: the name %s is reserved by the menu',
                     kind, name)
