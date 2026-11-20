"""色表 (GEO / MW / PO 共用)

色表選單定義在設定檔 [colormap] 的各個列表，順序即選單編號 (0 = 全部)：
  GEO、PO：ir / wv (單位 °C)；MW：tb89、tb37、tb150、tb、pct89、pct37、pct (單位 K)
  file : colormap 目錄中的 .cmap 檔名 (不含副檔名)，以 mpkit.gpf 讀取；
         顯示範圍預設取色表檔 levels 的最小/最大值
  mpl  : matplotlib 內建色表名稱 (例如 gray_r)，需搭配 vmin / vmax
"""
import os

import numpy as np

from .i18n import tr

_CACHE = {}
_CHUNK = 1 << 21            # 上色時每次處理的像素數 (控制暫存記憶體)


class CmapEntry:

    def __init__(self, item, cmap_dir, default=(-80., 50.)):
        self.default = default      # mpl 色表沒寫 vmin / vmax 時的範圍
        self.file = str(item.get('file', '') or '')
        self.mpl = str(item.get('mpl', '') or '')
        self.name = str(item.get('name') or self.file or self.mpl)
        self.vmin = item.get('vmin')
        self.vmax = item.get('vmax')
        self.cmap_dir = cmap_dir

    def __repr__(self):
        return 'CmapEntry(%s)' % self.name

    def load(self):
        """回傳 (cmap, vmin, vmax)。"""
        if self.file:
            path = os.path.join(self.cmap_dir, self.file)
            if path not in _CACHE:
                if not os.path.isfile(path + '.cmap'):
                    raise FileNotFoundError(tr('找不到色表檔：%s.cmap', 'Colormap file not found: %s.cmap') % path)
                import mpkit.gpf as gpf
                # gpf 以 os.path.join(內建目錄, 名稱) 找檔，傳入絕對路徑即直接使用該檔
                _CACHE[path] = gpf.cmap(path)
            d = _CACHE[path]
            levels = d['levels']
            vmin = float(np.min(levels)) if self.vmin is None else float(self.vmin)
            vmax = float(np.max(levels)) if self.vmax is None else float(self.vmax)
            return d['cmap'], vmin, vmax
        vmin = float(self.default[0]) if self.vmin is None else float(self.vmin)
        vmax = float(self.default[1]) if self.vmax is None else float(self.vmax)
        return self.mpl, vmin, vmax


def entries(cfg, group, default=(-80., 50.)):
    """group：[colormap] 的分組名稱 (ir、wv、tb89、pct37 ...)。"""
    return [CmapEntry(item, cfg.general.colormap_dir, default) for item in cfg.colormap[group]]


def ask_colormaps(cfg, group, prompt=None):
    """印出色表選單 (0.ALL 1.BW 2.BD ...)，回傳所選 CmapEntry 列表。"""
    from . import ui
    ents = entries(cfg, group)
    options = {i + 1: e.name for i, e in enumerate(ents)}
    keys = ui.ask_menu('', options, prompt or tr('選擇色表 (可複選，0 = 全部)：',
                                                 'Select colormap (several allowed, 0 = all): '))
    return [ents[k - 1] for k in keys]


def colorize(data, cmap, vmin, vmax, out=None):
    """單波段 -> 8 位元 RGBA：色階 = floor((值 - vmin) / (vmax - vmin) x N) (matplotlib 的定義，以 float64 計算)。
    NaN 為色表的 bad 顏色 (預設透明，顯示底色)；超出範圍為色表兩端 (或 under / over 顏色)。
    matplotlib 自己上色時先以 float32 縮放再還原，極細的色表 (例如 BD 約 5 萬階) 偶有像素差一階。"""
    import matplotlib
    cm = matplotlib.colormaps[cmap] if isinstance(cmap, str) else cmap
    data = np.asarray(data)
    if out is None:
        out = np.empty(data.shape + (4,), dtype=np.uint8)
    src, dst = data.reshape(-1), out.reshape(-1, 4)
    for s in range(0, src.size, _CHUNK):
        x = (src[s:s + _CHUNK].astype(np.float64) - float(vmin)) / (float(vmax) - float(vmin))
        dst[s:s + _CHUNK] = cm(np.ma.masked_invalid(x), bytes=True)
    return out


def rgb8(rgb):
    """0~1 的浮點 RGB -> 8 位元 (與 matplotlib 顯示浮點 RGB 的換算相同：截斷而非四捨五入)。"""
    return (np.clip(rgb, 0., 1.) * 255).astype(np.uint8)
