"""亮溫 / PCT 色表 (單位 K)：共用的色表模組 (satima/core/colormaps.py) 加上 MW 的分組

色表選單定義在設定檔 [colormap] 的各個列表 (tb89、tb37、tb150、tb、pct89、pct37、pct)，順序即選單編號。
"""
from ..core.colormaps import CmapEntry, colorize, rgb8      # noqa: F401
from ..core.colormaps import entries as _entries


def group_of(kind, pct=False):
    """頻率類別 (reader.kind_of) -> [colormap] 的分組。"""
    if pct:
        return {'tb89': 'pct89', 'tb37': 'pct37'}.get(kind, 'pct')
    return {'tb89': 'tb89', 'tb37': 'tb37', 'tb150': 'tb150', 'tb183': 'tb150'}.get(kind, 'tb')


def entries(cfg, group):
    """group：[colormap] 的分組名稱 (tb89、pct37 ...)。mpl 色表沒寫範圍時為 150~300 K。"""
    return _entries(cfg, group, (150., 300.))
