"""GEO 的畫布：共用畫布 (satima/core/canvas.py) 加上 GEO 的輸出格點 (解析式導航取樣，resample.py)。"""
from ..core.canvas import Canvas as _Canvas
from .resample import OutputGrid


class Canvas(_Canvas):
    GRID = OutputGrid
