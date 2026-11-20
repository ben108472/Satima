"""MW 的畫布：共用畫布 (satima/core/canvas.py) 加上 MW 的輸出格點 (grid.py；掃描帶取樣見 swath.py)。"""
from ..core.canvas import Canvas as _Canvas
from .grid import OutputGrid


class Canvas(_Canvas):
    GRID = OutputGrid
