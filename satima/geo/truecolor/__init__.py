"""CIRA 風格真實色 (SHAC True Color) 與可見光反射率校正工具

    from satima.truecolor import make_truecolor          # 完整真實色流程 (整張)
    from satima.truecolor import TrueColor               # 可分條計算 (節省記憶體模式)
    from satima.truecolor import correction, astronomy   # 個別校正 / 幾何計算
"""
from . import astronomy, correction
from .composite import TrueColor, make_truecolor

__all__ = ['astronomy', 'correction', 'TrueColor', 'make_truecolor']
