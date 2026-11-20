"""Satima 的極軌成像儀 (PO) 部分

  config     PO 自己的設定 (預設值與檢查；共同設定見 satima/core/config.py)
  reader     格式處理腳本的基底類別 SatReader 與 Segment (一段連續的觀測)
  family     流程：搜尋檔案並接成段、資料列表、範圍、頻道與色表、繪圖存檔
  canvas     畫布：掃描帶取樣到共用地圖的格點 (最近鄰)，再交給共用畫布畫線、標註、存檔
  swath      經度展開、範圍裁切、座標補值、颱風眼溫度
  truecolor  真實色共用演算法 (CREFL Rayleigh 校正、太陽天頂角校正、增強曲線)
  obstime    觀測時間的換算
  archive    Sentinel-3 SAFE 產品讀取 (不解壓縮)
  tiepoint   規則 tie-point 網格內插 (經緯度於地心直角座標內插)
各格式的處理腳本放在 sensors/po/ (檔名形如 VIIRS-HDF5.py)。
互動輸入 (ui)、色表 (colormaps)、檔名 (naming) 是三類衛星共用的 (satima/core)，這裡轉出同名的模組給處理腳本使用。
"""
from ..core import colormaps, naming, ui                          # noqa: F401
from . import archive, obstime, swath, tiepoint, truecolor      # noqa: F401
from .reader import SatReader, Segment                          # noqa: F401

__all__ = ['SatReader', 'Segment', 'ui', 'swath', 'truecolor', 'colormaps', 'naming', 'archive', 'tiepoint', 'obstime']
