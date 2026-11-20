"""Satima — 氣象衛星影像繪圖 (靜止衛星 GEO、微波 MW、極軌成像儀 PO 的整合版)

  core/   共用：設定檔、互動輸入、語言、檔名、PNG 輸出、主流程
  geo/    靜止衛星 (解析式導航)
  mw/     微波輻射計 (掃描帶)
  po/     極軌成像儀 (掃描帶)
各衛星讀取器放在 sensors/geo、sensors/mw、sensors/po (檔名形如 AHI-HSD.py)。
"""
