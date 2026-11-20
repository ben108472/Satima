"""Satima 的靜止衛星 (GEO) 部分

  config.py     GEO 自己的設定 (預設值與檢查；共同設定見 satima/core/config.py)
  family.py     流程：搜尋場景、資料列表、投影與範圍、產品、繪圖存檔
  nav.py        導航介面：GeosNav (CGMS 標準靜止衛星投影)、ProjGridNav (規則投影格點)、LonLatNav (經緯度陣列)
  reader.py     各衛星讀取器的基底類別 GeoReader 與 Channel / Scene
  resample.py   輸出格點與取樣 (經緯度 -> 原始行列號，直接取值)
  canvas.py     畫布 (共用畫布 + GEO 的輸出格點)
  products.py   產品：單波段、灰階底圖、色表強化、真實色
  landsea.py    灰階底圖：CWA 海陸、LIGHTS 夜光
  bluemarble.py 灰階底圖：NASA / NICT Blue Marble
  truecolor/    CIRA 風格真實色與 GeoColor
互動輸入、地圖 (含全圓盤 FD)、色表、檔名、存檔是三類衛星共用的 (satima/core)。
各衛星讀取器放在 sensors/geo/ (檔名形如 AHI-HSD.py)，新增衛星請見 docs/GEO.md。
"""
