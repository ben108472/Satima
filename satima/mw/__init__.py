"""Satima 的微波 (MW) 部分：被動微波輻射計的掃描帶影像

  config.py     MW 自己的設定 (預設值與檢查；共同設定見 satima/core/config.py)
  family.py     流程：搜尋場景、資料列表、多檔合併、範圍、產品、繪圖存檔
  reader.py     各儀器讀取器的基底類別 SwathReader 與 Channel / Scene / Geo
  swath.py      掃描帶取樣：輸出像素 -> 掃描帶小數座標 -> 高斯距離加權或索引空間內插 (含 A / B 掃描合併)
  grid.py       輸出格點
  canvas.py     畫布 (共用畫布 + MW 的輸出格點)
  colormaps.py  亮溫 / PCT 色表的分組
  products.py   產品：單頻道亮溫、PCT、RGB 合成
互動輸入、地圖、檔名、存檔是三類衛星共用的 (satima/core)。
各儀器讀取器放在 sensors/mw/ (檔名形如 PPS-L1C.py)，新增格式請見 docs/MW.md。
"""
