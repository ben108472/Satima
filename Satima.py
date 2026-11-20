"""
Satima — 氣象衛星影像繪圖主程式 (靜止衛星 GEO、微波 MW、極軌成像儀 PO)

用法
  python Satima.py                  在設定檔的 data_dir (預設為 Satima.py 旁的 data) 搜尋資料
  python Satima.py -d D:\\data       指定資料目錄
  python Satima.py -o D:\\png        指定圖片輸出目錄
  python Satima.py -t geo           只列某一類：geo / mw / po
  python Satima.py -f AHI-HSD       只搜尋某個格式
  python Satima.py -c my.toml       指定共同設定檔
  python Satima.py --lang en        介面語言 zh / en
  python Satima.py --mode dpi       圖片模式 native (原始解析度) / dpi
  python Satima.py --list           列出可用的格式

流程：選資料 -> 投影與地圖範圍 -> 產品與色表 -> 圖片寬度 (原始解析度模式) -> 標題 -> 繪圖存檔
設定檔：Satima.toml (共同) 與 Satima-GEO.toml / Satima-MW.toml / Satima-PO.toml (各類)。
各衛星讀取器放在 sensors/geo、sensors/mw、sensors/po (檔名形如 AHI-HSD.py)，新增衛星請見 README.md。
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from satima.core.session import main

if __name__ == '__main__':
    main()
