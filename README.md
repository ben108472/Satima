# Satima — 氣象衛星影像繪圖

Satima-GEO (靜止衛星)、Satima-MW (微波)、Satima-PO (極軌成像儀) 的整合版：一個主程式、一套操作介面、
分層的設定檔。English version: [README.en.md](README.en.md)

| 類別 | 資料 | 讀取器 | 詳細說明 |
|---|---|---|---|
| **GEO** 靜止衛星 | Himawari、GOES、GK-2A、MTG、MTSAT、COMS、GMS、Meteosat… | `sensors/geo/` (18 種格式) | [docs/GEO.md](docs/GEO.md) |
| **MW** 微波輻射計 | AMSR2/3、GMI、SSMIS、SSM/I、ATMS、AMSU-A、MWRI、NASA PPS L1C… | `sensors/mw/` (12 種格式) | [docs/MW.md](docs/MW.md) |
| **PO** 極軌成像儀 | VIIRS、MODIS、AVHRR、METimage、OLCI、SLSTR、SGLI | `sensors/po/` (8 種格式) | [docs/PO.md](docs/PO.md) |

## 執行

```
python Satima.py                  在設定檔的 data_dir (預設為 Satima.py 旁的 data) 搜尋資料
python Satima.py -d D:\data       指定資料目錄
python Satima.py -o D:\png        指定圖片輸出目錄
python Satima.py -t geo           只列某一類：geo / mw / po
python Satima.py -f AHI-HSD       只搜尋某個格式
python Satima.py -c my.toml       指定共同設定檔 (類別檔在同一個資料夾找)
python Satima.py --lang en        介面語言 zh / en
python Satima.py --mode dpi       圖片模式 native (原始解析度) / dpi
python Satima.py --list           列出可用的格式
```

流程三類相同：**選資料 → 投影 → 地圖範圍 → 產品 → 色表 / 樣式 → 圖片寬度 → 標題 → 繪圖存檔**。
只有選單的內容依類別不同 (全圓盤投影只有 GEO 有，PCT 與 RGB 只有 MW 有…)。

### 資料列表

三類的資料可以放在同一個資料夾。列表依類別分段、編號連續；只有一類有資料時不印標題：

```
[GEO]
  1   MTSAT-1R     IMAGER   HRIT   2013/11/07 11:30Z   FD DK01   5 bands
  2   Himawari-9   AHI      HSD    2026/10/20 21:00Z   FLDK      1 bands
[MW]
  3   F16       SSMIS   SDR   2017/09/06 07:25-09:20Z   CLASS   19 ch
  4   GCOM-W1   AMSR2   L1B   2025/09/22 04:34-05:24Z   JAXA    14 ch
  5   NOAA-20   ATMS    SDR   2025/09/28 05:50-06:01Z   CLASS   22 ch   21 檔
[PO]
  6   S-NPP     VIIRS    HDF5   2023/05/26 04:03-04:10Z   NOAA       2 bands   4 檔 (缺 1)
  7   METOP-B   AVHRR3   NAT    2026/04/12 09:31-11:13Z   EUMETSAT   IR/VIS
輸入資料編號 (多筆可用 (3,5) 或 [1,3,5])：
```

* 各欄：衛星、儀器、格式、時間、觀測範圍 (GEO) 或資料來源 (MW、PO)、波段 / 頻道數、**備註**。
* 備註欄：多個檔接成一筆時顯示檔數，中間有缺檔才加 `(缺 N)`；GEO 的區域觀測顯示影像的經緯度範圍。
* `3` = 單一筆；`(3,5)` = 第 3 到第 5 筆；`[1,3,5]` = 指定多筆。一次只能選同一種格式，選到不同格式會請你重選。
* GEO 選多筆 = 多個時次各畫一張；MW、PO 選多筆時，時間相接的可以接成一張 (見「極軌的接續判定」)。

### 地圖範圍

`[map] range_input` 決定先問哪一種，但兩種輸入法隨時都能用：

* 輸入 **4 個數字** = 南 北 西 東，例如 `20 27 117 124` (跨換日線用 0\~360，例如 `10 30 170 200`)
* 輸入 **2 個數字** = 中心緯度 經度，接著問範圍 (中心往外的緯度、經度各幾度；只輸入一個 = 兩者相同)

極軌 (PO) 的中心時間與颱風眼溫度以中心點為準：輸入四邊時為範圍中點，想指定颱風中心就輸入 2 個數字。
GEO 的區域觀測 (目標區、mesoscale) 不詢問，直接用整個影像；全圓盤 (FD) 不需要範圍。

投影：MERC (麥卡托)、CYL (等距圓柱)、LCC (蘭伯特正形圓錐)，GEO 另有 FD (全圓盤)。
麥卡托與 LCC 的地球模型一律是 WGS84 橢球，所以**同一範圍、同一大小的 GEO、MW、PO 圖逐像素對齊**，可以直接疊合。

## 設定檔

| 檔案 | 內容 |
|---|---|
| `Satima.toml` | 共同設定：語言、目錄、地圖、海岸線、圖片模式與大小、標題、標註、檔名 |
| `Satima-GEO.toml` | 靜止衛星：真實色、GeoColor、底圖、紅外色表… |
| `Satima-MW.toml` | 微波：內插、PCT、RGB 配方、色表… |
| `Satima-PO.toml` | 極軌成像儀：颱風眼溫度、色表… |

* 優先順序：程式內建預設值 < `Satima.toml` < `Satima-<類別>.toml` < `[readers.<格式>.<節>]` < 命令列。
* **每個設定只需要寫在一處**。類別檔可以覆寫共同檔的任何一節 (只寫想改的鍵)；同一個鍵兩邊都有時，類別檔為準。
  例如共同檔的海岸線是白色 0.2，`Satima-PO.toml` 寫 `[coastline] color = "#00BB00"` 就只有極軌用綠色。
* `[readers.AHI-HSD.general] data_dir = "D:/data/himawari"` 這種寫法可以再針對單一格式覆寫，寫在共同檔或類別檔都可以。
* 目錄可寫絕對或相對路徑；**相對路徑一律以 `Satima.py` 所在的資料夾為準**，與從哪裡執行無關
  (命令列 `-d` / `-o` 的相對路徑則依一般慣例以目前目錄為準)。
* 各鍵的說明都寫在設定檔的註解裡。英文註解的版本在 `config-en/` (鍵與值相同)；**程式只讀 `Satima.py` 旁邊的那四個檔**，
  `config-en/` 不會被讀到，想用英文註解的版本就把它複製出來蓋掉原本的。`python tests/check_config_sync.py` 可以檢查兩份的鍵與值是否一致。

### 語言

`[general] language = "zh"` (繁體中文) 或 `"en"` (English)，命令列 `--lang` 可臨時覆寫，一次只顯示一種。
只影響畫面上的訊息；產品與波段名稱、色表名稱、檔名、圖上的標註兩種語言都相同，所以兩種語言產生的檔案完全一樣。

### 圖片模式 (`[image] mode`)

三類衛星都適用，命令列 `--mode` 可臨時覆寫，也可以在類別檔或 `[readers.<格式>.image]` 個別指定。

| | `mode = "native"` (原始解析度) | `mode = "dpi"` |
|---|---|---|
| 圖片大小 | 地圖中央一個像素 = 這個產品的原始解析度 (極軌、微波為原始取樣間距) | 像素數 = 英吋 x dpi |
| 張數 | 一張 | 每個 `[[output.image]]` 一張 (預設 1000 dpi 與 500 dpi) |
| 執行時 | 詢問圖片寬度：Enter = 原始解析度、`x3` = 3 倍、或直接輸入像素數 | 不詢問 |
| 相關設定 | `[image] size / max_size / min_size / min_scale / suffix` | `[[output.image]] dpi / suffix / annotation` |

* dpi 模式的英吋數照 matplotlib 預設 (6.4 x 4.8 英吋的圖、座標軸約 4.96 x 3.70 英吋)：地圖等比例放進座標軸的範圍後裁掉白邊。
  寬的範圍 1000 dpi 約 4960 像素寬，方形範圍約 3700 像素高。預設值讀自 matplotlib 的設定，改 `matplotlibrc` 會跟著變。
* 字級與線寬的單位是點：原始解析度模式以 5 吋寬的圖為準，dpi 模式以實際的英吋數為準，所以不論輸出多少像素，
  字與線占圖片的比例都一樣。
* **兩種模式的影像都是先取樣到輸出地圖的格點，不經 matplotlib 重新取樣**，也不再使用 pcolormesh
  (PO 原本的 dpi 模式是 pcolormesh)。GEO 用解析式導航取樣，MW 用掃描帶內插 (`[swath] interpolation`)，
  PO 取最近的有效取樣點。想看到「一個取樣點一塊」的樣子，GEO 設 `[image] interpolation = "nearest"`、
  MW 設 `[swath] interpolation = "nearest"`。
* dpi 模式輸出多張時用檔名的 `{suffix}` 或 `{dpi}` 區分；樣板裡沒有這兩個欄位時，有 suffix 的圖會自動在檔名最後加上 `.suffix`。
* GEO、MW 每個 dpi 各取樣一次，輸出兩張的時間大約是一張的兩倍。

## 極軌的接續判定

極軌的檔案長短差很多：ATMS 一檔 32 秒、METimage 約 1 分鐘、VIIRS 約 85 秒、MODIS 5 分鐘、JAXA L1B 半圈、
SSMIS 整圈。程式把「相鄰」和「同一段」分成兩件事 (`satima/core/segments.py`)：

* **相鄰** (用來數缺檔)：下一檔起始 − 前一檔結束 ≤ 相鄰門檻。門檻跟著片段長度走：該組檔案長度中位數的一半，
  限制在 2\~120 秒 (ATMS 16 秒、VIIRS 約 43 秒、5 分鐘以上的檔 120 秒)。缺一個檔的間隔等於一檔的長度，一定超過門檻。
* **同一段** (資料列表的一列)：同衛星、同儀器、同格式，間隔不超過 10 分鐘就接成一段，中間缺檔照樣接，
  列表的備註欄標出 `(缺 N)`。同一地點的兩次通過至少隔一圈 (約 100 分鐘)，10 分鐘內不可能是另一次通過；
  缺口處取樣時會留空，不會被內插抹平。
* 軌道號不用來切段：軌道號在升交點過赤道時跳號，會剛好把低緯度颱風的通過切成兩半。
* **短片段** (一檔不到 15 分鐘：VIIRS、ATMS、METimage、MODIS…) 自動接成一筆，不詢問。
  **整圈 / 半圈的長檔**各自一列；一次選到時間相接的多列時依 `[swath] merge` 詢問要不要接成一張 (`"ask"` / `true` / `false`)。
* 同一個片段有重複的檔 (重新下載、不同產製時間) 時只用最新的。
* VIIRS：範圍確定後只讀與地圖相交的片段，前後各多帶一個。
* 圖上與檔名的時間是掃描帶**通過地圖範圍**的時間，不是檔案的起始時間。

目前的涵蓋範圍：ATMS (MW) 與全部 PO 讀取器使用這套判定。其他 MW 讀取器 (一個檔 = 整圈或半圈) 維持原本的
`[swath] merge` / `merge_gap`：選到時間相接的多個檔時詢問是否合併。

### VIIRS、ATMS 的檔案

* **NOAA CLASS 的合併檔** (定位與資料包在同一個檔，檔名前段列出所含的產品，例如 `GMTCO-SVM03-SVM04-SVM05_j01_…`、
  `GATMO-SATMS-TATMS_j01_…`) 可以讀，分開的檔與合併檔可以混用；幾個片段併成一個檔的聚合檔照檔名的起訖時間當成一個較長的片段。
* **VIIRS 定位檔**：優先用地形校正版 `GMTCO` / `GITCO`，找不到時改用未做地形校正的 `GMODO` / `GIMGO`，畫面上會註明。
  海面上兩者沒有差別；高山在掃描帶邊緣會偏移 (約為地形高度 × tan(視角)，3 km 的山在 60° 視角約 5 km)。
* **VIIRS DNB**：`GDNBO` 有地形校正的經緯度 (`Latitude_TC` / `Longitude_TC`) 時優先使用。
* **ATMS 定位檔**：`GATMO` 有各頻段的波束經緯度 (`BeamLatitude`) 就用；舊格式沒有時改用單一組的 `Latitude` / `Longitude`。

> **VIIRS 已用真實檔案驗證的部分** (S-NPP、NOAA-20、NOAA-21 的 I1、I5、DNB、M3 / M4 / M5，2026/10)：I5 色表、I1 反射率、
> DNB (ERF)、M5、TrueColor-M、TrueColor-I 的中心時間與颱風眼溫度和原本的 Satima-PO 數值相同，圖目視一致；片段自動接續、只讀相交的片段、
> 一次選兩顆衛星各一段、缺定位檔時的提示、DNB 的 `Latitude_TC` 都以真實檔案跑過。合併檔以真實資料組出來的檔
> (`tests/pack_viirs.py`) 測試，輸出與分開的檔逐像素相同。
>
> **尚未以真實檔案驗證的部分**：CLASS 下載的原始合併檔與聚合檔、VIIRS 的 `GMODO` / `GIMGO`、M 頻道的亮溫 (M12\~M16)、
> ATMS 的合併檔與舊格式定位檔 (其中 VIIRS 的幾項只用合成的測試檔 `tests/fake_viirs.py` 確認流程走得通)；
> MODIS、METimage、OLCI、SLSTR、SGLI、AVHRR-NC 沒有樣本，只確認程式可以載入。
> 這些用真實資料跑第一次時請留意，有問題的話把錯誤訊息與檔名留下來。

## 與原本三個程式的差異

**三類共同**
* 一個主程式 `Satima.py`；設定檔分成共同 + 三個類別檔；訊息可切換中文 / 英文。
* `[image] mode` (native / dpi) 三類都能用；dpi 模式的大小規則三類相同。
* 資料列表分段顯示，多一欄備註。
* 檔名樣板的預設值多了 `{suffix}` (空白時省略，原始解析度模式下檔名與以前相同)。

**GEO、MW**：原始解析度模式下的輸出與原程式逐像素相同 (見「驗證」)。GEO 多了中心 + 範圍的輸入法。

**PO** (改動最大)
* 拿掉 pcolormesh：兩種模式都是取樣到地圖格點 (原本的「原始解析度模式」的做法)。
* 麥卡托由圓球改為 WGS84 橢球：同一範圍的圖高會差一點 (例：6918 列變 6896 列)，與以前的圖不會逐像素相同。
* 多了 CYL、LCC 投影；範圍可以輸入四邊；訊息有中文。
* 相對路徑改以 `Satima.py` 所在的資料夾為準 (原本是執行時的目前目錄)，圖預設存到 `plot/` (原本存在資料旁邊)。
* 標題 (Picname) 除了填入檔名，也會畫在圖片左下角 (同 GEO、MW)；圖下方的標註改用 `[annotation] text` 樣板。
* 列表由「一個檔一列」改為「一段連續的觀測一列」；VIIRS 不再是一個波段一列，選了資料之後才從選單選頻道或真實色。
* 設定檔的鍵：`[output] facecolor` → `[image] background`、`[output] ask_title` / `title` → `[title] ask` / `default`、
  `[coastline] resolution` → `[map] resolution` (舊的鍵會提示改寫)。

## 目錄

```
Satima.py                 主程式
Satima.toml               共同設定檔
Satima-GEO.toml  Satima-MW.toml  Satima-PO.toml     各類的設定檔
config-en/                英文註解的設定檔範本 (程式不讀)
satima/
  core/                   三類共用
    session.py            命令列、讀取器載入、分段的資料列表、選資料
    config.py             分層設定檔的讀取、共同的預設值與檢查
    i18n.py               介面語言 tr(中文, English)
    ui.py                 互動輸入、選單、列表
    maps.py               地圖投影 (MERC / CYL / LCC / FD) 與底圖 (basemap / cartopy)
    canvas.py             圖的大小 (native / dpi)、海岸線、標註、存檔
    naming.py             輸出檔名樣板
    colormaps.py          色表
    segments.py           極軌的接續判定
    pngio.py  parallel.py 多核心 PNG 輸出
  geo/                    靜止衛星：reader、nav、resample、products、truecolor/、landsea、bluemarble
  mw/                     微波：reader、swath (內插)、grid、products
  po/                     極軌成像儀：reader、canvas (最近鄰取樣)、swath、truecolor、obstime、tiepoint、archive
  <類別>/family.py         該類的流程 (搜尋資料、列表的各列、範圍、產品、繪圖)
  <類別>/config.py         該類自己的預設值與檢查
sensors/geo/  sensors/mw/  sensors/po/    各格式的讀取器 (檔名含 '-' 的 .py 會自動載入)
colormap/                 色表 (.cmap)
conversion_table/         校正表
bluemarble/ landsea/ nightlights/   GEO 的底圖資料
data/  plot/              預設的資料與輸出目錄
docs/                     各類的詳細說明
tests/                    回歸測試與檢查工具
```

需要的套件見 `requirements.txt`。Python 3.11 以上 (使用內建的 tomllib)。缺少某個套件時只會略過用到它的格式。

### 新增讀取器

在 `sensors/<類別>/` 放一支 `<儀器>-<格式>.py`，最後一行 `READER = 讀取器類別`。各類的基底類別與寫法：

| 類別 | 基底類別 | 說明 |
|---|---|---|
| GEO | `satima.geo.reader.GeoReader` | [docs/GEO.md](docs/GEO.md)「新增衛星」 |
| MW | `satima.mw.reader.SwathReader` | [docs/MW.md](docs/MW.md)「新增格式」 |
| PO | `satima.po.SatReader` | [docs/PO.md](docs/PO.md)「新增衛星」 |

訊息請寫成 `tr('中文', 'English')` (`from satima.core.i18n import tr`)；讀取器的說明用 `DESCRIPTION` (中文) 與 `DESCRIPTION_EN`。

## 驗證 (2026/11)

`tests/regress.py` 以相同的輸入分別執行原本的三個程式與整合版，逐像素比對輸出：

```
python tests/regress.py run orig      執行原版 (需要原本的三個專案與樣本資料；路徑寫在腳本開頭)
python tests/regress.py run new       執行整合版
python tests/regress.py cmp orig new  比對
```

| 個案 | 內容 | 結果 |
|---|---|---|
| geo_ahi_merc | Himawari-9 AHI 麥卡托：BAND-13 色表、TrueColor | 2 張與原版逐像素相同 |
| geo_ahi_fd | AHI 全圓盤：BAND-13 灰階、CWA 海陸底圖 | 2 張相同 |
| geo_mtsat_lcc | MTSAT-1R LCC：VIS、IR1 色表 | 2 張相同 |
| mw_amsr2 | AMSR2 L1B：89V 色表、PCT89、NRL-37 | 3 張相同 |
| mw_atms | ATMS SDR 21 個片段自動接成一筆：88.2 GHz | 1 張相同 |
| mw_ssmis_l1c | PPS L1C (TMI)：10.65V | 1 張相同 |
| po_avhrr | MetOp-B AVHRR EPS native：CH4 色表 | 可出圖；因改用 WGS84 與共用畫布，與原版不逐像素相同 (目視一致) |
| po_viirs_i5 | NOAA-20 VIIRS I5 色表，3 個片段自動接成一筆 | 中心時間、颱風眼溫度與原版相同；圖 2897 x 2985 (原版 2894 x 3000)，目視一致 |
| po_viirs_dnb | NOAA-21 VIIRS DNB (ERF)，7 個片段中只讀相交的 5 個 | 中心時間與原版相同，目視一致 (新版用 `Latitude_TC`) |
| po_viirs_tc | NOAA-21 VIIRS TrueColor-I (M3 / M4 / M5 + I1) | 中心時間與原版相同，目視一致 |

dpi 模式：三類各以一個樣本輸出 1000 / 500 dpi 兩張。AVHRR 的 1000 dpi 圖為 2830 x 3696 像素，原版 PO 的 pcolormesh
為 2818 x 3696 (高度相同；寬度差 12 像素來自圓球改橢球)。
