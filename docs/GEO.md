# Satima：靜止衛星 (GEO)

> 這份是 **GEO (靜止衛星)** 的詳細說明，內容來自原本的 Satima-GEO。整合版的執行方式、資料列表、設定檔分層
> (`Satima.toml` + `Satima-GEO.toml`)、介面語言、圖片模式 (native / dpi) 與極軌的接續判定請看 [README.md](../README.md)。
> 以下提到的設定，共同的幾節 (`[general]`、`[map]`、`[coastline]`、`[image]`、`[title]`、`[annotation]`、`[output]`…)
> 寫在 `Satima.toml`，GEO 才有的寫在 `Satima-GEO.toml`；任何一節都可以在 `Satima-GEO.toml` 覆寫。
> 文中「原始解析度」的行為指 `[image] mode = "native"`；`mode = "dpi"` 時圖片大小改由 dpi 決定，取樣方式不變。


支援 Himawari-8/9 (AHI)、GOES-16\~19 (ABI)、GK-2A (AMI)、MTG-I (FCI)、MTSAT-1R/2 (Imager)、COMS-1 (MI)、
GMS-1\~5 / GOES-9 (VISSR)、GOES-8\~15 (Imager，GVAR / AREA)、SMS-1/2、GOES-1\~7 (VISSR / VAS，McIDAS AREA 與 Mode A / AA / AAA 原始資料)、
Meteosat-8\~11 (MSG SEVIRI，HRIT / native / NetCDF)、Meteosat-2\~7 (MFG MVIRI，FCDR / NetCDF)：
單波段、紅外/水氣色表強化、CIRA 真實色 (True Color)，以及加上夜間雲與城市燈光的 GeoColor；
MTG 的 FDHSI 與 HRFI 分開列出，HRFI 的 TrueColor 以 0.5 km 紅光銳化。

## 執行

```
python Satima.py                  # 搜尋設定檔中 data_dir 的資料
python Satima.py -d D:\data\h9    # 指定資料目錄
python Satima.py -f AHI-HSD       # 只找某個格式
python Satima.py -o D:\png        # 指定輸出目錄
python Satima.py --list           # 列出可用格式
```

流程：選資料 (可多個時次) → 地圖投影 (Enter = MERC) → 地圖範圍 (區域觀測自動帶入；FD 不需要) →
產品 (可複選) → 色表 → 圖片寬度 (Enter = 原始解析度；`x3` = 原始解析度的 3 倍)。

設定在 `Satima.toml` (共同) 與 `Satima-GEO.toml` (GEO 才有的與要覆寫的；UTF-8，可寫中文)；沒寫的項目用 `satima/core/config.py` 與 `satima/geo/config.py` 的預設值。

設定檔中的目錄 (`data_dir`、`output_dir`、`colormap_dir`、`conversion_table_dir`) 都可寫絕對路徑或相對路徑；
相對路徑一律以 Satima.py 所在目錄為準 (預設 `data/` 放資料、`plot/` 放圖)，與從哪裡執行無關。
命令列 `-d` / `-o` 的相對路徑則依一般慣例以目前目錄為準。

地圖套件由 `[map] backend` 選擇 `basemap` 或 `cartopy`。衛星影像的取樣完全相同，
只有底圖線條的資料來源不同：basemap 為 GSHHS，cartopy 為 Natural Earth。
`[map] resolution` 寫 `c/l/i/h/f` 或 `110m/50m/10m` 都可以，會自動對應。

### 地圖投影

執行時選擇 (直接 Enter = `[map] projection`，預設 merc；`ask_projection = false` 則不詢問)：

| 選單 | 投影 | 說明 |
|---|---|---|
| MERC | 麥卡托 | 預設，與舊版相同 (檔名不加投影) |
| CYL | 等距圓柱 | 經緯度等間距 |
| LCC | 蘭伯特正形圓錐 | 標準緯線 `lcc_parallels` 預設取南北範圍的 1/6 與 5/6 處，中央經線 `lcc_central_longitude` 預設為範圍中央；圖框為涵蓋整個經緯度範圍的最小矩形 |
| FD | 全圓盤 (衛星視角) | 不需輸入經緯度範圍；星下點、衛星高度取自資料 |

MERC 與 LCC 的地球模型一律為 WGS84 橢球 (同 PROJ `ellps=WGS84`，衛星導航的經緯度本來就是 WGS84 大地緯度)：
* basemap 預設是半徑 6370997 m 的圓球，Satima-GEO 改成與 cartopy 相同的橢球。
* 兩種底圖的格點、影像因此逐像素相同 (麥卡托、LCC、CYL 各實測，經緯度差 < 1e-12°)，只有海岸線資料不同。
* 舊版 (basemap 圓球) 的圖：同一範圍的麥卡托圖高不同 (例如北緯 5\~35° 為 1692 -> 1685 列)，位置最多差約半個像素 (0.9 km)。

FD 的圖框依各產品的波段決定：
* 標準靜止衛星投影的資料 (Himawari、GOES-R、GK-2A、MTG、MTSAT、MSG、MFG FCDR…)，全圓盤觀測時圖框就是原始影像
  (例如 MTSAT VIS 11000x11000、IR 2750x2750)；原始解析度輸出時，輸出像素與原始像素一一對應，等於直接畫原始影像
  (已逐像素驗證)。此時行列號由掃描角線性換算，不經過經緯度，大圖也很快。
* 區域觀測 (半圓盤、RSS、mesoscale 等) 仍畫整個圓盤，圖框為地球邊緣外加 1%。
* 自旋掃描的舊衛星 (GMS VISSR、GOES-1\~15 等) 原始影像不是理想的靜止衛星投影 (姿態、東西向取樣間距不同)，
  圖框為地球邊緣外加 1% 的理想投影，每個輸出像素仍以該衛星的精確導航取樣。
* basemap 不支援 GOES-R 的 sweep = x，FD 改用 sweep = y (仍以精確導航取樣，只是不再與原始像素一一對應)；
  cartopy 則支援。

檔名欄位 `{projection}`：MERC 為空 (與舊版檔名相同)，其餘為 CYL、LCC、FD，例如
`Himawari-9.AHI.BAND-13.BD.FD.202611080400Z.png`。

## 目錄

```
Satima.py                     進入點 (三類共用)
Satima.toml、Satima-GEO.toml   共同設定檔、GEO 的設定檔
satima/core/              三類共用的核心 (與衛星無關)
  config.py                 分層設定檔的讀取、共同的預設值、[readers.<格式>] 個別覆寫
  session.py                主流程：搜尋資料、分段的資料列表、選資料
  ui.py                     互動輸入
  maps.py                   地圖投影 (MERC / CYL / LCC / FD) 與底圖：basemap 或 cartopy ([map] backend)
  canvas.py                 圖的大小 (native / dpi)、標題與標註、存檔
  colormaps.py              色表 (.cmap)
  naming.py                 檔名樣板
  parallel.py               多核心 (執行緒平行，[memory] workers)
  pngio.py                  多核心 PNG 寫入 (分段平行壓縮)
satima/geo/               GEO 的程式
  config.py                 GEO 自己的預設值與檢查
  family.py                 GEO 的流程：列表的各列、投影與範圍、產品選單、逐張繪製存檔
  reader.py                 讀取器基底類別 GeoReader、Channel、Scene
  nav.py                    導航：GeosNav (靜止衛星投影)、ProjGridNav (規則投影格點)、LonLatNav (逐點經緯度)
  resample.py               輸出格點與取樣 (只讀涵蓋地圖範圍的視窗)
  products.py               產品選單與單波段/色表繪製
  landsea.py                灰階底圖：海陸遮罩 (GSHHS 柵格化) 與雲 + CWA 海陸 / LIGHTS 夜光底圖的合成
  bluemarble.py             灰階底圖 NASA / NICT：Blue Marble 2004 / 2002 (下載、JPEG 分塊 GeoTIFF) 與白雲半透明的合成
  _resample_nb.py           取樣的 numba 核心 (區塊平均、取值、全圓盤行列號與經緯度)
  truecolor/                真實色 (Rayleigh 校正、太陽天頂角校正、混合綠光、銳化)
    kernel.py                 逐像素合成的 numba 核心 ([truecolor] engine)
    geocolor.py               GeoColor (夜間紅外高雲、低雲、城市燈光與日夜混合)
    background.py             城市燈光 (內建分塊 GeoTIFF 讀取器，不需 GDAL；單波段 DEFLATE、RGB JPEG) 與海陸遮罩 (basemap 內建)
sensors/geo/              各衛星讀取器，檔名 <感測器>-<格式>.py
  AHI-HSD.py                Himawari-8/9 Himawari Standard Data (.DAT / .DAT.bz2)
  ABI-NC.py                 GOES-16~19 L1b Radiances NetCDF
  AMI-NC.py                 GK-2A L1B NetCDF
  FCI-NC.py                 MTG-I FCI L1c FDHSI / HRFI (zip 直接讀取)
  MTSAT-HRIT.py             MTSAT-1R / MTSAT-2 Imager JMA HRIT (.gz，或整包 .tar 直接讀取)
  MI-HDF.py                 COMS-1 MI L1B HDF-EOS5 (.he5) — 依舊版 COMS_HDF.py 改寫，尚未以實際資料驗證
  GMS-VISSR.py              GMS-1~5、GOES-9 JMA VISSR (VS/IR*.Z?? 與 VISSR_*.A.IMG，gzip，可放在 tar 中)
  GMS-TAPE.py               GMS-1 VISSR 磁帶映像 (FGGE 1978~1979；36-bit word，檔名為卷號，例如 VS2303)
  GMS-AREA.py               GMS McIDAS AREA (GMSX 導航，例如 GMS_IR_1_1979_007_1533.area)
  vissr_common.py           (輔助) JMA VISSR 參數區的欄位定義、校正表、導航參數與地球邊緣修正，GMS-VISSR、GMS-TAPE 共用
  vissr_nav.py              (輔助) JMA VISSR 導航 MGIVSR 與地球邊緣修正，GMS-VISSR、GMS-TAPE、GMS-AREA 共用
  vissr_nav_cy.pyx          (輔助) 同上的 Cython 版；build_vissr_nav_cy.py 編譯成 .pyd
  GOES-GVAR.py              GOES-8~15 Imager 原始 GVAR (NOAA CLASS goes??.yyyy.ddd.hhmmss、NSMC GOS*_GVR_*.dat)
  GOES-AREA.py              GOES-8~15 Imager McIDAS AREA (CLASS GVAR_IMG *.BAND_0?)
  goes_imager.py            (輔助) GOES I-M O&A 導航 (同 McIDAS nvxgvar)、Gould / BCD 數值、校正
  goes_nav_cy.pyx           (輔助) GOES 導航的 Cython 版
  gvar_format.py、gvar_area_format.py  (輔助) 舊版 GVAR_format.py、GVAR_AREA_format.py 原樣複製的格式定義
  GOES-PREGVAR.py           SMS-1/2、GOES-1~7 VISSR / VAS McIDAS AREA (goesNN.yyyy.ddd.hhmm.AREA_IR / AREA_VIS)
  mcidas_goes_nav.py        (輔助) McIDAS 'GOES' 導航 (自旋掃描，同 nvxgoes / VisAD GOESnav) 與精確反算
  mcgoes_nav_cy.pyx         (輔助) 同上的 Cython 版
  pregvar_area_format.py    (輔助) 舊版 AREA_format.py 原樣複製的格式定義
  GOES-STRETCHED.py         SMS-1/2、GOES-1~7 Mode A / AA / AAA 原始 stretched 資料 (NCEI 2024：*.A.data、*.AA-MSI.data、*.AAA-DS.data…)
  stretched_format.py       (輔助) Mode A / AA / AAA 的記錄配置、文件欄位、O&A、校正，以及 VAS 波段表與亮溫係數 (GOES-PREGVAR 共用)
  SEVIRI-HRIT.py            MSG SEVIRI EUMETSAT HRIT (H-000-MSG*，含 -C_ 壓縮、.bz2)
  SEVIRI-NAT.py             MSG SEVIRI Level 1.5 native (MSG?-SEVI-MSG15-*.nat，或 Data Store 的 .zip 不必解開)
  SEVIRI-NC.py              MSG SEVIRI Level 1.5 NetCDF (W_XX-EUMETSAT-Darmstadt,...SEVIRI_C_EUMG_*.nc) — 尚無樣本驗證
  seviri_common.py          (輔助) SEVIRI 頻道、定標 (nominal / GSICS / Meirink)、導航 (含 HRV 上下視窗)
  seviri_native_hdr.py      (輔助) 舊版 prologue.py (satpy 的 native 檔頭結構) 複製，np.bool 改 np.bool_
  MVIRI-FCDR.py             MFG MVIRI FCDR (MVIRI_FCDR-EASY/FULL_*.nc 或 .zip、FIDUCEO_FCDR_*.nc)
  MVIRI-NC.py               MFG MVIRI EUMETSAT NetCDF (...MVIRI_C_EUMS_*.nc，舊版 MFG.py 格式) — 尚無樣本驗證
colormap/                 色表檔
conversion_table/         校正表 (GK-2A 等)
nightlights/              GeoColor 的城市燈光：VIIRS 2024 年平均夜光 500 m (Cloud-Optimized GeoTIFF，61 MB)
landsea/                  CWA 底圖的海陸遮罩 (GSHHS 全解析度 30 角秒，不存在時自動產生)
bluemarble/               NASA / NICT 底圖的 Blue Marble 2004 (約 1 km) / 2002 (約 1.85 km) 分塊 GeoTIFF (不存在時自動下載)
```

## 處理方式

每個輸出像素的經緯度 → 該波段的導航 `nav.to_rowcol()` → 原始影像行列號 → 取值。

* 只讀取涵蓋地圖範圍的視窗 (邊界多留 2 個像素)，不需要整張影像的經緯度，也不需要 KD-tree。
  GOES 小範圍時讀取量比舊版少數十到數百倍。
* 輸出像素對應的原始像素存成整數索引 (最近的像素中心)，經緯度分段計算，大圖 (全圓盤 VIS 等) 不需整張的經緯度陣列。
  與投影無關：MERC / CYL / LCC / FD 都是同一套取樣，海岸線由同一個投影畫出，兩者一定吻合。
* 上色由 Satima-GEO 以 float64 依色階定義計算 (colormaps.colorize)，再把 8 位元影像交給 matplotlib 疊海岸線與標註。
  matplotlib 自己上色時會先以 float32 縮放，極細的色表 (例如 BD 約 5 萬階) 偶有少數像素差一階；灰階與真實色完全相同。
* 輸出比原始解析度粗時先做區塊平均；區塊起點對齊倍數，所以 0.5 km 紅光的 2x2 區塊與 1 km 藍光像素自然對齊 (真實色不會色散)。
* 行列號一律依 CGMS 規範由 1 起算 (第 1 列第 1 行 = 左上角)。

### 放大時的取值 (`[image] interpolation`，預設 bilinear)

輸出像素比原始像素細或相近時 (區塊平均倍數 = 1)，取值方式有兩種：

* **最近鄰 (nearest)**：直接取最近的原始像素。放大倍率通常不是整數，而且隨緯度 (麥卡托) 與衛星斜視角改變，
  所以部分原始像素被重複畫成 2 格、部分只有 1 格。畫面因此出現不規則的階梯與拼貼感，色表強化 (BD 等) 的色階邊緣最明顯。
  * 例：GOES-19 C13 (標稱 2 km)，北緯 5\~35° 的麥卡托圖，原始解析度寬 2612 像素。
  * 風眼附近 (北緯 17°) 約 8% 的原始像素被重複。
  * 北緯 30° 附近輸出像素約 1.84 km，原始像素約 2.1\~2.4 km，約 40% 的原始像素被重複 2\~3 次。
* **雙線性 (bilinear)**：以相鄰 4 個原始像素內插 (偏移量化為 1/254 像素，存成 int8；最近像素索引不變)，邊緣平順。
  * 鄰點有缺值 (地球邊緣、太空) 時改用最近像素，所以邊界與最近鄰相同。
  * 單一像素的極值會被鄰點沖淡一些。最高/最低雲頂溫度仍以原始像素計算，與最近鄰相同。
  * 輸出像素與原始像素一一對應時 (Himawari、MTSAT、SEVIRI 全圓盤原始解析度)，與最近鄰逐像素相同。

適用於單波段、色表強化與 GeoColor 的夜間紅外。真實色的白天波段維持最近鄰：1 km 波段以 0.5 km 紅光做比值銳化，
細節由紅光補回。分條、多核心、numba 與 numpy 版的結果都逐像素相同。

#### 呈現原始像素：原始解析度的倍數 (`[image] size = "x4"`，執行時也可輸入 `x3`)

舊版的 `m.pcolormesh` 與最近鄰相同：同樣的圖寬下，91\~99% 的輸出像素取到同一個原始像素，其餘只在兩個原始像素的交界差 1 格。
舊版看起來沒有拼貼感，是因為固定 2000 像素寬 (縮小)。以 GOES-19 北緯 5\~35° 為例，風眼附近約 36% 的原始像素沒有畫出來。

衛星像素投影到麥卡托等地圖後大小不一，均勻的輸出格點不可能剛好一格對一個原始像素：不是重複、丟失，就是內插。
要看到所有原始像素且數值不變，可以把圖寬放大成原始解析度的 3\~4 倍並用最近鄰 (`interpolation = "nearest"`)。
每個原始像素會畫成約 3x3\~4x4 的方塊 (即高解析度 pcolormesh 的效果)：

| 圖寬 (GOES-19 C13，北緯 5\~35°、西經 100\~50°) | 風眼附近沒畫到的原始像素 | 方塊大小不均 (風眼 / 北緯 30°) |
|---|---|---|
| 2000 (舊版) | 約 36% | 0% / 22% |
| 2612 (原始解析度) | 約 0 | 24% / 44% |
| 10448 (x4) | 約 0 | 15% / 13% |

方塊大小不均 = 每個原始像素所佔輸出像素數的標準差 / 平均。

* `size` 可寫 `"native"`、`"x2"`\~`"x8"` (執行時直接 Enter 的寬度) 或固定寬度。圖寬提示輸入 `x3` (或 `3x`) 則只用於這一次。
* 原始解析度一律放大成**整數倍**，每個原始像素的大小才一致：
  * 寬度不足 `min_size` 時，取達到下限的最小倍數，例如 357 像素 -> x6 = 2142。以前直接放大到 2000 (5.6 倍)，原始像素會不規則地畫成 5 或 6 格。
  * 寬度不足 `min_size` 時至少放大 `min_scale` 倍 (預設 4)，例如 1500 像素 -> x4 = 6000，而不是 x2 = 3000；
    `min_scale = 1` 則只取達到下限的最小倍數。原始寬度已達 `min_size` 時不受影響。
  * 長邊超過 `max_size` 時，改用放得下的最大倍數。原始解析度也放不下，才縮小到上限。
* basemap 與 cartopy 的取樣完全相同 (衛星影像由 Satima-GEO 取樣，地圖套件只畫海岸線)，所以麥卡托等地圖的重複情形兩者一樣：
  北緯 30° 附近都有 40\~50% 的原始像素被重複。唯一的差別在 GOES-R 全圓盤：basemap 不支援 sweep = x，沒有與原始像素對齊
  (原始解析度 5489 寬)；cartopy 則與原始影像 5424x5424 一一對應。
  * 以前 basemap 的麥卡托是圓球、cartopy 是 WGS84 橢球，格點錯開約半格，重複的列落在不同位置，所以同一處看起來略有不同。
    現在統一為 WGS84 (見「地圖投影」)，兩者逐像素相同。
* 全圓盤與原始像素對齊時 (Himawari、GK-2A、MTG、MTSAT、SEVIRI；GOES-R 需 cartopy)，原始解析度已經一格對一個原始像素，不套用倍數。
  放大只會把每個原始像素複製成 N x N；0.5 km、1 km 波段還會碰到長邊上限，變成非整數倍 (例如 22000 x4 -> 25000，1.14 倍)，
  又出現不規則的重複。GOES-R 用 basemap、GMS、GOES-1\~15 等舊衛星的全圓盤沒有對齊，仍會放大。
* 原始像素不是方的衛星 (GMS、SMS / GOES-1\~15 等掃描式儀器)：東西向的取樣間距比南北向的掃描線間距密，
  例如 GOES-8\~15 紅外是東西 2.29 km × 南北 4.01 km (可見光 0.57 × 1.00 km)，GMS-5 紅外是 3.43 × 5.02 km。
  * 輸出像素是方的。「原始解析度」以**較細的方向**為準 (GOES 紅外 2.29 km)，所有原始像素都保留；較粗方向的每個原始像素畫成多格
    (GOES 約 1.75 格、GMS-5 約 1.47 格)。執行時會印出兩個方向的大小與採用的解析度。
  * 以前以南北向的標稱解析度 (GOES 紅外 4 km、GMS 5 km) 為準，東西向被壓縮：全圓盤只用到 57% 的原始像素 (GOES)，現在是 100%，
    圖寬由 2744 變成 4793 (GOES 紅外)、2195 變成 3204 (GMS-5 紅外)、10977 變成 19172 (GOES 可見光)。
  * 像素大小由導航量測 (影像中央附近相鄰像素的掃描角間距，見 `maps.source_pixel_km`)，不必在各讀取器另外設定；
    McIDAS AREA 等檔案若經過抽樣，量到的就是檔案實際的間距。標準靜止衛星投影 (Himawari、GOES-R、MTSAT、SEVIRI 等) 的像素是方的，不受影響。
  * 所有投影 (麥卡托、LCC 等) 的原始解析度都以此為準。南北向重複的列數不是整數，要每個原始像素大小一致可再用 `x4`
    (GOES：每個原始像素約 4 × 7 格)。
* 像素數是倍數的平方。上例 x4 為 10448x6770，耗時 10.6 秒，PNG 約 11 MB。
* 海岸線與文字以圖寬為準等比例放大 (圖固定 5 吋寬)，相對於衛星像素的粗細不變。放大檢視時覺得太粗，可調小 `[coastline] linewidth`。

### 真實色 (`[truecolor]`)

依 CIRA GeoColor 的白天真實色 (SHAC True Color；Miller et al. 2016 BAMS、2020 JTECH)：

1. 太陽天頂角校正。
2. Rayleigh 散射校正 (`rayleigh_correction`)。藍、綠、紅、近紅外都會校正，並依 10.3 µm 亮溫削減雲頂上方的路徑反射 (283 K 不減，233 K 減 70%)。
   * `"full"` (預設)：pyspectral 的路徑反射扣除後，再除以上下行穿透率並考慮大氣球面反照率。這是 SeaDAS / CREFL 的作法，CIRA 也改自 SeaDAS。
   * `"subtract"`：舊版作法，只扣除路徑反射。雲會偏黃，晴空偏暗，晨昏線附近偏橘。
   * `cloud_rayleigh_transmittance = true` (預設)：高雲的穿透率與球面反照率也用同一個削減係數 (雲頂上方只剩這個比例的大氣)。
     只削減路徑反射、穿透率仍用整層大氣時，高雲的藍光被除以過小的穿透率 (太陽天頂角 80° 時藍光約 0.60、紅光約 0.85)，雲頂偏藍。
   * `ozone = 0.30` (預設，atm-cm；0 或 false 關閉)：臭氧 Chappuis 帶的吸收校正，反射率除以臭氧穿透率。
     * 臭氧吸收紅光最多、藍光最少。吸收係數以 Leckner (1978) 的表對各波段光譜響應平均，AHI 的 0.47 / 0.51 / 0.64 µm 為
       0.011 / 0.040 / 0.082 每 atm-cm；斜程倍數用臭氧層高度 22 km 的球面幾何 (SPCTRAL2)。
     * 臭氧量是固定值。實際約 0.24\~0.45，熱帶偏低、中高緯度春季偏高，高緯度會有殘差。
   * 這兩項 CIRA 論文沒有，是為了日出日落時雲頂偏藍而加的。深對流雲 (10.4 µm < −55°C) 校正後的 藍/紅 比值 (Himawari-9，未乘增益)：

     | 太陽天頂角 | < 50° | 50\~70° | 70\~78° | 78\~82° | 82\~86° |
     |---|---|---|---|---|---|
     | 兩項關閉 | 1.16 | 1.17\~1.26 | 1.35 | 1.42 | 1.53 |
     | 兩項開啟 | 1.06 | 1.06\~1.07 | 1.09 | 1.10 | 1.12 |

     兩項都關閉時，輸出與加入之前逐像素相同。下方與 CIRA SLIDER 的比對表是兩項關閉時的結果。
3. 綠光
   * AHI / AMI / FCI：0.91 × 0.51 µm + 0.09 × 0.86 µm (hybrid green，`hybrid_green_fraction`)。
     CIRA 論文寫的是約 7%。和 CIRA 的實際圖片比對，7% 時四顆衛星的植被都比 CIRA 偏紅；9% 時色調最接近，荒漠只變綠 1\~2 個色階。
   * ABI 沒有綠光，`abi_green = "ahi"` 照 CIRA 的作法處理：先以 AHI 合成 0.51 µm 綠光，再同樣混入 0.86 µm。
     * CIRA 的查表沒有公開，這裡用 Himawari-9 兩個時次 (北緯 0\~50 度) 擬合出 G0.51 = 0.976 × 藍^0.668 × 紅^0.347。和 AHI 實際綠光的誤差中位數小於 1 個色階。
     * 也可以改設 `[藍, 紅, 近紅外]` 線性權重，例如舊版 CIMSS 的 `[0.45, 0.45, 0.10]`。用這組權重陸地會明顯偏綠。
4. 0.5 km 紅光比例銳化、晨昏線淡出、log10 拉伸 (反射率 0.025\~1.20 → [-1.6, 0.176]，與 CIRA 論文相同)。
5. 色調曲線 (`tone`) 與飽和度 (`saturation`)，預設 `tone = "cira"`、`saturation = 1.0`，輸出與加入這兩項之前逐位元相同。
   * CIRA 原式有兩個先天限制：
     * 反射率截在 1.2，白點卻在 10^0.176 = 1.5，輸出最亮只到 241；厚雲 (反射率 0.75 左右) 約 212，看起來偏灰，
       反射率 0.6\~1.0 只佔 32 個色階，颱風雲頂的層次很少。
     * 反射率低於 0.025 一律為 0。Rayleigh 校正後晴空海面的紅光約 0.02，五\~六成的海面像素紅色是 0，海面呈塊狀的藍黑色。
   * `tone = "enhanced"`：同一條 log 曲線，中間調幾乎不變，只改兩端。
     * 白點降到 `white_point` (預設 1.25)，顯示值超過 `shoulder` (0.95) 的部分平滑壓縮到 1，最亮的雲頂不會整片死白。
     * 暗部以平滑下限 (`shadow_toe`，預設 0.006) 取代硬截斷，低於 0.025 的反射率仍有層次。
   * `saturation`：拉伸後以亮度為中心放大色差。顯示亮度 0.50 以下為設定的倍數，往上遞減，0.85 以上不放大，雲不會染色。
   * `tone = "enhanced"`、`saturation = 1.25` 與 CIRA 原式比較 (8 位元色階中位數)：

     | 衛星 (時間 UTC) | 曲線 | 晴空海面 | 紅 = 0 的海面 | 植被 | 亮雲 | 全圖最大值 |
     |---|---|---|---|---|---|---|
     | Himawari-9 (2026-10-21 01:00) | cira | 0, 28, 57 | 52% | 66, 85, 64 | 212, 213, 218 | 241 |
     | | enhanced | 9, 33, 68 | 0% | 66, 91, 65 | 222, 224, 228 | 254 |
     | GOES-16 (2019-09-01 16:00) | cira | 0, 37, 67 | 57% | 70, 94, 76 | 209, 214, 218 | 241 |
     | | enhanced | 6, 43, 80 | 5% | 69, 101, 77 | 219, 224, 228 | 251 |

   * 下方與 CIRA SLIDER 的逐像素比對都是 `tone = "cira"` 的結果；改用 enhanced 後雲與海面會比 CIRA 亮。

**Himawari-9 的波段相對增益 (`band_gain`)**：Himawari-9 的雲在真實色中偏青綠 (修正前雲的色階 綠−紅 約 +5、藍−紅 約 +9)。
以深對流雲 (10.4 µm < −50°C、紅光反射率 > 0.8) 檢查各波段，這種雲在可見光幾乎是白色：

| 衛星 (時間 UTC) | 藍 | 綠 | 紅 | 近紅外 | 藍/紅 | 綠/紅 |
|---|---|---|---|---|---|---|
| Himawari-8 (2020-09-04 03:00) | 0.927 | 0.902 | 0.897 | 0.919 | 1.035 | 1.008 |
| Himawari-9 (2026-10-21 01:00) | 0.992 | 0.934 | 0.886 | 0.910 | 1.124 | 1.057 |
| GOES-19 (2026-10-28 13:00) | 0.885 | — | 0.882 | 0.926 | 1.006 | — |

(反射率為太陽與衛星天頂角 < 40° 的中位數。) 三顆衛星的紅光、近紅外一致，只有 Himawari-9 的 B01、B02 偏亮。
Himawari-9 另外兩個時次 (2023-07-01 03:00、2026-11-09 05:00) 的 藍/紅 也是 1.11\~1.12、綠/紅 1.04\~1.05，不是偶發的。
改用檔頭的更新定標係數也無法消除這個差異。
(上表是 TOA 反射率，未做臭氧與 Rayleigh 校正。)

加入臭氧與高雲穿透率校正後，同樣的深對流雲校正後的比值 (太陽天頂角 < 50°、未乘增益)：

| 衛星 (時間 UTC) | 藍/紅 | 綠/紅 |
|---|---|---|
| Himawari-8 (2020-09-04 03:00) | 0.98 | 0.98 |
| Himawari-9 (2026-10-20 22:00、10-21 01:00) | 1.06、1.07 | 1.02、1.02 |
| GOES-16 (2019-09-01 16:00) | 1.06 | — |
| GOES-19 (2025-10-28、2026-10-28 13:00) | 0.95、0.95 | — |

各衛星的定標差異仍在，所以保留 `band_gain`。真實色 / GeoColor 預設把 Himawari-9 的 B01 乘 0.94、B02 乘 0.98，
修正後兩個比值在太陽天頂角 < 70° 都是 1.00，80° 時為 1.03、1.01 (兩項校正關閉時的舊值為 0.92、0.95)。
`band_gain` 的寫法是 `{ 衛星 = { 波段 = 倍數 } }`，其他衛星也可以加；單波段影像不受影響；設 `band_gain = false` 可關閉。

與 CIRA SLIDER GeoColor 逐像素比對，數字為色階差 (0\~255) 的絕對值中位數，以「新版 / 舊版」表示：

| 衛星 (時間 UTC) | 紅 | 綠 | 藍 | 植被色調 (綠−紅) |
|---|---|---|---|---|
| GOES-19 CONUS (2026-11-08 18:01) | 4 / 3 | 4 / 8 | 4 / 11 | 0 / +8 |
| Himawari-9 全圓盤 (2026-11-09 05:00) | 2 / 4 | 3 / 5 | 5 / 7 | 0 / −6 |
| GK-2A 全圓盤 (2026-11-08 04:00) | 3 / 4 | 4 / 5 | 8 / 5 | +2 / −5 |
| MTG-I1 全圓盤 (2026-11-09 13:20) | 2 / 4 | 3 / 10 | 5 / 20 | −4 / −10 |

* 植被色調一欄是「我們的 (綠−紅) − CIRA 的 (綠−紅)」：負值表示比 CIRA 偏紅、不夠綠，正值表示偏綠。
* 雲和陸地的色調一致：雲是白的，美國西部、撒哈拉、澳洲的乾燥地表都和 CIRA 一樣是棕褐色。
* 剩下的差異：
  * CIRA 各衛星的產品彼此並不完全一致。例如 GK-2A 晴空的藍色比本程式暗約 12 個色階，MTG 的則亮約 10 個色階。
  * 淺海處 CIRA 比較青綠。CIRA 的說明書也提到，它的查表在淺海可能有雜訊或不正確。
  * 這些差異應該來自 CIRA 未公開的 SeaDAS 參數與查表。

### GeoColor (`[geocolor]`，第 2 種真實色風格)

依 CIRA GeoColor V2.0 (Miller et al. 2020 JTECH 第 5 節)。資料有 3.9 µm 與 10.4 µm 波段時 (AHI、ABI、GK-2A、MTG)，
產品選單的 RGB 會在 TrueColor 之後多一個 GeoColor：

* **白天**：與 TrueColor 相同，包括 SZA 80°\~90° 的晨昏線淡出 (CIRA 的白天圖層在晨昏帶同樣會變暗)。
* **夜間**：由上到下疊三層。
  1. 高雲：10.4 µm 亮溫越冷越白。
  2. 低雲：10.4 − 3.9 µm 亮溫差畫成淺藍；陸地與水面用不同門檻，海陸遮罩取自 basemap 內建資料。
  3. 地表：城市燈光畫成鈉燈的橘黃色；其他陸地為深紫，水面為黑色。
* **晨昏線**：白天與夜間依 cos(太陽天頂角) 混合。
  * 日出、日落時 3.9 µm 含太陽反射，低雲會短暫變淡；這段時間白天圖層佔多數，CIRA 也有同樣情形。
* **城市燈光**：`nightlights/` 內的 VIIRS 2024 年平均夜光，500 m 全解析度，依輸出像素大小自動選用檔內的縮圖層。
  * 資料來源：EOG VNL v2，取自 Zenodo 鏡像 (Hengl, OpenGeoHub, doi:10.5281/zenodo.17294744，CC BY 4.0)。
  * 用內建的分塊讀取器讀取，不需要 GDAL。basemap 自帶的 GEOS 與 GDAL 版本衝突，basemap 先載入時 GDAL 無法匯入。
  * 設 `nightlights = ""` 則不畫燈光。
* **ABI 合成綠光的藍色下限 (0.3 倍紅色)**：合成綠光是用 AHI 擬合的。晨昏線附近藍光幾乎被 Rayleigh 校正扣光，沒有下限時綠色也會趨近 0，晨昏帶變成紅色。
  這項也會影響 TrueColor。

預設參數取自 CIRA 目前 (2025 年更新後) SLIDER 上的圖片，與論文有三處不同：
* 紅外範圍 `ir_range = [178, 301]` K (論文 [200, 280])。冷端都隨緯度增加 20 K。
* 日夜混合 `twilight_cos = [0.0, 0.18]`、`twilight_power = 1.75`，SZA 約 80°\~90° 漸變 (論文 [0.1, 0.3]、1.5，約 72°\~84°)。
  CIRA 的說明是「縮窄晨昏帶、白天延伸最多約一小時」。
* 水面背景為黑色。CIRA 夜間晴空海面的灰色來自紅外圖層。

這些參數是用 GOES-19 (2026-11-12 23:01 UTC，日落橫跨美國) 與 MTG-I1 (2026-11-09 13:20 UTC) 逐像素反推的，兩顆衛星結果一致。
GOES-19 夜間各圖層和 CIRA 的色階差 (中位數)：

| 區域 | 本程式 RGB | CIRA RGB |
|---|---|---|
| 晴空陸地 | 54, 47, 64 | 60, 51, 65 |
| 晴空海面 | 14, 14, 14 | 13, 13, 13 |
| 低雲 (陸地) | 176, 211, 251 | 177, 212, 252 |
| 高雲 (< 230 K) | 165, 165, 166 | 166, 166, 166 |
| 城市燈光 | 201, 176, 144 | 191, 167, 135 |
| 晨昏 SZA 70°\~84° | 108, 99, 98 | 104, 104, 105 |

城市燈光比 CIRA 亮約 10 個色階，可能是 CIRA 用的夜光資料年份較舊 (論文為 2015 年)。
想要論文原版的外觀，把上面三個參數改回括號內的值即可。

### 灰階底圖 (GRAY 的第 2 層)

色表選了 `1.GRAY` 之後會再問底圖 (可複選，0 = 全部)：

```
BAND-13 colormaps: 0.ALL 1.GRAY 2.BD 3.COLOR ...
  選擇色表 (可複選，0 = 全部)：1
BAND-13 GRAY basemaps: 0.ALL 1.NONE 2.CWA 3.LIGHTS 4.NASA 5.NICT
  選擇底圖 (可複選，0 = 全部)：
```

* `NONE`：一般的灰階圖 (無底圖)。
* `CWA`：海陸底圖 (`[landsea]`，見下節)。紅外與可見光 / 近紅外都有 (可見光只有灰階，直接問底圖)。
* `LIGHTS`：夜光底圖 (`[nightmap]`)，只有紅外。
* `NASA`、`NICT`：Blue Marble 底圖 (`[bluemarble]`、`[nictmap]`)，紅外與可見光 / 近紅外都有
  (可見光的選單為 `1.NONE 2.CWA 3.NASA 4.NICT`)。
* 水氣看不到地表，沒有底圖，不會詢問。
* 名稱即檔名的 `{colorbar}`，例如 `BAND-13-CWA`、`BAND-13-LIGHTS`、`BAND-13-NASA`、`BAND-13-NICT`；可用各節的 `name` 修改。

### 海陸底圖 CWA (`[landsea]`)

灰階雲疊在單色的陸地 (橄欖綠) 與海洋 (深藍) 上，仿 CWA 官網的「紅外線彩色雲圖」。

* 合成：輸出 = a × 灰階 + (1 − a) × 底色。
  * 紅外：灰階為亮溫的分段線性對照 (`ir_gray`，+80°C 黑 \~ −95°C 白)；a = clip((26.1 − T) / 70.8, 0, 1)^0.8，
    約 +26°C 以上全透明、−45°C 以下不透明。
  * 可見光：灰階 = 反射率；a 由反射率 0.10 (全透明) 漸變到 0.50 (不透明)。
* 底色只有兩色，海岸處依陸地比例混合。官網陸地上的地形紋理其實是紅外本身 (夜間山區地表較冷，透明度與灰階隨之改變)，
  不是另外的山影圖。
* 預設值是以官網 2026-11-13 11:30 UTC 的 Himawari-9 全圓盤「紅外線彩色」與「紅外線灰階」圖，對照同時次的 B13 亮溫
  逐像素反推的 (透明度擬合誤差 0.009)。同時次以本程式出圖，與官網逐像素的色階差中位數為 R 3、G 3、B 5
  (官網為 JPEG，且有黃色海岸線與經緯線)。
* 海陸遮罩 (`landsea/gshhs_f_30s.tif`)：basemap 內建的 GSHHS 全解析度海岸線 (湖泊算水、湖中島算陸地)，
  以 numba 掃描線演算法柵格化成 30 角秒 (約 0.93 km) 的陸地比例 (垂直 4 倍取樣、水平精確覆蓋比例，海岸平滑無鋸齒)，
  存成分塊 GeoTIFF，含 2、4、8… 倍的縮圖層。取樣時選不比輸出像素粗的最粗一層並雙線性內插，跨換日線只讀用到的分塊。
  檔案不存在時第一次使用會自動產生 (numba 約 30 秒；沒有 numba 時改用 matplotlib，約 12 分鐘)。
* 預設不畫海岸線、國界與經緯線 (`lines = false`)；改成 true 則依 `[coastline]` 等設定畫線。

### 夜光底圖 LIGHTS (`[nightmap]`)

雲疊在 VIIRS 年平均夜光 (城市燈光) 上，做法同 GeoColor 的夜間圖層，只是不含需要 3.9 µm 的低雲圖層：

* 底色：`[geocolor]` 的 `land_color` (深紫) / `water_color` (黑)，海岸處依 `[landsea] mask` 的陸地比例混合；
  城市燈光處改為鈉燈的橘黃色 (與 GeoColor 相同)。燈光檔為 `[geocolor] nightlights`，找不到時略過此產品。
* 雲：輸出 = N_IR × 白 + (1 − N_IR) × 底色，N_IR = 1 − N(亮溫)[IRmin, IRmax]，即 GeoColor 高雲圖層的透明度
  (`[geocolor] ir_range`，預設 178\~301 K，冷端在緯度 30°\~60° 之間增加 20 K)。
* 調整 `[geocolor]` 的參數會同時影響 GeoColor 與 LIGHTS。`[nightmap]` 只有 `name` 與 `lines` 兩項。

### Blue Marble 底圖 NASA / NICT (`[bluemarble]`、`[nictmap]`)

白雲依亮溫 / 反射率半透明，疊在 NASA Blue Marble 上，仿 NICT「Himawari-8 Real-time Web」
(https://himawari8.nict.go.jp/) 各波段疊 Blue Marble 的畫面。兩種只差在底圖，雲的合成相同 (設定都在 `[bluemarble]`)：

| 選單 | 設定段 | 底圖 | 解析度 | 外觀 |
|---|---|---|---|---|
| `NASA` | `[bluemarble]` | Blue Marble Next Generation 2004 (NASA GIBS) | 約 0.98 km | 較暗：海洋近黑，青藏高原有積雪 |
| `NICT` | `[nictmap]` | Blue Marble 2002 `land_shallow_topo` (NICT 網站所用的版本) | 約 1.85 km | 較亮：海洋為均勻的藍紫色 (RGB 10, 10, 51)，淺海為青綠色，高原無雪 |

* NICT 網站的這種畫面是兩層：每個時次的雲圖是**純白色 + 透明度**的 PNG，疊在固定的 Blue Marble 上。
  所以合成為：輸出 = a × 白 + (1 − a) × Blue Marble。無資料的像素只有底圖，地球以外透明。
* 透明度：以 NICT 網站的雲圖 (10d 層，5500 像素全圓盤) 對照同時次的 Himawari-9 HSD 逐像素反推：
  * 紅外 (`ir_alpha`)：B13，2026-10-20 21:00 UTC 全圓盤與 2026-11-09 05:00 UTC 北半部。
    不透明度與 **10.4 µm 的輻射量**成線性 (不是亮溫線性)：10.5°C 以上全透明、−80°C 為 230.5 (約 90%)、−84°C 約 235。
    兩個時次與 NICT 的差都只有 0.4 個色階 (8 位元量化誤差)；若改成亮溫線性，最多會差 28 個色階。
    所有紅外波段都用同一條亮溫 → 不透明度曲線。
  * 可見光 (`vis_alpha`)：B03，2026-11-09 05:00 UTC。不透明度 = 反射率 × 255 (NICT 的 2 km 是 0.5 km 的 4×4 平均，
    沒有太陽天頂角校正)，差 0.3 個色階。這正好是 `[visible]` 的預設值 (不校正、倍數 1、無 gamma)；
    改了 `[visible]` 也會套用到這裡。
* 底圖 (`image`) 都是 NASA 公有領域的圖 (引用時註明 NASA Earth Observatory / Goddard Space Flight Center)，
  存成分塊 JPEG GeoTIFF (含縮圖層；內容相同的大洋分塊共用同一份資料；GDAL 可讀)。
  取樣時選不比輸出像素粗的最粗一層，雙線性內插 (同海陸遮罩)。
  檔案不存在時第一次使用會自動下載；下載失敗時略過此產品並印出原因，同一次執行不再重試。
  * **NASA**：Blue Marble Next Generation (MODIS 2004 年合成)，取自 NASA GIBS 的 WMTS 圖層 `BlueMarble_NextGeneration`
    (EPSG:4326)。第 6 層約 0.98 km，40960×20480，4273 個圖磚，存檔 67.5 MB，實測約 6 分鐘。
    GIBS 的 512×512 JPEG 圖磚原封不動當成 GeoTIFF 的分塊 (不重新壓縮)，縮圖層 (約 2、4、8… km) 也直接取 GIBS 較粗的層級。
    要更細的底圖 (第 7 層約 0.5 km，估計約 250 MB、下載約半小時) 可執行
    `python -c "from satima import bluemarble; bluemarble.build('bluemarble/BlueMarble_NextGeneration_GIBS_500m.tif', 7)"`，
    再把 `image` 指到這個檔案。
  * **NICT**：Blue Marble 2002 `land_shallow_topo` (陸地、淺海與地形陰影)，取自 Wikimedia Commons 的
    [Land shallow topo 2048.jpg](https://commons.wikimedia.org/wiki/File:Land_shallow_topo_2048.jpg)
    (檔名雖為 2048，實際為 21600×10800，1/60°，57.6 MB)。切成 512×512 分塊以 JPEG (品質 92) 重新壓縮，存檔 44.9 MB；
    轉檔時暫時需要約 1.5 GB 記憶體。
    Wikimedia 有時會限流 (HTTP 429 Too many requests)：可稍後再試，或先用瀏覽器下載該檔，再執行
    `python -c "from satima import bluemarble; bluemarble.build_bm2002('bluemarble/BlueMarble_2002_land_shallow_topo_2km.tif', r'下載的檔案.jpg')"`。
* **與 NICT 網站底圖的比對** (NICT 的 `FULL_24h/BlueMarble`，全圓盤 4d 層與東亞 10d 層)：
  * NICT 用的是 Blue Marble 2002：陸地與海洋的平均色差都在 1 個色階內 (陸地 −0.2、海洋 R +0.9 / B −1.0)。
    NASA (2004 版) 則陸地暗 11\~26 個色階、海洋藍色暗 33 個色階。
  * 逐像素的紋理仍有差異 (陸地 rms 約 15 個色階)：NICT 的底圖相對真實位置偏了約 2\~3 個像素 (東亞 4\~6 km，往西北)，
    平移後 rms 降到 9。Satima-GEO 的 Blue Marble 2002 與 GSHHS 海岸線對齊在半個像素 (約 1 km) 內，所以不跟著 NICT 偏移；
    雲則兩者都是依衛星導航放在正確位置。
  * NASA 存放原始 Blue Marble 檔案的 eoimages.gsfc.nasa.gov 測試時連不上，所以改從 GIBS 與 Wikimedia Commons 取得。
* 預設不畫海岸線、國界與經緯線 (`lines = false`，`[bluemarble]`、`[nictmap]` 各自設定)。

### 節省記憶體 (`[memory] low_memory`，預設開啟)

輸出影像分條處理：每一條各自算行列號、區塊平均、取樣、真實色合成、上色，處理完即釋放，
只保留最後的 8 位元影像。每個波段的原始資料仍只讀一次 (讀取範圍以稀疏取樣估計，某一條超出時另外補讀)，
所以讀取器不必改，結果與整張處理逐像素相同 (已比對 MERC / CYL / LCC / FD、真實色、色表、各衛星)。

Himawari 全圓盤 (Python 程序峰值，16 GB 電腦實測)：

| 產品 | 整張處理 | 分條處理 |
|---|---|---|
| 真實色 8000x8000 | 9.4 GB | 4.3 GB |
| B03 8000x8000 | 6.2 GB | 2.9 GB |
| B01 原始解析度 11000x11000 | 6.3 GB | 2.5 GB |

全圓盤仍需讀入整張原始影像 (B03 22000x22000 的 float32 約 1.9 GB)，這部分不受分條影響。

### 多核心 (`[memory] workers`，預設 auto) 與 numba 加速 (`[truecolor] engine`)

* **各條平行**：輸出影像分條後，由多個執行緒同時處理 (satima/core/parallel.py)。numpy、numba、zlib、bz2 計算時會釋放 GIL，
  所以用執行緒就能多核心，不必複製大陣列。結果與單核心逐像素相同。
* **numba 核心**：真實色的逐像素合成融合成單一迴圈 (truecolor/kernel.py)，包括太陽/衛星角度、Rayleigh 查表內插與校正、
  綠光、銳化、淡出、拉伸。區塊平均、取值、全圓盤行列號與經緯度也改用 numba 核心 (_resample_nb.py)。
  * numpy 版受記憶體頻寬與 GIL 限制，8 核心只快約 2 倍；numba 版幾乎隨核心數線性加速。
  * 與 numpy 版的差異只來自浮點運算順序：色階差最多 1\~2，只佔 0.001%\~0.08% 的像素。
  * 第一次執行要編譯 (約 10\~20 秒，結果快取在 `__pycache__`)。沒有 numba 時自動改用 numpy。
* **Himawari 分段平行解壓**：同一波段的 10 個 HSD 分段同時解壓 (bz2) 與定標。
* **存檔**：
  * matplotlib 只在透明畫布上畫一次海岸線與標註，再與影像疊合 (basemap、cartopy 都是)。原本 `savefig(bbox_inches='tight')`
    會畫兩次，每次都把整張影像重新取樣。
  * PNG 改以多核心壓縮 (satima/core/pngio.py，同 pigz 的分段壓縮)。
  * 與原本的輸出只在線條重疊處偶有 1\~3 個色階差。
  * **cartopy 原本由 matplotlib 繪製影像，會偏移像素**：GeoAxes 的座標轉換是非仿射的，matplotlib 畫影像時會再重新取樣一次。
    即使畫布與影像一樣大，GOES-19 麥卡托也有 54% 的像素被換成上下左右相鄰的像素 (偏移約半格，取整後 0 或 1 格)。
    * 色階邊緣因此像被打散，看起來比較不像階梯，但像素位置最多錯 1 格 (約 2 km)。
    * 現在 cartopy 也走疊合：影像逐像素不變，與 basemap 相同。線條與 matplotlib 繪製的最多差 1 個色階 (0.03\~0.06% 的像素)。
    * 全圓盤地球外的影像本來就是透明的 (背景色)；地球邊緣依 Satima-GEO 逐像素的可見性判斷，不再被 cartopy 的圓盤多邊形多裁掉幾個像素。
* **圖高改為格點的整數列數**：以前圖高是浮點數，matplotlib 會把影像縮放約 0.4 像素，最近鄰取樣時漏掉一列。
  現在影像、海岸線、格點三者一對一對齊。
  * 英吋 x dpi 的浮點誤差 (例如 1691.9999999999998) 會被 Agg 截成少一列，所以圖的寬高各多給 0.001 像素。

Himawari-9 真實色 500 m 原始解析度全圓盤 (22000x22000，手上只有 4/10 段，約全圓盤的 40%)，
i7-9700 (8 核) / 16 GB 實測：

| | 改版前 | 現在 |
|---|---|---|
| 讀取與計算 | 467.5 秒 | 41.3 秒 |
| 畫圖與存檔 | 103.5 秒 | 9.2 秒 |
| 合計 | 571 秒 (9.5 分) | 50.6 秒 |
| 記憶體峰值 | 9.8 GB | 8.3 GB |

完整全圓盤 (10 段) 約為 2.5 倍，改版前約 24 分鐘，現在約 2 分鐘。
每多一個執行緒約多用數百 MB 記憶體；記憶體不足時可把 `workers` 調小。

## 新增衛星

在 `sensors/geo/` 新增 `<感測器>-<格式>.py` (檔名要有 `-`；沒有 `-` 的 .py 視為輔助模組，不會被載入)，
寫一個繼承 `GeoReader` 的類別，最後一行 `READER = 類別名稱`。只需要寫和衛星有關的部分：

| 項目 | 內容 |
|---|---|
| `NAME` `SENSOR` `DESCRIPTION` | 格式名稱 (= 檔名)、感測器名稱、`--list` 的說明。資料列表依序顯示衛星、`SENSOR`、格式 (`NAME` 中 '-' 之後的部分，例如 HSD、GVAR；場景可用 `Scene.format` 自訂，例如 NSMC GVAR)、時間、範圍 |
| `PATTERNS` | 檔案搜尋樣式 (glob) |
| `CHANNELS` | `{資料中的代號: Channel(代號, 中心波長 µm, 星下點解析度 km, 'vis'/'nir'/'swir'/'wv'/'ir', 顯示名稱, 功能代稱)}`；顯示名稱用於選單、標註、檔名 (例如 `BAND-13`)，功能代稱見下表 |
| `TRUECOLOR` | (可省略) `{'TrueColor': dict(red=, green=, blue=, nir=, ir=, lo_blur=)}`；沒有綠光時 `green=None` (ABI 合成綠光) |
| `scenes(folder)` | 把檔案分組成 `Scene` (同時次、同範圍的所有波段) |
| `nav(scene, channel)` | 回傳導航物件 (見下) |
| `read(scene, channel, (r0, r1, c0, c1))` | 讀取視窗並校正：反射率 0\~1 / 亮溫 °C，無效值 NaN |
| `satellite_position(scene)` | (星下點經度, 衛星高度 km)，真實色用 |
| `platform(scene)` | (可省略) pyspectral 的衛星名稱，預設為 `scene.satellite` |

### 波段名稱

* 選單、圖片標註、檔名一律 `BAND-01`\~`BAND-16` (各機構的頻道編號；MTG 依波長順序)。
  MTSAT 官方沒有頻道編號，沿用 JMA 的 `VIS`、`IR1`\~`IR4` (IR1 10.8、IR2 12.0、IR3 6.75 水氣、IR4 3.75 µm)；
  COMS 沿用 KMA 的 `VIS`、`IR1`、`IR2`、`WV`、`SWIR`；
  GMS-1\~4 為 `VIS`、`IR`，GMS-5 / GOES-9 為 `VIS`、`IR1`、`IR2`、`IR3` (水氣) (JMA 名稱)。
  MTG 的 FDHSI 與 HRFI 分成兩筆資料，標註與檔名再接上資料集標記：`BAND-14 HRFI` / 檔名 `BAND-14_HRFI`、
  `TrueColor FDHSI` / 檔名 `TrueColor_FDHSI`。HRFI 波段沿用同波段的編號 (BAND-03、BAND-08、BAND-09、BAND-14)。
* 注意：MTG 多了 0.9 µm (BAND-05)，之後的編號與其他三顆錯開，例如 MTG 的 BAND-13 是 9.7 µm 臭氧、
  BAND-14 才是 10.5 µm 窗區 (其他三顆的 BAND-13 是窗區)。選單上都有波長，可對照下表。
* 內部 (讀檔、pyspectral) 仍用資料中的代號：B13、C13、IR105、ir_105_hr。
* 功能代稱 (選單用，英文) 沿用 NOAA ABI 官方代稱，其他衛星依功能對應，新增衛星請沿用同一組寫法：

| 代稱 | GOES ABI | Himawari AHI | GK-2A AMI | MTG FCI |
|---|---|---|---|---|
| Blue | BAND-01 0.47 | BAND-01 0.47 | BAND-01 0.47 | BAND-01 0.444 |
| Green | — | BAND-02 0.51 | BAND-02 0.51 | BAND-02 0.51 |
| Red | BAND-02 0.64 | BAND-03 0.64 | BAND-03 0.64 | BAND-03 0.64 (HRFI 有) |
| Veggie | BAND-03 0.86 | BAND-04 0.86 | BAND-04 0.86 | BAND-04 0.865 |
| Near-IR WV | — | — | — | BAND-05 0.914 |
| Cirrus | BAND-04 1.37 | — | BAND-05 1.38 | BAND-06 1.38 |
| Snow/Ice | BAND-05 1.6 | BAND-05 1.6 | BAND-06 1.61 | BAND-07 1.61 |
| Cloud Particle Size | BAND-06 2.2 | BAND-06 2.3 | — | BAND-08 2.25 (HRFI 有) |
| Shortwave Window | BAND-07 3.9 | BAND-07 3.9 | BAND-07 3.83 | BAND-09 3.8 (HRFI 有) |
| Upper-level WV | BAND-08 6.2 | BAND-08 6.2 | BAND-08 6.24 | BAND-10 6.3 |
| Mid-level WV | BAND-09 6.9 | BAND-09 6.9 | BAND-09 6.95 | — |
| Lower-level WV | BAND-10 7.3 | BAND-10 7.3 | BAND-10 7.34 | BAND-11 7.35 |
| Cloud-top Phase | BAND-11 8.4 | BAND-11 8.6 | BAND-11 8.59 | BAND-12 8.7 |
| Ozone | BAND-12 9.6 | BAND-12 9.6 | BAND-12 9.63 | BAND-13 9.66 |
| Clean LW Window | BAND-13 10.3 | BAND-13 10.4 | BAND-13 10.4 | BAND-14 10.5 (HRFI 有) |
| LW Window | BAND-14 11.2 | BAND-14 11.2 | BAND-14 11.21 | — |
| Dirty LW Window | BAND-15 12.3 | BAND-15 12.4 | BAND-15 12.36 | BAND-15 12.3 |
| CO2 LW | BAND-16 13.3 | BAND-16 13.3 | BAND-16 13.31 | BAND-16 13.3 |

波長 (µm) 取各機構官方表：NOAA、JMA、KMA GK2A Fact Sheet、EUMETSAT。

### 導航

導航依投影方式選擇，不必把所有投影塞進同一套程式：

* 靜止衛星投影：`GeosNav(...)`，或由 CGMS 參數 `GeosNav.from_cgms(shape, cfac, lfac, coff, loff, ...)`
* 其他規則格點投影 (蘭伯特、麥卡托、等經緯度…)：`ProjGridNav(shape, proj 參數, x0, dx, y0, dy)`
* 只有逐點經緯度陣列的舊衛星：`LonLatNav(lon, lat)`
* 都不適用時繼承 `Nav`，實作 `to_rowcol()`、`to_lonlat()`、`key()` 即可
  (例如 GMS 的自旋掃描導航 `sensors/geo/vissr_nav.py` 的 `VissrNav`)

### GMS (VISSR) 導航

* 依 JMA MSC Technical Note No.23 的 MGIVSR，每一點依自己的觀測時刻由姿態、軌道預報計算 (不內插)；
  反算直接解掃描角 (原程式的反算以 VMIS 小角度近似)，觀測時刻反覆修正到收斂。
  與舊版 Cython (gms_nav.pyx) 的正算差 0.0000 km、與檔頭簡易座標表的星下點行列號差 0.01。
* `[navigation] engine`：`numba` (預設，不需編譯)、`cython` (`python sensors/geo/build_vissr_nav_cy.py` 編譯一次)、
  `numpy`；三者結果相同 (差 < 1e-10 像素)。2000x2000 的輸出約 0.4\~0.5 秒。
* 行號用 LCW 的行號 (由 1 起算，同 JMA)；satpy 多加了 1 行。GOES-9 紅外的中心像素號要減去檔頭的
  「Pixel difference of VISSR center」(1.5)，以地球邊緣驗證 (不減時東西偏 -1.54 像素)。
* `[navigation] earth_edge_correction` (預設開)：MSC Technical Note No.26 的地球邊緣修正，以每條掃描線 LCW
  的實際地球東西邊緣與預測邊緣比較，最小平方求南北、東西偏移與大氣寬度。西澳晴空海岸線比對：
  修正前偏 2\~5 km，修正後約 1.5\~3 km。GMS-AREA 沒有逐線資料，不做此修正。
* 完整檔頭：`reader.header(scene, 'IR1')` 回傳所有紀錄 (名稱與欄位同舊版 _VISSR_*.py)，含每條線的 Line_Control_Word。

### GMS-1 磁帶映像 (GMS-TAPE)

FGGE 時期 (1978~1979) 的 GMS-1 VISSR 磁帶映像，一卷磁帶一個檔案 (檔名為卷號：`V` 開頭共 6 個字、沒有副檔名，
例如 `VS2303`、`V10585`)，可見光與紅外都在同一個檔案。

* 內容就是 JMA VISSR 格式 (同 GMS-VISSR)，只是每個 4-byte word 存成一個 36-bit word：
  整數 36 位元；實數為 27 位元尾數 + 9 位元指數 (倍精確度 63 + 9)；文字為 9 位元的 EBCDIC；
  紅外每像素 9 位元、可見光 6 位元連續存放。紀錄前有 4 bytes 的長度，開頭與結尾是 EBCDIC 的磁帶標籤。
  解碼方式取自隨磁帶提供的 vissr.c / gms_area.c，欄位位置與 JMA「FORMAT OF VISSR ARCHIVE DATA」一一對應。
* 定標與導航都用磁帶自己的參數：紅外亮溫表、可見光各感測器的反照率表、座標轉換參數、姿態與軌道預報；
  導航、地球邊緣修正與 GMS-VISSR 是同一套程式 (`sensors/geo/vissr_common.py`)，行號用 LCW 的行號。
* 只取影像掃描線 (Line name = 1)；黑體、測試、註記、灰階的紀錄不用。
* 可見光缺行：這批磁帶用備援感測器組，其中可見光第 2 感測器故障 (Mode Block 的 `VIS redundant 2ch` = 0)，
  每 4 條掃描線缺 1 條。`[visible] fill_missing_lines = true` (預設) 以上下兩條的反照率平均補上；
  `false` 則保留為無資料 (圖上是底色的橫線)。補的是各自定標後的反照率，不是第 2 感測器的觀測值。
* 完整檔頭：`reader.header(scene, 'IR')`，名稱與欄位同 GMS-VISSR；數值保留磁帶的精確度 (整數 int64、實數 float64)，
  另有 `Labels` (磁帶標籤)。
* 第 3 個參數區塊：位置對應 JMA 文件的 ⑧ (狀態表的 Data segment 8) 與 ⑨ β 角取樣，但文件對這兩段的內容都寫「Omitted」，
  沒有欄位定義。磁帶上只有 ⑧ 有內容 (運用模式、行號範圍、感測器故障表、畫面常數等 S/DB 的畫面設定，多數與 Mode Block 重複)，
  定標與導航都用不到，以原始 word 保留 (`Data_segment_8`)；⑨ 全為 0。

### GOES-8\~15 (GVAR) 讀取與導航

* 波段名稱 `BAND-01`\~`BAND-06` (= GVAR 頻道 1\~6；GOES-8\~11 有 1\~5，GOES-12\~15 有 1、2、3、4、6)。
* 原始 GVAR 的容錯 (CRC、Data Valid 旗標錯誤但影像正常的情形很多，都不當作丟棄的依據)：
  區塊頭 3 份逐 byte 多數決 (第 1 份壞也能辨識)；字數錯誤時從上一個區塊重新尋找下一個區塊頭；
  NSMC 檔依固定 26624 bytes 的槽位讀取並自動偵測起點；Line Documentation 的頻道、偵測器等欄位
  以「區塊內第幾筆紀錄」的全檔眾數修正；每筆紀錄依自己的掃描序號 RISCT 定位行號，並與 Block 0 交叉驗證；
  導航參數取所有 Block 0 的 O&A 參數眾數。
* 導航同 McIDAS nvxgvar (與舊版 le2ll 差 0.0000 km)；IMC 關閉時 (例如 2003\~2005 在 155°E 的 GOES-9)
  依每條線的掃描時間計算 (舊版與 McIDAS 用單一時間，一張影像內最多差約 18 km)。
  紅外每條線的中心在其 4 條可見光線的中間 (+1.5 線；舊版與 McIDAS 未加，10.7 µm 海岸線比對證實 +1.5 正確)。
* 校正：可見光 = bias + gain x 計數值 (+ 二次項)，再乘 albedo 係數 (舊版公式為 (計數值 + bias) x gain，太空偏亮)；
  紅外依 conversion_table/GOES 的表。
  GOES-14/15 的 GVAR 沒有可見光係數，用 NOAA OSPO 公布的各偵測器發射前係數 (Weinreb & Han 2009；GOES-14 不再借用 GOES-15)。
  `[visible] post_launch_calibration` (預設開) 再乘以 NOAA 的發射後校正 (感測器衰減)：GOES-11/12/13/15 依 STAR 每月修正值
  (conversion_table/GOES/GOES_Imager_Vis_OpCal.txt，次月第二個星期二起生效)，GOES-8/10 依 OSPO 線性衰減式，
  GOES-14 固定 x 1.0695 (2009/12 發射後測試以 MODIS 比對，只見於 PLT 週報摘要)；GOES-9 沒有官方值。
  GVAR Block 0 的可見光係數一直是發射前值 (GOES-13 2013\~2017 不變)，所以不會重複修正。
  紅外表 (NOAA OSPO) 與 OSPO 公布的最新係數 (GOES-13 CH6 為 STAR 修正版、GOES-14/15 為 RevH + STAR) 相差 < 0.02 K；
  表只列 180\~340 K，高溫端 (3.9 µm 可到 343 K) 以表本身擬合的 T = a + b x Teff(n, R) 延伸，低溫端 (太空) 不延伸。
  偵測器依 GVAR 規格：LIDET 奇數 = 北 (表的 DET1)、偶數 = 南 (DET2)；LSIDE 1023 = side 2。
  GOES-12 CH6 (13.3 µm) 的 NOAA 係數沒有光譜響應修正，GSICS 與 IASI 比對約有 3 K 偏差，目前未修正。
* 海岸線比對 (晴空沙漠海岸)：可見光 1 km 以內；紅外與可見光的東西對準約差 0.4 個紅外像素 (依衛星而異)。

### SMS、GOES-1\~7 (GVAR 以前) 讀取與導航

* 一筆資料 = 同一衛星、時間的 AREA_VIS 與 AREA_IR；選單的儀器欄顯示 VISSR (SMS、GOES-1\~3) 或 VAS (GOES-4\~7)。
  波段：`VIS`、`IR` (VISSR 模式的紅外)、VAS 的 `BAND-01`\~`BAND-12` (依檔案的 filter map，例如 7、8、10)。
* 影像來源：VISR (每像素 1 byte 的 McIDAS 亮度值)、VAS (mode AA，1987/03 以前) 與 AAA (GOES-7)，2 bytes 最多 3 波段交錯，
  每行依波段表找該波段的位置。
* 校正：可見光為平方根編碼的 6 位元值 x 4，反射率 = (B/252)^2 (同型 GMS VISSR 的反照率表即此關係，rms 0.002；
  SMS / GOES 無絕對校正。舊版直接以 B/255 顯示，等同已做 gamma，要相同觀感可開 `[visible] gamma`)；
  紅外亮度值 B >= 176 時 418 - B，否則 330 - B/2 (K)；VAS mode AA 輻射 = (P - YSUBZ) x 2^(F + DELTAF - 15)；
  AAA 輻射 = (IAB2 x P/32 - IAB1) / 2^(15 - IFAB)；輻射 -> 亮溫 T = (FK2/ln(FK1/R + 1) - TC1)/TC2
  (舊版寫成 ((FK2 - TC1)/TC2)/ln(...)，偏暖約 0.3 K)。FK1\~TC2 沿用舊版的一組 VAS 係數 (GOES-4\~7 共用，沒有各衛星的值)。
* 導航同 McIDAS nvxgoes (VisAD GOESnav)：同一圈 (自旋) 的 8 條可見光線同一時刻，Kepler 軌道、自旋軸、姿態；
  反算直接解旋轉 (McIDAS 的反算以 beta 紀錄內插時間，與其正算並不互逆)，與正算互為反函數 (1e-7 線)。
  與舊版編譯的 nav.pyd 相差 0\~0.5 km (epoch 時間的秒數 McIDAS 四捨五入、舊版截斷)；舊版 Python 原始碼 nav.py
  的 ietimh、NUMSEN、掃描序號用了浮點除法，與 McIDAS 平均差約 11 km (舊讀取器實際載入的是 nav.pyd)。
* 紅外每圈取樣一次，AREA 紅外行號 = 該圈第一條可見光線，像素中心取 +3.5 線 (8 條可見光線的中央)。
  同時次可見光 / 紅外雲邊比對 27 區皆為 +3.25\~+8 線；東西向 1988 年 +7、1994 年 -3 像素，屬當年配準誤差，不修正。
* 這個年代的導航本身有數 km 誤差 (可見光海岸線：南北 +2\~+6 km、東西 -7\~+2 km，依影像而異)，照 McIDAS 公式不另修正。
* 導航區塊沒有軌道參數的檔案 (例如 1979 年的 GOES-1、GOES-3，McIDAS 也無法導航) 會列出原因後略過。
* 完整檔頭：`reader.header(scene, 'B08')` 回傳 AREA 目錄、導航區塊、AAA 校正區塊、每行前置資料 (_Data_Block) 與波段表。

### SMS、GOES-1\~7 原始資料 (Mode A / AA / AAA)

* 檔名 (NCEI 2024 重處理版，NOAA CLASS)：`goesNN.YYYY.DDD.HHMM.A.data`、`.AA-MSI.data`、`.AA-DS.data`、`.AAA-MSI.data`、
  `.AAA-DS.data`、`smsNN.YYYY.DDD.HHMM.A.data`。一個檔 = 一筆資料 (VIS 與 IR 在同一個檔)；資料列表的格式欄為 MODE A、MODE AA-MSI…
  (MSI = 多頻道影像、DS = dwell sounding)。
* 格式依 NOAA 的格式文件：Mode A 為 NCC 1978 的 VISSR Digital Archive Users' Guide 與 1982/10 的 Modified IR Documentation，
  Mode AA 為 Output Format Rev. B / C，Mode AAA 為 SFP 002 v2.5A (1987/09)。檔案中一筆記錄 = 一圈 (0.6 秒) 去掉同步碼後的資料：

  | | 記錄 | IR | VIS |
  |---|---|---|---|
  | Mode A | 99244 bytes | 文件 128 字 + 3818 個樣本，9 位元字 (1 個旗標位元 + 8 位元) | 8 段 x (文件 512 字 + 15288 個樣本)，6 位元 |
  | Mode AA | 109536 bytes | IR1、IR2 各：文件 16 + 影像 3822 + EOLC 8 + 共同文件 2048 字，10 位元 | 同上 |
  | Mode AAA | 109516 bytes | 兩個 IR 區塊各：文件 16 + 影像 3822 + EOLC 8 + 共同文件 512 + 格網 1024 + Mode A 文件 128…，10 位元 | 同上 |

  記錄偶有不完整 (AA-DS 樣本中一筆缺 VIS 第 1 段、之後的區塊整體前移；AAA-MSI 中一筆沒有 IR1 的損壞記錄；檔尾截斷)：
  有同名的 `.indx` (CLASS 附的區塊索引) 時依它逐區塊定位，可讀到不完整的記錄；沒有時以 VIS 段碼定位，只用完整的記錄
  (6 個樣本兩種方式的區塊位置完全相同，`.indx` 多讀回 3 筆不完整記錄)。`.indx`、`.meta` 都不是必要的。
  `.meta` (JSON) 的導航來源 (`navigation_source`，樣本皆為 native = 原始資料流的 O&A)、O&A 眾數 (與讀取器組出的相同)
  只放進 `header()`。
  Mode A 的 IR 第 3819\~3822 個位置是行尾標記，為 NaN。
* 波段：Mode A 為 `VIS`、`IR`；Mode AA / AAA 為 `VIS` 與 VAS 的 `BAND-01`\~`BAND-12`。每圈的濾鏡輪位置 (預測標頭) 決定該圈的波段：
  MSI 例如 8、7、8、7、8、10、8、10 輪替 (band 8 每兩圈、band 7 / 10 每四圈)；DS 在部分掃描步駐留、依序取各波段，
  同一波段、同一步有多圈時輻射取平均。校正圈 (submode 4) 與 frame code 不是 picture 的圈不使用。
* VIS 線號 (由 1 起算，掃描 N 的 8 條線 = (N-1) x 8 + 1 \~ + 8，同 McIDAS AREA)：
  AAA 的每段 VIS 文件 (第 9\~12 字) 直接記錄掃描數 N，段 j = 第 j+1 條線；AA 依 Rev. C 表 2：段 0\~3 = 掃描 N 的 V5\~V8、
  段 4\~7 = 掃描 N+1 的 V1\~V4；Mode A：同時次的 Mode A 與 Mode AA (GOES-6 1985/09/24 20:31) 逐段比對，
  Mode A 掃描 s 的第 k 段 = AA 掃描 s-1 的第 k 段 (400 條線全部吻合)，所以與 AA 相同 (線號 (s-2) x 8 + 5 + k)。
* IR 線位置 (掃描步，一步 = 8 條 VIS 線)：Mode A 為記錄的掃描數。VAS 每圈有上 (IR1)、下 (IR2) 兩個偵測器：
  小偵測器 (HgCdTe 小) 下 = 上 + 1 步；大偵測器 (HgCdTe 大、InSb) 各涵蓋 2 步，下 = 上 + 4 步。
  AAA 的 IR 文件直接記錄已修正上下偵測器位移的 adjusted earth count；AA 沒有這個欄位，以影像相關性量得相同的排列。
  DS 的小偵測器波段只涵蓋一半的掃描步，空缺以最近的線補上 (最多 2 步，同 McIDAS AREA 以最近一圈重複填線)。
* 校正：
  * VIS：反射率 = (計數值 / 63)^2 (同 GOES-PREGVAR)。
  * Mode A IR：NCC 1978 附錄 II，計數值 < 176 為 330 - c/2，否則 418 - c (K)。
  * Mode AA：輻射 = (32 x V - YZQ) x 2^(F + DELTAF - 15)；YZQ = 共同文件第 3 段的 Y_Z + Q，DELTAF 由 DELTA SCALE FACTOR 查表，
    逐圈、上下偵測器各一組 (即 McIDAS AREA 的 YSUBZ、RAWDELTAF；AREA 的 P = 32 x V)。
  * Mode AAA：輻射 = (AB_k(1) x V - AB_k(0)) / 2^(15 - FAB_k)；AB、FAB 由第 3 欄第 63\~112 字組成 (20 圈輪替)，
    k = 偵測器 / 濾鏡組合 1\~38 (表 5C1)。
  * 輻射 -> 亮溫同 GOES-PREGVAR。
* 導航：McIDAS 'GOES' 導航 (同 GOES-PREGVAR)，參數取自每圈文件中輪替的 O&A (20 圈一組，重複出現時取眾數)：
  VIRGS 的 Kepler 軌道 (EPY、EPH、SMA、ECC、INC、MA、AP、RAN，與同日 McIDAS AREA 導航區塊的值相同)、
  自旋軸 = SPRA / SPDC 在 13 小時區間內插後的反向、pitch = -ZETA、yaw = ETA (同日 AREA 的 pitch / yaw 即此值)。
  沒有 Kepler 元素的檔案 (例如 1977/11 以前的舊 O&A 格式) 列出原因後略過。
  照檔案內的導航參數計算、不另外修正；IR 照資料本身的掃描編號，比照 McIDAS 慣例 (線中心 = 該步第一條 VIS 線 + 3.5)，
  不修正配準 (原始 AA 與同日 McIDAS AREA 的 IR 地球邊緣位置相同，差 0\~1 個樣本)。
* 驗證 (PREGVAR 的 6 個樣本，2026/11)：
  * 校正：同時次的 GOES-6 Mode AA band 8 與 Mode A IR (兩條獨立的校正路徑) 亮溫差中位數 1.3 K (冷雲 0.3 K、暖地表 1.5 K；
    Mode A 溫度表本身的精度為 2\~4 K)。GOES-7 AAA 組出的 AB / FAB 與 McIDAS AREA 校正區塊的 IAB / IFAB 逐值相同
    (band 8：R = 0.24 x V - 10，同文件附錄 F)。
  * 導航 (可見光與海岸線比對，影像相對導航的位移；線、像素約 1 km)：GOES-6 Mode A / AA 偏北約 6 線、偏東 3 像素；
    GOES-7 AAA 偏南 2\~7 線、偏東 5 像素；SMS-2 偏北約 14 線、偏西 3\~14 像素 (該檔 O&A 本身的誤差，四個區域一致)。
    原始 O&A 導航與 1 小時前 AREA 的導航相差約 6.5 km (兩者的軌道解不同)。
  * IR 與 VIS 的配準 (未修正)：IR 相對同一標稱位置的 VIS，南北 Mode A、AA 偏北約 10 線，AAA 偏南約 2.5 線；
    東西偏西 6\~22 像素 (Mode A、AA)、偏東 10 像素 (AAA)，依衛星與格式而異，所以不以樣本的值修正。
* 完整檔頭：`reader.header(scene)` 回傳 {'mode', 'records' (各記錄的文件欄位陣列), 'offsets'、'ir_pos'、'vis_pos' (區塊位置)、
  'index_file'、'meta' (`.meta` 的內容)、'oa' (O&A), 'nav' (導航參數), 'calibration'}。

### MSG SEVIRI (Meteosat-8\~11)

* 波段 `BAND-01`\~`BAND-12` (= EUMETSAT 頻道 1\~12：VIS006 … IR_134、HRV；選單中 BAND-12 標示「HRV (High Resolution
  Visible)」，標註與檔名仍為 BAND-12)；資料列表顯示 FD / RSS / ROI 與 IODC。
* 三種格式結果相同 (數值與 satpy 差 < 1e-4，經緯度 0.0000 km；HRIT 以 satpy 的 Zenodo 公開樣本
  MSG4 2018/02/28 15:00 驗證，native 以 Meteosat-9 IODC 全圓盤與 Meteosat-11 RSS 驗證)。
* 影像第 0 列最南、第 0 行最東；2017/12 以前 (地球模型 1) 的資料依 EUMETSAT 說明修正 1.5 km 偏移。
* HRV (1 km) 全圓盤分下、上兩個視窗，各以自己的東界定位 (資料不補成 11136 行寬)。
  上視窗的排列實際從 epilogue 的 LowerNorthLineActual 那一行開始 (比 UpperSouthLineActual 早一行；
  以影像連續性驗證，同 satpy)。IODC 的 HRV 在兩視窗交界處有一段斜向的無資料區，原始資料即如此。
* 定標：`[seviri] calibration` = nominal (預設) / gsics (紅外) / meirink (可見光)；反射率含日地距離修正，
  亮溫依輻射種類 (有效 / 光譜) 用 EUMETSAT 的 α、β 或二次式。輻射為 0 的像素為 NaN (satpy 會給 -276 °C)。
* 壓縮 HRIT (-C_) 需 pyPublicDecompWT，或 `[seviri] xrit_decompress` 指定 xRITDecompress；只在記憶體 / 暫存目錄解壓。
* native 的 zip：列表只解壓檔頭；第一次讀影像時把 .nat 解壓到記憶體 (約 270 MB，只保留一個)。
* 舊版 prologue.py 的 np.bool 在 numpy 1.24 以後已移除，舊版 MSG_HRIT 讀 prologue 會出錯；新版已修正。

### MFG MVIRI (Meteosat-2\~7)

* 波段 `VIS` (5000x5000)、`IR` 11.5 µm、`WV` 6.4 µm (2500x2500)。
* FCDR：導航為 18 度視場的標準靜止衛星投影 (即舊版 regfeo)，投影經度取自檔名 (E5700 = 57.0)；
  與 satpy 相同 (經緯度 0.0000 km、亮溫差 < 1e-4)，影像方向另以檔內的衛星天頂角驗證 (南到北 0.2 度、反向 8 度)。
  VIS 的 FCDR 值是雙向反射因子 (已除以 cos 太陽天頂角)，乘回 cos(太陽天頂角) 成為與其他衛星一致的 TOA 反射率。
  部分檔案的時間變數比檔名晚一天，以檔名為準。
* MVIRI-NC (舊版 MFG.py 的格式) 依舊版移植，以檔案的經緯度定位、亮溫用 EUMETSAT 的 A、B 係數；尚無樣本驗證。

## 驗證 (2026/11)

* 新讀取器與舊版 (已修正 1-based 的) 讀取器比較：AHI、AMI、FCI 導航 0.0000 像素、數值一致；
  ABI 與 satpy 的導航差距 < 2 m (舊版 GOES16_NC 因 float32 計算有 0.01\~0.04 像素誤差)。
* 真實色各波段在輸出格點上的區塊對應 99.98\~100% 一致，聯合擬合位移為 0。
* GMS-VISSR / GMS-AREA：同一時次 (GMS-1 1979/10/11 02:33) 兩種格式的亮溫差 < 0.26 K；
  GMS-AREA 重畫 IR197901071533.png 結果相同。
* GMS-TAPE (VS2283、VS2288、VS2303、V10585 四卷)：
  * 與先前轉成 JMA 格式再用 GMS-VISSR 讀的同一卷 (VS2303) 比較：可見光逐像素相同、紅外差 < 0.00003 K、
    經緯度差 < 0.00007° (轉檔時實數縮成 32 位元的捨入)。
  * 與同時次的 McIDAS AREA (GMS_VIS_1_1979_282_2333.area，VS2288) 比較：有觀測的線 AREA 亮度 = 計數值 x 4 (全部相同)，
    缺的線 AREA 也是上下兩條的平均；AREA 的 GMSX 導航參數 (中心行 5064.5、中心像素 6688.5、行號起點) 與磁帶相同，
    同一行列號的經緯度平均差 0.55 km (GMSX 的步進角只存 4 位有效數字)。
  * 地球邊緣修正量：南北 -0.2~+1.3 行、東西 +0.2~+1.6 像素 (紅外)；澳洲海岸線比對吻合。
  * GMS-VISSR 抽出共用程式後，GMS-1、GMS-3、GMS-5、GOES-9 共 11 個波段的檔頭、導航、影像值與修改前逐位元相同。
* 尚未有樣本驗證：GK-2A 蘭伯特投影 (lc) 的區域資料。
