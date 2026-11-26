# Satima：極軌成像儀 (PO)

> 這份是 **PO (極軌成像儀)** 的詳細說明，改寫自原本的 Satima-PO。整合版的執行方式、資料列表、設定檔分層
> (`Satima.toml` + `Satima-PO.toml`)、介面語言、圖片模式 (native / dpi) 與極軌的接續判定請看 [README.md](../README.md)。
> PO 是三類之中改動最大的：與原本 Satima-PO 的差異列在 README 的「與原本三個程式的差異」。

VIIRS、MODIS、AVHRR、METimage、OLCI、SLSTR、SGLI、FY-1C/1D MVISR 的單頻道 (反射率、亮溫配色表、DNB) 與真實色。
8 支原本各自獨立的腳本先整合成 Satima-PO，再併入 Satima；各格式的校正與演算法都留在自己的腳本裡。
FY-1C/1D MVISR (`MVISR-L1B`) 由 `SatimaV3-FY1D-NSMC.py`、NSMC 格式的 NOAA AVHRR (`AVHRR-NSMC`) 由 `SatimaV3-NOAA-NSMC_accnumpy.py` 改寫，
2026/11 加入 (見「FY-1C/1D MVISR」、「NOAA AVHRR/3 (NSMC 格式)」)。

## 使用方式

```
python Satima.py -t po               只列極軌成像儀的資料
python Satima.py -f VIIRS-HDF5       只搜尋某個格式
```

```
1   S-NPP     VIIRS    HDF5   2023/05/26 04:03-04:10Z   NOAA       22 bands   5 檔
2   METOP-B   AVHRR3   NAT    2026/04/12 09:31-11:13Z   EUMETSAT   IR/VIS
輸入資料編號 (多筆可用 (3,5) 或 [1,3,5])：1
Projection: 1.MERC  2.CYL  3.LCC
輸入投影編號 (Enter = MERC)：
地圖中心 (緯度 經度)：18 121
範圍 (中心往外的緯度、經度各幾度；只輸入一個 = 兩者相同)：5 8
(各格式的頻道 / 產品選單)
選擇色表 (可複選，0 = 全部)：0
圖片寬度 (…；Enter = 原始取樣間距 …)：
圖片標題：Haikui
```

* 列表的一列是**一段連續的觀測**：VIIRS (約 85 秒一檔)、METimage (約 1 分鐘)、MODIS (5 分鐘) 這類短片段，同衛星、
  時間相接的自動接成一筆，備註欄顯示檔數與缺檔數。AVHRR 這類整圈的檔各自一列。判定的規則見 README「極軌的接續判定」。
* VIIRS 不再是一個波段一列：選了資料之後，選單列出這一段每個片段都有的頻道，以及 M3 / M4 / M5 (加 I1) 都齊全時的真實色。
* 地圖範圍：`Satima-PO.toml` 預設 `[map] range_input = "center"` (先問中心、再問範圍，同原本的 Satima-PO)；
  輸入 4 個數字則當成 南 北 西 東。中心經度可用 -180\~180 或 0\~360，跨換日線的範圍也能直接畫。
* **中心時間**：圖上與檔名的時間是最接近中心點的像元的實際觀測時間。中心點 = 你輸入的中心；輸入四邊時為範圍中點。
* 紅外線頻道會印出中心點附近的最高亮溫 (`[eyetemp]`)，所以畫颱風時把中心設在颱風眼。
* 色表、真實色樣式、DNB 縮放都可以複選，`0` = 全部。
* 圖片標註與檔名的頻道名稱同 GEO：`BAND-04`、`BAND-03A` (2026/11 以前是 `B04`、`B03A`)；VIIRS (`I05`)、SGLI (`VN11`)、METimage 的頻道本名不變。
* 標題 (原本的 Picname) 預設會詢問 (`Satima-PO.toml` 的 `[title] ask = true`)，可留空；會填入檔名的 `{title}`。PO 預設 `[title] draw = false`，只用於檔名、不畫在圖上；設成 `true` 則也畫在圖片左下角 (共同預設，GEO、MW 同)。

## 設定檔重點

| 項目 | 設定 | 檔案 |
|---|---|---|
| 資料 / 輸出目錄 | `[general] data_dir`、`output_dir`；相對路徑以 `Satima.py` 所在的資料夾為準；輸出目錄可帶欄位，如 `plot/{sensor}/{time:%Y%m%d}` | 共同 |
| 海岸線 | `[coastline] draw / color / linewidth`；解析度在 `[map] resolution` | 共同 (PO 覆寫成綠色 0.5) |
| 標註 | `[annotation] draw = false` 關閉；或個別 `[[output.image]] annotation = false`；文字樣板 `[annotation] text` | 共同 |
| 檔名 | `[output] filename`，空欄位自動省略；PO 另有 `{res}` (MODIS、SGLI) | 共同 |
| 圖片模式 | `[image] mode = "native"` (原始解析度) 或 `"dpi"` | 共同 |
| dpi | 每個 `[[output.image]]` 一張圖：`dpi`、`suffix` (只用於 `mode = "dpi"`) | 共同 |
| 底色 | `[image] background` (掃描帶以外)；OLCI、SGLI 在 `[readers.<格式>.image]` 設成白色 | 共同 / PO |
| 範圍提醒 | `[general] range_warning`：半寬超過此值時提醒記憶體用量 | PO |
| 長檔的合併 | `[swath] merge`、`merge_gap` | PO |
| 颱風眼溫度 | `[eyetemp] enable / radius_deg / legacy_counts` | PO |
| 色表選單 | `[colormap] ir / wv`，可增減或調整順序 | PO |
| 各格式個別設定 | `[readers.<格式>.<節>]`，例如 `[readers.MODIS-HDF4.general] data_dir = "D:/modis"` | 任一 |

預設檔名範例：`Haikui.JPSS-1.VIIRS.I5.BD.202608190028Z.png`、dpi 模式的第二張 `S-NPP.VIIRS.I1.202305260403Z.500DPI.png`

## 取樣與圖片大小

掃描帶的每個取樣點自帶經緯度。兩種圖片模式的畫法相同 (`satima/po/canvas.py`)，只有圖片大小的決定方式不同：

| | `mode = "native"` | `mode = "dpi"` |
|---|---|---|
| 圖片大小 | 地圖中央一個像素 = 原始取樣間距 (可用 `[image] size` 放大整數倍或指定寬度) | `[[output.image]]` 的 dpi x matplotlib 預設的英吋數 (方形範圍 1000 dpi 約 3700 像素高) |
| 輸出張數 | 一張 (`[image] suffix` 填入檔名的 `{suffix}`) | 每個 `[[output.image]]` 一張 |
| 取值 | 最近鄰：每個像素取最近的有效取樣點，數值不內插 | 同左 (不再使用 pcolormesh) |

- **原始取樣間距**由這次所畫資料的經緯度量出 (範圍內沿掃描、沿軌間距較細者)，不需各格式另外設定；真實色銳化後的
  250 m / 375 m 也自動適用。地圖範圍不含星下點時，取到的是範圍內最細的間距 (掃描帶邊緣的像元本來就比較大)。
- `[image] size`：`"native"`、`"x2"`\~`"x8"` (整數倍) 或固定寬度 (像素)。寬度不足 `min_size` (預設 2000) 時自動放大成
  整數倍且至少 `min_scale` 倍 (預設 4)，長邊超過 `max_size` 時改用放得下的最大整數倍；`min_size = 1` 可關閉自動放大。
- 輸出格點是共用地圖 (MERC / CYL / LCC，WGS84) 的規則格點。於地圖平面上以 KD-tree 找最近的取樣點：麥卡托與 LCC 是
  正形投影，平面上最近即實際上最近；等距圓柱先把東西向乘上中央緯度的 cos，近似等向 (緯度範圍很大時邊緣會略有偏差)。
- 離最近取樣點超過該處約半個取樣格的像素視為掃描帶外，顯示 `[image] background`。資料無效的取樣點不參與取樣，所以
  bow-tie 重疊、VIIRS bow-tie deletion 處由相鄰掃描的有效像元補上；星下點附近掃描與掃描之間的空隙
  (軌道降低後的 Aqua MODIS) 也會補上。缺檔的缺口超過距離上限，留空不補。
- 輸出像素比取樣粗時 (dpi 模式畫大範圍、或指定較小的寬度) 等於抽樣。
- 選「0 = 全部色表」時取樣只算一次；dpi 模式輸出多張時每張各取樣一次，讀取器的讀檔與校正只做一次。
- LCC 的圖框是涵蓋整個經緯度範圍的最小矩形，四角會超出輸入的範圍；讀取器裁切資料時用的是圖框實際涵蓋的經緯度範圍。

## 格式一覽

| 腳本 | 衛星 / 感測器 | 檔案 | 產品 |
|---|---|---|---|
| VIIRS-HDF5 | S-NPP、NOAA-20/21 VIIRS | `SV*_ops.h5` + GMTCO / GITCO / GDNBO (或 GMODO / GIMGO)；CLASS 的合併檔也可讀 | 單頻道、DNB (4 種縮放)、TrueColor-M / TrueColor-I |
| MODIS-HDF4 | Terra、Aqua MODIS | `MOD02*/MYD02*` + MOD03/MYD03 | 熱紅外、可見光、真實色 (RGB / CORR / WV) |
| AVHRR-NAT | Metop-A/B/C AVHRR/3 | `AVHR_*_1B_*.zip` (EPS native) | CH1–CH5 |
| AVHRR-NC | Metop-A/B/C AVHRR/3 | `W_XX-EUMETSAT-...METOP*.nc` | CH1–CH5 |
| METimage-NC | Metop-SG METimage (VII) | `*VII-1B*.zip / .nc` | 20 頻道、真實色 (General / CIRA / Worldview) |
| OLCI-SAFE | Sentinel-3 OLCI | `S3*_OL_1_*.SEN3` (zip / tar / 資料夾) | Oa01–Oa21、真實色 (RGB / CORR / WV) |
| SLSTR-SAFE | Sentinel-3 SLSTR | `S3*_SL_1_RBT*.SEN3` (zip / tar / 資料夾) | S1–S6 反射率、S7–S9 亮溫 |
| SGLI-HDF5 | GCOM-C SGLI | `GC1SG1_*_1BSG_*.h5` | VN01–VN11、SW01–SW04 反射率、TI01–TI02 亮溫、真實色 (RGB / CORR / WV) |
| MVISR-L1B | FY-1C/1D MVISR | `FY1?_AVHRR_HRPT_L1_*.L1B` (NSMC HRPT 1B) | CH1–CH10：CH3–CH5 亮溫，其餘反射率；真實色 (RGB / CORR / WV) |
| AVHRR-NSMC | NOAA-15~19 AVHRR/3 | `NA??_AVHRR_HRPT_L1_*.L1B` (NSMC HRPT 1B) | CH1–CH5 (CH3A / CH3B 依掃描線切換) |
| AVHRR-L1B | TIROS-N~NOAA-19、MetOp AVHRR/1、/2、/3 | NOAA Level 1b (POD / KLM)：CLASS 的 `NSS.GHRR/LHRR/HRPT/FRAC.*`、千葉大學 HRPT；**依內容辨識，不看檔名** | CH1–CH5 (AVHRR/3 的 CH3A / CH3B 依掃描線切換；AVHRR/1 沒有 CH5) |
| AVHRR-L0 | TIROS-N~NOAA-19 AVHRR/1、/2、/3 | HRPT 原始框 (Level 0)：接收站收錄的 `.raw` / `.raw.xz` 等，可直接放 `.tar`；**依內容辨識**；定位需要 TLE | 同 AVHRR-L1B |

SGLI 的 VN07 / VN10 是高增益 (海色) 頻道，雲上約在反射率 0.16 / 0.13 就飽和，VN08 約 0.49 飽和；畫雲請用 VN11 (868.5 nm)。

各格式的起訖時間：VIIRS、METimage、AVHRR-NAT、OLCI、SLSTR 取自檔名；MODIS 檔名只有起始時間，以標稱的 5 分鐘計；
AVHRR-NC、SGLI 的檔名沒有結束時間，列表只顯示起始時間 (間隔 10 分鐘內的檔仍會接成一段)。AVHRR-NAT、AVHRR-NC 一次只讀一個檔。
MVISR、AVHRR-NSMC、AVHRR-L1B 的起訖時間取自檔案的表頭 (精確到毫秒)。AVHRR-L0 的起訖時間取自框的時間碼與框數。

### VIIRS 的檔案

* 檔名 `<產品[-產品…]>_<衛星>_d<日期>_t<起>_e<訖>_b<軌道>_c<產製時間>_<來源>_ops.h5`。一般的檔只有一個產品 (`SVM15`、`GMTCO`…)；
  NOAA CLASS 的合併檔在檔名前段列出所含的產品 (`GMTCO-SVM03-SVM04-SVM05_…`、`GDNBO-SVDNB_…`)。兩種可以混用：
  程式登記每個檔含有哪些頻道與定位，各頻道與定位各自指向含有它的檔。HDF5 內部的群組名稱兩種檔相同。
* 定位：M 頻道 `GMTCO` → `GMODO`，I 頻道 `GITCO` → `GIMGO`，DNB `GDNBO`。箭頭左邊是地形校正版，找不到才用右邊，
  並在畫面上註明。真實色的 Rayleigh 校正在定位檔沒有高度資料時用海平面。
* DNB：`GDNBO` 有 `Latitude_TC` / `Longitude_TC` (地形校正) 時優先使用。
* 同一個片段的同一個產品有多個檔時用產製時間 (`c…`) 最新的。
* 範圍確定後只讀與地圖相交的片段 (每隔 64 個像元取一點判斷，範圍外加 2 度)，前後各多帶一個。
* 一段裡只有部分片段有的頻道也會列在選單 (會提示)；選了之後沒有該頻道的片段略過。真實色需要的頻道缺任何一個片段時會停下來並指出是哪一個。
* 有資料但找不到定位檔的頻道不列入選單；資料列表的頻道欄會標出缺哪一種 (例如 `0 bands，缺 GITCO`)。
* 只有一個產品可選時不詢問。
* DNB 的 ERF + 飽和修正 (縮放 3) 會統計載入範圍內的飽和比例，所以結果與載入了哪些片段有關 (現在是與地圖相交的片段加前後各一個)。
* **驗證狀況** (2026/11，S-NPP、NOAA-20、NOAA-21 的真實檔案)：I5 色表、I1、M5、DNB (ERF)、TrueColor-M、TrueColor-I 的中心時間與颱風眼溫度
  和原本的 Satima-PO 相同，圖目視一致；片段接續、只讀相交片段、一次選多段、缺定位檔的提示、`Latitude_TC` 都以真實檔案跑過。
  合併檔以真實資料組出來的檔 (`tests/pack_viirs.py`) 測試，輸出與分開的檔逐像素相同。
  **還沒有真實檔案可測**：CLASS 下載的原始合併檔與聚合檔、`GMODO` / `GIMGO`、M 頻道亮溫 (M12\~M16)、重複檔去除與缺檔標示
  (後幾項只用 `tests/fake_viirs.py` 的合成檔確認流程走得通)。

### FY-1C/1D MVISR

* 檔案：NSMC 的 HRPT 1B，每筆 28400 bytes 的定長紀錄 (big-endian)，格式見國家衛星氣象中心《風雲一號 (02) 批極軌氣象衛星數據格式》§3.1。
  檔名寫的是 AVHRR，儀器其實是 10 頻道的 MVISR。一個檔是一次接收的過境 (約 9\~20 分鐘)，選到時間相接的多個檔時依 `[swath] merge` 合併。
* 範圍確定後先以定位點挑出與地圖相交的掃描線 (範圍外多帶 3 度)，只解碼、內插這些掃描線。
* **定標**：每條掃描線都有 10 個頻道的定標係數。格式說明寫「(2,10)，(1,10) 存放斜率、(2,10) 存放截距」，這是 Fortran 的下標，
  實際存放順序是 斜率1、截距1、斜率2、截距2…。原腳本把前 10 個當斜率、後 10 個當截距，配對全錯，這是原本「沒有好的定標方法」的主因。
  正確配對後的意義與 NOAA pre-KLM 1b 相同：
  - 可見光 / 近紅外 (CH1、CH2、CH6\~CH10)：反照率 (%) = 斜率 / 2^30 × 計數 + 截距 / 2^22；反射率 = 反照率 / 100，未除 cos(太陽天頂角) (同 AVHRR-NAT)。
    零點落在太空計數 (12\~29)。
  - 紅外 (CH3\~CH5)：輻射 (mW/(m² sr cm⁻¹)) = 斜率 / 2^30 × 計數 + 截距 / 2^22，斜率為負，每條掃描線依星上黑體更新。
    亮溫 T\* = c2·ν / ln(1 + c1·ν³ / L)，T = (T\* − A) / B；中心波數 ν 與頻帶修正 A、B 取自 RTTOV 係數檔的 FILTER_FUNCTIONS
    (濾波函數由 NSMC 提供；`rtcoef_fy1_4_mvisr.dat`、`rtcoef_fy1_3_mvisr.dat`)：

    | 頻道 | FY-1D ν (cm⁻¹) | A | B | FY-1C ν (cm⁻¹) | A | B |
    |---|---|---|---|---|---|---|
    | CH3 | 2699.586 | 3.3894 | 0.995275 | 2691.460 | 3.6385 | 0.994675 |
    | CH4 | 925.553 | 0.7552 | 0.997388 | 912.074 | 0.8458 | 0.996996 |
    | CH5 | 850.420 | 0.7725 | 0.997133 | 829.868 | 0.7262 | 0.997266 |

  - 檢驗 (FY-1D 2003\~2006 年 5 個檔)：亮溫約 190\~302 K；分裂窗差 CH4 − CH5 在暖洋面 (CH4 > 290 K) 中位數 1.5\~2.3 K、
    冷雲頂 (CH4 < 215 K) 中位數約 1 K，都在正常範圍。中心波數改用 1/波長 (12.0 µm = 833 cm⁻¹) 時冷雲頂的差約 3.5 K，不合理。
  - 限制：沒有做非線性修正 (格式說明與檔案內都沒有相關係數)。可見光的係數在 2003\~2006 年的檔案中完全相同 (應為發射前係數)，
    沒有反映感測器衰減，反射率的絕對值僅供參考 (畫圖不受影響)。
* **定位**：51 個定位點位於第 25、65、…、2025 像元 (1 起算，同 NOAA LAC 1b)。以衛星天頂角反推的掃描角驗證
  (第 1 點約 −54.0°、間隔 2.16° = 40 像元 × 0.0541°)。原腳本把 51 點平均分布在整條掃描線上，掃描帶邊緣約偏 24 個像元。
  沿掃描以地心角為座標、在地心直角座標線性內插，兩端各 24 個像元外插；留一法估計內部約 0.4 km、邊緣約 1\~4 km
  (定位點本身以 1/128° 量化，約 0.9 km)。日本關東位於掃描帶東緣的夜間紅外影像，海岸線與陸海溫差的交界吻合。
* 經緯度是有號整數；原腳本讀成無號整數，南半球、西半球的掃描線會被略過。
* 各頻道邊緣有幾個像元是配準時補上的 0 或固定值 (CH5/9/10 第 1\~2、CH7 第 1\~4、CH3/CH8 最後 2 個像元)，設為無效。
* 每條掃描線最後一個 4-byte 字依序是 CH9、CH10、0 (實際資料；格式說明寫成 0、CH9、CH10)。資料頭的掃描線數與實際筆數不符，以檔案大小計算。
* **真實色** (頻道選單的 0)：R = CH1 (0.58–0.68 µm)、G = CH9 (0.53–0.58)、B = CH7 (0.43–0.48)，模式同 SGLI：1.RGB、2.CORR、3.WV。
  - 感測器衰減：2003\~2006 年檔案的可見光係數完全相同，2011 年的檔案已更新 (CH1 斜率 0.0893 → 0.1292)。厚冰雲的 TOA 反射率在可見光應接近中性
    (以 CREFL 的正向模式估算，雲頂上方有 200\~400 hPa 的 Rayleigh 時 CH7/CH1 仍在 1.00 ± 1%)，實測 CH7/CH1 卻逐年下降：

    | 日期 | 2003/09 | 2004/10 | 2005/07 | 2005/09\~10 (3 軌) | 2006/04 | 2006/10 | 2011/05 (新係數) |
    |---|---|---|---|---|---|---|---|
    | CH9/CH1 | 1.00 | 0.99 | 0.97 | 0.95 | 0.92 | 0.94 | 0.93 |
    | CH7/CH1 | 0.95 | 0.94 | 0.84 | 0.82\~0.83 | 0.77 | 0.72 | 0.85 |

    直接合成會偏黃，年份越晚越嚴重。因此真實色先乘上**頻道間增益**：以整軌的厚冰雲 (CH4 < −30°C、CH1/cos(θs) > 0.5、θs < 75°、
    衛星天頂角 < 50°) 的 TOA 比值中位數，令 CH9、CH7 與 CH1 一致；同期的三軌彼此相差 < 2%。厚冰雲不足 100 點時 (晨昏)，
    改用同一組定標係數的實測值依日期內插。執行時會印出增益與依據。單頻道輸出不受影響；CH1 本身的絕對衰減無法由此得知。
  - Rayleigh / 臭氧：沒有可見光頻道的光譜響應函數，以格式說明的波段範圍 (方波) 並以太陽光譜加權平均：
    τR = 0.0577 / 0.0945 / 0.2133、aO3 = 0.0912 / 0.0911 / 0.0054 (CH1 / CH9 / CH7)。同樣方法重算 MODIS B1/B3/B4、VIIRS M3/M4/M5/I1 的
    CREFL 係數，τR 誤差 < 4%、臭氧多在 6% 以內。臭氧 0.319 cm-atm，海平面，不做水氣修正。
  - 幾何：檔案內 51 個定位點的太陽天頂角、衛星天頂角、相對方位角，與定位相同以地心角內插。與 pyorbital 的天文計算比對，太陽天頂角相差 < 0.02°；
    相對方位角 = 太陽方位角 − 衛星方位角 (0° 為背向散射，同 CREFL 的約定)。
  - 已知的限制：厚雲校正後略帶米黃 (B/R 約 0.9，CREFL 以整層大氣校正高雲，MODIS 也一樣)；太陽天頂角大時太陽反光偏卡其色
    (CREFL 假設朗伯地表)；2003\~2006 年的晴空海面藍光偏低、海面偏橄欖色 (對 CH7 的定標與 τR 很敏感；
    τR 改為 0.18 時這三個場景的海面正常，但 2011 年的新係數用目前的 τR 已經偏藍，所以不以調整 τR 處理)。
* **驗證狀況** (2026/11)：FY-1D 13 個檔 (2003\~2011)、10 個頻道都可讀，CH4 色表、CH2、真實色 (2003/09 杜鵑、2005/09 龍王颱風、2011/05 黃海) 出圖目視正常。
  FY-1C 沒有樣本 (係數已放入，未實測)。

### NOAA AVHRR/3 (NSMC 格式)

* 檔案：國家衛星氣象中心的 NOAA-15\~19 HRPT 1B (`NA16_AVHRR_HRPT_L1_ORB_MLT_NUL_…L1B`)，每筆 22016 bytes，格式見 NSMC
  《NOAA 系列極軌氣象衛星數據格式》§3.2。與 NOAA KLM 1b 的 16-bit 格式相近，但幾個地方不同：可見光係數的縮放 (斜率 ×10¹⁰、截距 ×10⁷)、
  紅外係數的順序 (二次項、一次項、常數項)、衛星代碼 (NOAA-16/17/18 = 3/11/13)；表頭的輻射換算常數填錯 (NOAA-18 的檔是 NOAA-16 / 17 的值)，不使用。
* 定標：可見光為雙增益的反照率 (交點欄位在 2008 年以前是 0，由兩條直線計算；NOAA-17 CH3A 的斜率2 溢位 int32，由連續條件反推)，
  反射率 = 反照率 / 100。紅外的線性輻射已含太空輻射 Ns (以太空計數約 990 驗證)，CH4 / CH5 再做 NOAA KLM 的非線性修正；
  亮溫 T = (T\* − A) / B，ν、A、B 取自 NOAA KLM User's Guide (與 RTTOV 的 NOAA-18 係數相同)。
* 通道 3A / 3B：掃描線 bit 域的通道 3 旗標在 2008 年以前與 2010 年以後的意義相反，同一檔內也有不一致，改由 3B 定標係數的線性項正負判斷
  (負 = 3B；3A 時該欄位存的是 3A 的反照率係數)，已以計數分布驗證。
* 定標係數的品質：個別掃描線的紅外係數會跳動 (約 1\~4 K)，與 51 條掃描線的移動中位數相差超過門檻時改用中位數。
  NOAA-17 2010/10/28 的紅外係數整檔無效 (AVHRR 故障期間)。
* 資料受損：品質旗標沒有標出來。整條錯位 (CH4 − CH5 不合理的像元 > 10%，或與上下掃描線不連續) 時整條遮蔽，
  連續不超過 2 條的以上下內插補上；單點位元錯誤只遮蔽 CH4 / CH5 中分裂窗差不合理的像元 (避免被當成颱風眼的最高亮溫)；
  可見光的單條線段雜訊以上下平均取代。NOAA-16 2005/09/28 整軌、2005/07/16 北緯 32 度以北是錯位的亂碼，會遮蔽 (範圍內沒有資料)。
* 定位：51 個定位點在第 25\~2025 像元 (NSMC 格式說明)，以地心角內插 (同 MVISR-L1B)。台灣 (掃描角約 34°) 與海南島 (接近掃描帶邊緣) 的
  CH2 影像海岸線吻合。
* **驗證狀況** (2026/11)：NOAA-16 / 17 / 18 共 23 個檔 (2004\~2016) 可讀；暖洋面 CH4 − CH5 約 0.8\~3.5 K。莫蘭蒂 (2016/09/13) CH4、
  尤特 (2013/08/11) CH1、NOAA-17 CH3A 出圖目視正常。NOAA-15、NOAA-19 沒有樣本。

### NOAA AVHRR Level 1b (CLASS、千葉大學)

* 檔案：NOAA Level 1b 的兩種格式，格式說明是 NCEI 的 POD Guide (TIROS-N\~NOAA-14) 與 NOAA KLM User's Guide (NOAA-15 以後)，
  `D:\github\HRPT\doc` 有下載的 PDF。
  - POD：GAC 每筆 3220 bytes (409 像元，資料從第 3 筆起)，LAC/HRPT 14800 bytes (2048 像元)。
  - KLM：GAC 4608、LAC/HRPT/FRAC 15872 bytes (10-bit)。
  - 資料夾中的每個檔依內容判斷：檔案大小與記錄長度、表頭欄位、前置檔頭 (0 / 122 bytes TBM / 512 bytes)；改過檔名也能讀。
  - 衛星以檔內的資料集名稱 (`NSS.HRPT.NN.…`) 為準，其次是表頭的衛星 ID (POD 與 KLM 的編號不同)。
  - 千葉大學的檔 (處理區塊 `B0000000`) 在列表的資料來源顯示為 Chiba。它把 NOAA-18 存成 POD 格式，表頭 ID 7 在 POD 裡是 NOAA-9。
* 頻道：AVHRR/3 同 AVHRR-NSMC (CH3A / CH3B 依掃描線 bit 域切換；KLM 衛星存成 POD 格式時依頻道 3 定標斜率的正負)。
  AVHRR/2 (NOAA-7、9、11、12、13、14) 為 CH1\~CH5；AVHRR/1 (TIROS-N、NOAA-6、8、10) 為 CH1\~CH4。
* 可見光定標：反射率 = 反照率 / 100，未除 cos(太陽天頂角)。
  - KLM：記錄中的業務係數 (雙增益，NOAA 逐次更新)。
  - POD：記錄中的單一斜率 (多半是發射前值；NOAA-14 1996/11 起逐月更新)。
  - KLM 衛星存成 POD 格式 (或 KLM 檔沒有可用的係數)：改用 NOAA STAR 的在軌係數 (見下面 AVHRR-L0 的可見光)，
    沒有時用附錄 D 的發射前雙增益係數。
* 紅外定標：亮溫 T = (T\* − A) / B。
  - KLM：N = a0 + a1·C + a2·C²。記錄的係數已含太空輻射與非線性修正；a2 的縮放在格式版本 <= 2 為 10⁶，>= 3 的 CH4 / CH5 為 10⁷。
  - POD：1b 只是太空與黑體兩點的線性定標，太空點的輻射依衛星與年代不同 (0 或負值；NOAA-9 在 1986\~87 年改過幾次)。
    所以改用 KLM 的模型重算：由遙測的太空、黑體計數求線性位置，黑體點輻射取 1b 係數，再套 Walton 的非線性修正。
    NOAA-14 與指南的二次式完全相同；其他衛星與指南的 Ns = 0 修正表相差約 0.1\~0.4 K。
  - 接收站產品或紅外係數整軌固定 (千葉 NOAA-17)、以及 KLM 衛星存成 POD 格式 (千葉 NOAA-18 的紅外斜率錯了約 3.5 倍)：
    依 KLM Guide §7.1.2.4 由 PRT、黑體與太空計數自行定標。以 CLASS 的 NOAA-15/18/19 驗證，與業務係數相差 ±0.05 K 以內。
    PRT 每 5 條掃描線有一條參考讀數，指南寫 0，實際常是 1\~5 (POD 衛星都不是 0)：三個讀數的中位數 < 20 即為參考。
    接收站的資料常有個別讀數的位元錯誤，所以每條取三個讀數的中位數，四顆 PRT 各自去掉離群值、平滑後再平均
    (千葉 NOAA-17 / 18 的黑體溫度原本會跳動數 K，現在整軌平順)。
* 定位：51 個定位點 (LAC/HRPT 第 25\~2025 像元、GAC 第 5\~405 像元)，以地心角內插。定位無效或與前後掃描線不連續的掃描線遮蔽。
* 時間：時間碼受損的掃描線 (千葉的檔開頭常見；CLASS 也有個別錯誤，例如日數錯成 19) 改用推算值。
  資料缺口的前後不會被誤判；掃描線序號在缺口處不可靠，不使用。
* 不同來源的合併：同衛星的 LAC、HRPT、FRAC 解析度相同，在列表中一次選取 (例如 `[10,11]`) 就合併成一張
  (GAC 另成一組)。同一次過境的不同來源 (CLASS 的 LAC 與接收站的 HRPT、CLASS 與千葉) 依下列規則合併：
  - 掃描線依時間對齊，再以頻道 4 計數驗證：同一筆資料應幾乎逐位元相同，時間差了整數條掃描線時依計數校正。
  - 每條掃描線選用「無資料與尖峰比例」最低的來源。以重疊段整體最好的來源為主，其他來源明顯較好時才換。
  - 選中的掃描線裡無資料或尖峰的像元，以其他來源同一位置的好值補上。
  - 例：米契 (1998/10/26 NOAA-14) 的 LAC (20:28–20:34) 與 HRPT (20:31–20:45) 重疊 903 條，計數 99.94% 相同；
    合併後涵蓋 1\~58N，HRPT 開頭的條紋由 LAC 取代。人造測試：時間偏移 1 秒可自動對齊，
    主要來源的 300 個尖峰有 287 個由另一來源補回。
* **驗證狀況** (2026/11)：30 個檔全部可讀 (NOAA-9\~19，GAC / LAC / HRPT，CLASS 與千葉)。
  - 暖洋面 CH4 − CH5 為 1.3\~3.4 K，冷雲 0\~1.9 K。
  - 定位與颱風眼位置：海燕 (2013/11/07 GAC)、米契 (1998/10/26 NOAA-14)、尤特 (2013/08/11 千葉 NOAA-18) 的眼與最佳路徑相符，
    海岸線吻合。
  - 接收站資料的位元錯誤雜點 (千葉 NOAA-17 約 10%、NOAA-18 約 2% 的像元分裂窗差不合理；CLASS < 0.1%)
    與過境頭尾的條紋，目前照原值畫出。
  - 沒有樣本：TIROS-N、NOAA-6/7/8/13、MetOp (NOAA 1b)、FRAC。

### NOAA AVHRR HRPT 原始框 (Level 0)

* 檔案：接收站直接收錄的 HRPT 小框，每條掃描線一框 (11090 個 10-bit 字，每字存成 16 bits，小端或大端自動判斷)。
  - 未壓縮或 `.xz` / `.gz` / `.bz2`，也可以直接放未壓縮的 `.tar` (列表讀出其中的每個檔，不必先解開)；資料目錄往下兩層的子資料夾也會搜尋。
    同一個檔同時在 `.tar` 與資料夾中時只列一次。
  - 依框同步字辨識，不看檔名。衛星由 ID 字的衛星位址判斷 (AAPP 的對照表，TIROS-N\~NOAA-19；同一位址依年份區分)，
    路徑中的 `noaa18`、`n14` 等字樣優先。
  - 原始框沒有年份：取路徑中的日期 (`2003-02-04/`、`20130811.063140.noaa18.raw.xz`)，沒有時取檔案 (或 tar 成員) 的修改時間。
  - 樣本：`D:\github\L0` 的三個 tar (2003/02/04 NOAA-14\~17、2013/08/11 與 2013/11/07 NOAA-15/16/18/19，共 42 次過境)。
* TLE：原始框沒有軌道資料，以 TLE (SGP4，需要 `sgp4` 或 `pyorbital` 套件) 推算衛星位置。
  - 讀 `Satima-PO.toml` 的 `[tle] dir` (預設 `tle/`，含子資料夾與 .zip) 與資料目錄中的 `.tle` / `.txt`，取 epoch 最接近觀測的一組；
    相差超過 `[tle] max_age` 天 (預設 5) 時警告。最接近的三組推算的位置相差超過 20 km 時也會警告 (歷史檔偶有互相矛盾的根數，
    例如 CelesTrak 1988 年的 NOAA-9 有兩組交錯、相差 160 km)。
  - `tle/` 附 CelesTrak 的歷史檔 (TIROS-N\~NOAA-17，1980\~2004)。2005 年以後用 `python tools/fetch_spacetrack_tle.py` 從 Space-Track
    下載 (免費帳號；帳號密碼放在家目錄的 `.spacetrack`，或環境變數 `SPACETRACK_IDENTITY` / `SPACETRACK_PASSWORD`)，
    存到 `tle/spacetrack/` (依 Space-Track 的使用條款不進版控)。NOAA-12/14/15/16/17/18/19 全部約 85 次查詢、5 分鐘、11 MB。
* 定位：與 AAPP 相同的幾何 (大地天底、掃描面垂直於慣性速度、每像元 25 µs、與 WGS84 橢球求交點)，在 51 個定位點計算後同
  AVHRR-L1B 內插；加上各衛星的 roll、yaw 與步進角 (`avhrr_hrpt.NAV_PARAMS`)。
  - 參數以海岸線擬合：紅外影像 (冬夜陸地比海冷、夏日相反) 與 Satima 的 30 角秒陸海遮罩，分 128×128 的區塊求最吻合的位移，
    再以最小平方解出各衛星的 roll、yaw、步進角 (每次過境另有一個時間差，吸收時鐘與 TLE 的誤差)。
    資料：2003/02/04 與 2013/11/07 日本附近 14 次升降軌過境；殘差 RMS 由 0.9\~3.3 km 降到 0.3\~0.5 km。
  - 結果 (本程式的符號：roll > 0 往第 2048 像元那側)：NOAA-18 roll −0.18°、NOAA-19 yaw +0.33°、NOAA-15 yaw +0.18°、
    NOAA-16 yaw −0.13° 且步進角 0.05396° (AAPP 為 0.053955°)，其他衛星 0.05406\~0.05409° (AAPP 由取樣頻率算出 0.0540865°)；
    各衛星共有約 −0.03° 的 roll。不修正時掃描邊緣差到 5\~10 km。與 NOAA 1b 表頭的固定姿態修正 (NOAA-18 roll 0.15°、
    NOAA-19 yaw 0.25°) 與 AAPP 的預設偏差同量級。NOAA-14 只有 2 次受損的降軌可用 (roll −0.076°、yaw −0.056°)；
    沒有擬合的衛星用 AAPP 的步進角與 roll −0.03°。
  - 天底：KLM Guide 與 AAPP 都是橢球法線 (大地天底)。海岸線比對的南北誤差：大地天底約 −1 km，地心天底約 −3.5 km
    (升降軌反號，與時鐘誤差分得開)。另外，2013 年以後的 NOAA 業務定位也吻合大地天底，1998 年的 CLASS 1b 則較接近地心天底。
  - 獨立驗證 (擬合沒有用到的 2013/08/11 過境)：多數區塊的誤差由 1\~5 km 降到 1 km 內。
  - 不做時鐘修正：原始時間碼沒有 NOAA 的時鐘修正，會留下固定的沿軌差。以 2013/11/07 的 CLASS GAC 比對，NOAA-18 快 0.225 秒，
    約沿軌 1.5 km (2013/08/11 的海岸線也顯示約 −1.8 km)；NOAA-19 0.05 秒、NOAA-15 0.003 秒。
* 時間：同 AVHRR-L1B 修復受損的時間碼。以下的框為無資料：框同步字對不到 4 個、KLM 的 ID 字顯示為測試訊號、
  與前一框重複 (時間碼與小框編號相同；同 AAPP)。
  接收中斷的時段 (例如 2003/02/04 NOAA-17 在 00:59:59\~01:00:14 少了約 90 框) 在圖上是一條無資料帶。
* 定標：
  - 紅外一律由 PRT 黑體溫度、黑體與太空計數自行定標 (同上)。PRT 係數：NOAA-15 以後為 KLM 附錄 D；更早的衛星為 pygac 的值，
    與 AAPP 的 avhcal.txt 相同，以 CLASS 1b 驗證 (NOAA-9/10/11/14)，換算的黑體溫度與 1b 係數在黑體計數處的溫度相差中位數 ±0.012 K 以內。
    NOAA-14 的頻帶修正 (ν、A、B) pygac 與 AAPP 的值不同，換算的亮溫只差 0.01 K 以內。
  - 端對端驗證：2013/11/07 的原始框平均成 GAC 像元 (GAC 像元 n = LAC 樣本 5n\~5n+3，計數相差 0.25) 後與 CLASS GAC 的業務定標比較，
    CH4 亮溫差中位數 NOAA-18 −0.020 K、NOAA-19 −0.033 K、NOAA-15 −0.022 K (第 5\~95 百分位約 ±0.15 K)。
  - 可見光：NOAA-16\~19 用 NOAA STAR 的在軌係數 (F. Wu 以利比亞沙漠做的替代定標，每月一組；AAPP 附的檔，放在
    `conversion_table/AVHRR/`)，取 Center 日期最接近觀測者。涵蓋期間：NOAA-16 2003/05\~2014、NOAA-17 2004/02\~2010、
    NOAA-18 2005/10\~2017、NOAA-19 2009\~2017。早於第一組 120 天以上時，在發射日 (發射前係數) 與第一組之間內插；
    相差超過 400 天時改用發射前係數。差異很大：2013 年 NOAA-18 的 CH2 比發射前高 24.5%、NOAA-16 CH2 高 32%。
    NOAA-15 與更早的衛星沒有 STAR 的係數，用發射前係數 (NOAA-15 以後為附錄 D 的雙增益，更早為 POD Guide 表 3.3.2-1)，反射率會偏低。
  - 3A / 3B：同 AAPP，頻道 3 的太空計數 (10 個樣本中至少 5 個 < 100 為 3A、> 900 為 3B) 且與 ID 字的 bit 0 (1 = 3A) 一致，
    不一致或無法判斷的掃描線以前後內插，再取 11 條的中位數。
* AAPP 8.15 的對照 (程式碼在 `D:\github\L0\AAPP_8.15.tgz`)：上述的大地天底、步進角、位址表、重複框、3A/3B、STAR 係數都取自 AAPP。
  AAPP 沒有處理 AVHRR 掃描馬達抖動造成的錯位與重複片段 (下一點)，也只用 TLE 的時間與標稱的 1/6 秒推算各條掃描線。
* **資料本身的損壞** (2026/11 樣本)：2003/02/04 的 NOAA-14 約 22% 像元、NOAA-15/16 約 2\~4% 是錯位或重複的片段
  (同一段數值原封不動地出現在相鄰幾框且位置逐框漂移，呈階梯狀條帶)；NOAA-17 與 2013 年的檔 < 0.5% (2013 年的 NOAA-15 約 3%)。
  型態符合 AVHRR 掃描馬達抖動時 MIRP 重送舊的掃描線，不是解碼的問題；目前照原值畫出，颱風眼溫度也可能被這些值拉高。

## 真實色校正檢查

| | 太陽天頂角校正 | Rayleigh (CREFL) | O3 | H2O | 地形高度 | 銳化 | 增強 |
|---|---|---|---|---|---|---|---|
| VIIRS | SDR 已內含 | ✔ | 氣候值 | ✔ | ✔ GMTCO Height | I1 比值銳化 | Simple(sRGB) / WorldView / CIRA |
| MODIS | ✔ (88° 起漸減) | ✔ | 氣候值 | ✔ | ✔ MOD03 Height | B1 比值銳化 | CIRA / WorldView |
| METimage | ✔ | ✔ (ISRF 加權 τ) | 氣候值 | ✔ | — 海平面 | — | CIRA / WorldView |
| OLCI | ✔ | ✔ (含海平面氣壓) | ECMWF | — | ✔ altitude | unsharp mask | CIRA (白點正規化) / WorldView |
| SGLI | ✔ | ✔ (JAXA ATBD 的 τR) | 氣候值 300 DU (JAXA Koz) | — | — 海平面 | — (皆 250 m) | CIRA / WorldView |
| MVISR | ✔ | ✔ (波段範圍估算的 τR) | 氣候值 0.319 cm-atm | — | — 海平面 | — | CIRA / WorldView；另有厚冰雲的頻道間增益 |

SGLI 真實色採 JAXA 的 L1B RGB 組合 R = VN08、G = VN05、B = VN03。VN08 / VN05 在亮雲上會飽和 (約 0.49 / 0.67)，飽和像元以
VN11 (868.5 nm 低增益) 補值，避免雲偏青色；此補值為 Satima 的處理，非 JAXA 官方作法。

CREFL 的 `chand`、`csalbr`、合成公式、太陽天頂角校正與增強曲線四支腳本原本各自複製一份，數學上完全相同，已集中到
`satima/po/truecolor.py`。已用相同輸入比對原腳本與新版，結果一致 (差異在 1e-6 以內)。各衛星特有的係數、水氣公式、
地形修正、銳化方式與增強微調仍留在各自的腳本。

## 修正的原腳本錯誤 (Satima-PO 時期)

- VIIRS：頻道名稱 `M10` 會被讀成 `M1`；M8–M11、I3 被當成亮溫而讀不到資料；I4、M12–M14、M16 沒有色表選單 (會當掉)；單頻道的資料列與經緯度列分別去除填充列，兩者列數可能不一致
- AVHRR-NAT：可見光頻道會發生 NameError；CH3A/3B 的輻射值縮放應為 1e4；CH3B 中心波數的縮放應為 1e2 (原本 CH3B 亮溫約為 -212°C，整張圖空白)；未依掃描線區分 3A/3B；VEADR 筆數、dummy MDR 寫死；每條掃描線都 print 一次輻射值
- AVHRR-NAT / AVHRR-NC：亮溫的頻帶修正公式反了。官方 (AVHRR L1B Product Guide §6.3.4) 為 T = A + B·T\*，原本寫成 (T\* − A)/B，亮溫偏高：CH3B 約 1.5–2.5 K、CH4 約 0.1–0.5 K、CH5 約 0–0.35 K，雲頂愈冷差愈多。已用 AVHRR/3 光譜響應函數積分驗證 (以 SRF 擬合出的 A、B 與檔案中的係數一致，官方公式誤差 < 0.003 K)
- METimage：頻道 11 (2.25 µm 反射率) 會出現色表選單並輸出 13 張相同的圖；選單最後一個頻道誤印為 18；tar 分支缺少 import
- SLSTR：可見光由 πL/1500 的近似改為 πL/(F0·cosθs) (F0 依偵測器取自產品，θs 由 tie-point 網格內插，已用天文計算驗證)；舊版產品才補乘 ESA 替代校正係數 (新版已內含)；移除會把暗海面畫成白色的處理；兩種 dpi 的檔名頻道格式不一致 (B8 / B08)
- SGLI：缺值 (Error_DN) 被當成最大值 (掃描範圍外約 14% 畫成亮灰色)；SW01–SW04 誤用熱紅外亮溫轉換表；1 km 頻道 (SW01/02/04) 定位錯誤；經緯度改為直角座標樣條內插 (誤差約原本的 1/14、快約 10 倍)
- 中心時間 (MODIS 以外的 7 支)：原本以檔案起訖時間依裁切範圍的中間列線性推算，改為取最接近中心點的像元的實際觀測時間 (MODIS 原本的作法)，各格式的時間尺度 (TAI / GPS / UTC) 統一換成 UTC (`satima/po/obstime.py`)。實測舊作法：裁切範圍被檔案邊界截斷時差 24 s (VIIRS)；整軌 NAT 檔在高緯度頭尾重疊時差 50 分鐘；SLSTR 錐狀掃描邊緣差 50 s
- 共用：色表選單 `10.NHC` 應為 `10.AVN`；選檔範圍檢查用的是檔名長度 `len(h5f)`；裁切範圍少一列一行；範圍內沒有資料時會當掉；`eval(input())` 改為安全的 `ast.literal_eval`
- MVISR (`SatimaV3-FY1D-NSMC.py`，2026/11 加入時)：定標係數的斜率 / 截距配對錯誤；可見光把反照率當成輻射再除以假設的太陽輻照度；
  亮溫的中心波數用 1/波長、沒有頻帶修正；定位點的像元位置錯誤 (邊緣偏約 24 像元)；經緯度讀成無號整數 (南半球、西半球的掃描線被略過)；
  資料頭的結束時間欄位型態錯誤；以資料頭的掃描線數決定讀幾筆 (與實際筆數不符)；每次都印出整個資料陣列
- AVHRR-NSMC (`SatimaV3-NOAA-NSMC_accnumpy.py`，2026/11 加入時)：亮溫的頻帶修正寫成 T = A + B·T\* (NOAA KLM 為 (T\* − A) / B；NOAA-18 CH4
  在 200 K 偏高約 0.3 K、CH3B 在 300 K 偏高約 1.6 K)；可見光把反照率乘 100π / F 再除 100，反射率放大約 2.4 倍 (反照率 41% 以上就飽和)；
  選單有通道 3 但選了會中斷；定位點平均分布在 0\~2047 像元 (應為 25\~2025)，且使用 scipy 1.14 已移除的 interp2d；QC 以第 2 / 98 百分位剔除係數
  (只套用在 CH1、CH4，每個檔固定改掉 4% 的掃描線)；掃描線時間都用第一條的日期；範圍內沒有資料時會中斷。
  可見光係數多除的 1000 / 10、紅外係數的順序、非線性修正都是對的 (NSMC 格式與 NOAA KLM 不同)

## 新增衛星

在 `sensors/po/` 放一支 `<感測器>-<格式>.py`，結構照其他腳本：單一檔案類別命名為 `<感測器>File`、讀取器命名為
`<感測器>Reader` (例如 `AvhrrFile` / `AvhrrReader`)。

```python
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr
from satima.po import SatReader, ui

class XxxFile:                    # 單一檔案
    def __init__(self, fullname):
        self.fullname, self.satellite = fullname, 'XXX'
        self.time = ...           # 起始時間 (UTC datetime)
        self.end = ...            # 結束時間；檔名沒有的話不設 (或在讀取器設 GRANULE_SECONDS = 標稱長度)
    def extract(self, band): ...  # 反射率 0–1 / 亮溫 °C
    def geocoord(self): ...       # lats, lons (與資料同大小)
    def row_times(self, nrows):   # 每一列的觀測時間 (UTC datetime64，換算見 satima/po/obstime.py)；
        ...                       # 檔案沒有逐列時間時可用 obstime.linear(start, end, nrows)

class XxxReader(SatReader):       # 讀取器
    NAME = 'XXX-FMT'
    SENSOR = 'XXX'
    DESCRIPTION = '中文說明'
    DESCRIPTION_EN = 'English description'
    PATTERNS = ('XXX_*.h5',)
    SOURCE = 'XXX'                # 資料列表的來源欄

    def make_file(self, name):
        return XxxFile(name)
    def select_band(self):
        self.band = ui.ask_int(tr('輸入頻道：', 'Input the band: '), (1, 2))
    def select_variants(self):    # 紅外線回傳色表列表，其餘 [None]
        return ask_colormaps(self.cfg, 'ir') if self.band == 2 else [None]
    def prepare(self):            # self.files = 這一段的檔案 (依時間排序)；多個時沿軌串接 (self.concat)
        f = self.files[0]
        data = f.extract(self.band)
        lats, lons = f.geocoord()
        if self.slice_swath(lats, lons) is None:    # 設定 self.lats / self.lons / self.box
            return False
        self.data = self.crop(data)
        self.time = self.center_time(f.row_times(lats.shape[0]))   # 最接近中心點的像元的觀測時間
        return True
    def plot(self, canvas, variant):
        if variant is None:
            canvas.mesh(self.lons, self.lats, self.data, cmap='gray', vmin=0., vmax=1.)
            cm = ''
        else:
            cm = canvas.bt_mesh(self.lons, self.lats, self.data, variant)
        return dict(band='B%02d' % self.band, colorbar=cm)

READER = XxxReader
```

* 資料列表、接續判定、選資料、投影與範圍、色表選單、海岸線、標註、檔名、存檔都由主流程處理，不需要另外寫。
* `canvas.mesh()` / `bt_mesh()` / `rgb()` 只是把這次要畫的資料交給畫布；取樣在主流程決定圖片大小之後才做，
  所以 `plot()` 每個色表 / 樣式只會被呼叫一次。
* `plot()` 回傳的 `band`、`colorbar` 會填入檔名與圖下方的標註。
* 一個波段一個檔的格式 (像 VIIRS) 要覆寫 `segments()` 把檔案整理成片段、覆寫 `begin()` 依所選的頻道挑檔，寫法見 `sensors/po/VIIRS-HDF5.py`。
* 只能一次讀一個檔的格式設 `MULTI_FILE = False`。
