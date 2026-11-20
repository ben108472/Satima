# Satima：極軌成像儀 (PO)

> 這份是 **PO (極軌成像儀)** 的詳細說明，改寫自原本的 Satima-PO。整合版的執行方式、資料列表、設定檔分層
> (`Satima.toml` + `Satima-PO.toml`)、介面語言、圖片模式 (native / dpi) 與極軌的接續判定請看 [README.md](../README.md)。
> PO 是三類之中改動最大的：與原本 Satima-PO 的差異列在 README 的「與原本三個程式的差異」。

VIIRS、MODIS、AVHRR、METimage、OLCI、SLSTR、SGLI 的單頻道 (反射率、亮溫配色表、DNB) 與真實色。
8 支原本各自獨立的腳本先整合成 Satima-PO，再併入 Satima；各格式的校正與演算法都留在自己的腳本裡。

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
* 標題 (原本的 Picname) 預設會詢問 (`Satima-PO.toml` 的 `[title] ask = true`)，可留空；會填入檔名並畫在圖片左下角。

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

SGLI 的 VN07 / VN10 是高增益 (海色) 頻道，雲上約在反射率 0.16 / 0.13 就飽和，VN08 約 0.49 飽和；畫雲請用 VN11 (868.5 nm)。

各格式的起訖時間：VIIRS、METimage、AVHRR-NAT、OLCI、SLSTR 取自檔名；MODIS 檔名只有起始時間，以標稱的 5 分鐘計；
AVHRR-NC、SGLI 的檔名沒有結束時間，列表只顯示起始時間 (間隔 10 分鐘內的檔仍會接成一段)。AVHRR 兩種格式一次只讀一個檔。

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

## 真實色校正檢查

| | 太陽天頂角校正 | Rayleigh (CREFL) | O3 | H2O | 地形高度 | 銳化 | 增強 |
|---|---|---|---|---|---|---|---|
| VIIRS | SDR 已內含 | ✔ | 氣候值 | ✔ | ✔ GMTCO Height | I1 比值銳化 | Simple(sRGB) / WorldView / CIRA |
| MODIS | ✔ (88° 起漸減) | ✔ | 氣候值 | ✔ | ✔ MOD03 Height | B1 比值銳化 | CIRA / WorldView |
| METimage | ✔ | ✔ (ISRF 加權 τ) | 氣候值 | ✔ | — 海平面 | — | CIRA / WorldView |
| OLCI | ✔ | ✔ (含海平面氣壓) | ECMWF | — | ✔ altitude | unsharp mask | CIRA (白點正規化) / WorldView |
| SGLI | ✔ | ✔ (JAXA ATBD 的 τR) | 氣候值 300 DU (JAXA Koz) | — | — 海平面 | — (皆 250 m) | CIRA / WorldView |

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
