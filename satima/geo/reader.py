"""各衛星讀取器的基底類別

每支讀取器 (sensors/<感測器>-<格式>.py) 只寫和衛星有關的部分：
  CHANNELS   波段表 {資料中的代號: Channel}，順序即選單順序；Channel.band 為顯示/檔名用名稱 (BAND-13 等)，
             Channel.alias 為英文功能代稱 (Blue、Clean LW Window 等，見 README)
  TRUECOLOR  真實色配置 {產品名稱: recipe}，例如 {'TrueColor': dict(red='B03', green='B02', ...)}；
             recipe 再有 ir (10.4 µm) 與 swir (3.9 µm) 時，選單自動多一個 GeoColor (夜間疊紅外雲與城市燈光)
  scenes()   把檔案分組成「場景」(同一時次、同一觀測範圍的所有波段)
  nav()      該波段的導航 (satima/nav.py 的 GeosNav 等)
  read()     讀取指定視窗並校正：反射率 0~1 (可見光、近紅外) / 亮溫 °C (紅外、水氣)；無效值為 NaN
  satellite_position()  星下點經度與衛星離地高度 (真實色計算觀測幾何用)
讀取器不保留「上一次讀取」的狀態：同一波段的 nav() / read() 不論呼叫順序結果都相同。
檔案列表、範圍、產品選單、取樣、繪圖、存檔都由主流程處理。
"""
import datetime
from dataclasses import dataclass, field
from glob import glob

import numpy as np

from ..core.i18n import tr


@dataclass(frozen=True)
class Channel:
    name: str           # 資料中的官方代號 (讀檔、pyspectral 查 RSR 用)，例如 B13、C02、IR105、ir_105_hr
    wavelength: float   # 標稱中心波長 (µm)，依各機構官方表
    res: float          # 星下點解析度 (km)
    kind: str           # vis / nir / swir / wv / ir
    band: str = ''      # 選單、圖片標註、檔名 {band} 顯示的名稱，例如 BAND-13 (空白 = name)
    alias: str = ''     # 功能代稱 (英文，選單用)，沿用 NOAA ABI 官方代稱，其他衛星依功能對應

    @property
    def reflective(self):
        """反射率波段 (可見光、近紅外)；其餘為亮溫。"""
        return self.kind in ('vis', 'nir')

    @property
    def display(self):
        return self.band or self.name

    @property
    def wl_text(self):
        """例如 10.4µm (標註與選單用)。"""
        text = '%g' % self.wavelength
        return (text if '.' in text else text + '.0') + 'µm'      # 12.0µm 而不是 12µm

    @property
    def label(self):
        return ' '.join(t for t in (self.display, self.wl_text, self.alias) if t)


@dataclass
class Scene:
    satellite: str                  # 例如 Himawari-9
    time: datetime.datetime         # 觀測時間 (UTC)
    area: str                       # 觀測範圍，例如 FLDK、R301、CONUS、FD
    files: dict = field(default_factory=dict)     # 波段名稱 -> 檔案 (或檔案列表)
    auto_range: bool = False        # 區域觀測 (目標區、mesoscale 等)：地圖範圍直接用整個影像
    extra: dict = field(default_factory=dict)     # 讀取器自用
    region: str = ''                # 觀測範圍的名稱 (資料列表用)，例如 TARGET、JAPAN (同舊版)
    tag: str = ''                   # 接在 {band} 後面的資料集標記 (例如 MTG 的 FDHSI / HRFI)，
                                    # 用來區分同一時次不同資料集的產品 (標註 BAND-14 HRFI、檔名 BAND-14_HRFI)
    format: str = ''                # 資料列表「格式」欄 (例如 NSMC GVAR)；空白時用讀取器名稱 '-' 之後的部分 (HSD、GVAR、AREA…)
    sensor: str = ''                # 這筆資料的儀器 (同一讀取器有多種儀器時，例如 VISSR / VAS)；空白時用讀取器的 SENSOR

    @property
    def area_label(self):
        """資料列表顯示的觀測範圍，例如 TARGET R301、JAPAN JP01、FLDK (名稱與代號相同時只顯示一次)。"""
        if not self.region or self.region == self.area:
            return self.area
        return '%s %s' % (self.region, self.area)


class GeoReader:
    NAME = ''           # 格式名稱 (= 腳本檔名，例如 'AHI-HSD')
    SENSOR = ''         # 感測器名稱 (檔名 {sensor} 與標註)
    DESCRIPTION = ''    # 格式選單中的說明
    PATTERNS = ()       # 檔案搜尋樣式 (glob，於資料目錄中搜尋)
    CHANNELS = {}       # 波段表
    TRUECOLOR = {}      # 真實色配置

    READ_CACHE = 4      # 最近讀取的視窗數 (同一波段畫多種色表、真實色重複使用時免重讀)

    def __init__(self, cfg):
        self.cfg = cfg
        self._read_cache = {}

    # ---- 檔案 --------------------------------------------------------------
    @classmethod
    def find(cls, folder='.'):
        names = set()
        for pattern in cls.PATTERNS:
            names.update(glob(pattern, root_dir=folder))
        return sorted(names)

    def scenes(self, folder='.'):
        """回傳 Scene 列表 (依時間排序)。"""
        raise NotImplementedError

    def channels(self, scene):
        """這個場景有的波段 (依 CHANNELS 的順序)。"""
        return [ch for name, ch in self.CHANNELS.items() if name in scene.files]

    def format_label(self, scene):
        """資料列表的格式欄：場景自訂的 format，否則為讀取器名稱 '-' 之後的部分 (AHI-HSD -> HSD)。"""
        return scene.format or self.NAME.split('-', 1)[-1]

    def sensor_label(self, scene):
        """資料列表的儀器欄與檔名 {sensor}：場景自訂的 sensor，否則為讀取器的 SENSOR。"""
        return scene.sensor or self.SENSOR

    def truecolor_products(self, scene):
        """這個場景可做的真實色產品 {名稱: recipe}。"""
        out = {}
        for name, recipe in self.TRUECOLOR.items():
            need = [recipe[k] for k in ('red', 'green', 'blue', 'nir') if recipe.get(k)]
            if all(b in scene.files for b in need):
                out[name] = recipe
        return self.with_geocolor(out, scene)

    @staticmethod
    def with_geocolor(products, scene):
        """有 10.4 µm (ir) 與 3.9 µm (swir) 的真實色，另加一個 GeoColor (名稱的 TrueColor 換成 GeoColor)。"""
        out = dict(products)
        for name, recipe in products.items():
            if all(recipe.get(k) and recipe[k] in scene.files for k in ('ir', 'swir')):
                out[name.replace('TrueColor', 'GeoColor')] = dict(recipe, style='geocolor')
        return out

    # ---- 資料 --------------------------------------------------------------
    def nav(self, scene, channel):
        raise NotImplementedError

    def read(self, scene, channel, window):
        """window = (r0, r1, c0, c1)：陣列索引範圍 (含 r0、不含 r1)。回傳 float32 陣列 (r1-r0, c1-c0)。"""
        raise NotImplementedError

    def read_cached(self, scene, channel, window):
        key = (id(scene), channel, tuple(window))
        if key in self._read_cache:
            self._read_cache[key] = self._read_cache.pop(key)     # 移到最新
        else:
            self._read_cache[key] = self.read(scene, channel, window)
            while len(self._read_cache) > self.READ_CACHE:
                self._read_cache.pop(next(iter(self._read_cache)))
        return self._read_cache[key]

    def satellite_position(self, scene):
        """(星下點經度 (度), 衛星離地高度 (km))。"""
        raise NotImplementedError

    def platform(self, scene):
        """pyspectral 的衛星名稱 (真實色 Rayleigh 校正用)。"""
        return scene.satellite

    def auto_georange(self, scene):
        """區域觀測：以影像邊緣與內部取樣點的經緯度範圍作為地圖範圍 (結果存在 scene.extra，資料列表與繪圖共用)。"""
        if 'georange' not in scene.extra:
            scene.extra['georange'] = self._image_range(scene)
        return scene.extra['georange']

    def _image_range(self, scene):
        name = next(iter(scene.files))
        nav = self.nav(scene, name)
        ny, nx = nav.shape
        rows, cols = np.meshgrid(np.linspace(0, ny - 1, 41), np.linspace(0, nx - 1, 41), indexing='ij')
        lon, lat = nav.to_lonlat(rows, cols)
        ok = np.isfinite(lon) & np.isfinite(lat)
        if not ok.any():
            raise SystemExit(tr('無法計算 %s 的影像範圍', 'Cannot determine the image range of %s') % scene.area)
        lon = lon[ok]
        if lon.max() - lon.min() > 180:     # 跨換日線：改用 0~360
            lon = lon % 360.
        return float(lat[ok].min()), float(lat[ok].max()), float(lon.min()), float(lon.max())

    # ---- 共用工具 ----------------------------------------------------------
    @staticmethod
    def full_window(nav):
        return 0, nav.shape[0], 0, nav.shape[1]


class FileCache:
    """最近使用的開啟檔案 (LRU)：超過上限時關閉最久沒用的檔案。
    多時次、多波段 (例如上百個 GOES mesoscale 檔) 時避免同時開太多檔案 (Windows 約 512 個上限)。"""

    def __init__(self, opener, maxsize=32):
        self._opener = opener
        self.maxsize = maxsize
        self._items = {}

    def get(self, key):
        if key in self._items:
            self._items[key] = self._items.pop(key)        # 移到最新
        else:
            self._items[key] = self._opener(key)
            while len(self._items) > self.maxsize:
                self._items.pop(next(iter(self._items))).close()
        return self._items[key]

    def clear(self):
        while self._items:
            self._items.popitem()[1].close()


def to_float(arr, dtype=np.float32):
    """masked array / 整數 -> 浮點陣列，遮罩處為 NaN。"""
    return np.ma.filled(np.ma.asarray(arr, dtype=dtype), np.nan)
