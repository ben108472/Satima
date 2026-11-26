"""格式處理腳本的基底類別

每支格式處理腳本 (sensors/po/<感測器>-<格式>.py) 的結構相同：
  1. 衛星特有的常數與演算法 (保留原函式)
  2. 單一檔案類別 <感測器>File (例如 ViirsFile)：extract / geocoord / row_times / plot
     屬性 fullname、satellite、time (起始時間)；有的話再加 end (結束時間)。
     row_times 回傳每一列的觀測時間 (UTC datetime64)
  3. 讀取器類別 <感測器>Reader (例如 ViirsReader，繼承 SatReader)：只寫和衛星有關的部分，
     資料列表、選資料、範圍、海岸線、圖說、存檔都由主流程 (satima/po/family.py、satima/core) 負責
  4. READER = 讀取器類別

主流程呼叫順序：
  find() -> make_file()                列出檔案
  segments(files)                      把檔案接成「一段連續的觀測」(資料列表的一列；見 satima/core/segments.py)
  (使用者選資料、輸入範圍)               -> begin(段) 設定 self.files；self.georange、self.center
  select_band()                        頻道/產品選單 (第一段問一次，其餘各段沿用)
  select_variants()                    色表/樣式選單，回傳要畫的版本列表
  prepare()                            讀檔、串接、裁切；範圍內無資料回傳 False
                                       中心時間：self.time = self.center_time(串接後的逐列時間)
  plot(canvas, variant)                每個版本畫一次，回傳 dict(band=, colorbar=, label=)
                                       canvas 見 satima/po/canvas.py (mesh / bt_mesh / rgb)
"""
import datetime
from dataclasses import dataclass, field
from glob import glob

import numpy as np

from . import obstime, swath
from ..core import segments as seg
from ..core.i18n import tr


@dataclass
class Segment:
    """一段連續的觀測 (資料列表的一列)：一個或多個時間相接的檔案。"""
    satellite: str
    time: datetime.datetime             # 起始時間 (UTC)
    end: datetime.datetime              # 結束時間；檔名沒有結束時間、也沒有標稱長度時 = 最後一檔的起始時間
    files: list = field(default_factory=list)     # 檔案物件 (依時間排序；VIIRS 等一個波段一個檔的格式含各波段的檔)
    slots: int = 1                      # 片段數 (一個時間一個；備註欄的「N 檔」)
    missing: int = 0                    # 段內缺的片段數
    product: str = ''                   # 列表的波段 / 產品欄
    short: bool = False                 # 短片段接成的 (多段時間相接時自動合併，不詢問)
    extra: dict = field(default_factory=dict)


class SatReader:
    NAME = ''           # 格式名稱 (= 腳本檔名，例如 'VIIRS-HDF5')
    SENSOR = ''         # 感測器名稱 (檔名 {sensor} 與圖說)
    DESCRIPTION = ''    # 格式選單中的說明
    PATTERNS = ()       # 檔案搜尋樣式 (glob，於資料目錄中搜尋)
    MULTI_FILE = True   # 是否可把多個檔案沿軌串接
    SOURCE = ''         # 資料來源 (列表用：CLASS、EUMETSAT、LAADS、JAXA…)
    GRANULE_SECONDS = 0     # 檔名沒有結束時間時，一個檔的標稱長度 (秒；例如 MODIS 300)。0 = 不明

    def __init__(self, cfg):
        self.cfg = cfg
        self.all_files = []
        self.files = []
        self.georange = None    # (latmin, latmax, lonmin, lonmax)：圖框涵蓋的經緯度範圍
        self.center = None      # (clat, clon)
        self.time = None        # 影像中心時間 (datetime)
        self.res = ''           # 解析度 (檔名 {res})
        self.box = None         # 裁切範圍 (yi, yj, xi, xj)
        self.lats = None        # 裁切後的經緯度 (經度已依中心經度展開)
        self.lons = None

    # ---- 檔案列表 ----------------------------------------------------------
    @classmethod
    def find(cls, folder=None):
        """在 folder (預設為目前目錄) 搜尋檔案，回傳相對於 folder 的檔名。"""
        names = set()
        for pattern in cls.PATTERNS:
            names.update(glob(pattern, root_dir=folder))
        return sorted(names, key=cls.sort_key)

    @staticmethod
    def sort_key(name):
        return name

    def make_file(self, name):
        """由檔名建立單一檔案物件。"""
        raise NotImplementedError

    def listing(self, f):
        """(舊介面) 檔案列表的一列：(衛星, 產品/頻道, 時間, [其他])。現在只用第 2 欄當作列表的產品欄。"""
        return (getattr(f, 'satellite', ''), self.SENSOR, f.time)

    @property
    def satellite(self):
        return getattr(self.files[0], 'satellite', '') if self.files else ''

    # ---- 接續判定 (見 satima/core/segments.py) --------------------------------
    def span(self, f):
        """某個檔的 (起始, 結束) 時間。檔案物件有 end 就用；否則用標稱長度 GRANULE_SECONDS (不明時長度為 0)。"""
        end = getattr(f, 'end', None)
        if end is None:
            end = f.time + datetime.timedelta(seconds=float(self.GRANULE_SECONDS or 0))
        if end < f.time:                    # 檔名只有時分秒、跨日
            end += datetime.timedelta(days=1)
        return f.time, end

    def chain_key(self, f):
        """可以互相接續的檔案的分組依據 (同衛星、同產品)。"""
        return (getattr(f, 'satellite', ''), str(self.listing(f)[1]))

    def merge_key(self, segment):
        """一次選到多段時，可以合併成一張的分組依據 (見 family._merge)；預設同衛星、同產品。
        例如 AVHRR 的 LAC、HRPT、FRAC 解析度相同，可以合併 (同一次過境的不同來源由讀取器挑選較好的掃描線)。"""
        return (segment.satellite, segment.product)

    def segments(self, files):
        """把檔案接成幾段連續的觀測 [Segment, ...]。
        短片段 (一檔不到 15 分鐘) 間隔在 10 分鐘內就接成一段，中間缺檔照樣接、另外標示缺幾個；
        整圈 / 半圈的長檔各自一列 (一次選多列時才詢問要不要合併，見 family.run)。"""
        if not self.MULTI_FILE:
            return [self._segment([f], 0) for f in files]
        groups = {}
        for f in files:
            groups.setdefault(self.chain_key(f), []).append(f)
        out = []
        for items in groups.values():
            start = lambda f: self.span(f)[0]
            end = lambda f: self.span(f)[1]
            lengths = [(end(f) - start(f)).total_seconds() for f in items]
            if any(v >= seg.SHORT for v in lengths):
                out.extend(self._segment([f], 0) for f in items)
                continue
            for members, missing in seg.chain(items, start, end):
                out.append(self._segment(members, missing, short=True))
        return out

    def _segment(self, files, missing, short=False):
        files = sorted(files, key=lambda f: f.time)
        return Segment(getattr(files[0], 'satellite', ''), files[0].time, max(self.span(f)[1] for f in files),
                       files=files, slots=len(files), missing=missing, product=str(self.listing(files[0])[1]),
                       short=short)

    def begin(self, segment):
        """開始處理某一段：設定 self.files (一個波段一個檔的格式在這裡依所選的波段挑檔)。"""
        self.files = list(segment.files)

    def source_label(self, segment):
        return segment.extra.get('source') or self.SOURCE

    # ---- 互動選單 ----------------------------------------------------------
    def select_band(self):
        pass

    def select_variants(self):
        return [None]

    # ---- 資料 --------------------------------------------------------------
    def prepare(self):
        raise NotImplementedError

    def plot(self, canvas, variant):
        raise NotImplementedError

    # ---- 共用工具 ----------------------------------------------------------
    def slice_swath(self, lats, lons):
        """依 self.georange 計算裁切範圍並設定 self.box / self.lats / self.lons。
        回傳 [yi, yj, 總列數]，範圍內無資料回傳 None。"""
        lons = swath.wrap_lon(lons, self.center[1])
        box = swath.swath_box(self.georange, lats, lons)
        self.box = box
        if box is None:
            return None
        yi, yj, xi, xj = box
        self.lats = swath.as_float(lats[yi:yj, xi:xj])
        self.lons = lons[yi:yj, xi:xj]
        return [yi, yj, lats.shape[0]]

    def crop(self, arr):
        return swath.crop(self.box, arr)

    def center_time(self, rowtime):
        """影像中心時間 = 裁切範圍內最接近中心點 (clat, clon) 的像元的觀測時間 (UTC)。
        rowtime：裁切前的觀測時間 (datetime64，見 satima/po/obstime.py)，與傳給 slice_swath 的經緯度對應：
                 1-D = 每一列一個時間 (掃描式感測器)；2-D = 每個像元一個時間 (例如 SLSTR 錐狀掃描)。
        中心點不在資料範圍內時取最近的像元 (印出距離)；沒有任何有效時間時回傳 None。"""
        yi, yj, xi, xj = self.box
        rt = np.asarray(rowtime, dtype=obstime.US)
        rt = rt[yi:yj, xi:xj] if rt.ndim == 2 else rt[yi:yj, None]
        clat, clon = self.center
        dist = np.asarray(np.ma.filled(swath.angular_distance(self.lats, self.lons, clat, clon), np.nan),
                          dtype=np.float64)
        dist[~np.isfinite(dist)] = np.inf
        dist[np.broadcast_to(np.isnat(rt), dist.shape)] = np.inf
        r, c = np.unravel_index(np.argmin(dist), dist.shape)
        if not np.isfinite(dist[r, c]):
            print(tr('注意：範圍內沒有有效的觀測時間。', 'WARNING: no valid observation time in the selected range.'))
            return None
        t = obstime.to_datetime(rt[r, c if rt.shape[1] > 1 else 0])
        lon = (float(self.lons[r, c]) + 180.) % 360. - 180.
        print(tr('中心時間 %s UTC (最近的像元 %.3fN %.3fE，距中心 %.2f 度)',
                 'Center time %s UTC (nearest pixel %.3fN %.3fE, %.2f deg from center)')
              % (t.strftime('%Y-%m-%d %H:%M:%S'), self.lats[r, c], lon, dist[r, c]))
        return t

    def report_eyetemp(self, data):
        """紅外線頻道：印出中心點附近最高亮溫 (設定檔 [eyetemp])。"""
        e = self.cfg.eyetemp
        if e.enable:
            swath.eyetemp(data, self.lats, self.lons, self.center[0], self.center[1], e.radius_deg)
        if e.legacy_counts:
            swath.eyetemp_counts(data)

    @staticmethod
    def concat(arrays):
        """沿軌 (axis 0) 串接；有 masked array 時保留遮罩。"""
        if len(arrays) == 1:
            return arrays[0]
        if any(np.ma.isMaskedArray(a) for a in arrays):
            return np.ma.concatenate(arrays, axis=0)
        return np.concatenate(arrays, axis=0)
