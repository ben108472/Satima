"""各微波儀器讀取器的基底類別

微波輻射計的資料是「掃描帶」：每個取樣點自帶經緯度，沒有解析式的導航。每支讀取器
(sensors/<儀器>-<格式>.py) 只寫和檔案格式有關的部分：
  scenes()    把檔案整理成「場景」(一個檔案 = 一段軌道)
  channels()  這個場景有的頻道 (Channel)；Channel.group 指出它用哪一組經緯度
  geo()       某一組的地理定位 [Geo, ...]：一般只有一個；AMSR 的 89 GHz A / B 掃描合併時為兩個
  tb()        某頻道的亮溫 [陣列, ...] (K，無效值為 NaN)，個數與形狀對應 geo()
檔案列表、範圍、產品選單、內插 (swath.py)、繪圖、存檔都由主流程處理。
多檔合併 (同一衛星時間相接的幾個檔畫成一張) 也在這裡 (SwathReader.swath)：沿掃描方向串接，
以掃描時間去掉重疊的掃描 (JAXA L1B 每個檔前後各多 30 條)。
"""
import datetime
from dataclasses import dataclass, field
from glob import glob

import numpy as np

from ..core.i18n import tr

EPOCH = datetime.datetime(1970, 1, 1)

# 頻率 (GHz) -> 類別：色表分組、PCT / RGB 配對用
_KINDS = ((0., 8.5, 'tb06'), (8.5, 12.5, 'tb10'), (17., 20.5, 'tb19'), (20.5, 26., 'tb22'), (30., 40., 'tb37'),
          (45., 70., 'sounding'), (80., 95., 'tb89'), (140., 170., 'tb150'), (170., 200., 'tb183'))


def kind_of(freq):
    for lo, hi, kind in _KINDS:
        if lo <= freq < hi:
            return kind
    return 'tb'


def freq_text(freq):
    """19.35 -> '19.35'、89.0 -> '89'。"""
    return ('%.3f' % float(freq)).rstrip('0').rstrip('.')


@dataclass(frozen=True)
class Channel:
    name: str           # 選單、標註、檔名用，例如 89V、36.5H、183.31±7V (檔名中 ± 寫成 +-)
    freq: float         # 中心頻率 (GHz)
    pol: str            # 極化：V / H / QV / QH (橫向掃描儀器的準垂直 / 準水平)；不明為空
    group: str          # 地理定位組 (SwathReader.geo 的參數)
    key: object = None  # 讀取器自用 (資料集名稱、索引等)
    note: str = ''      # 選單附註，例如 TA (天線溫度)

    @property
    def kind(self):
        return kind_of(self.freq)

    @property
    def display(self):
        return self.name

    @property
    def freq_label(self):
        return freq_text(self.freq) + 'GHz'


@dataclass
class Scene:
    satellite: str                  # 例如 GCOM-W1、F18、GPM
    sensor: str                     # 例如 AMSR2、SSMIS、GMI
    time: datetime.datetime         # 檔案的觀測起始時間 (UTC)
    end: datetime.datetime          # 觀測結束時間
    file: str = ''                  # 檔案路徑
    format: str = ''                # 資料列表「格式」欄；空白時用讀取器名稱 '-' 之後的部分
    orbit: str = ''                 # 軌道 / granule 編號 (列表用)
    extra: dict = field(default_factory=dict)       # 讀取器自用
    parts: list = None              # 多檔合併時的各個場景 (依時間排序)

    @property
    def area(self):                 # 排序用 (同 Satima-GEO 的介面)
        return self.orbit


class Geo:
    """一組地理定位：lat / lon (掃描數, 每條掃描的點數)，time 為每條掃描的時間 (1970 起算的秒數，不明為 NaN)。
    key 用來辨識同一組幾何 (內插座標的快取)。"""

    def __init__(self, lat, lon, time, key):
        self.lat = np.asarray(lat)
        self.lon = np.asarray(lon)
        self.time = np.asarray(time, dtype=np.float64)
        self.key = key
        self.shape = self.lat.shape

    def spacing(self):
        """(沿掃描, 沿軌) 的取樣間距 km：掃描帶中央附近的中位數。"""
        ny, nx = self.shape
        j = nx // 2
        rows = slice(max(0, ny // 2 - 50), min(ny, ny // 2 + 50))

        def dist(lat1, lon1, lat2, lon2):
            p1, p2 = np.deg2rad(lat1.astype(np.float64)), np.deg2rad(lat2.astype(np.float64))
            dl = np.deg2rad(lon2.astype(np.float64) - lon1.astype(np.float64))
            h = np.sin((p2 - p1) / 2.) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.) ** 2
            return 2. * 6371. * np.arcsin(np.sqrt(h))

        with np.errstate(invalid='ignore'):
            across = dist(self.lat[rows, j - 1], self.lon[rows, j - 1], self.lat[rows, j], self.lon[rows, j])
            along = dist(self.lat[:-1, j], self.lon[:-1, j], self.lat[1:, j], self.lon[1:, j])
        across = across[np.isfinite(across)]
        along = along[np.isfinite(along) & (along > 0)]
        return (float(np.median(across)) if across.size else np.nan,
                float(np.median(along)) if along.size else np.nan)


class SwathReader:
    NAME = ''           # 格式名稱 (= 腳本檔名，例如 'PPS-L1C')
    SENSOR = ''         # 儀器名稱 (列表用；一支讀取器讀多種儀器時由 Scene.sensor 決定)
    DESCRIPTION = ''    # 格式選單中的說明
    PATTERNS = ()       # 檔案搜尋樣式 (glob，於資料目錄中搜尋)
    SOURCE = ''         # 資料來源 (列表用：PPS、JAXA、CLASS、NSMC、EUMETSAT…)；個別場景可用 Scene.extra['source'] 覆寫

    GEO_CACHE = 6       # 最近使用的地理定位組數

    def __init__(self, cfg):
        self.cfg = cfg
        self._geo_cache = {}

    # ---- 檔案 --------------------------------------------------------------
    @classmethod
    def find(cls, folder='.'):
        names = set()
        for pattern in cls.PATTERNS:
            names.update(glob(pattern, root_dir=folder))
        return sorted(names)

    def scenes(self, folder='.'):
        """回傳 Scene 列表。"""
        raise NotImplementedError

    def channels(self, scene):
        """這個場景的頻道列表 (選單順序)。"""
        raise NotImplementedError

    def format_label(self, scene):
        return scene.format or self.NAME.split('-', 1)[-1]

    def source_label(self, scene):
        return scene.extra.get('source') or self.SOURCE

    def sensor_label(self, scene):
        return scene.sensor or self.SENSOR

    def note(self, scene):
        """資料列表最後一欄的備註：多檔接成一筆時的檔數，有缺檔時另外標出 (Scene.extra 的 nfiles、missing)。"""
        from ..core import segments
        n = scene.extra.get('nfiles') or len(scene.extra.get('granules', ()))
        return segments.note(n, scene.extra.get('missing', 0))

    # ---- 資料 (各讀取器實作) ------------------------------------------------
    def geo(self, scene, group):
        """-> [Geo, ...]"""
        raise NotImplementedError

    def tb(self, scene, channel):
        """-> [亮溫陣列 (K，float32，無效為 NaN), ...]，對應 geo(scene, channel.group)。"""
        raise NotImplementedError

    # ---- 主流程使用的介面 ---------------------------------------------------
    def geo_cached(self, scene, group):
        key = (id(scene), group)
        if key in self._geo_cache:
            self._geo_cache[key] = self._geo_cache.pop(key)
        else:
            self._geo_cache[key] = self.geo(scene, group)
            while len(self._geo_cache) > self.GEO_CACHE:
                self._geo_cache.pop(next(iter(self._geo_cache)))
        return self._geo_cache[key]

    def swath(self, scene, channel):
        """某頻道的 [(Geo, 亮溫), ...]；多檔合併的場景沿掃描方向串接。"""
        parts = scene.parts or [scene]
        if len(parts) == 1:
            s = parts[0]
            return list(zip(self.geo_cached(s, channel.group), self.tb(s, channel)))
        per = []
        for s in parts:
            ch = self.channel_named(s, channel.name)
            if ch is None:
                raise ValueError(tr('%s 沒有 %s，無法合併', '%s has no %s; cannot merge') % (s.file, channel.name))
            per.append((self.geo_cached(s, ch.group), self.tb(s, ch)))
        out = []
        for k in range(len(per[0][0])):
            geos, tbs = [p[0][k] for p in per], [p[1][k] for p in per]
            ckey = ('merged', id(scene), channel.group, k)
            if ckey not in self._geo_cache:
                self._geo_cache[ckey] = _concat_geo(geos, ckey)
            geo, keep = self._geo_cache[ckey]
            out.append((geo, np.concatenate([t[m] for t, m in zip(tbs, keep)], axis=0)))
        return out

    def channel_named(self, scene, name):
        for ch in self.channels(scene):
            if ch.name == name:
                return ch
        return None

    def pairs(self, scene):
        """同頻率同時有 V、H 的頻道對 {類別: [(V 頻道, H 頻道), ...]} (PCT、RGB 用)。"""
        by = {}
        for ch in self.channels(scene):
            if ch.pol in ('V', 'H'):
                by.setdefault((ch.kind, ch.freq, ch.group), {})[ch.pol] = ch
        out = {}
        for (kind, _, _), d in by.items():
            if 'V' in d and 'H' in d:
                out.setdefault(kind, []).append((d['V'], d['H']))
        return out


def _concat_geo(geos, key):
    """沿掃描方向串接多個檔案的地理定位；以掃描時間去掉與前一檔重疊的掃描。
    回傳 (Geo, [各檔保留的掃描遮罩])。"""
    if len({g.shape[1] for g in geos}) != 1:
        raise ValueError(tr('各檔每條掃描的點數不同，無法合併',
                            'the files have different numbers of samples per scan; cannot merge'))
    keep, last, step = [], -np.inf, None
    for g in geos:
        t = g.time
        if step is None:
            d = np.diff(t[np.isfinite(t)])
            step = float(np.median(d)) if d.size else 1.
        m = np.isfinite(t) & (t > last + 0.25 * step)
        if not np.isfinite(t).any():        # 沒有掃描時間：整個檔照接
            m = np.ones(t.size, bool)
        elif m.any():
            last = float(t[m].max())
        keep.append(m)
    lat = np.concatenate([g.lat[m] for g, m in zip(geos, keep)], axis=0)
    lon = np.concatenate([g.lon[m] for g, m in zip(geos, keep)], axis=0)
    time = np.concatenate([g.time[m] for g, m in zip(geos, keep)], axis=0)
    return Geo(lat, lon, time, key), keep


def merge_scenes(scenes):
    """把多個場景合成一個 (依時間排序)。"""
    parts = sorted(scenes, key=lambda s: s.time)
    first = parts[0]
    return Scene(first.satellite, first.sensor, first.time, max(s.end for s in parts), file=first.file,
                 format=first.format, orbit=first.orbit, parts=parts)


# ---- 共用工具 --------------------------------------------------------------
def to_seconds(dt):
    """datetime -> 1970 起算的秒數。"""
    return (dt - EPOCH).total_seconds()


def from_seconds(sec):
    return EPOCH + datetime.timedelta(seconds=float(sec))


# TAI - UTC 增加 1 秒的日期 (1993-01-01 之後；當時 TAI - UTC = 27 秒)
_LEAPS = [to_seconds(datetime.datetime(*d)) for d in (
    (1993, 7, 1), (1994, 7, 1), (1996, 1, 1), (1997, 7, 1), (1999, 1, 1), (2006, 1, 1), (2009, 1, 1), (2012, 7, 1),
    (2015, 7, 1), (2017, 1, 1))]
_TAI93 = to_seconds(datetime.datetime(1993, 1, 1))


def tai93_to_unix(t):
    """TAI93 (1993-01-01 起算、含閏秒的秒數；JAXA AMSR 的 Scan Time) -> 1970 起算的 UTC 秒數。無效 (<= 0) 為 NaN。"""
    t = np.asarray(t, dtype=np.float64)
    u = t + _TAI93
    u = u - np.searchsorted(np.array(_LEAPS), u, side='right')
    return np.where(t > 0, u, np.nan)


_TAI58 = to_seconds(datetime.datetime(1958, 1, 1))


def iet_to_unix(t):
    """IET (IDPS Epoch Time：1958-01-01 起算、含閏秒的微秒數；JPSS 的 StartTime) -> 1970 起算的 UTC 秒數。
    無效 (<= 0) 為 NaN。只適用 1993 年以後 (TAI - UTC >= 27 秒)。"""
    t = np.asarray(t, dtype=np.float64)
    u = t * 1e-6 + _TAI58 - 27.
    u = u - np.searchsorted(np.array(_LEAPS), u, side='right')
    return np.where(t > 0, u, np.nan)


def clean_latlon(lat, lon):
    """無效的經緯度 (超出範圍、填補值) 改為 NaN；回傳 float32。"""
    lat = np.array(lat, dtype=np.float32)
    lon = np.array(lon, dtype=np.float32)
    bad = ~(np.isfinite(lat) & np.isfinite(lon) & (np.abs(lat) <= 90.) & (np.abs(lon) <= 360.))
    lat[bad] = np.nan
    lon[bad] = np.nan
    return lat, lon


def band_name(freq, pol, offset=None, suffix=''):
    """頻道名稱：89V、36.5H、183.31±7V。"""
    text = freq_text(freq)
    if offset:
        text += '±' + freq_text(offset)
    return text + pol + suffix
