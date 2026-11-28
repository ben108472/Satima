"""
MERSI-HDF5 — 風雲三號 MERSI 系列 L1 (NSMC，HDF5)

衛星 / 儀器
  FY-3D            MERSI-II   25 頻道 (1~4、24、25 有 250 m)
  FY-3F、FY-3H     MERSI-III  頻道配置同 MERSI-II
  FY-3E            MERSI-LL   1 微光 (LL，0.7 µm) + 6 紅外 (6、7 有 250 m)
  FY-3G            MERSI-RM   8 頻道，全部 500 m
  FY-3A/B/C        MERSI-1    20 頻道 (1~5 有 250 m；5 為 11.25 µm 紅外)  * 沒有樣本，照格式說明與 satpy 寫，未驗證

所需檔案 (同一目錄；一個片段 5 分鐘，時間相接的片段自動接成一筆)
  FY3X_MERSI_{GBAL|GRAN}_L1_YYYYMMDD_HHmm_<種類>_{MS|Vn}.HDF
  資料  1000M (1 km)、0250M (250 m)、0500M (FY-3G，500 m)
  定位  GEO1K (1 km，含角度)、GEOQK (250 m)、GEOHK (FY-3G 500 m，含角度)、GEOXX (FY-3A/B 250 m)
  * 經緯度：同解析度的定位檔 > (250 m) GEO1K 逐掃描帶內插 > 資料檔裡每 5 / 20 點一個的 tie point (未做地形校正)。
  * 真實色、微光縮放需要角度：取自 GEO1K / GEOHK (250 m 由 1 km 內插)；沒有時太陽角由時間與經緯度計算，
    衛星角由星下點推算 (球面幾何)，畫面上會註明。
  * 範圍確定後先用 tie point 判斷哪些掃描帶、哪幾行在地圖內，只讀那一塊 (250 m 的整個片段很大)。

定標 (NSMC 各衛星的「一級產品使用說明」，2026/11 自 satellite.nsmc.org.cn 的文件頁取得)
  反射通道  dn = DN x Slope + Intercept；Ref = (k0 + k1 dn + k2 dn^2) / 100 (VIS_Cal_Coeff，單位 %)
            MERSI-RM：Ref = k0 + k1 dn (RSB_Cal_Coeff，已是 0~1)
            Ref 為 NSMC 的定標反射率 (= π L / E0，未除 cos 太陽天頂角)；單頻道照這個值畫 (灰階 0~1)
  紅外通道  RAD = DN x Slope + Intercept  (mW m-2 sr-1 cm)；Te = Planck^-1(RAD, 等效中心波數)；Tbb = A Te + B
            等效中心波長取自檔案 (Effect_Center_WaveLength)，A、B 取自 TBB_Trans_Coefficient(_A / _B)；
            檔案裡沒有或數值不合理時改用文件的表格 (IR_TABLE)。(satpy 用 (Te - B) / A，與 NSMC 文件相反)
            MERSI-1 的第 5 頻道：Slope 誤植為 1 (應為 0.01)，A、B = 1.0047、-0.8549 (同 satpy)
  微光 (LL) EV_1KM_LL 為輻亮度 W m-2 sr-1；ERF 縮放換成 W cm-2 sr-1 後與 VIIRS DNB 相同 (satima/po/dnb.py)
  時間      Calibration/EV_start_time：2000-01-01 12:00 UTC 起算的秒 (FY-3C 為當日的十進位小時)

產品
  單頻道  反射率 (灰階 0~1) / 亮溫 (紅外、水氣色表；MERSI-II/III 的 22、MERSI-LL 的 4 用水氣色表) / 微光 LL
  真實色  1~3 頻道 (MERSI-1/II/III)，250 m 或 1 km：1.General (只做太陽天頂角校正) 2.CIRA 3.NASA Worldview
          CREFL 沒有 FY-3 的官方係數：O3 / H2O 借用 MODIS 相近波段，Rayleigh 光學厚度由等效中心波長計算
"""
import datetime
import os
import re
from dataclasses import dataclass, field

import h5py
import numpy as np

import fy3_common as fc
from satima.core import segments as seg
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr
from satima.po import SatReader, Segment, ui, dnb, swath

_NAME = re.compile(r'^(FY3[A-Z])_MERSI_([A-Z]{4})_L1_(\d{8})_(\d{4})_([0-9A-Z]{5})_(MS|V\d+)', re.I)
DATA_KINDS = ('0250M', '0500M', '1000M')
RES_LABEL = {'0250M': '250m', '0500M': '500m', '1000M': '1km'}
NATIVE_GEO = {'0250M': ('GEOQK', 'GEOXX'), '0500M': ('GEOHK',), '1000M': ('GEO1K',)}
SAT_SENSOR = {'A': 'MERSI-1', 'B': 'MERSI-1', 'C': 'MERSI-1', 'D': 'MERSI-II', 'E': 'MERSI-LL',
              'F': 'MERSI-III', 'G': 'MERSI-RM', 'H': 'MERSI-III'}
SAT_HEIGHT = {'FY-3G': 407.}        # km；其餘約 836 km
GEO_GROUPS = ('Geolocation', '', 'GEO')
WV_BANDS = {'MERSI-II': (22,), 'MERSI-III': (22,), 'MERSI-LL': (4,)}
TC_BANDS = {'R': 3, 'G': 2, 'B': 1}                 # MERSI-1 / II / III
FAR_MARGIN = 0.3                    # 判斷掃描帶是否在地圖內時，範圍外加的度數

# 文件的紅外參數表：(等效中心波數 cm-1, A, B)，依紅外頻道順序。檔案裡沒有或不合理時才用
IR_TABLE = {
    'FY-3D': ((2634.359, 2471.654, 1382.621, 1168.182, 933.364, 836.941),
              (1.00103, 1.00085, 1.00125, 1.00030, 1.00133, 1.00065),
              (-0.4759, -0.3139, -0.2662, -0.0513, -0.0734, 0.0875)),
    'FY-3E': ((2623.369, 2466.214, 1384.461, 1164.837, 926.606, 837.013),
              (1.00090, 1.00057, 1.00110, 1.00026, 1.00121, 1.00102),
              (-0.501577, -0.313816, -0.246650, -0.083270, -0.300752, -0.250539)),
    'FY-3F': ((2641.6461, 2476.5705, 1384.4775, 1165.6627, 927.37106, 828.01898),
              (1.00078, 1.00050, 1.00105, 1.00034, 1.00162, 1.00121),
              (-0.566653, -0.346261, -0.434122, -0.120750, -0.4480995, -0.289908)),
    'FY-3G': ((2624.158, 929.837, 830.676), (1.00069, 1.00143, 1.00114), (-0.485743, -0.425257, -0.306088)),
    'FY-3H': ((2634.322236, 2468.717622, 1385.940533, 1169.168806, 931.601394, 835.885240),
              (1.001381, 1.000840, 1.001284, 1.000295, 1.001578, 1.001128),
              (-0.418879, -0.255056, -0.357935, -0.082549, -0.483007, -0.342837)),
}
MERSI1_IR = (1e4 / 11.25, 1.0047, -0.8549)

LL_MODES = {1: 'ERF Dynamic', 2: 'ERF Dynamic + saturation correction', 3: 'Linear (pi*L/E0, 0~1)'}
LL_NAMES = {1: 'ERF', 2: 'ERF-SAT', 3: 'Linear'}


@dataclass(frozen=True)
class Band:
    num: int
    wl: float               # 標稱中心波長 µm
    kind: str               # 'ref' / 'bt' / 'll'
    sources: tuple          # ((資料種類, 資料集, 波段索引或 None), ...)
    cal: int                # 反射通道：定標係數的列；紅外：紅外頻道的順序

    def source(self, kind):
        for k, name, idx in self.sources:
            if k == kind:
                return name, idx
        return None

    @property
    def name(self):
        return 'LL' if self.kind == 'll' else 'CH%02d' % self.num

    @property
    def text(self):
        return '%s(%s)' % (self.name, ('%.3f' % self.wl).rstrip('0').rstrip('.'))


def _bands_2_3():
    wl = (0.47, 0.55, 0.65, 0.865, 1.38, 1.64, 2.13, 0.412, 0.443, 0.49, 0.555, 0.67, 0.709, 0.746, 0.865, 0.905,
          0.936, 0.94, 1.03, 3.8, 4.05, 7.2, 8.55, 10.8, 12.0)
    out = []
    for b in range(1, 26):
        if b <= 4:
            src = (('0250M', 'EV_250_RefSB_b%d' % b, None), ('1000M', 'EV_250_Aggr.1KM_RefSB', b - 1))
            out.append(Band(b, wl[b - 1], 'ref', src, b - 1))
        elif b <= 19:
            out.append(Band(b, wl[b - 1], 'ref', (('1000M', 'EV_1KM_RefSB', b - 5),), b - 1))
        elif b <= 23:
            out.append(Band(b, wl[b - 1], 'bt', (('1000M', 'EV_1KM_Emissive', b - 20),), b - 20))
        else:
            src = (('0250M', 'EV_250_Emissive_b%d' % b, None), ('1000M', 'EV_250_Aggr.1KM_Emissive', b - 24))
            out.append(Band(b, wl[b - 1], 'bt', src, b - 20))
    return tuple(out)


def _bands_ll():
    wl = (0.7, 3.8, 4.05, 7.2, 8.55, 10.8, 12.0)
    out = [Band(1, wl[0], 'll', (('1000M', 'EV_1KM_LL', None),), 0)]
    for b in range(2, 6):
        out.append(Band(b, wl[b - 1], 'bt', (('1000M', 'EV_1KM_Emissive', b - 2),), b - 2))
    for b in (6, 7):
        src = (('0250M', 'EV_250_Emissive_b%d' % b, None), ('1000M', 'EV_250_Aggr.1KM_Emissive', b - 6))
        out.append(Band(b, wl[b - 1], 'bt', src, b - 2))
    return tuple(out)


def _bands_rm():
    wl = (0.65, 0.865, 0.94, 1.38, 1.64, 3.8, 10.8, 12.0)
    out = [Band(b, wl[b - 1], 'ref', (('0500M', 'EV_Reflectance', b - 1),), b - 1) for b in range(1, 6)]
    out += [Band(b, wl[b - 1], 'bt', (('0500M', 'EV_Emissive', b - 6),), b - 6) for b in range(6, 9)]
    return tuple(out)


def _bands_1():
    wl = (0.47, 0.55, 0.65, 0.865, 11.25, 1.64, 2.13, 0.412, 0.443, 0.49, 0.52, 0.565, 0.65, 0.685, 0.765, 0.865,
          0.905, 0.94, 0.98, 1.03)
    out = []
    for b in range(1, 21):
        if b <= 4:
            src = (('0250M', 'EV_250_RefSB_b%d' % b, None), ('1000M', 'EV_250_Aggr.1KM_RefSB', b - 1))
            out.append(Band(b, wl[b - 1], 'ref', src, b - 1))
        elif b == 5:
            src = (('0250M', 'EV_250_Emissive', None), ('1000M', 'EV_250_Aggr.1KM_Emissive', None))
            out.append(Band(b, wl[b - 1], 'bt', src, 0))
        else:
            out.append(Band(b, wl[b - 1], 'ref', (('1000M', 'EV_1KM_RefSB', b - 6),), b - 2))
    return tuple(out)


BANDS = {'MERSI-II': _bands_2_3(), 'MERSI-III': _bands_2_3(), 'MERSI-LL': _bands_ll(), 'MERSI-RM': _bands_rm(),
         'MERSI-1': _bands_1()}


def _sat_letter(sat):
    return sat[-1].upper()


# ===== 檔案與片段 =====
class MersiPart:
    """一個實體檔案。"""

    def __init__(self, name):
        m = _NAME.match(os.path.basename(name))
        if not m:
            raise ValueError(tr('不認得的 MERSI 檔名：%s', 'Unexpected MERSI file name: %s') % name)
        self.fullname = name
        self.satellite = 'FY-3' + m.group(1)[-1].upper()
        self.kind = m.group(5).upper()
        self.time = datetime.datetime.strptime(m.group(3) + m.group(4), '%Y%m%d%H%M')


class Granule:
    """一個 5 分鐘片段：同衛星、同起始時間的各種檔案 {種類: 檔名}。"""

    def __init__(self, satellite, time):
        self.satellite = satellite
        self.time = time
        self.end = time + datetime.timedelta(minutes=5)
        self.files = {}

    @property
    def sensor(self):
        return SAT_SENSOR.get(_sat_letter(self.satellite), 'MERSI-II')

    def __lt__(self, other):
        return self.time < other.time


@dataclass
class Window:
    """一個片段裡要讀的範圍：掃描帶 s0~s1 (不含)、行 c0~c1 (不含)。"""
    granule: Granule
    kind: str
    s0: int
    s1: int
    c0: int
    c1: int
    shape: tuple            # 整個片段的 (列, 行)
    rps: int                # 每個掃描帶的列數
    tie: object = None      # (lat, lon, 列間隔, 行間隔)：資料檔的 tie point
    notes: set = field(default_factory=set)

    @property
    def rows(self):
        return slice(self.s0 * self.rps, self.s1 * self.rps)

    @property
    def cols(self):
        return slice(self.c0, self.c1)


def _ds_shape(f, kind, sensor):
    """資料檔中任一個波段資料集的 (列, 行)。"""
    for b in BANDS[sensor]:
        src = b.source(kind)
        if src:
            ds = fc.find_ds(f, src[0])
            if ds is not None:
                return ds.shape[-2:]
    raise ValueError(tr('%s 裡找不到資料集', 'no data set found in %s') % f.filename)


def _nscans(f, nrows):
    n = fc.attr_array(f, 'Number Of Scans')
    if n is not None and n.size and n[0] > 0 and nrows % int(n[0]) == 0:
        return int(n[0])
    t = fc.find_ds(f, 'EV_start_time', ('Calibration', 'Data', ''))
    if t is not None and t.shape and nrows % t.shape[0] == 0:
        return t.shape[0]
    return nrows // 10


def _step(ds, key, n_full, n_tie):
    text = fc.attr_text(ds, key)
    nums = [int(v) for v in re.findall(r'\d+', text)[:2]]
    if len(nums) == 2 and nums[1] > nums[0]:
        return nums[1] - nums[0]
    return max(1, int(round(n_full / float(n_tie))))


def _tie_points(f, shape):
    """資料檔裡的 tie point 經緯度 (FY-3D 起：Geolocation/ 或 GEO/ 的 Latitude，每 5 / 20 點一個)。"""
    lat = fc.find_ds(f, 'Latitude', GEO_GROUPS)
    lon = fc.find_ds(f, 'Longitude', GEO_GROUPS)
    if lat is None or lon is None or lat.ndim != 2:
        return None
    if lat.shape == tuple(shape):                   # FY-3A/B 的 1000M：檔內就是逐點的經緯度
        return None
    rstep = _step(lat, 'Line_number', shape[0], lat.shape[0])
    cstep = _step(lat, 'Pixel_number', shape[1], lat.shape[1])
    la, lo = fc.clean_latlon(lat[:], lon[:])
    return la, lo, rstep, cstep


def _geo_ds(f, shape):
    """定位檔 (或資料檔) 中逐點的 Latitude / Longitude 資料集；形狀不符時 None。"""
    lat = fc.find_ds(f, 'Latitude', GEO_GROUPS)
    lon = fc.find_ds(f, 'Longitude', GEO_GROUPS)
    if lat is None or lon is None or lat.shape != tuple(shape):
        return None
    return lat, lon


def _scaled_read(ds, rows, cols, fill_mask=True):
    raw = ds[rows, cols]
    slope, inter = fc.scaled(ds)
    out = raw.astype(np.float64) * slope + inter
    if fill_mask:
        fv = fc.attr_array(ds, 'FillValue')
        if fv is not None and fv.size:
            out[raw == fv[0]] = np.nan
        vr = fc.attr_array(ds, 'valid_range')
        if vr is not None and vr.size == 2 and vr[0] < vr[1]:
            out[(raw < vr[0]) | (raw > vr[1])] = np.nan
    return out


# ===== 讀取器 =====
class MersiReader(SatReader):
    NAME = 'MERSI-HDF5'
    SENSOR = 'MERSI'
    SOURCE = 'NSMC'
    DESCRIPTION = '風雲三號 MERSI-1/II/III/LL/RM L1 (NSMC，HDF5)'
    DESCRIPTION_EN = 'Fengyun-3 MERSI-1/II/III/LL/RM L1 (NSMC, HDF5)'
    PATTERNS = ('FY3?_MERSI_*_L1_*.HDF',)
    GRANULE_SECONDS = 300
    MULTI_FILE = True

    def make_file(self, name):
        return MersiPart(name)

    # ---- 片段與接續 ----------------------------------------------------------
    def segments(self, parts):
        grans = {}
        for p in sorted(parts, key=lambda p: p.fullname):       # 同一片段同一種檔有多個版本時留最後一個 (版本號大者)
            g = grans.setdefault((p.satellite, p.time), Granule(p.satellite, p.time))
            g.files[p.kind] = p.fullname
        by_sat = {}
        for g in grans.values():
            if any(k in g.files for k in DATA_KINDS):           # 只有定位檔的片段不列
                by_sat.setdefault(g.satellite, []).append(g)
        out = []
        for sat, items in by_sat.items():
            for members, missing in seg.chain(items, lambda g: g.time, lambda g: g.end):
                kinds = [RES_LABEL[k] for k in DATA_KINDS if any(k in m.files for m in members)]
                out.append(Segment(sat, members[0].time, max(m.end for m in members), files=members,
                                   slots=len(members), missing=missing, product='/'.join(kinds), short=True,
                                   extra={'sensor': members[0].sensor}))
        return out

    def sensor_label(self, segment):
        return segment.extra.get('sensor') or self.SENSOR

    def merge_key(self, segment):
        return (segment.satellite,)

    def begin(self, segment):
        self.granules = sorted(segment.files)
        self.files = self.granules
        self.sensor = self.granules[0].sensor

    # ---- 選單 --------------------------------------------------------------
    def _datasets(self, kind):
        """某種資料檔 (取第一個有的片段) 裡有的資料集名稱。"""
        for g in self.granules:
            if kind in g.files:
                with h5py.File(g.files[kind], 'r') as f:
                    names = set()
                    f.visit(lambda n: names.add(n.split('/')[-1]))
                return names
        return set()

    def select_band(self):
        bands = BANDS[self.sensor]
        kinds = [k for k in DATA_KINDS if any(k in g.files for g in self.granules)]
        options, groups = [], []
        for kind in kinds:
            have = self._datasets(kind)
            for kname, title in (('ref', 'Reflectance'), ('bt', 'Brightness temperature'), ('ll', 'Low light')):
                ks = []
                for b in bands:
                    src = b.source(kind)
                    if b.kind == kname and src and src[0] in have:
                        options.append(('band', kind, b))
                        ks.append(len(options))
                if ks:
                    groups.append(('%s %s' % (RES_LABEL[kind], title), ks))
        tc_ks = []
        if self.sensor in ('MERSI-1', 'MERSI-II', 'MERSI-III'):
            for kind in kinds:
                have = self._datasets(kind)
                if all(b.source(kind) and b.source(kind)[0] in have for b in bands if b.num in TC_BANDS.values()):
                    options.append(('tc', kind, None))
                    tc_ks.append(len(options))
        if tc_ks:
            groups.append(('RGB', tc_ks))
        if not options:
            raise SystemExit(tr('這一段沒有可用的頻道。', 'No usable channel in this segment.'))

        def text(opt):
            typ, kind, b = opt
            return 'TrueColor (%s)' % RES_LABEL[kind] if typ == 'tc' else b.text

        print('Products:')
        ui.print_columns([(title, ['%2d.%s' % (k, text(options[k - 1])) for k in ks]) for title, ks in groups])
        k = 1 if len(options) == 1 else ui.ask_int(tr('輸入產品編號：', 'Input the product: '),
                                                   range(1, len(options) + 1))
        typ, kind, b = options[k - 1]
        self.option = (typ, kind, b)
        self.option_sensor = self.sensor
        self.res = RES_LABEL[kind]

    def select_variants(self):
        typ, kind, b = self.option
        if typ == 'tc':
            return ui.ask_menu('Style:', fc.TC_STYLE_DESC, tr('選擇樣式 (可複選，0 = 全部)：', 'Select style: '))
        if b.kind == 'll':
            return ui.ask_menu('Low light scaling:', LL_MODES, tr('選擇微光縮放 (可複選，0 = 全部)：',
                                                                  'Select low light scaling: '))
        if b.kind == 'bt':
            return ask_colormaps(self.cfg, 'wv' if b.num in WV_BANDS.get(self.sensor, ()) else 'ir')
        return [None]

    def _band_for(self, sensor, ref):
        """所選的頻道在這一段的儀器上對應哪一個 (不同的 MERSI 依種類與波長對應)。"""
        if sensor == self.option_sensor:
            return ref
        for b in BANDS[sensor]:
            if b.kind == ref.kind and abs(b.wl - ref.wl) <= 0.03 * ref.wl:
                return b
        return None

    # ---- 範圍 --------------------------------------------------------------
    def _window(self, g, kind):
        """依 tie point (或定位檔抽樣) 判斷這個片段在地圖內的掃描帶與行；不在地圖內回傳 None。"""
        with h5py.File(g.files[kind], 'r') as f:
            shape = _ds_shape(f, kind, g.sensor)
            rps = shape[0] // _nscans(f, shape[0])
            tie = _tie_points(f, shape)
            coarse = tie
            if coarse is None:
                geo = _geo_ds(f, shape)
                if geo is not None:
                    coarse = fc.clean_latlon(geo[0][::8, ::8], geo[1][::8, ::8]) + (8, 8)
        if coarse is None:
            coarse = self._coarse_from_geo(g, kind, shape)
        if coarse is None:
            raise IOError(tr('%s %s 找不到定位資料 (定位檔或資料檔內的 tie point)',
                             '%s %s: no geolocation (geolocation file or tie points in the data file)')
                          % (g.satellite, g.time.strftime('%Y/%m/%d %H:%M')))
        lat, lon, rstep, cstep = coarse
        latmin, latmax, lonmin, lonmax = self.georange
        m = FAR_MARGIN
        lonw = swath.wrap_lon(lon, self.center[1])
        with np.errstate(invalid='ignore'):
            inside = (lat > latmin - m) & (lat < latmax + m) & (lonw > lonmin - m) & (lonw < lonmax + m)
        rows = np.flatnonzero(inside.any(axis=1))
        if rows.size == 0:
            return None
        cols = np.flatnonzero(inside.any(axis=0))
        nscan = shape[0] // rps
        r0, r1 = (rows[0] - 1) * rstep, (rows[-1] + 2) * rstep
        s0, s1 = max(0, r0 // rps), min(nscan, -(-r1 // rps))
        c0, c1 = max(0, (cols[0] - 1) * cstep), min(shape[1], (cols[-1] + 2) * cstep)
        return Window(g, kind, int(s0), int(s1), int(c0), int(c1), tuple(shape), int(rps), tie)

    @staticmethod
    def _coarse_from_geo(g, kind, shape):
        """資料檔沒有 tie point 時，由定位檔抽樣 (同解析度每 8 點、1 km 對 250 m 每 2 點)。"""
        for gk in NATIVE_GEO[kind]:
            if gk in g.files:
                with h5py.File(g.files[gk], 'r') as f:
                    geo = _geo_ds(f, shape)
                    if geo is not None:
                        return fc.clean_latlon(geo[0][::8, ::8], geo[1][::8, ::8]) + (8, 8)
        if kind == '0250M' and 'GEO1K' in g.files:
            with h5py.File(g.files['GEO1K'], 'r') as f:
                geo = _geo_ds(f, (shape[0] // 4, shape[1] // 4))
                if geo is not None:
                    return fc.clean_latlon(geo[0][::2, ::2], geo[1][::2, ::2]) + (8, 8)
        return None

    # ---- 經緯度、時間、角度 ---------------------------------------------------
    def _geo(self, w):
        """這個範圍的經緯度 (float32)。"""
        g, kind = w.granule, w.kind
        for gk in NATIVE_GEO[kind] + (kind,):
            if gk in g.files:
                with h5py.File(g.files[gk], 'r') as f:
                    geo = _geo_ds(f, w.shape)
                    if geo is not None:
                        return fc.clean_latlon(geo[0][w.rows, w.cols], geo[1][w.rows, w.cols])
        if kind == '0250M' and 'GEO1K' in g.files:
            with h5py.File(g.files['GEO1K'], 'r') as f:
                geo = _geo_ds(f, (w.shape[0] // 4, w.shape[1] // 4))
                if geo is not None:
                    w.notes.add(tr('250 m 經緯度由 GEO1K 逐掃描帶內插', '250 m lat/lon interpolated from GEO1K'))
                    r1 = w.rps // 4
                    j0, j1 = max(0, w.c0 // 4 - 1), min(w.shape[1] // 4, w.c1 // 4 + 2)
                    rows = slice(w.s0 * r1, w.s1 * r1)
                    lat, lon = fc.clean_latlon(geo[0][rows, j0:j1], geo[1][rows, j0:j1])
                    return fc.geo_interp(lat, lon, 4 * np.arange(r1) + 1.5, w.rps, 4 * np.arange(j0, j1) + 1.5,
                                         np.arange(w.c0, w.c1))
        if w.tie is not None:
            w.notes.add(tr('經緯度由資料檔內的 tie point 內插 (未做地形校正)',
                           'lat/lon interpolated from the tie points in the data file (no terrain correction)'))
            lat, lon, rstep, cstep = w.tie
            tps = max(1, w.rps // rstep)
            rows = slice(w.s0 * tps, w.s1 * tps)
            return fc.geo_interp(lat[rows], lon[rows], rstep * np.arange(tps), w.rps, cstep * np.arange(lat.shape[1]),
                                 np.arange(w.c0, w.c1))
        raise IOError(tr('%s 找不到定位資料', '%s: no geolocation') % g.files[kind])

    def _times(self, w):
        """每一列的觀測時間 (datetime64[us])。"""
        g = w.granule
        with h5py.File(g.files[w.kind], 'r') as f:
            ds = fc.find_ds(f, 'EV_start_time', ('Calibration', 'Data', ''))
            t = ds[:].astype(np.float64).ravel() if ds is not None else None
            begin = fc.attr_text(f, 'Observing Beginning Date')
        nscan = w.shape[0] // w.rps
        out = None
        if t is not None and t.size == nscan:
            good = t[np.isfinite(t) & (t > 0)]
            if good.size and np.median(good) > 1e6:                  # 2000-01-01 12:00 起的秒數
                out = fc.nsmc_times(t)
            elif good.size and good.max() <= 24.5:                   # FY-3C：當日的十進位小時
                day = np.datetime64(begin or g.time.strftime('%Y-%m-%d'), 'D').astype('datetime64[us]')
                hrs = np.where(np.isfinite(t) & (t >= 0), t, np.nan)
                with np.errstate(invalid='ignore'):
                    hrs = np.where(hrs < g.time.hour - 12, hrs + 24., hrs)     # 片段跨午夜
                out = day + np.round(np.nan_to_num(hrs) * 3.6e9).astype(np.int64).astype('timedelta64[us]')
                out[~np.isfinite(hrs)] = np.datetime64('NaT', 'us')
        if out is None:
            for gk in ('GEO1K', 'GEOQK', 'GEOHK', 'GEOXX'):
                if gk in g.files:
                    with h5py.File(g.files[gk], 'r') as f:
                        dc = fc.find_ds(f, 'Day_Count', ('Timedata', ''))
                        ms = fc.find_ds(f, 'Millisecond_Count', ('Timedata', ''))
                        if dc is not None and ms is not None and dc.shape[0] == nscan:
                            out = fc.day_ms_times(dc[:], ms[:].astype(np.float64) * fc.scaled(ms)[0])
                            break
        if out is None:                         # 備援：起訖時間平均分配
            start = np.datetime64(g.time, 'us')
            out = start + (np.arange(nscan) * (300e6 / nscan)).astype(np.int64).astype('timedelta64[us]')
        return np.repeat(out[w.s0:w.s1], w.rps)

    def _angle_source(self, g, kind, shape):
        """(檔名, 解析度倍數)：有 SolarZenith 的定位檔；倍數 1 = 同解析度、4 = 1 km 對 250 m。"""
        for gk in NATIVE_GEO[kind] + (kind,):           # FY-3A/B 的 1000M：角度在資料檔內
            if gk in g.files:
                with h5py.File(g.files[gk], 'r') as f:
                    ds = fc.find_ds(f, 'SolarZenith', GEO_GROUPS)
                    if ds is not None and ds.shape == tuple(shape):
                        return g.files[gk], 1
        if kind == '0250M' and 'GEO1K' in g.files:
            with h5py.File(g.files['GEO1K'], 'r') as f:
                if fc.find_ds(f, 'SolarZenith', GEO_GROUPS) is not None:
                    return g.files['GEO1K'], 4
        return None, 0

    def _angles(self, w, lat, lon, times, need=('sza', 'saz', 'vza', 'vaz', 'height')):
        """角度 (度) 與地形高度 (m)：dict。定位檔沒有角度時以計算值代替。"""
        g = w.granule
        path, factor = self._angle_source(g, w.kind, w.shape)
        out = {}
        names = {'sza': 'SolarZenith', 'saz': 'SolarAzimuth', 'vza': 'SensorZenith', 'vaz': 'SensorAzimuth',
                 'lza': 'MoonZenith', 'height': ('Altitude', 'DEM')}
        if path:
            with h5py.File(path, 'r') as f:
                for key in need:
                    cands = names[key] if isinstance(names[key], tuple) else (names[key],)
                    ds = next((d for d in (fc.find_ds(f, n, GEO_GROUPS) for n in cands) if d is not None), None)
                    if ds is None:
                        continue
                    if factor == 1:
                        out[key] = _scaled_read(ds, w.rows, w.cols)
                    else:
                        r1 = w.rps // 4
                        j0, j1 = max(0, w.c0 // 4 - 1), min(w.shape[1] // 4, w.c1 // 4 + 2)
                        src = _scaled_read(ds, slice(w.s0 * r1, w.s1 * r1), slice(j0, j1))
                        src = fc.angle_interp(src, *self._up_axes(w, r1, j0, j1)) if key in ('saz', 'vaz') else \
                            fc.scan_interp(src, *self._up_axes(w, r1, j0, j1))
                        out[key] = src
        if 'height' not in out:
            for gk in NATIVE_GEO[w.kind]:           # GEOQK 有 250 m 的 Altitude
                if gk in g.files:
                    with h5py.File(g.files[gk], 'r') as f:
                        ds = fc.find_ds(f, 'Altitude', GEO_GROUPS)
                        if ds is not None and ds.shape == tuple(w.shape):
                            out['height'] = _scaled_read(ds, w.rows, w.cols)
        if ('sza' in need and 'sza' not in out) or ('saz' in need and 'saz' not in out):
            w.notes.add(tr('太陽角由時間與經緯度計算 (定位檔沒有角度)',
                           'solar angles computed from time and lat/lon (no angles in the geolocation files)'))
            out['sza'], out['saz'] = fc.solar_angles(times, lat, lon)
        if ('vza' in need and 'vza' not in out) or ('vaz' in need and 'vaz' not in out):
            w.notes.add(tr('衛星角由星下點推算', 'satellite angles estimated from the sub-satellite track'))
            out['vza'], out['vaz'] = self._view_angles(w, lat, lon)
        if 'height' in out:
            out['height'] = np.nan_to_num(np.clip(out['height'], 0., None), nan=0.)
        return out

    @staticmethod
    def _up_axes(w, r1, j0, j1):
        return 4 * np.arange(r1) + 1.5, w.rps, 4 * np.arange(j0, j1) + 1.5, np.arange(w.c0, w.c1)

    def _view_angles(self, w, lat, lon):
        g = w.granule
        tie = w.tie
        if tie is None:
            coarse = self._coarse_from_geo(g, w.kind, w.shape)
            if coarse is None:
                return np.zeros_like(lat, dtype=np.float64), np.zeros_like(lat, dtype=np.float64)
            tie = coarse
        tlat, tlon, rstep, cstep = tie
        nlat, nlon = fc.nadir_track(tlat, tlon)
        # 星下點換到每一列 (直角座標內插)
        rows = np.arange(w.s0 * w.rps, w.s1 * w.rps)
        tpos = np.arange(tlat.shape[0]) * rstep
        la, lo = np.deg2rad(nlat), np.deg2rad(nlon)
        xyz = [np.interp(rows, tpos, v) for v in (np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la))]
        rl = np.rad2deg(np.arctan2(xyz[2], np.hypot(xyz[0], xyz[1])))
        rn = np.rad2deg(np.arctan2(xyz[1], xyz[0]))
        return fc.view_angles(lat, lon, rl, rn, SAT_HEIGHT.get(g.satellite, 836.))

    # ---- 定標 --------------------------------------------------------------
    @staticmethod
    def _vis_coeffs(f):
        ds = fc.find_ds(f, 'VIS_Cal_Coeff', ('Calibration', ''))
        if ds is None:
            ds = fc.find_ds(f, 'RSB_Cal_Coeff', ('Calibration', ''))
        if ds is not None:
            c = ds[:].astype(np.float64)
            slope, inter = fc.scaled(ds)
            return (c * slope + inter).reshape(-1, c.shape[-1])
        for key in ('VIR_Cal_Coeff', 'VIS_Cal_Coeff'):           # FY-3A/B：全域屬性 (19 x 3)
            v = fc.attr_array(f, key)
            if v is not None and v.size % 3 == 0:
                return v.reshape(-1, 3)
        raise ValueError(tr('%s 沒有反射通道的定標係數', '%s has no reflective-band calibration coefficients')
                         % f.filename)

    def _ir_params(self, f, g, b):
        """紅外頻道的 (等效中心波數, A, B)。"""
        sensor, sat = g.sensor, g.satellite
        if sensor == 'MERSI-1':
            wn, a, bb = MERSI1_IR
        else:
            tab = IR_TABLE.get(sat)
            n_ir = sum(1 for x in BANDS[sensor] if x.kind == 'bt')
            nb = len(BANDS[sensor])
            ecw = fc.attr_array(f, 'Effect_Center_WaveLength')
            if ecw is None:
                ds = fc.find_ds(f, 'Effect_Center_WaveLength', ('Calibration',))
                if ds is None:
                    ds = fc.find_ds(f, 'Effect_Center_Wave_Length', ('Calibration',))
                ecw = ds[:].astype(np.float64).ravel() if ds is not None else None
            wl = ecw[b.num - 1] if ecw is not None and ecw.size >= nb else np.nan
            if abs(wl - b.wl) < 0.1 * b.wl:
                wn = 1e4 / wl
            else:
                wn = tab[0][b.cal] if tab else 1e4 / b.wl
            coef = None
            va, vb = fc.attr_array(f, 'TBB_Trans_Coefficient_A'), fc.attr_array(f, 'TBB_Trans_Coefficient_B')
            if va is not None and vb is not None and va.size >= n_ir and vb.size >= n_ir:
                coef = va[b.cal], vb[b.cal]
            else:
                v = fc.attr_array(f, 'TBB_Trans_Coefficient')
                if v is not None and v.size >= 2 * n_ir:
                    coef = v[b.cal], v[n_ir + b.cal]
            if coef is None or not (0.95 < coef[0] < 1.05 and abs(coef[1]) < 5.):
                coef = (tab[1][b.cal], tab[2][b.cal]) if tab else (1., 0.)
            a, bb = coef
        return float(wn), float(a), float(bb)

    def _read_band(self, w, b):
        """這個範圍的某頻道：反射率 (0~1)、亮溫 (°C) 或微光輻亮度 (W m-2 sr-1)。"""
        g = w.granule
        name, idx = b.source(w.kind)
        with h5py.File(g.files[w.kind], 'r') as f:
            ds = fc.find_ds(f, name)
            if ds is None:
                raise IOError(tr('%s 沒有 %s', '%s has no %s') % (g.files[w.kind], name))
            raw = ds[idx, w.rows, w.cols] if idx is not None else ds[w.rows, w.cols]
            slope, inter = fc.scaled(ds, idx)
            vr = fc.attr_array(ds, 'valid_range')
            if b.kind == 'll':
                out = raw.astype(np.float64) * slope + inter
                hi = vr[1] if vr is not None and vr.size == 2 else 1e4
                out[~((raw >= 0) & (raw <= hi))] = np.nan
                return out
            if b.kind == 'bt':
                if g.sensor == 'MERSI-1' and slope == 1.:
                    slope = 0.01                        # MERSI-1 第 5 頻道的 Slope 誤植 (同 satpy)
                wn, a, bb = self._ir_params(f, g, b)
                rad = raw.astype(np.float64) * slope + inter
                tb = a * fc.planck_bt(rad, wn) + bb
                tb[(raw >= 65533) | (raw <= 0)] = np.nan
                return tb - 273.15
            coef = self._vis_coeffs(f)[b.cal]
            vmax = min(vr[1], 65532.) if vr is not None and vr.size == 2 else 4095.
            dn = raw.astype(np.float64)
            satur = raw == 65534
            bad = (raw >= 65533) & ~satur
            dn[satur] = vmax
            dn = dn * slope + inter
            ref = coef[0] + coef[1] * dn + (coef[2] * dn * dn if coef.size > 2 else 0.)
            if g.sensor != 'MERSI-RM':
                ref = ref / 100.
            ref[bad] = np.nan
            return ref

    @staticmethod
    def _effective_wl(f, b):
        ecw = fc.attr_array(f, 'Effect_Center_WaveLength')
        if ecw is not None and ecw.size >= b.num and abs(ecw[b.num - 1] - b.wl) < 0.05 * b.wl:
            return float(ecw[b.num - 1])
        return b.wl

    # ---- 資料 --------------------------------------------------------------
    def prepare(self):
        typ, kind, ref = self.option
        b = self._band_for(self.sensor, ref) if typ == 'band' else None
        if typ == 'band' and (b is None or not b.source(kind)):
            raise ValueError(tr('%s 沒有 %s (%s)', '%s has no %s (%s)') % (self.sensor, ref.text, RES_LABEL[kind]))
        if typ == 'tc' and self.sensor not in ('MERSI-1', 'MERSI-II', 'MERSI-III'):
            raise ValueError(tr('%s 沒有真實色需要的頻道', '%s has no true color channels') % self.sensor)
        grans = [g for g in self.granules if kind in g.files]
        if not grans:
            raise ValueError(tr('這一段沒有 %s 的資料檔', 'this segment has no %s data file') % kind)
        if len(grans) < len(self.granules):
            print(tr('注意：%d 個片段沒有 %s 的檔，略過。', 'NOTE: %d granule(s) have no %s file; skipped.')
                  % (len(self.granules) - len(grans), kind))
        wins = [w for w in (self._window(g, kind) for g in grans) if w is not None]
        if not wins:
            return False
        c0, c1 = min(w.c0 for w in wins), max(w.c1 for w in wins)
        for w in wins:
            w.c0, w.c1 = c0, c1
        if len(wins) < len(grans) or len(grans) > 1:
            print(tr('讀取 %d / %d 個片段的 %d 個掃描帶、第 %d~%d 行',
                     'Reading %d of %d granules: %d scans, columns %d-%d')
                  % (len(wins), len(grans), sum(w.s1 - w.s0 for w in wins), c0, c1 - 1))

        lats, lons, rowtime = [], [], []
        for w in wins:
            la, lo = self._geo(w)
            lats.append(la)
            lons.append(lo)
            rowtime.append(self._times(w))
        lats, lons, rowtime = np.concatenate(lats), np.concatenate(lons), np.concatenate(rowtime)
        notes = set().union(*(w.notes for w in wins))
        if self.slice_swath(lats, lons) is None:
            return False
        self.time = self.center_time(rowtime)
        self.wins = wins
        self.typ, self.band = typ, b
        self._tc = None
        self._ll = {}

        if typ == 'tc':
            refl, wls = {}, {}
            bands = {x.num: x for x in BANDS[self.sensor]}
            for c, num in TC_BANDS.items():
                refl[c] = self.crop(self.concat([self._read_band(w, bands[num]) for w in wins]))
                with h5py.File(wins[0].granule.files[kind], 'r') as f:
                    wls[c] = self._effective_wl(f, bands[num])
            ang = [self._angles(w, la, lo, t) for w, la, lo, t in self._split(wins, lats, lons, rowtime)]
            angles = {k: self.crop(np.concatenate([a[k] for a in ang])) for k in ('sza', 'saz', 'vza', 'vaz')}
            height = self.crop(np.concatenate([a['height'] for a in ang])) if all('height' in a for a in ang) else None
            with h5py.File(wins[0].granule.files[kind], 'r') as f:
                des = fc.attr_array(f, 'EarthSun Distance Ratio')
            des = float(des[0]) if des is not None and 0.9 < des[0] < 1.1 else 1.
            self._tc = fc.TrueColor(refl, wls, des, (angles['sza'], angles['saz'], angles['vza'], angles['vaz']),
                                    height)
        else:
            self.data = self.crop(self.concat([self._read_band(w, b) for w in wins]))
            if b.kind == 'll':
                ang = [self._angles(w, la, lo, t, need=('sza', 'lza'))
                       for w, la, lo, t in self._split(wins, lats, lons, rowtime)]
                self.sza = self.crop(np.concatenate([a['sza'] for a in ang]))
                self.lza = self.crop(np.concatenate([a['lza'] for a in ang])) if all('lza' in a for a in ang) else None
                with h5py.File(wins[0].granule.files[kind], 'r') as f:
                    ds = fc.find_ds(f, 'Solar_Irradiance_LL', ('Calibration', ''))
                    self.e0 = float(ds[:].ravel()[0]) if ds is not None else 429.
            elif b.kind == 'bt':
                self.report_eyetemp(self.data)
        notes |= set().union(*(w.notes for w in wins))
        for n in sorted(notes):
            print(tr('注意：%s', 'NOTE: %s') % n)
        return True

    @staticmethod
    def _split(wins, lats, lons, rowtime):
        """把串接後的經緯度、時間切回各片段。"""
        k = 0
        for w in wins:
            n = (w.s1 - w.s0) * w.rps
            yield w, lats[k:k + n], lons[k:k + n], rowtime[k:k + n]
            k += n

    def ll_scaled(self, mode):
        if mode in self._ll:
            return self._ll[mode]
        rad = self.data
        if mode == 3:
            out = np.pi * rad / self.e0
        else:
            lza = self.lza
            if lza is None:
                print(tr('注意：定位檔沒有月亮天頂角，ERF 以 90 度計算。',
                         'NOTE: no lunar zenith angle in the geolocation file; ERF assumes 90 degrees.'))
                lza = np.full(rad.shape, 90.)
            phase, frac = dnb.moon_phase_angle(self.time)
            print(tr('LL ERF：月亮相位角 %.1f 度 (照明 %.0f%%)', 'LL ERF: lunar phase angle %.1f deg (%.0f%% lit)')
                  % (phase, frac * 100))
            out = dnb.erf_dnb(rad * 1e-4, self.sza, lza, frac, saturation=(mode == 2), moon_phase=phase)
        self._ll[mode] = out
        return out

    def plot(self, canvas, variant):
        if self.typ == 'tc':
            canvas.rgb(self.lons, self.lats, self._tc.rgb(variant))
            return dict(band='TrueColor', colorbar=fc.TC_STYLES[variant], label='TrueColor')
        b = self.band
        if b.kind == 'll':
            canvas.mesh(self.lons, self.lats, self.ll_scaled(variant), cmap='gray', vmin=0., vmax=1.)
            return dict(band=b.name, colorbar=LL_NAMES[variant], label=b.name)
        if b.kind == 'ref':
            canvas.mesh(self.lons, self.lats, np.ma.masked_greater(np.ma.masked_invalid(self.data), 1.5),
                        cmap='gray', vmin=0., vmax=1.)
            return dict(band=b.name, colorbar='', label=b.name)
        cm = canvas.bt_mesh(self.lons, self.lats, self.data, variant)
        return dict(band=b.name, colorbar=cm, label=b.name)


READER = MersiReader
