"""
AVHRR-NAT — Metop-A/B/C AVHRR/3 L1B (EUMETSAT EPS native .nat，zip 壓縮)
原始腳本：SatimaV3METOP-NAT-loop.py   紀錄格式：avhrr_l1b_nat_format.py

所需檔案
  AVHR_xxx_1B_M0n_<開始>Z_<結束>Z_N_O_<處理時間>Z.zip  (內含同名 .nat)

頻道
  1.CH1 0.630  2.CH2 0.862  3.CH3A 1.61  (反射率 πL/F，未除 cos(太陽天頂角)，同原程式)
  4.CH3B 3.74  5.CH4 10.8  6.CH5 12.0    (亮溫)

與原腳本的差異 (原腳本的錯誤)
  - 可見光頻道原本呼叫 VISConvert(nc_obj, raw, band)，變數不存在會 NameError；改為 GIADR 與輻射值
  - SCENE_RADIANCES 的縮放：CH3A/3B 為 1e4 (其餘 1e2)，原本一律除 100，CH3B 亮溫會錯
  - GIADR 的 CH3B 中心波數縮放為 1e2 (CH4/CH5 為 1e3)，原本一律除 1000，CH3B 亮溫約 -212°C 而全被遮蔽
  - 亮溫的頻帶修正為 T = A + B·T* (Product Guide §6.3.4)，原本寫成 (T* - A) / B，
    亮溫偏高 (CH3B 約 1.6-2.4 K，CH4 約 0.1-0.6 K，CH5 約 0-0.3 K，愈冷差愈多)
  - CH3A/3B 共用同一個資料欄位、依掃描線切換；依 FRAME_INDICATOR 第 16 位元遮蔽另一模式的掃描線
  - 依各紀錄的 GRH (RECORD_CLASS / RECORD_SIZE) 讀檔，不再假設 VEADR 只有 1 筆、MDR 中沒有 dummy 紀錄
  - 移除每條掃描線都 print 一次 SCENE_RADIANCES 的除錯輸出
  - 中心時間：取最接近輸入中心點的像元，用該掃描線 GRH 的紀錄時間
    (原本以檔名起訖時間依裁切範圍的中間列線性推算，軌道中有資料缺口 (dummy MDR) 時會偏移)
"""
import zipfile
from datetime import datetime

import numpy as np

from avhrr_l1b_nat_format import (_Generic_Record_Header, _Main_Product_Header_Record,
                                  _Secondary_Product_Header_Record,
                                  _Global_Internal_Auxiliary_Data_Record_Subclass_1,
                                  _Measurement_Data_Record)
from satima.po import SatReader, ui
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr

BAND_MENU = '1.BAND-1(0.630) 2.BAND-2(0.862) 3.BAND-3A 4.BAND-3B 5.BAND-4(10.8) 6.BAND5(12)'
BAND_LABEL = {1: '01', 2: '02', 3: '03A', 4: '03B', 5: '04', 6: '05'}   # 檔名/圖說的頻道名稱
BT_BANDS = (4, 5, 6)
# 頻道 -> SCENE_RADIANCES 的索引 (CH3A / CH3B 共用索引 2，依掃描線切換)
RADIANCE_INDEX = {1: 0, 2: 1, 3: 2, 4: 2, 5: 3, 6: 4}
# SCENE_RADIANCES 的縮放 (EPS AVHRR/3 L1B 格式：10^2, 10^2, 10^4, 10^2, 10^2)
RADIANCE_SCALE = (1e2, 1e2, 1e4, 1e2, 1e2)

# EPS 紀錄類別 (GRH RECORD_CLASS)
MPHR, SPHR, IPR, GEADR, GIADR, VEADR, VIADR, MDR = range(1, 9)

# ---------------------------------------------------------------------------
# Tie-point -> pixel navigation (pure numpy, no geotiepoints)
# EUMETSAT AVHRR L1B Product Guide (EUM/OPS-EPS/MAN/04/0029) p.36, p.72-73:
#   EARTH_LOCATION_FIRST     -> earth view 1     (0-based column 0)
#   EARTH_LOCATIONS (NP=103) -> every 20th point, points 5..2045 (0-based 4, 24, ..., 2044)
#   EARTH_LOCATION_LAST      -> earth view 2048  (0-based column 2047)
# p.17: pixels are geolocated by linear or Lagrangian interpolation.
# ---------------------------------------------------------------------------
AVHRR_TIEPOINT_COLS = np.r_[0, np.arange(4, 2048, 20), 2047]


def _lagrange_weights(nodes, targets, order):
    """Local Lagrange weights for (possibly non-uniform) nodes.

    For each target, take `order` consecutive nodes around the bracketing
    segment (shifted inwards at the edges). Returns (node index, weight),
    both shaped (len(targets), order).
    """
    nodes = np.asarray(nodes, dtype=np.float64)
    targets = np.asarray(targets, dtype=np.float64)
    n = nodes.size
    seg = np.clip(np.searchsorted(nodes, targets, side='right') - 1, 0, n - 2)
    start = np.clip(seg - (order // 2 - 1), 0, n - order)
    idx = start[:, None] + np.arange(order)
    xn = nodes[idx]
    w = np.ones_like(xn)
    for m in range(order):
        for k in range(order):
            if k != m:
                w[:, m] *= (targets - xn[:, k]) / (xn[:, m] - xn[:, k])
    return idx, w


def interpolate_avhrr_lonlat_lagrange(lons, lats, nav_sample_rate=20,
                                      pixels_per_scanline=2048, order=4,
                                      chunk_size=512):
    """Expand EPS AVHRR tie-point navigation to all 2048 pixels; return (lons, lats).

    order=4: Lagrangian cubic (default), order=2: linear.
    Interpolation is done on unit-sphere Cartesian coordinates (x, y, z) so it
    is safe across the dateline and near the poles; lon returned in [-180, 180].
    Inputs are degrees (already divided by 1e4), shaped (scanlines, 105).
    """
    lons = np.asarray(lons, dtype=np.float64)
    lats = np.asarray(lats, dtype=np.float64)
    if nav_sample_rate != 20 or pixels_per_scanline != 2048:
        raise ValueError(tr('只支援 NAV_SAMPLE_RATE=20、每條掃描線 2048 個像元',
                            'Only NAV_SAMPLE_RATE=20 and 2048 pixels are supported'))
    if lons.ndim != 2 or lons.shape != lats.shape or lons.shape[1] != AVHRR_TIEPOINT_COLS.size:
        raise ValueError(tr('經度與緯度的形狀都必須是 (掃描線數, 105)',
                            'Longitude and latitude must both have shape (scanlines, 105)'))
    if order < 2 or order > AVHRR_TIEPOINT_COLS.size:
        raise ValueError(tr('order 須 >= 2 (2 = 線性，4 = 三次 Lagrange)',
                            'order must be >= 2 (2 = linear, 4 = cubic Lagrange)'))
    if not isinstance(chunk_size, (int, np.integer)) or chunk_size < 1:
        raise ValueError(tr('chunk_size 須為正整數', 'chunk_size must be a positive integer'))
    if (not np.all(np.isfinite(lons)) or not np.all(np.isfinite(lats))
            or np.any(np.abs(lons) > 180) or np.any(np.abs(lats) > 90)):
        raise ValueError(tr('定位的 tie point 無效；請檢查缺值與 1e-4 的縮放',
                            'Invalid navigation tiepoints; check missing values and 1e-4 scaling'))

    idx, w = _lagrange_weights(AVHRR_TIEPOINT_COLS, np.arange(pixels_per_scanline), order)

    nlines = lons.shape[0]
    out_lons = np.empty((nlines, pixels_per_scanline), dtype=np.float64)
    out_lats = np.empty_like(out_lons)
    for start in range(0, nlines, chunk_size):
        stop = min(start + chunk_size, nlines)
        lon_r = np.deg2rad(lons[start:stop])
        lat_r = np.deg2rad(lats[start:stop])
        coslat = np.cos(lat_r)
        x = (coslat * np.cos(lon_r))[:, idx]   # (lines, 2048, order)
        y = (coslat * np.sin(lon_r))[:, idx]
        z = np.sin(lat_r)[:, idx]
        x = np.einsum('lpk,pk->lp', x, w)
        y = np.einsum('lpk,pk->lp', y, w)
        z = np.einsum('lpk,pk->lp', z, w)
        out_lons[start:stop] = np.rad2deg(np.arctan2(y, x))
        out_lats[start:stop] = np.rad2deg(np.arctan2(z, np.hypot(x, y)))
    return out_lons, out_lats


# ===== 單一檔案 =====
class AvhrrNatFile:

    def __init__(self, fullname, channel=None, time=None):
        self.fullname = fullname

        spacecrafts = {'M01': 'METOP-B', 'M02': 'METOP-A', 'M03': 'METOP-C'}

        self.satellite = spacecrafts.get(self.fullname[12:15], self.fullname[12:15])
        self.channel = 'IR/VIS'
        if time:
            self.time = datetime.strptime(time, '%Y%m%d%H%M%S')
        else:
            self.time = datetime.strptime(self.fullname[16:30], '%Y%m%d%H%M%S')
        try:                    # 檔名第 2 個時間是結束時間：AVHR_xxx_1B_M01_<起>Z_<訖>Z_…
            self.end = datetime.strptime(self.fullname[32:46], '%Y%m%d%H%M%S')
        except ValueError:
            self.end = None

    def __lt__(self, other):
        return self.time < other.time

    @staticmethod
    def grh_time(raw):
        """GRH 的 RECORD_START_TIME (CDS short：2000-01-01 起的日數 uint16 + 當日毫秒 uint32，UTC) -> datetime64。"""
        b = np.frombuffer(np.ascontiguousarray(raw).tobytes(), dtype=np.uint8).reshape(-1, 6).astype(np.int64)
        day = (b[:, 0] << 8) | b[:, 1]
        ms = (b[:, 2] << 24) | (b[:, 3] << 16) | (b[:, 4] << 8) | b[:, 5]
        return np.datetime64('2000-01-01', 'us') + (day * 86_400_000 + ms).astype('timedelta64[ms]')

    @staticmethod
    def read_records(buf):
        """依 GRH 逐筆讀取 EPS 紀錄；回傳 dict(mphr, sphr, giadr1, mdr, mdr_time)。
        mdr 為 _Measurement_Data_Record 結構陣列 (每列一條掃描線)，大小不符的 dummy MDR 略過；
        mdr_time 為各 MDR 的 GRH 紀錄開始時間 (UTC datetime64)。"""
        grh_size = _Generic_Record_Header.itemsize
        mdr_size = _Measurement_Data_Record.itemsize
        out = dict(mphr=None, sphr=None, giadr1=None)
        offsets, skipped = [], 0
        pos, n = 0, len(buf)
        while pos + grh_size <= n:
            grh = np.frombuffer(buf, dtype=_Generic_Record_Header, count=1, offset=pos)[0]
            size = int(grh['RECORD_SIZE'])
            if size < grh_size or pos + size > n:
                print(tr('注意：第 %d 位元組處的紀錄不完整，停止讀取。',
                         'WARNING: truncated record at byte %d; stop reading.') % pos)
                break
            cls_, sub = int(grh['RECORD_CLASS']), int(grh['RECORD_SUBCLASS'])
            body = pos + grh_size
            if cls_ == MPHR:
                out['mphr'] = np.frombuffer(buf, dtype=_Main_Product_Header_Record, count=1, offset=body)
            elif cls_ == SPHR:
                out['sphr'] = np.frombuffer(buf, dtype=_Secondary_Product_Header_Record, count=1, offset=body)
            elif cls_ == GIADR and sub == 1:
                out['giadr1'] = np.frombuffer(buf, dtype=_Global_Internal_Auxiliary_Data_Record_Subclass_1,
                                              count=1, offset=body)
            elif cls_ == MDR:
                if size - grh_size == mdr_size:
                    offsets.append(body)
                else:
                    skipped += 1
            pos += size
        if skipped:
            print(tr('略過 %d 筆 dummy / 不明的 MDR 紀錄', 'Skipped %d dummy/unknown MDR records') % skipped)
        if not offsets:
            raise ValueError(tr('找不到 AVHRR 的 MDR 紀錄', 'No AVHRR MDR records found'))
        offsets = np.asarray(offsets)
        if offsets.size > 1 and np.all(np.diff(offsets) == mdr_size + grh_size):
            # 連續存放：一次讀入 (GRH + MDR) 結構陣列
            rec = np.dtype([('grh', _Generic_Record_Header), ('mdr', _Measurement_Data_Record)])
            recs = np.frombuffer(buf, dtype=rec, count=offsets.size, offset=int(offsets[0]) - grh_size)
            out['mdr'] = recs['mdr']
            grh = recs['grh']
        else:
            out['mdr'] = np.concatenate([np.frombuffer(buf, dtype=_Measurement_Data_Record, count=1, offset=int(o))
                                         for o in offsets])
            grh = np.concatenate([np.frombuffer(buf, dtype=_Generic_Record_Header, count=1, offset=int(o) - grh_size)
                                  for o in offsets])
        out['mdr_time'] = AvhrrNatFile.grh_time(grh['RECORD_START_TIME'])
        return out

    def extract(self, band):

        with zipfile.ZipFile(self.fullname, 'r') as zf:
            names = zf.namelist()
            member = self.fullname[:-4] + '.nat'
            if member not in names:
                nat = [s for s in names if s.lower().endswith('.nat')]
                if not nat:
                    raise FileNotFoundError(tr('%s 裡面找不到 .nat 檔', 'No .nat file inside %s') % self.fullname)
                member = nat[0]
            buf = zf.read(member)
        rec = self.read_records(buf)

        # MPHR
        mphr_blk = rec['mphr']
        for key in ('TOTAL_RECORDS', 'TOTAL_MPHR', 'TOTAL_SPHR', 'TOTAL_IPR', 'TOTAL_GEADR',
                    'TOTAL_GIADR', 'TOTAL_VEADR', 'TOTAL_VIADR', 'TOTAL_MDR'):
            print('%s: %s' % (key, mphr_blk[key][0].decode('utf-8').strip()[-6:].strip()))

        # SPHR
        sphr_blk = rec['sphr']
        self.nav_sample_rate = int(sphr_blk['NAV_SAMPLE_RATE'][0].decode().split('=')[1])
        self.pixels_per_scanline = int(sphr_blk['EARTH_VIEWS_PER_SCANLINE'][0].decode().split('=')[1])
        if self.nav_sample_rate != 20 or self.pixels_per_scanline != 2048:
            raise ValueError(tr('這個讀取器只能讀 NAV_SAMPLE_RATE=20、每條掃描線 2048 個像元的檔案',
                                'This reader requires NAV_SAMPLE_RATE=20 and 2048 pixels'))

        # GIADR#1 (輻射 -> 反射率 / 亮溫係數)
        giadr_1_blk = rec['giadr1']
        if giadr_1_blk is None:
            raise ValueError(tr('找不到 GIADR (輻射換算) 紀錄', 'GIADR (radiance conversion) record not found'))

        # MDR
        mdr = rec['mdr']
        total_mdr = mdr.shape[0]
        print('MDR read: %d' % total_mdr)
        self._rowtime = rec['mdr_time']
        if np.any(mdr['NUM_NAVIGATION_POINTS'] != 103):
            raise ValueError(tr('定位點的數目不是預期的 103', 'Unexpected navigation point count'))
        self.lats = np.empty((total_mdr, 105), dtype=np.float64)
        self.lons = np.empty((total_mdr, 105), dtype=np.float64)
        self.lats[:, 0], self.lons[:, 0] = mdr['EARTH_LOCATION_FIRST'][:, 0], mdr['EARTH_LOCATION_FIRST'][:, 1]
        self.lats[:, 1:-1] = mdr['EARTH_LOCATIONS'][:, 0::2]
        self.lons[:, 1:-1] = mdr['EARTH_LOCATIONS'][:, 1::2]
        self.lats[:, -1], self.lons[:, -1] = mdr['EARTH_LOCATION_LAST'][:, 0], mdr['EARTH_LOCATION_LAST'][:, 1]

        idx = RADIANCE_INDEX[band]
        data = mdr['SCENE_RADIANCES'][:, idx, :].astype(np.float64)
        data /= RADIANCE_SCALE[idx]

        if band in (3, 4):
            # FRAME_INDICATOR (32 位元，此處存為 2 個 int16) 第 16 位元：1 = CH3A，0 = CH3B
            is3a = (mdr['FRAME_INDICATOR'][:, 0].astype(np.int32) & 1) == 1
            other = ~is3a if band == 3 else is3a
            print(tr('CH3A 掃描線：%d，CH3B 掃描線：%d', 'CH3A scanlines: %d, CH3B scanlines: %d')
                  % (is3a.sum(), (~is3a).sum()))
            if other.all():
                print(tr('注意：這一軌沒有 %s 的掃描線。', 'WARNING: this orbit has no %s scanlines.')
                      % ('CH3A' if band == 3 else 'CH3B'))
            data[other] = np.nan

        if band <= 3:
            out = self.VISConvert(giadr_1_blk, data, band)
        else:
            with np.errstate(invalid='ignore', divide='ignore'):
                out = self.BTConvert(giadr_1_blk, data, band)

        return out

    def BTConvert(self, blk, R, band):

        channel = {4: 'CH3B', 5: 'CH4', 6: 'CH5'}

        C1 = 1.191044E-5
        C2 = 1.438769

        # 中心波數縮放：CH3B 為 1e2 (例 266432 -> 2664.32 cm-1)，CH4/CH5 為 1e3 (原本一律除 1000，CH3B 會錯)
        wnc = blk[f'{channel[band]}_CENTRAL_WAVENUMBER'][0] / (100 if band == 4 else 1000)
        A = blk[f'{channel[band]}_CONSTANT1'][0] / 100000
        B = blk[f'{channel[band]}_CONSTANT2_SLOPE'][0] / 1000000

        # AVHRR L1B Product Guide §6.3.4：T* = C2·γ / ln(1 + C1·γ³/R)，T = A + B·T*
        # (原本寫成 (T* - A) / B，CH3B 高估約 1.6-2.4 K、CH4 0.1-0.6 K、CH5 0-0.3 K)
        bt = A + B * (C2 * wnc) / np.log(1. + (C1 * wnc ** 3) / R) - 273.15

        return bt

    def VISConvert(self, blk, R, band):

        channel = {1: 'CH1', 2: 'CH2', 3: 'CH3A'}

        f_sol = blk[f'{channel[band]}_SOLAR_FILTERED_IRRADIANCE'][0] / 10

        refl = (R * np.pi) / f_sol

        return refl

    def geocoord(self, band, order=4):
        # Tie points (scanlines, 105) are scaled by 1e4 -> degrees.
        # Expand to (scanlines, 2048) per AVHRR L1B Product Guide
        # (order=4: Lagrangian cubic, order=2: linear). Pure numpy, no geotiepoints.
        lats = self.lats / 10000.
        lons = self.lons / 10000.
        lons, lats = interpolate_avhrr_lonlat_lagrange(
            lons, lats, self.nav_sample_rate, self.pixels_per_scanline, order=order)
        return lats, lons

    def plot(self, canvas, lats, lons, data, variant, band):
        if data.size == 0:
            return ''
        if band in (1, 2, 3):
            canvas.mesh(lons, lats, data, cmap='gray', vmin=0., vmax=1.)
            return ''
        return canvas.bt_mesh(lons, lats, data, variant)

    def row_times(self):
        """每條掃描線的觀測時間 (UTC)：各 MDR 的 GRH 紀錄開始時間 (須先呼叫 extract)。"""
        return self._rowtime


# ===== 讀取器 =====
class AvhrrNatReader(SatReader):
    NAME = 'AVHRR-NAT'
    SENSOR = 'AVHRR3'
    DESCRIPTION = 'Metop AVHRR/3 L1B (EPS 原生格式，zip)'
    DESCRIPTION_EN = 'Metop AVHRR/3 L1B (EPS native, zip)'
    SOURCE = 'EUMETSAT'
    PATTERNS = ('AVHR_*_1B*.zip',)
    MULTI_FILE = False

    def make_file(self, name):
        # AVHR_xxx_1B_M01_20260412093103Z_20260412111303Z_N_O_20260412102021Z.zip
        return AvhrrNatFile(name, time=name[16:30])

    def listing(self, f):
        return (f.satellite, f.channel, f.time)

    def select_band(self):
        self.file = self.files[0]
        print(self.file.fullname)
        print(BAND_MENU)
        self.band = ui.ask_int(tr('輸入頻道編號：', 'Input the band: '), range(1, 7))

    def select_variants(self):
        if self.band in BT_BANDS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare(self):
        f = self.file
        data = f.extract(self.band)
        lats, lons = f.geocoord(self.band)
        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        self.data = self.crop(data)
        self.time = self.center_time(f.row_times())
        if self.band in BT_BANDS:
            self.report_eyetemp(self.data)
        return True

    def slice(self, georange, lats, lons):
        return self.slice_swath(lats, lons)

    def plot(self, canvas, variant):
        cm = self.file.plot(canvas, self.lats, self.lons, self.data, variant, self.band)
        channel = BAND_LABEL[self.band]
        return dict(band='B%s' % channel, colorbar=cm, label='BAND%s' % channel)


READER = AvhrrNatReader
