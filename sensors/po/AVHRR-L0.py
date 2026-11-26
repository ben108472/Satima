"""
AVHRR-L0 — TIROS-N ~ NOAA-19 AVHRR 的 HRPT 原始框 (Level 0：接收站直接收錄的 HRPT 小框，例如千葉大學 CEReS 的 .raw.xz)
框結構、檔案辨識與定位見 avhrr_hrpt.py；定標、時間碼修復、多來源合併沿用 AVHRR-L1B。

所需檔案
  每條掃描線一框 (11090 個 10-bit 字，每字存成 16 bits，小端或大端)；未壓縮或 .xz / .gz / .bz2，也可以直接放
  .tar 封存檔 (列表會讀出其中的每個檔，不必先解開)。資料目錄往下兩層的子資料夾也會搜尋。不依檔名辨識：
  依框同步字判斷，衛星由 ID 字的衛星位址判斷 (AAPP 的對照表；同一位址依年份區分，路徑中的 noaa12、n12 等字樣優先)。
  原始框沒有年份：取路徑中的日期 (例如 2003-02-04/、20130811.063140.noaa18.raw)，沒有時取檔案的修改時間。
  同一個檔同時在 .tar 與資料夾中時只列一次。
  千葉大學 CEReS 的 L0 (20141007163202.n19.l0.bz2、n1297121720.gi.bz2)：開頭 240 bytes 的 ANSI 磁帶標籤 (vol1/hdr1/hdr2)
  辨識時跳過，年份取自這兩種檔名 (原本的日期寫法都對不上時才用)；細節見 avhrr_hrpt.py。

定位 (需要 TLE 與 sgp4 或 pyorbital 套件)
  - 原始框沒有軌道資料：以 [tle] dir (預設 Satima.py 旁的 tle/) 與資料目錄中的 TLE 推算衛星位置 (SGP4)，取 epoch 最接近
    觀測時間的一組；相差超過 [tle] max_age 天時警告 (TLE 推算的誤差約每天 1 km 以上)。TLE 可以是任何副檔名的
    文字檔 (2 行或 3 行格式，一個檔可含多顆衛星、多個時間) 或 .zip。
  - TLE 來源：1980~2004 年 CelesTrak 的歷史檔 (https://celestrak.org/NORAD/archives/，tle/ 已附 TIROS-N~NOAA-17)；
    2005 年以後需從 Space-Track (https://www.space-track.org，免費帳號，gp_history) 或 CelesTrak 的申請表取得。
  - 掃描幾何：大地天底 (KLM Guide §2.2、AAPP)、掃描面垂直於慣性速度、每像元 25 µs，與 WGS84 橢球求交點；在 51 個
    定位點計算後同 AVHRR-L1B 內插。各衛星的 roll、yaw 與步進角 (avhrr_hrpt.NAV_PARAMS) 以海岸線擬合：NOAA-18 的 roll
    約 −0.18°、NOAA-19 的 yaw 約 0.33°、NOAA-16 的步進角 0.05396° (其餘 0.05406~0.05409°)，不修正時掃描邊緣差到 5~10 km。
    沒有擬合的衛星用 AAPP 的步進角 0.0540865° 與共同的 roll −0.03°。不做時鐘修正 (NOAA 的星上時鐘誤差通常在 1 秒內，
    1 秒 ≈ 沿軌 6.6 km；樣本的擬合結果在 ±0.6 秒內)。
  - 驗證：海岸線擬合後的殘差 RMS 0.3~0.5 km (每顆衛星 2~4 次過境，涵蓋 2003 與 2013 年)。
  - 歷史 TLE 檔偶有互相矛盾的根數 (例如 CelesTrak 1988 年的 NOAA-9 有兩組 epoch 交錯、推算位置差 160 km)：
    最接近的三組推算到觀測時間的位置相差超過 TLE_SPREAD_KM 時警告。
時間
  時間碼 (年中日數 + 當日毫秒) 的修復同 AVHRR-L1B (殘差的中位數濾波)；框同步字對不到 4 個以上、KLM 的 ID 字
  顯示為測試訊號、或與前一框重複 (時間碼與小框編號相同) 的框為無資料。
定標
  - 紅外：一律由 PRT 黑體溫度、黑體與太空計數自行定標 (KLM Guide §7.1.2.4)，非線性修正、Planck 常數同 AVHRR-L1B；
    PRT 係數：NOAA-15 以後為 KLM 附錄 D，更早的衛星為 pygac 的值 (以 CLASS 1b 驗證，見 avhrr_l1b_noaa.PRT)。
  - 可見光：NOAA-16~19 用 NOAA STAR 的在軌係數 (conversion_table/AVHRR，每月一組，取 Center 最接近觀測者；見
    avhrr_l1b_noaa.star_vis)。NOAA-15 與 STAR 涵蓋期間以外用發射前係數 (NOAA-15 以後為附錄 D 的雙增益，更早為
    POD Guide 表 3.3.2-1)，沒有在軌衰減修正，反射率會偏低 (例如 2013 年的 NOAA-18 CH2 約低 20%)。
  - 3A / 3B：同 AAPP，頻道 3 的太空計數 (3A 約 40、3B 約 990) 且與 ID 字的 bit 0 (1 = 3A) 一致，再取 11 條的中位數。
"""
import importlib.util
import os
import sys

import numpy as np
from scipy.ndimage import median_filter

import avhrr_hrpt as hrpt
import avhrr_l1b_noaa as fmt
from satima.core.i18n import tr


def _load_l1b():
    """AVHRR-L1B 讀取器模組 (同一個 Satima 讀取器載入方式的模組名稱；已載入時直接用)。"""
    name = 'satima_sensor_AVHRR_L1B'
    mod = sys.modules.get(name)
    if mod is None:
        path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'AVHRR-L1B.py')
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
    return mod


l1b = _load_l1b()

TIE_COLS = 24. + 40. * np.arange(51)        # 定位點 (同 LAC/HRPT 1b)
IMAGE = slice(750, 10990)                   # 影像字 751~10990
TLE_EXT = {'.tle', '.txt', '.2le', '.3le'}  # 資料目錄中當作 TLE 讀的副檔名 (TLE 目錄則讀所有檔)
TLE_SPREAD_KM = 20.                         # 最接近的幾組 TLE 推算的位置相差超過此值時警告
_SNIFF = {}


def sniff(path):
    """hrpt.sniff 的快取 (資料夾掃描與建立檔案物件各做一次)。"""
    archive, member = hrpt.split_member(path)
    try:
        st = os.stat(archive)
    except OSError:
        return None
    key = (os.path.abspath(archive), member, st.st_size, st.st_mtime)
    if key not in _SNIFF:
        _SNIFF[key] = hrpt.sniff(path)
    return _SNIFF[key]


def tle_table(cfg):
    """[tle] dir (遞迴) 與資料目錄 (只有 TLE_EXT 的檔) 中的 TLE。"""
    folder = str(cfg.tle.dir or '').strip()
    if folder and not os.path.isabs(folder):
        folder = os.path.join(cfg.get('_home', ''), folder)
    tables = [hrpt.load_tles([folder])] if folder else []
    here = os.getcwd()
    try:
        names = [n for n in os.listdir(here) if os.path.splitext(n)[1].lower() in TLE_EXT]
    except OSError:
        names = []
    if names:
        tables.append(hrpt.load_tle_files([os.path.join(here, n) for n in names]))
    out = {}
    for t in tables:
        for num, items in t.items():
            out.setdefault(num, []).extend(items)
    return {num: sorted(set(items)) for num, items in out.items()}


# ===== 單一檔案 =====
class AvhrrL0File(l1b.AvhrrL1bFile):

    def __init__(self, fullname, cfg):
        s = sniff(fullname)
        if s is None:
            raise ValueError(tr('不是 HRPT 原始框', 'not HRPT raw frames'))
        s.kind, s.layout, s.version, s.station, s.header = 'LAC', 'L0', 0, True, b''
        self.fullname = fullname
        self.s = s
        self.cfg = cfg
        self.satellite = s.satellite
        self.generation = s.generation
        self.time = s.time
        self.end = s.end
        self.channel = 'HRPT'
        self.station = True
        self.nlines = s.nframes or 0
        self.npix = 2048
        self._nav = None
        self._cache = {}
        self._w = None
        self.label = os.path.basename(hrpt.split_member(fullname)[1] or fullname)
        self.table_dir = cfg.general.conversion_table_dir if cfg.get('general') else None

    # ---- 原始框 --------------------------------------------------------------------------
    def _frames(self):
        """(框數, 11090) uint16；第一次用到時解壓讀入。"""
        if self._w is None:
            raw = hrpt.read_bytes(self.fullname)
            self._w = hrpt.frame_words(raw, self.s.order, self.s.offset)
            self.nlines = self._w.shape[0]
        return self._w

    def release(self):
        """釋放原始框 (定位、遙測、3A/3B 已留下)。"""
        self._w = None

    # ---- 時間與定位 ------------------------------------------------------------------------
    def navigation(self):
        if self._nav is not None:
            return self._nav
        w = self._frames()
        good = hrpt.sync_count(w) >= hrpt.SYNC_MIN
        if self.generation == 3:
            good &= (w[:, 6] >> 1) & 1 == 1                 # KLM：0 = 測試用的偽隨機碼
        day, ms = hrpt.time_code(w)
        # 重複的框 (時間碼與小框編號都和前一框相同；接收站重送)，同 AAPP 的作法略過
        frame_no = (w[:, 6] >> 7) & 3
        dup = np.r_[False, (day[1:] == day[:-1]) & (ms[1:] == ms[:-1]) & (frame_no[1:] == frame_no[:-1])]
        good &= ~dup
        year = np.full(day.shape, self.s.year, dtype=np.int64)
        year[day < self.s.doy - 300] += 1                   # 跨年
        year[day > self.s.doy + 300] -= 1
        day = np.where(good, day, 0)                        # 同步不良的框：時間碼視為受損
        times, fixed = self._times(year, day, ms, None)
        tles = tle_table(self.cfg)
        near = hrpt.nearest_tles(tles, self.satellite, self.time)
        tle = near[0] if near else None
        if tle is None:
            raise ValueError(tr('%s：找不到 %s (NORAD %d) 的 TLE，無法定位。請把 TLE 檔放到 %s 或資料目錄 (1980~2004 年：'
                                'CelesTrak 歷史檔；2005 年以後：Space-Track)',
                                '%s: no TLE for %s (NORAD %d); cannot navigate. Put TLE files in %s or the data folder '
                                '(1980-2004: CelesTrak archives; 2005 on: Space-Track)')
                             % (self.label, self.satellite, hrpt.NORAD[self.satellite], self.cfg.tle.dir))
        epoch, line1, line2 = tle
        age = (self.time - epoch).total_seconds() / 86400.
        if age >= 0:
            msg = tr('%s：TLE epoch %s (觀測前 %.1f 天)', '%s: TLE epoch %s (%.1f days before the observation)')
        else:
            msg = tr('%s：TLE epoch %s (觀測後 %.1f 天)', '%s: TLE epoch %s (%.1f days after the observation)')
        print(msg % (self.label, epoch.strftime('%Y-%m-%d %H:%M'), abs(age)))
        if abs(age) > float(self.cfg.tle.max_age):
            print(tr('注意：TLE 與觀測相差 %.1f 天 (超過 [tle] max_age = %g)，定位誤差可能達數公里以上',
                     'NOTE: the TLE is %.1f days from the observation (over [tle] max_age = %g); navigation errors may '
                     'reach several km or more') % (abs(age), float(self.cfg.tle.max_age)))
        spread = hrpt.tle_spread(near, self.time)
        if spread > TLE_SPREAD_KM:
            print(tr('注意：最接近的 %d 組 TLE (epoch %s) 推算的衛星位置相差達 %.0f km，TLE 檔中可能有錯誤的根數，定位可能不準',
                     'NOTE: the %d nearest TLEs (epochs %s) disagree by up to %.0f km in satellite position; the TLE files '
                     'may contain bad elements and the navigation may be off')
                  % (len(near), ', '.join(e.strftime('%m-%d %H:%M') for e, _, _ in near), spread))
        roll, pitch, yaw, step = hrpt.NAV_PARAMS.get(self.satellite, hrpt.NAV_DEFAULT)
        lat, lon = hrpt.geolocate(line1, line2, times, TIE_COLS, attitude=(roll, pitch, yaw), step_deg=step)
        good &= np.isfinite(lat).all(axis=1)
        lat[~good], lon[~good] = np.nan, np.nan
        self._nav = dict(lat=lat, lon=lon, time=times, good=good)
        bad = int((~good).sum())
        if bad or fixed:
            print(tr('%s：%d 框無效 (同步不良或重複)、%d 條的時間碼受損 (改用推算值)',
                     '%s: %d invalid frames (bad sync or duplicates), %d with damaged time codes (replaced by the expected time)')
                  % (self.label, bad, fixed))
        return self._nav

    # ---- 計數值與遙測 ----------------------------------------------------------------------
    def counts(self, band, r0, r1):
        w = self._frames()[r0:r1, IMAGE].reshape(-1, self.npix, 5)
        return np.array(w[:, :, l1b.COUNT_INDEX[band]])

    def telemetry(self):
        if 'tm' not in self._cache:
            w = self._frames()
            prt, ict, space = w[:, 17:20], w[:, 22:52], w[:, 52:102]     # 字 18~20、23~52、53~102
            self._cache['tm'] = dict(prt=prt.astype(np.int64),
                                     ict={c: ict[:, k::3].astype(np.int64) for k, c in enumerate((3, 4, 5))},
                                     space={c: space[:, c - 1::5].astype(np.int64) for c in range(1, 6)})
        return self._cache['tm']

    def ch3_mode(self):
        """KLM：同 AAPP，以頻道 3 的太空計數判斷 (10 個樣本中至少 5 個 < 100 為 3A、> 900 為 3B)，並要與 ID 字 bit 0
        (1 = 3A) 一致；不一致、同步不良或無法判斷的掃描線以前後內插，再做 MODE_WINDOW 條的中位數濾波。"""
        if 'mode' not in self._cache:
            w = self._frames()
            n = w.shape[0]
            if self.generation < 3:
                is3a = np.zeros(n, dtype=bool)
            else:
                space3 = w[:, 52:102][:, 2::5]
                low = (space3 < 100).sum(axis=1) >= 5
                high = (space3 > 900).sum(axis=1) >= 5
                bit = (w[:, 6] & 1) == 1
                ok = (hrpt.sync_count(w) >= hrpt.SYNC_MIN) & (low != high) & (bit == low)
                is3a = (median_filter(l1b.fill_rows(low.astype(np.float64), ok), size=l1b.MODE_WINDOW,
                                      mode='nearest') > 0.5) if ok.any() else np.zeros(n, dtype=bool)
            self._cache['mode'] = (is3a, np.zeros(n, dtype=bool))
        return self._cache['mode']

    # ---- 定標 ------------------------------------------------------------------------------
    def _vis_coefficients(self, ch):
        """NOAA-15 以後：NOAA STAR 的在軌係數 (沒有時發射前雙增益係數，見 AvhrrL1bFile._vis_fallback)；
        更早的衛星：POD Guide 的發射前係數。"""
        key = ('vis', ch)
        if key not in self._cache:
            if self.generation == 3:
                out = self._vis_fallback(ch, tr('原始框', 'raw frames'))
            else:
                p = fmt.VIS_PRELAUNCH_POD.get(self.satellite, {}).get(ch)
                out = tuple(np.full(self.nlines, v, dtype=np.float64) for v in p) if p else None
                if p:
                    print(tr('%s CH%d：發射前的可見光係數 (沒有在軌衰減修正，反射率偏低)',
                             '%s CH%d: pre-launch visible coefficients (no in-orbit degradation correction; '
                             'reflectances are low)') % (self.label, ch))
            self._cache[key] = out
        return self._cache[key]

    def _ir_radiance(self, ch, counts, r0, r1):
        if ('rawmsg', ch) not in self._cache:
            self._cache[('rawmsg', ch)] = True
            print(tr('%s CH%d：由 PRT、黑體與太空計數定標 (KLM Guide §7.1.2.4)',
                     '%s CH%d: calibrated from the PRT, blackbody and space counts (KLM Guide 7.1.2.4)')
                  % (self.label, ch))
        return self._ir_raw(ch, counts.astype(np.float64), r0, r1)


# ===== 讀取器 =====
class AvhrrL0Reader(l1b.AvhrrL1bReader):
    NAME = 'AVHRR-L0'
    SENSOR = 'AVHRR'
    DESCRIPTION = ('TIROS-N~NOAA-19 AVHRR HRPT 原始框 (Level 0；.raw / .raw.xz、千葉 .l0.bz2 / .gi.bz2，可直接讀 .tar；'
                   '依內容辨識，定位需要 TLE)')
    DESCRIPTION_EN = ('TIROS-N~NOAA-19 AVHRR HRPT raw frames (Level 0; .raw / .raw.xz, Chiba .l0.bz2 / .gi.bz2, also inside '
                      '.tar; detected by content; navigation needs TLEs)')
    SOURCE = 'Station'
    PATTERNS = ('*',)

    @classmethod
    def find(cls, folder=None):
        """資料夾 (往下兩層) 中依內容判斷為 HRPT 原始框的檔；同一個檔同時在 .tar 與資料夾中時只留資料夾中的。"""
        root = folder or '.'
        seen = {}
        for name in hrpt.candidates(root):
            s = sniff(os.path.join(root, name) if folder else name)
            if s is None:
                continue
            key = (s.satellite, s.time.replace(microsecond=0), s.nframes)
            old = seen.get(key)
            if old is None or (hrpt.TAR_SEP in old and hrpt.TAR_SEP not in name):
                seen[key] = name
        return sorted(seen.values(), key=cls.sort_key)

    def make_file(self, name):
        return AvhrrL0File(name, self.cfg)

    def source_label(self, segment):
        return self.SOURCE

    def select_band(self):
        tles = tle_table(self.cfg)
        for f in self.files:
            if hrpt.pick_tle(tles, f.satellite, f.time) is None:
                print(tr('注意：%s 沒有 %s 的 TLE，將無法定位 (見 [tle] dir)', 'NOTE: %s has no TLE for %s and cannot be '
                         'navigated (see [tle] dir)') % (f.label, f.satellite))
        super().select_band()

    def prepare(self):
        try:
            return super().prepare()
        finally:
            for f in self.files:
                f.release()


READER = AvhrrL0Reader
