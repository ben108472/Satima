"""MSG SEVIRI (Meteosat-8~11) — EUMETSAT Level 1.5 native 格式 (.nat，EUMETSAT Data Store / UMARF)
檔名：MSG4-SEVI-MSG15-0100-NA-20210302124244.185000000Z-NA.nat (MSG1~4 = Meteosat-8~11)；
  Data Store 下載的同名 .zip 可直接放 (不必解開；讀影像時在記憶體中解壓)。
格式 (MSG Level 1.5 Native Format File Definition；結構同舊版 prologue.py，見 sensors/seviri_native_hdr.py)：
  檔頭 450400 bytes = 主檔頭 3674 (ASCII 名稱 : 值) + 次檔頭 1440 (選取的頻道 SelectedBandIDs、範圍、行列數)
         + GP 封包頭 38 + 15_DATA_HEADER 445248 (衛星、星下點經度、參考格點、定標、地球模型…)；
         沒有主 / 次檔頭的檔案 (開頭不是 "FormatName ... NATIVE") 檔頭為 445286 bytes，視為 12 頻道全圓盤。
  影像：每條掃描線一組 = 選取的可見光 / 紅外頻道各一筆 + HRV 3 筆 (有選 HRV 時)；
         每筆 = 行頭 65 bytes (封包頭、衛星、行號、頻道、取樣時間、有效 / 輻射 / 幾何品質) + 10 位元資料。
         可見光 / 紅外每行補到 4 的倍數行；HRV 全圓盤每行 5568 行 (11136 的一半，上下視窗)。
  檔尾 380363 bytes：實際掃描時間、實際範圍 (HRV 上下視窗)。
影像方向：第 0 列最南 (SouthLineSelectedRectangle)、第 0 行最東 (EastColumnSelectedRectangle)；
  導航以參考格點的格距 (km) 計算，2017/12 以前 (地球模型 1) 的 1.5 km 偏移見 seviri_common.py。
完整檔頭：reader.header(scene, 'IR_108') 回傳 {'header', 'trailer', 'line_headers' (該頻道每行的行頭)}。
"""
import datetime
import os
import re
import zipfile

import numpy as np

import seviri_common as sc
import seviri_native_hdr as hdr
from satima.core.i18n import tr
from satima.geo.reader import GeoReader, Scene

_NAME = re.compile(r'^(MSG\d)-SEVI-MSG15-\d{4}-\w+-(\d{14})\.\d+Z-.*\.nat$', re.I)
_HEADER = hdr.Msg15NativeHeaderRecord().get()                      # 450400 bytes
_DATA_HEADER = np.dtype([('GP_PK_HEADER', hdr.GSDTRecords.gp_pk_header), ('GP_PK_SH1', hdr.GSDTRecords.gp_pk_sh1),
                         ('15_DATA_HEADER', hdr.L15DataHeaderRecord().get())]).newbyteorder('>')   # 445286 bytes
_TRAILER = hdr.Msg15NativeTrailerRecord().get()                    # 380363 bytes
_ARCHIVE = b'FormatName                  : NATIVE'
_CDS = [('Days', '>u2'), ('Milliseconds', '>u4')]


def _line_header():
    return np.dtype([('GP_PK_HEADER', hdr.GSDTRecords.gp_pk_header), ('GP_PK_SH1', hdr.GSDTRecords.gp_pk_sh1),
                     ('Version', 'u1'), ('SatId', '>u2'), ('Time', '>u2', (5,)), ('LineNo', '>u4'), ('ChanId', 'u1'),
                     ('AcqTime', _CDS), ('LineValidity', 'u1'), ('LineRQuality', 'u1'),
                     ('LineGQuality', 'u1')]).newbyteorder('>')


_LINE = _line_header()                                             # 65 bytes


def _value(rec):
    """主 / 次檔頭的 ASCII 紀錄 -> 值字串。"""
    return bytes(rec['Value']).decode('ascii', 'replace').split(':')[0].strip()


class _Source:
    """.nat 的內容：一般檔案以 seek 讀取；EUMETSAT Data Store 的 zip (deflate) 列表時只解壓檔頭，
    讀影像時把整個 .nat 解壓到記憶體 (約 270 MB，只保留最近使用的一個)。"""
    _loaded = {}

    def __init__(self, path):
        self.path, self.member = path, None
        if path.lower().endswith('.zip'):
            with zipfile.ZipFile(path) as z:
                self.member = next((n for n in z.namelist() if n.lower().endswith('.nat')), None)
            if self.member is None:
                raise ValueError(tr('zip 中沒有 .nat 檔', 'no .nat file in the zip'))

    @property
    def name(self):
        return os.path.basename(self.member or self.path)

    def head(self, n):
        if self.member is None:
            with open(self.path, 'rb') as f:
                return f.read(n)
        with zipfile.ZipFile(self.path) as z, z.open(self.member) as f:
            return f.read(n)

    def read(self, offset, n):
        if self.member is None:
            with open(self.path, 'rb') as f:
                f.seek(offset)
                return f.read(n)
        data = _Source._loaded.get(self.path)
        if data is None:
            print(tr('  解壓 %s ...', '  Extracting %s ...') % self.name)
            _Source._loaded.clear()
            with zipfile.ZipFile(self.path) as z:
                data = z.read(self.member)
            _Source._loaded[self.path] = data
        return data[offset:offset + n]


class NatFile:
    def __init__(self, path):
        self.path = path
        self.src = _Source(path)
        raw = self.src.head(_HEADER.itemsize)
        if raw[:len(_ARCHIVE)] == _ARCHIVE:
            h = np.frombuffer(raw, _HEADER, count=1)[0]
            self.size = _HEADER.itemsize
            sec = h['15_SECONDARY_PRODUCT_HEADER']
            sel = _value(sec['SelectedBandIDs'])
            south, north = int(_value(sec['SouthLineSelectedRectangle'])), int(_value(sec['NorthLineSelectedRectangle']))
            east, west = int(_value(sec['EastColumnSelectedRectangle'])), int(_value(sec['WestColumnSelectedRectangle']))
            nl_visir, nl_hrv = int(_value(sec['NumberLinesVISIR'])), int(_value(sec['NumberLinesHRV']))
            nc_hrv = int(_value(sec['NumberColumnsHRV']))
        else:
            h = np.frombuffer(raw[:_DATA_HEADER.itemsize], _DATA_HEADER, count=1)[0]
            self.size = _DATA_HEADER.itemsize
            sel, south, north, east, west = 'X' * 12, 1, 3712, 1, 3712
            nl_visir, nl_hrv, nc_hrv = 3712, 11136, 11136
        self.header = h
        self.data_header = h['15_DATA_HEADER']
        self.channels = [i + 1 for i, ch in enumerate(sel[:12]) if ch.upper() == 'X']
        self.visir = [n for n in self.channels if n != 12]
        self.south, self.north, self.east, self.west = south, north, east, west
        self.nlines, self.ncols = north - south + 1, west - east + 1
        self.full_disk = self.nlines >= 3712 and self.ncols >= 3712
        # 快速掃描 (RSS)：北端 1392 行、全寬 (同 satpy 的判斷；不必讀檔尾)
        self.rss = (not self.full_disk) and self.ncols == 3712 and self.nlines == 1392 and north == 3712
        self.cols_visir = self.ncols + (4 - self.ncols % 4) % 4              # UMARF 補到 4 的倍數
        self.cols_hrv = nc_hrv if self.ncols < 3712 else nc_hrv // 2
        self.lines_hrv = nl_hrv
        self.rec_visir = _LINE.itemsize + self.cols_visir * 5 // 4
        self.rec_hrv = _LINE.itemsize + self.cols_hrv * 5 // 4
        self.group = len(self.visir) * self.rec_visir + (3 * self.rec_hrv if 12 in self.channels else 0)
        self.ngroups = nl_visir
        self._trailer = None

    @property
    def trailer(self):
        if self._trailer is None:
            raw = self.src.read(self.size + self.group * self.ngroups, _TRAILER.itemsize)
            self._trailer = np.frombuffer(raw, _TRAILER, count=1)[0] if len(raw) == _TRAILER.itemsize else None
        return self._trailer

    def records(self, n, r0, r1):
        """頻道 n 第 r0~r1-1 列的 (行頭, 10 位元資料 bytes) (HRV 的列為 HRV 列號)。"""
        if n == 12:
            g0, g1 = r0 // 3, (r1 + 2) // 3
            base = len(self.visir) * self.rec_visir
            size = self.rec_hrv
        else:
            g0, g1 = r0, r1
            base = self.visir.index(n) * self.rec_visir
            size = self.rec_visir
        g0, g1 = max(0, g0), min(self.ngroups, g1)
        if g1 <= g0:
            return np.zeros(0, _LINE), []
        buf = self.src.read(self.size + g0 * self.group, (g1 - g0) * self.group)
        k = len(buf) // self.group
        groups = np.frombuffer(buf[:k * self.group], np.uint8).reshape(k, self.group)
        if n == 12:
            recs = groups[:, base:base + 3 * size].reshape(k * 3, size)
            first = g0 * 3
            recs = recs[r0 - first:r1 - first]
        else:
            recs = groups[:, base:base + size]
        heads = np.frombuffer(np.ascontiguousarray(recs[:, :_LINE.itemsize]).tobytes(), _LINE)
        return heads, recs[:, _LINE.itemsize:]


class SeviriNatReader(GeoReader):
    NAME = 'SEVIRI-NAT'
    SENSOR = 'SEVIRI'
    DESCRIPTION = 'MSG (Meteosat-8~11) SEVIRI Level 1.5 native (MSG?-SEVI-MSG15-*.nat，或 Data Store 的 .zip)'
    DESCRIPTION_EN = 'MSG (Meteosat-8~11) SEVIRI Level 1.5 native (MSG?-SEVI-MSG15-*.nat, or the Data Store .zip)'
    PATTERNS = ('MSG?-SEVI-MSG15-*.nat', 'MSG?-SEVI-MSG15-*.zip')
    CHANNELS = {n: sc.make_channel(i + 1) for i, n in enumerate(sc.NAMES)}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._nats = {}
        self._cache = {}

    def _nat(self, path):
        if path not in self._nats:
            self._nats[path] = NatFile(path)
        return self._nats[path]

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            path = os.path.join(folder, name)
            try:
                nat = self._nat(path)
            except (OSError, ValueError) as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                continue
            dh = nat.data_header
            satid = int(dh['SatelliteStatus']['SatelliteDefinition']['SatelliteId'])
            ts = dh['ImageAcquisition']['PlannedAcquisitionTime']['TrueRepeatCycleStart']
            m = _NAME.match(nat.src.name)
            start = sc.cds(ts['Days'], ts['Milliseconds']) or (
                datetime.datetime.strptime(m.group(2), '%Y%m%d%H%M%S') if m else None)
            if start is None:
                print(tr('  [略過] %s：沒有觀測時間', '  [skip] %s: no observation time') % name)
                continue
            ssp = float(dh['ImageDescription']['ProjectionDescription']['LongitudeOfSSP'])
            reduced = nat.rss
            region = 'IODC' if ssp > 30. else ''
            area = 'RSS' if reduced else ('FD' if nat.full_disk else 'ROI')
            s = Scene(sc.SATELLITES.get(satid, 'Meteosat-?'), sc.nominal_time(start, reduced), area,
                      files={sc.NAMES[n - 1]: path for n in nat.channels}, region=region, format='NAT',
                      extra={'satid': satid, 'start': start})
            if area == 'ROI':
                try:
                    s.extra['georange'] = self._image_range(s)
                    s.auto_range = True
                except (SystemExit, ValueError):
                    pass
            out.append(s)
        return sorted(out, key=lambda s: (s.time, s.satellite))

    def channels(self, scene):
        return [self.CHANNELS[n] for n in sc.NAMES if n in scene.files]

    def header(self, scene, channel):
        nat = self._nat(scene.files[channel])
        n = sc.NUMBER[channel]
        rows = nat.lines_hrv if n == 12 else nat.ngroups
        heads, _ = nat.records(n, 0, rows)
        trailer = nat.trailer
        return {'header': nat.header, 'trailer': trailer, 'line_headers': heads}

    def _geo(self, nat):
        dh = nat.data_header
        em = dh['GeometricProcessing']['EarthModel']
        return (float(dh['ImageDescription']['ProjectionDescription']['LongitudeOfSSP']), float(em['EquatorialRadius']) * 1000.,
                (float(em['NorthPolarRadius']) + float(em['SouthPolarRadius'])) / 2. * 1000., int(em['TypeOfEarthModel']))

    def nav(self, scene, channel):
        path = scene.files[channel]
        key = ('nav', path, channel == 'HRV')
        if key in self._cache:
            return self._cache[key]
        nat = self._nat(path)
        dh = nat.data_header
        ssp, a, b, model = self._geo(nat)
        if channel != 'HRV':
            g = dh['ImageDescription']['ReferenceGridVIS_IR']
            if int(g['GridOrigin']) != 2:
                raise ValueError(tr('參考格點原點 GridOrigin = %d (只支援 2：東南角)',
                                    'reference grid origin GridOrigin = %d (only 2, the south-east corner, is supported)')
                                 % int(g['GridOrigin']))
            nav = sc.geos_nav((nat.nlines, nat.ncols), nat.south, nat.east, float(g['LineDirGridStep']) * 1000. / sc.H,
                              float(g['ColumnDirGridStep']) * 1000. / sc.H, ssp, a, b, model)
        else:
            g = dh['ImageDescription']['ReferenceGridHRV']
            ls, cs = float(g['LineDirGridStep']) * 1000. / sc.H, float(g['ColumnDirGridStep']) * 1000. / sc.H
            shape = (nat.lines_hrv, nat.cols_hrv)
            mk = lambda south, east: sc.geos_nav(shape, south, east, ls, cs, ssp, a, b, model, hrv=True)
            trailer = nat.trailer
            # 全圓盤與標準快速掃描 (RSS) 依檔尾的實際範圍 (下 / 上視窗)；自選範圍 (ROI) 依選取範圍換算
            if (nat.full_disk or nat.rss) and trailer is not None:
                cov = trailer['15TRAILER']['ImageProductionStats']['ActualL15CoverageHRV']
                ls_, ln_, le_ = int(cov['LowerSouthLineActual']), int(cov['LowerNorthLineActual']), int(cov['LowerEastColumnActual'])
                us_, un_, ue_ = int(cov['UpperSouthLineActual']), int(cov['UpperNorthLineActual']), int(cov['UpperEastColumnActual'])
                navs, starts = [mk(ls_, le_)], [-1e9]
                if un_ > us_ > 0:
                    navs.append(mk(ls_, ue_))
                    # 上視窗的排列自 LowerNorthLineActual 這一行開始 (見 SEVIRI-HRIT.py)
                    starts.append(ln_ - ls_)
                nav = sc.SplitNav(navs, starts, shape) if len(navs) > 1 else navs[0]
            else:
                nav = mk(3 * nat.south - 2, 3 * nat.east - 2)
        self._cache[key] = nav
        return nav

    def satellite_position(self, scene):
        return self._geo(self._nat(next(iter(scene.files.values()))))[0], sc.H / 1000.

    def _calibration(self, scene, channel):
        key = ('cal', scene.files[channel], channel)
        if key not in self._cache:
            nat = self._nat(scene.files[channel])
            dh = nat.data_header
            n = sc.NUMBER[channel]
            rp = dh['RadiometricProcessing']
            cal, mp = rp['Level15ImageCalibration'][n - 1], rp['MPEFCalFeedback'][n - 1]
            c = sc.Calibration(scene.extra['satid'], n, cal['CalSlope'], cal['CalOffset'],
                               rad_type=int(dh['ImageDescription']['Level15ImageProduction']['PlannedChanProcessing'][n - 1]),
                               gsics=(float(mp['GSICSCalCoeff']), float(mp['GSICSOffsetCount'])),
                               mode=self.cfg.seviri.calibration, time=scene.extra['start'])
            if c.mode != 'nominal':
                print(tr('  %s 定標：%s', '  %s calibration: %s') % (channel, c.mode))
            self._cache[key] = c
        return self._cache[key]

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        nat = self._nat(scene.files[channel])
        n = sc.NUMBER[channel]
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        heads, data = nat.records(n, r0, r1)
        if len(heads) == 0:
            return out
        ncol = nat.cols_hrv if n == 12 else nat.cols_visir
        counts = sc.unpack10(np.ascontiguousarray(data).tobytes()).reshape(len(heads), ncol)
        width = nat.cols_hrv if n == 12 else nat.ncols
        counts = counts[:, c0:min(c1, width)]
        vals = self._calibration(scene, channel)(counts)
        bad = sc.bad_lines(heads['LineValidity'], heads['LineRQuality'], heads['LineGQuality'])
        vals[bad] = np.nan
        out[:len(heads), :vals.shape[1]] = vals
        return out


READER = SeviriNatReader
