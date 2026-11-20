"""MSG SEVIRI (Meteosat-8~11) — EUMETSAT HRIT
檔名：H-000-MSG4__-MSG4________-IR_108___-000001___-201802281500-__ (結尾 C_ = 小波壓縮；可再加 .bz2)
  服務：MSGn________ = 0 度全圓盤、MSGn_IODC___ = 印度洋、MSGn_RSS____ = 快速掃描 (歐洲)
  每個頻道 8 個分段 (每段 464 行，第 1 段最南)，HRV 24 段；同一時次另有 PRO (prologue) 與 EPI (epilogue)。
  prologue：衛星、星下點經度、地球模型、定標係數 (結構同舊版 prologue.py，見 sensors/seviri_native_hdr.py)；
  epilogue：實際掃描時間、HRV 上下視窗的範圍。沒有 prologue 的時次無法定標，會列出原因後略過。
檔頭紀錄 (CGMS HRIT + EUMETSAT MSG)：0 主檔頭、1 影像結構、2 導航 (CFAC/LFAC/COFF/LOFF，CFAC、LFAC 為負)、
  3 資料函數、4 註記、5 時間、128 分段識別 (頻道、分段序號、計畫起訖分段)、129 逐行品質 (行號、時間、有效 / 品質旗標)。
影像：10 位元，第 0 列最南、第 0 行最東 (導航與定標見 seviri_common.py)。
壓縮檔 (-C_)：EUMETSAT 小波壓縮，需安裝 pyPublicDecompWT (pip install pyPublicDecompWT)，
  或在 Satima-GEO.toml 的 [seviri] xrit_decompress 指定 EUMETSAT xRITDecompress 執行檔 (可寫 "wsl /路徑/xRITDecompress")；
  解壓在記憶體 / 暫存目錄中進行，不會寫入資料目錄。
完整檔頭：reader.header(scene, 'IR_108') 回傳 {'prologue', 'epilogue', 'segments': {分段: 檔頭紀錄}}。
"""
import bz2
import collections
import datetime
import os
import re
import shlex
import shutil
import subprocess
import tempfile

import numpy as np

import seviri_common as sc
import seviri_native_hdr as hdr
from satima.core.i18n import tr
from satima.geo.reader import GeoReader, Scene

_NAME = re.compile(r'^([HL])-000-(MSG\d)_*-(MSG\d)_*([A-Z]*)_*-([A-Z0-9_]{9})-([A-Z0-9_]{9})-(\d{12})-(__|C_)(\.bz2)?$')
_SERVICE = {'': 'FD', 'IODC': 'IODC', 'RSS': 'RSS'}
_PROLOGUE = hdr.HritPrologue().get()
_EPILOGUE = np.dtype(hdr.Msg15NativeTrailerRecord().seviri_l15_trailer).newbyteorder('>')
_CDS = [('Days', '>u2'), ('Milliseconds', '>u4')]
_RECORDS = {
    0: ('Primary_Header', np.dtype([('Header_Type', 'u1'), ('Header_Record_Length', '>u2'), ('File_Type_Code', 'u1'),
                                    ('Total_Header_Length', '>u4'), ('Data_Field_Length', '>u8')])),
    1: ('Image_Structure', np.dtype([('Header_Type', 'u1'), ('Header_Record_Length', '>u2'), ('NB', 'u1'),
                                     ('NC', '>u2'), ('NL', '>u2'), ('Compression_Flag', 'u1')])),
    2: ('Image_Navigation', np.dtype([('Header_Type', 'u1'), ('Header_Record_Length', '>u2'),
                                      ('Projection_Name', 'S32'), ('CFAC', '>i4'), ('LFAC', '>i4'), ('COFF', '>i4'),
                                      ('LOFF', '>i4')])),
    5: ('Time_Stamp', np.dtype([('Header_Type', 'u1'), ('Header_Record_Length', '>u2'), ('CDS_P_Field', 'u1')] + _CDS)),
    7: ('Key_Header', np.dtype([('Header_Type', 'u1'), ('Header_Record_Length', '>u2'), ('Key_Number', 'u1'),
                                ('Seed', '>f8')])),
    128: ('Segment_Identification', np.dtype([('Header_Type', 'u1'), ('Header_Record_Length', '>u2'),
                                              ('GP_SC_ID', '>i2'), ('Spectral_Channel_ID', 'i1'),
                                              ('Segment_Sequence_Number', '>u2'),
                                              ('Planned_Start_Segment_Number', '>u2'),
                                              ('Planned_End_Segment_Number', '>u2'),
                                              ('Data_Field_Representation', 'i1')])),
}
_TEXT = {3: 'Image_Data_Function', 4: 'Annotation', 6: 'Ancillary_Text'}
_LINE_QUALITY = np.dtype([('Line_Number_in_Grid', '>i4'), ('Line_Mean_Acquisition', _CDS), ('Line_Validity', 'u1'),
                          ('Line_Radiometric_Quality', 'u1'), ('Line_Geometric_Quality', 'u1')])


def read_header(raw):
    """HRIT 檔頭 -> {紀錄名稱: numpy 結構 / 文字}；129 (逐行品質) 為結構陣列。"""
    out, pos = {}, 0
    total = None
    while pos + 3 <= len(raw) and (total is None or pos < total):
        rtype = raw[pos]
        length = int(np.frombuffer(raw[pos + 1:pos + 3], '>u2')[0])
        if length < 3:
            break
        rec = raw[pos:pos + length]
        if rtype in _RECORDS:
            name, dt = _RECORDS[rtype]
            out[name] = np.frombuffer(rec[:dt.itemsize].ljust(dt.itemsize, b'\0'), dt, count=1)[0]
            if rtype == 0:
                total = int(out[name]['Total_Header_Length'])
        elif rtype in _TEXT:
            out[_TEXT[rtype]] = rec[3:].decode('ascii', 'replace')
        elif rtype == 129:
            n = (length - 3) // _LINE_QUALITY.itemsize
            out['Line_Quality'] = np.frombuffer(rec[3:3 + n * _LINE_QUALITY.itemsize], _LINE_QUALITY)
        pos += length
    return out


def _wsl_path(p):
    p = os.path.abspath(p)
    return '/mnt/%s%s' % (p[0].lower(), p[2:].replace('\\', '/')) if p[1:2] == ':' else p.replace('\\', '/')


class _Files:
    """讀取 (含 .bz2 與小波壓縮) 的整個檔案內容，最近使用的保留在記憶體。"""

    def __init__(self, cfg, maxsize=40):
        self.cfg, self.maxsize = cfg, maxsize
        self._items = collections.OrderedDict()

    def get(self, path):
        if path in self._items:
            self._items.move_to_end(path)
            return self._items[path]
        with open(path, 'rb') as f:
            raw = f.read()
        if path.lower().endswith('.bz2'):
            raw = bz2.decompress(raw)
        h = read_header(raw)
        if int(h['Image_Structure']['Compression_Flag']) == 1 if 'Image_Structure' in h else False:
            raw = self._decompress(path, raw)
        self._items[path] = raw
        while len(self._items) > self.maxsize:
            self._items.popitem(last=False)
        return raw

    def _decompress(self, path, raw):
        try:
            from pyPublicDecompWT import xRITDecompress
            d = xRITDecompress()
            d.decompress(raw)
            return bytes(d.data())
        except ImportError:
            pass
        cmd = str(self.cfg.seviri.get('xrit_decompress', '') or '').strip()
        if not cmd:
            raise ValueError(tr('%s 為小波壓縮的 HRIT：請安裝 pyPublicDecompWT (pip install pyPublicDecompWT)，'
                                '或在 Satima-GEO.toml 的 [seviri] xrit_decompress 指定 xRITDecompress 執行檔',
                                '%s is a wavelet-compressed HRIT file: install pyPublicDecompWT '
                                '(pip install pyPublicDecompWT), or set [seviri] xrit_decompress in Satima-GEO.toml '
                                'to the xRITDecompress executable') % os.path.basename(path))
        tmp = tempfile.mkdtemp(prefix='satima_xrit_')
        try:
            src = os.path.join(tmp, os.path.basename(path).replace('.bz2', ''))
            with open(src, 'wb') as f:
                f.write(raw)
            args = shlex.split(cmd, posix=False)
            wsl = args[0].lower() == 'wsl'
            subprocess.run(args + [_wsl_path(src) if wsl else src], cwd=tmp, check=True, capture_output=True,
                           timeout=300)
            outs = [n for n in os.listdir(tmp) if n != os.path.basename(src)]
            if not outs:
                raise ValueError(tr('xRITDecompress 沒有產生輸出 (%s)', 'xRITDecompress produced no output (%s)')
                                 % os.path.basename(path))
            with open(os.path.join(tmp, outs[0]), 'rb') as f:
                return f.read()
        finally:
            shutil.rmtree(tmp, ignore_errors=True)


class SeviriHritReader(GeoReader):
    NAME = 'SEVIRI-HRIT'
    SENSOR = 'SEVIRI'
    DESCRIPTION = 'MSG (Meteosat-8~11) SEVIRI EUMETSAT HRIT (H-000-MSG*，含 -C_ 壓縮與 .bz2)'
    DESCRIPTION_EN = 'MSG (Meteosat-8~11) SEVIRI EUMETSAT HRIT (H-000-MSG*, including -C_ compressed and .bz2)'
    PATTERNS = ('H-000-MSG*', 'L-000-MSG*')
    CHANNELS = {n: sc.make_channel(i + 1) for i, n in enumerate(sc.NAMES)}
    TRUECOLOR = {}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._files = _Files(cfg)
        self._info = {}

    # ---- 檔案 --------------------------------------------------------------
    def scenes(self, folder='.'):
        groups = {}
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            platform, service = m.group(3), _SERVICE.get(m.group(4), m.group(4) or 'FD')
            chan, seg = m.group(5).rstrip('_'), m.group(6).rstrip('_')
            t = datetime.datetime.strptime(m.group(7), '%Y%m%d%H%M')
            key = (platform, service, t)
            g = groups.setdefault(key, {'files': {}, 'pro': None, 'epi': None, 'compressed': False})
            path = os.path.join(folder, name)
            if seg == 'PRO':
                g['pro'] = path
            elif seg == 'EPI':
                g['epi'] = path
            elif chan in sc.NUMBER and seg.isdigit():
                g['files'].setdefault(chan, {})[int(seg)] = path
                g['compressed'] |= m.group(8) == 'C_'
        out = []
        for (platform, service, t), g in sorted(groups.items(), key=lambda kv: (kv[0][2], kv[0][0])):
            if not g['files']:
                continue
            if g['pro'] is None:
                print(tr('  [略過] %s %s %s：缺 prologue (PRO) 檔，無法定標',
                         '  [skip] %s %s %s: the prologue (PRO) file is missing; cannot calibrate')
                      % (platform, service, t.strftime('%Y%m%d %H%M')))
                continue
            satid = sc.SATID.get(platform, 324)
            s = Scene(sc.SATELLITES[satid], t, 'RSS' if service == 'RSS' else 'FD',
                      files=dict(g['files']), extra={'pro': g['pro'], 'epi': g['epi'], 'satid': satid},
                      region='IODC' if service == 'IODC' else '', format='HRIT (C)' if g['compressed'] else 'HRIT')
            out.append(s)
        return out

    def channels(self, scene):
        return [self.CHANNELS[n] for n in sc.NAMES if n in scene.files]

    def _meta(self, scene):
        """prologue / epilogue -> 定標、導航所需的資訊。"""
        key = scene.extra['pro']
        if key in self._info:
            return self._info[key]
        raw = self._files.get(scene.extra['pro'])
        h = read_header(raw)
        thl = int(h['Primary_Header']['Total_Header_Length'])
        pro = np.frombuffer(raw[thl:thl + _PROLOGUE.itemsize], _PROLOGUE, count=1)[0]
        epi = None
        if scene.extra.get('epi'):
            eraw = self._files.get(scene.extra['epi'])
            eh = read_header(eraw)
            ethl = int(eh['Primary_Header']['Total_Header_Length'])
            if len(eraw) >= ethl + _EPILOGUE.itemsize:
                epi = np.frombuffer(eraw[ethl:ethl + _EPILOGUE.itemsize], _EPILOGUE, count=1)[0]
        em = pro['GeometricProcessing']['EarthModel']
        start = sc.cds(*pro['ImageAcquisition']['PlannedAcquisitionTime']['TrueRepeatCycleStart'][['Days', 'Milliseconds']].tolist())
        if epi is not None:
            fs = epi['ImageProductionStats']['ActualScanningSummary']['ForwardScanStart']
            start = sc.cds(fs['Days'], fs['Milliseconds']) or start
        info = dict(prologue=pro, epilogue=epi, satid=int(pro['SatelliteStatus']['SatelliteDefinition']['SatelliteId']),
                    sub_lon=float(pro['ImageDescription']['ProjectionDescription']['LongitudeOfSSP']),
                    earth_model=int(em['TypeOfEarthModel']), a=float(em['EquatorialRadius']) * 1000.,
                    b=(float(em['NorthPolarRadius']) + float(em['SouthPolarRadius'])) / 2. * 1000.,
                    start=start or scene.time)
        self._info[key] = info
        return info

    def header(self, scene, channel=None):
        info = self._meta(scene)
        segs = {}
        for c in ([channel] if channel else scene.files):
            for seg, path in sorted(scene.files[c].items()):
                segs[(c, seg)] = read_header(self._files.get(path))
        return {'prologue': info['prologue'], 'epilogue': info['epilogue'], 'segments': segs}

    def _layout(self, scene, channel):
        """頻道的分段配置：(各分段 (序號, 路徑, 檔頭), 每段行數, 行數 (寬), 計畫起始分段, 計畫分段數)。"""
        key = ('layout', scene.extra['pro'], channel)
        if key in self._info:
            return self._info[key]
        segs = []
        for seg, path in sorted(scene.files[channel].items()):
            h = read_header(self._files.get(path))
            segs.append((seg, path, h))
        h0 = segs[0][2]
        sid = h0['Segment_Identification']
        nl, nc = int(h0['Image_Structure']['NL']), int(h0['Image_Structure']['NC'])
        p0, p1 = int(sid['Planned_Start_Segment_Number']), int(sid['Planned_End_Segment_Number'])
        out = (segs, nl, nc, p0, p1 - p0 + 1)
        self._info[key] = out
        return out

    # ---- 導航 --------------------------------------------------------------
    def nav(self, scene, channel):
        key = ('nav', scene.extra['pro'], channel)
        if key in self._info:
            return self._info[key]
        info = self._meta(scene)
        segs, nl, nc, p0, nseg = self._layout(scene, channel)
        seg, _, h = segs[0]
        nv = h['Image_Navigation']
        cstep = np.deg2rad(2. ** 16 / abs(int(nv['CFAC'])))
        lstep = np.deg2rad(2. ** 16 / abs(int(nv['LFAC'])))
        shape = (nl * nseg, nc)
        # 該分段的 LOFF 對應其第一 (最南) 行：行號 L = center - LOFF + 1 + 分段內列號
        hrv = channel == 'HRV'
        center = sc.CENTER_HRV if hrv else sc.CENTER_VISIR
        first_line = center - int(nv['LOFF']) + 1 - (seg - p0) * nl
        if not hrv:
            first_col = center - int(nv['COFF']) + 1
            nav = sc.geos_nav(shape, first_line, first_col, lstep, cstep, info['sub_lon'], info['a'], info['b'],
                              info['earth_model'])
        else:
            nav = self._hrv_nav(info, shape, first_line, lstep, cstep)
        self._info[key] = nav
        return nav

    def _hrv_nav(self, info, shape, first_line, lstep, cstep):
        """HRV：依 epilogue 的實際範圍分成下 (南) 與上 (北) 視窗，各以自己的東界為第 0 行。"""
        epi = info['epilogue']
        mk = lambda east: sc.geos_nav(shape, first_line, east, lstep, cstep, info['sub_lon'], info['a'], info['b'],
                                      info['earth_model'], hrv=True)
        if epi is None:
            print(tr('  [HRV] 缺 epilogue (EPI)，以 prologue 的計畫範圍定位',
                     '  [HRV] the epilogue (EPI) is missing; geolocating with the planned coverage in the prologue'))
            cov = info['prologue']['ImageDescription']['PlannedCoverageHRV']
            lower = (int(cov['LowerSouthLinePlanned']), int(cov['LowerNorthLinePlanned']), int(cov['LowerEastColumnPlanned']))
            upper = (int(cov['UpperSouthLinePlanned']), int(cov['UpperNorthLinePlanned']), int(cov['UpperEastColumnPlanned']))
        else:
            cov = epi['ImageProductionStats']['ActualL15CoverageHRV']
            lower = (int(cov['LowerSouthLineActual']), int(cov['LowerNorthLineActual']), int(cov['LowerEastColumnActual']))
            upper = (int(cov['UpperSouthLineActual']), int(cov['UpperNorthLineActual']), int(cov['UpperEastColumnActual']))
        navs, starts = [mk(lower[2])], [lower[0] - first_line]
        if upper[1] > upper[0] > 0 and upper[0] > lower[1] - 1:
            navs.append(mk(upper[2]))
            # 實際資料從 LowerNorthLineActual 這一行起就是上視窗的排列 (Zenodo 2018/02/28 樣本：行號 8256 的資料
            # 與上方各行連續 (平均差 0.003)、與下視窗不連續 (0.054)；同 satpy)，比 UpperSouthLineActual 早一行
            starts.append(lower[1] - first_line)
        starts[0] = -1e9
        return sc.SplitNav(navs, starts, shape) if len(navs) > 1 else navs[0]

    def satellite_position(self, scene):
        return self._meta(scene)['sub_lon'], sc.H / 1000.

    # ---- 讀取與定標 --------------------------------------------------------
    def _calibration(self, scene, channel):
        key = ('cal', scene.extra['pro'], channel)
        if key not in self._info:
            info = self._meta(scene)
            pro = info['prologue']
            n = sc.NUMBER[channel]
            rp = pro['RadiometricProcessing']
            cal = rp['Level15ImageCalibration'][n - 1]
            mp = rp['MPEFCalFeedback'][n - 1]
            c = sc.Calibration(info['satid'], n, cal['CalSlope'], cal['CalOffset'],
                               rad_type=int(pro['ImageDescription']['Level15ImageProduction']['PlannedChanProcessing'][n - 1]),
                               gsics=(float(mp['GSICSCalCoeff']), float(mp['GSICSOffsetCount'])),
                               mode=self.cfg.seviri.calibration, time=info['start'])
            if c.mode != 'nominal':
                print(tr('  %s 定標：%s', '  %s calibration: %s') % (channel, c.mode))
            self._info[key] = c
        return self._info[key]

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        segs, nl, nc, p0, nseg = self._layout(scene, channel)
        cal = self._calibration(scene, channel)
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        rowbytes = nc * 10 // 8
        mask_bad = True
        for seg, path, _ in segs:
            s0 = (seg - p0) * nl                      # 此分段第 0 列在整張影像的列號
            a, b = max(r0, s0), min(r1, s0 + nl)
            if a >= b:
                continue
            raw = self._files.get(path)
            h = read_header(raw)
            thl = int(h['Primary_Header']['Total_Header_Length'])
            buf = raw[thl + (a - s0) * rowbytes:thl + (b - s0) * rowbytes]
            n = len(buf) // rowbytes
            if n == 0:
                continue
            counts = sc.unpack10(buf[:n * rowbytes]).reshape(n, nc)[:, c0:c1]
            vals = cal(counts)
            q = h.get('Line_Quality')
            if mask_bad and q is not None and len(q) == nl:
                qq = q[a - s0:a - s0 + n]
                bad = sc.bad_lines(qq['Line_Validity'], qq['Line_Radiometric_Quality'], qq['Line_Geometric_Quality'])
                vals[bad] = np.nan
            out[a - r0:a - r0 + n, :vals.shape[1]] = vals
        return out


READER = SeviriHritReader
