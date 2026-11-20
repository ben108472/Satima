"""MTSAT-1R / MTSAT-2 Imager — JMA HRIT (gzip，可放在 tar 中)

檔名：HRIT_MTSAT2_20150707_0101_DK02IR1.gz (可直接放，或整包 HRIT_MTSAT2_201507070115.tar 不必解開)
  衛星 MTSAT1 = MTSAT-1R、MTSAT2 = MTSAT-2；觀測範圍 DK01 全圓盤、DK02 北半球 (NH)、DK03 南半球 (SH)
  波段 VIS (1 km)、IR1 10.8 µm、IR2 12.0 µm、IR3 6.75 µm (水氣)、IR4 3.75 µm (4 km)
  tar 內的 HLM_*.gz 只有檔頭、沒有影像，不使用。
  MTSAT 官方沒有頻道編號，名稱沿用 JMA 的 VIS、IR1~IR4。感測器沒有專屬名稱 (JMA 稱 Imager)，
  所以本格式以衛星系列命名。
格式：HRIT (CGMS LRIT/HRIT Global Specification)：
  導航 = 標準靜止衛星投影 (CFAC/LFAC/COFF/LOFF，行列號由 1 起算)，衛星經度在投影名稱 GEOS(145.00) 中
  校正 = 檔頭 Data Function 內的「計數值:=亮溫 (K) / 反照率 (%)」對照點，點與點之間線性內插 (同 satpy)。
         與 conversion_table/MTSAT 的外部表相同，但不必再依衛星與日期 (MTSAT-2 於 2010/03/15 換表) 選表。
資料以 gzip 壓縮，只能從頭依序解壓：讀取時只解壓到需要的最後一列為止。
"""
import datetime
import gzip
import os
import re
import tarfile

import numpy as np

from satima.core.i18n import tr
from satima.geo.nav import GeosNav
from satima.geo.reader import Channel, FileCache, GeoReader, Scene

_NAME = re.compile(r'HRIT_MTSAT([12])_(\d{8})_(\d{4})_(DK\d\d)(VIS|IR[1-4])\.gz$')
_SATELLITE = {'1': 'MTSAT-1R', '2': 'MTSAT-2'}
_REGION = {'DK01': 'FD', 'DK02': 'NH', 'DK03': 'SH'}        # 全圓盤 / 北半球 / 南半球 (同舊版 Full / NH / SH)
# (名稱, 中心波長 µm, 解析度 km, 種類, 功能代稱)
_BANDS = (('VIS', 0.68, 1.0, 'vis', 'Visible'), ('IR1', 10.8, 4.0, 'ir', 'LW Window'),
          ('IR2', 12.0, 4.0, 'ir', 'Dirty LW Window'), ('IR3', 6.75, 4.0, 'wv', 'Water Vapor'),
          ('IR4', 3.75, 4.0, 'swir', 'Shortwave Window'))
# CGMS 標準地球與衛星距離 (km)
_REQ, _RPOL, _DISTANCE = 6378.169, 6356.5838, 42164.0

# ---- HRIT 檔頭紀錄 (CGMS LRIT/HRIT Global Spec 與 JMA HRIT Mission Specific Implementation) ----
# 每筆紀錄開頭都是 Header_Type (1 byte) + Header_Record_Length (2 bytes，含開頭)，大端序。
# 依類型代碼解析 (不假設紀錄的順序與有無)；結構與欄位名稱同舊版 MTSAT_HRIT.py。
# Header Type 0 - Primary Header
_Primary_Header = np.dtype([('Header_Type', '>u1'), ('Header_Record_Length', '>u2'), ('File_Type_Code', '>u1'),
                            ('Total_Header_Length', '>u4'), ('Data_Field_Length', '>u8')])
# Header Type 1 - Image Structure (NB 每像素位元數、NC 行數、NL 列數)
_Image_Structure = np.dtype([('Header_Type', '>u1'), ('Header_Record_Length', '>u2'), ('NB', '>u1'),
                             ('NC', '>u2'), ('NL', '>u2'), ('Compression_Flag', '>u1')])
# Header Type 2 - Image Navigation (投影名稱含衛星經度，例如 GEOS(145.00))
_Image_Navigation = np.dtype([('Header_Type', '>u1'), ('Header_Record_Length', '>u2'), ('Projection_Name', 'S32'),
                              ('CFAC', '>i4'), ('LFAC', '>i4'), ('COFF', '>i4'), ('LOFF', '>i4')])
# Header Type 5 - Time Stamp (CCSDS 日數 (1958/01/01 起) + 毫秒)
_Time_Stamp = np.dtype([('Header_Type', '>u1'), ('Header_Record_Length', '>u2'), ('CDS_P_Field', '>u1'),
                        ('CDS_T_Field_1', '>u2'), ('CDS_T_Field_2', '>u4')])
# Header Type 128 - Image Segment Identification
_Image_Segment_Identification = np.dtype([('Header_Type', '>u1'), ('Header_Record_Length', '>u2'),
                                          ('Image_Segm_Seq_No', '>u1'), ('Total_No_Image_Segm', '>u1'),
                                          ('Line_No_Image_Segm', '>u2')])
_RECORDS = {0: ('Primary_Header', _Primary_Header), 1: ('Image_Structure', _Image_Structure),
            2: ('Image_Navigation', _Image_Navigation), 5: ('Time_Stamp', _Time_Stamp),
            128: ('Image_Segment_Identification', _Image_Segment_Identification)}
# 文字型紀錄 (開頭 3 bytes 之後為 ASCII 文字)
_TEXT_RECORDS = {3: 'Image_Data_Function',               # 校正表「計數值:=物理量」
                 4: 'Annotation',                        # 檔名
                 130: 'Image_Compensation_Information',  # 逐列的 COFF / LOFF 補償
                 131: 'Image_Observation_Time',          # 逐列的觀測時間 (MJD)
                 132: 'Image_Quality_Information'}


def read_header(raw):
    """HRIT 檔頭 bytes -> {紀錄名稱: numpy 結構 (固定格式紀錄) 或 文字 (文字型紀錄)}。"""
    hrit, pos = {}, 0
    while pos + 3 <= len(raw):
        rtype = raw[pos]
        length = int(np.frombuffer(raw[pos + 1:pos + 3], '>u2')[0])
        if length < 3:
            break
        rec = raw[pos:pos + length]
        if rtype in _RECORDS:
            name, dtype = _RECORDS[rtype]
            hrit[name] = np.frombuffer(rec[:dtype.itemsize], dtype, count=1)[0]
        elif rtype in _TEXT_RECORDS:
            hrit[_TEXT_RECORDS[rtype]] = rec[3:].decode('ascii', 'replace')
        pos += length
    return hrit


def cds_time(stamp):
    """Time Stamp 紀錄 -> datetime (同舊版 read_cds_time)。"""
    return datetime.datetime(1958, 1, 1) + datetime.timedelta(days=int(stamp['CDS_T_Field_1']),
                                                              milliseconds=int(stamp['CDS_T_Field_2']))


def _parse_header(raw):
    """讀取需要的資訊；完整的檔頭紀錄保留在 'hrit' (時間戳、分段、補償、逐列觀測時間等)。"""
    hrit = read_header(raw)
    s, n = hrit['Image_Structure'], hrit['Image_Navigation']
    m = re.search(r'\(([-\d.]+)\)', n['Projection_Name'].decode('ascii', 'replace'))
    return dict(total=len(raw), nc=int(s['NC']), nl=int(s['NL']), cfac=int(n['CFAC']), lfac=int(n['LFAC']),
                coff=int(n['COFF']), loff=int(n['LOFF']), sub_lon=float(m.group(1)) if m else None,
                lut=_calibration(hrit['Image_Data_Function']), hrit=hrit)


def _calibration(text):
    """Data Function 的對照點 -> 65536 項的查表 (°C 或 0~1 反照率)；空間 / 錯誤值為 NaN。"""
    items = [t.strip() for t in re.split(r'[\r\n]+', text) if t.strip()]
    unit = next((t.split(':=')[1] for t in items if t.startswith('_UNIT')), '')
    pts = np.array([[float(v) for v in t.split(':=')] for t in items if re.match(r'^\d+:=', t)])
    pts = pts[pts[:, 0] < 65535]
    last = int(pts[:, 0].max())
    lut = np.full(65536, np.nan, dtype=np.float32)
    lut[:last + 1] = np.interp(np.arange(last + 1), pts[:, 0], pts[:, 1])
    if unit.upper().startswith('KELVIN'):
        lut -= 273.15
        lut[last] = np.nan          # 最大計數值 (1023) 為空間 / 飽和 (約 130 K)，不是實際亮溫
    elif 'ALBEDO' in unit.upper():
        lut /= 100.                 # 反照率 (%) -> 0~1
    return lut


class MtsatHritReader(GeoReader):
    NAME = 'MTSAT-HRIT'
    SENSOR = 'IMAGER'
    DESCRIPTION = 'MTSAT-1R / MTSAT-2 Imager (JMA HRIT .gz / .tar)'
    PATTERNS = ('HRIT_MTSAT*_DK*.gz', 'HRIT_MTSAT*.tar')
    CHANNELS = {b[0]: Channel(*b[:4], alias=b[4]) for b in _BANDS}
    TRUECOLOR = {}                  # 只有一個寬頻可見光，無法合成真實色

    def __init__(self, cfg):
        super().__init__(cfg)
        self._tars = FileCache(tarfile.open, maxsize=8)
        self._headers = {}

    def scenes(self, folder='.'):
        groups = {}

        def add(base, source):
            m = _NAME.search(base)
            if not m:
                return
            sat, date, hhmm, area, band = m.groups()
            key = (sat, date, hhmm, area)
            if key not in groups:
                t = datetime.datetime.strptime(date + hhmm, '%Y%m%d%H%M')
                groups[key] = Scene(_SATELLITE[sat], t, area, region=_REGION.get(area, ''))
            # 同一個檔案同時有單獨的 .gz 與 tar 內的成員時，用單獨的 .gz (不必經過 tar)
            if band not in groups[key].files or source[1] is None:
                groups[key].files[band] = source

        for name in self.find(folder):
            path = os.path.join(folder, name)
            if name.endswith('.tar'):
                try:
                    members = self._tars.get(path).getnames()
                except (tarfile.TarError, OSError) as err:
                    print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                    continue
                for member in members:
                    add(os.path.basename(member), (path, member))
            else:
                add(name, (path, None))
        return sorted(groups.values(), key=lambda s: (s.time, s.area))

    # ---- 檔案 --------------------------------------------------------------
    def _open(self, source):
        """回傳解壓縮中的串流 (gzip 只能依序讀)。"""
        path, member = source
        if member is None:
            return gzip.open(path, 'rb')
        return gzip.GzipFile(fileobj=self._tars.get(path).extractfile(member))

    def _header(self, source):
        if source not in self._headers:
            with self._open(source) as g:
                head = g.read(_Primary_Header.itemsize)
                total = int(np.frombuffer(head, _Primary_Header, count=1)[0]['Total_Header_Length'])
                self._headers[source] = _parse_header(head + g.read(total - len(head)))
        return self._headers[source]

    # ---- 讀取 --------------------------------------------------------------
    def nav(self, scene, channel):
        h = self._header(scene.files[channel])
        return GeosNav.from_cgms((h['nl'], h['nc']), h['cfac'], h['lfac'], h['coff'], h['loff'], h['sub_lon'],
                                 (_DISTANCE - _REQ) * 1000., _REQ * 1000., _RPOL * 1000.)

    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        source = scene.files[channel]
        h = self._header(source)
        nbytes = h['nc'] * 2
        with self._open(source) as g:
            g.seek(h['total'] + r0 * nbytes)          # 只解壓到需要的列 (前面的列解壓後丟棄)
            buf = g.read((r1 - r0) * nbytes)
        counts = np.frombuffer(buf, '>u2').reshape(-1, h['nc'])[:, c0:c1]
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        out[:counts.shape[0]] = h['lut'][counts]
        if counts.shape[0] < r1 - r0:
            print(tr('  注意：%s 檔案不完整，缺少部分影像列',
                     '  NOTE: the %s file is incomplete; some image lines are missing') % channel)
        return out

    def satellite_position(self, scene):
        h = self._header(next(iter(scene.files.values())))
        return h['sub_lon'], _DISTANCE - _REQ


READER = MtsatHritReader
