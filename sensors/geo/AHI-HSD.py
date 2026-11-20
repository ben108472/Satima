"""Himawari-8/9 AHI — Himawari Standard Data (HSD)

檔名：HS_H09_20261108_0400_B13_FLDK_R20_S0310.DAT[.bz2]
  觀測範圍 FLDK (全圓盤，分 10 段)、JP01~04 (日本域)、R301~304 (目標區)；
  日本域與目標區每 10 分鐘觀測 4 次，時間為檔名時間 + 2.5 分 x (編號 - 1)。
格式：JMA Himawari Standard Data User's Guide v1.3；導航為 CGMS 標準靜止衛星投影 (行列號由 1 起算)。
"""
import bz2
import datetime
import os
import re

import numpy as np

from satima.core import parallel
from satima.core.i18n import tr
from satima.geo.nav import GeosNav
from satima.geo.reader import Channel, GeoReader, Scene

# ---- 檔頭區塊 (HSD User's Guide 表 4) ------------------------------------------
_BLOCK_01 = np.dtype([('HeaderBlockNumber', 'u1'), ('BlockLength', '<u2'), ('TotalNumberOfHeaderBlocks', '<u2'),
                      ('ByteOrder', 'u1'), ('SatelliteName', 'S16'), ('ProcessingCenterName', 'S16'),
                      ('ObservationArea', 'S4'), ('OtherObservationInformation', 'S2'), ('ObservationTimeline', '<u2'),
                      ('ObservationStartTime', '<f8'), ('ObservationEndTime', '<f8'), ('FileCreationTime', '<f8'),
                      ('TotalHeaderLength', '<u4'), ('TotalDataLength', '<u4'), ('QualityFlag1', 'u1'),
                      ('QualityFlag2', 'u1'), ('QualityFlag3', 'u1'), ('QualityFlag4', 'u1'),
                      ('FileFormatVersion', 'S32'), ('FileName', 'S128'), ('Spare', 'S40')])
_BLOCK_02 = np.dtype([('HeaderBlockNumber', 'u1'), ('BlockLength', '<u2'), ('NumberOfBitsPerPixel', '<u2'),
                      ('NumberOfColumns', '<u2'), ('NumberOfLines', '<u2'), ('CompressionFlag', 'u1'), ('Spare', 'S40')])
_BLOCK_03 = np.dtype([('HeaderBlockNumber', 'u1'), ('BlockLength', '<u2'), ('SubLon', '<f8'), ('CFAC', '<u4'),
                      ('LFAC', '<u4'), ('COFF', '<f4'), ('LOFF', '<f4'), ('Distance', '<f8'),
                      ('EarthEquatorialRadius', '<f8'), ('EarthPolarRadius', '<f8'), ('EarthConst1', '<f8'),
                      ('EarthConst2', '<f8'), ('EarthConst3', '<f8'), ('EarthConstStd', '<f8'),
                      ('ResamplingTypes', '<u2'), ('ResamplingSize', '<u2'), ('Spare', 'S40')])
_BLOCK_05 = np.dtype([('HeaderBlockNumber', 'u1'), ('BlockLength', '<u2'), ('BandNumber', '<u2'),
                      ('CentralWaveLength', '<f8'), ('ValidNumberOfBitsPerPixel', '<u2'),
                      ('CountValueOfErrorPixels', '<u2'), ('CountValueOfPixelsOutsideScanArea', '<u2'),
                      ('Gain', '<f8'), ('Constant', '<f8')])
_IR_BAND = np.dtype([('c0', '<f8'), ('c1', '<f8'), ('c2', '<f8'), ('C0', '<f8'), ('C1', '<f8'), ('C2', '<f8'),
                     ('c', '<f8'), ('h', '<f8'), ('k', '<f8'), ('Spare', 'S40')])
_VIS_BAND = np.dtype([('c*', '<f8'), ('Spare', 'S104')])
_BLOCK_07 = np.dtype([('HeaderBlockNumber', 'u1'), ('BlockLength', '<u2'), ('TotalNumberOfSegments', 'u1'),
                      ('SegmentSequenceNumber', 'u1'), ('FirstLineNumber', '<u2'), ('Spare', 'S40')])

_NAME = re.compile(r'HS_H(\d\d)_(\d{8})_(\d{4})_B(\d\d)_(FLDK|JP\d\d|R3\d\d)_R(\d\d)_S(\d\d)(\d\d)\.DAT(\.bz2)?$')


def _open(path):
    return bz2.open(path, 'rb') if path.endswith('.bz2') else open(path, 'rb')


def _parse_header(raw):
    """依區塊長度切出各區塊，回傳需要的欄位。"""
    blocks, pos = {}, 0
    while pos + 3 <= len(raw):
        number, length = raw[pos], int.from_bytes(raw[pos + 1:pos + 3], 'little')
        if length == 0:
            break
        blocks[number] = raw[pos:pos + length]
        pos += length

    def rec(num, dtype, offset=0):
        return np.frombuffer(blocks[num][offset:offset + dtype.itemsize], dtype, count=1)[0]

    hdr = dict(b1=rec(1, _BLOCK_01), b2=rec(2, _BLOCK_02), b3=rec(3, _BLOCK_03), b5=rec(5, _BLOCK_05), b7=rec(7, _BLOCK_07))
    hdr['cal'] = rec(5, _VIS_BAND if hdr['b5']['BandNumber'] <= 6 else _IR_BAND, _BLOCK_05.itemsize)
    return hdr


def read_hsd(path, header_only=False):
    """讀取一個 HSD 分段檔：回傳 (檔頭, 計數值陣列 或 None)。"""
    with _open(path) as f:
        head = f.read(_BLOCK_01.itemsize)
        total = int(np.frombuffer(head, _BLOCK_01, count=1)[0]['TotalHeaderLength'])
        hdr = _parse_header(head + f.read(total - len(head)))
        if header_only:
            return hdr, None
        lines, cols = int(hdr['b2']['NumberOfLines']), int(hdr['b2']['NumberOfColumns'])
        counts = np.frombuffer(f.read(lines * cols * 2), '<u2').reshape(lines, cols)
    return hdr, counts


def calibrate(hdr, counts):
    """計數值 -> 反射率 (B01~B06) / 亮溫 °C (B07~B16)；錯誤與掃描範圍外的像素為 NaN。"""
    b5, cal = hdr['b5'], hdr['cal']
    invalid = (counts == b5['CountValueOfErrorPixels']) | (counts == b5['CountValueOfPixelsOutsideScanArea'])
    rad = b5['Gain'] * counts.astype(np.float64) + b5['Constant']
    if b5['BandNumber'] <= 6:
        out = cal['c*'] * rad
    else:
        lam = b5['CentralWaveLength'] * 1e-6
        c, h, k = cal['c'], cal['h'], cal['k']
        with np.errstate(invalid='ignore', divide='ignore'):
            ebt = (h * c / (k * lam)) / np.log1p(2 * h * c ** 2 * lam ** -5 / (rad * 1e6))
        out = cal['c0'] + cal['c1'] * ebt + cal['c2'] * ebt ** 2 - 273.15
    out = out.astype(np.float32)
    out[invalid] = np.nan
    return out


_KINDS = {1: 'vis', 2: 'vis', 3: 'vis', 4: 'nir', 5: 'nir', 6: 'nir', 7: 'swir', 8: 'wv', 9: 'wv', 10: 'wv'}
_WAVE = (0.47, 0.51, 0.64, 0.86, 1.6, 2.3, 3.9, 6.2, 6.9, 7.3, 8.6, 9.6, 10.4, 11.2, 12.4, 13.3)   # JMA 標稱值
_RES = {1: 1.0, 2: 1.0, 3: 0.5, 4: 1.0}
_ALIAS = ('Blue', 'Green', 'Red', 'Veggie', 'Snow/Ice', 'Cloud Particle Size', 'Shortwave Window',
          'Upper-level WV', 'Mid-level WV', 'Lower-level WV', 'Cloud-top Phase', 'Ozone',
          'Clean LW Window', 'LW Window', 'Dirty LW Window', 'CO2 LW')


class AhiHsdReader(GeoReader):
    NAME = 'AHI-HSD'
    SENSOR = 'AHI'
    DESCRIPTION = 'Himawari-8/9 AHI (Himawari Standard Data)'
    PATTERNS = ('HS_H0*_B*_*_R*_S*.DAT', 'HS_H0*_B*_*_R*_S*.DAT.bz2')
    CHANNELS = {'B%02d' % b: Channel('B%02d' % b, _WAVE[b - 1], _RES.get(b, 2.0), _KINDS.get(b, 'ir'),
                                     'BAND-%02d' % b, _ALIAS[b - 1])
                for b in range(1, 17)}
    TRUECOLOR = {'TrueColor': dict(red='B03', green='B02', blue='B01', nir='B04', ir='B13', swir='B07',
                                   lo_blur=(0.4, 0.4))}

    def __init__(self, cfg):
        super().__init__(cfg)
        self._headers = {}

    def scenes(self, folder='.'):
        groups = {}
        for name in self.find(folder):
            m = _NAME.match(name)
            if not m:
                continue
            sat, date, hhmm, band, area, _, seg, _, _ = m.groups()
            key = (sat, date, hhmm, area)
            if key not in groups:
                t = datetime.datetime.strptime(date + hhmm, '%Y%m%d%H%M')
                if area != 'FLDK':
                    t += datetime.timedelta(minutes=2.5 * (int(area[2:]) - 1))
                region = {'R3': 'TARGET', 'JP': 'JAPAN'}.get(area[:2], '')          # 同舊版的名稱
                groups[key] = Scene('Himawari-%d' % int(sat), t, area, auto_range=(area != 'FLDK'), region=region)
            segs = groups[key].files.setdefault('B' + band, {})
            seg = int(seg)
            if seg not in segs or name.endswith('.DAT'):      # 同一段同時有 .DAT 與 .bz2 時用 .DAT
                segs[seg] = os.path.join(folder, name)
        scenes = []
        for s in groups.values():
            s.files = {b: [segs[k] for k in sorted(segs)] for b, segs in s.files.items()}
            scenes.append(s)
        return sorted(scenes, key=lambda s: (s.time, s.area))

    def _header(self, path):
        if path not in self._headers:
            self._headers[path] = read_hsd(path, header_only=True)[0]
        return self._headers[path]

    def nav(self, scene, channel):
        hdr = self._header(scene.files[channel][0])
        b2, b3, b7 = hdr['b2'], hdr['b3'], hdr['b7']
        shape = (int(b2['NumberOfLines']) * int(b7['TotalNumberOfSegments']), int(b2['NumberOfColumns']))
        req, rpol = float(b3['EarthEquatorialRadius']), float(b3['EarthPolarRadius'])
        return GeosNav.from_cgms(shape, b3['CFAC'], b3['LFAC'], b3['COFF'], b3['LOFF'], b3['SubLon'],
                                 (float(b3['Distance']) - req) * 1000., req * 1000., rpol * 1000.)

    def read(self, scene, channel, window):
        """各分段平行解壓與定標 ([memory] workers；bz2 解壓與 numpy 運算時會釋放 GIL)。"""
        r0, r1, c0, c1 = window
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        covered = np.zeros(r1 - r0, dtype=bool)
        jobs = []
        for path in scene.files[channel]:
            hdr = self._header(path)
            first = int(hdr['b7']['FirstLineNumber']) - 1          # 該段第一列 (陣列索引)
            last = first + int(hdr['b2']['NumberOfLines'])
            a, b = max(r0, first), min(r1, last)
            if a >= b:
                continue

            def job(path=path, hdr=hdr, first=first, a=a, b=b):
                _, counts = read_hsd(path)
                out[a - r0:b - r0] = calibrate(hdr, counts[a - first:b - first, c0:c1])

            jobs.append(job)
            covered[a - r0:b - r0] = True
        parallel.run(jobs, parallel.workers(self.cfg))
        if not covered.all():
            print(tr('  注意：%s 缺少部分分段檔，該部分為空白',
                     '  NOTE: some segment files of %s are missing; that part is left blank') % channel)
        return out

    def satellite_position(self, scene):
        b3 = self._header(next(iter(scene.files.values()))[0])['b3']
        return float(b3['SubLon']), float(b3['Distance']) - float(b3['EarthEquatorialRadius'])


READER = AhiHsdReader
