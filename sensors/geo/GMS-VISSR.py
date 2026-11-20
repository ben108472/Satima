"""GMS-1~5、GOES-9 VISSR — JMA/MSC VISSR 影像資料 (gzip，可放在 tar 中)

檔案 (tar 不必解開，也可直接放 tar 內的檔案)：
  GMS-1~4  VISSR_GMS3_198907150600.tar 內的 VS890715.Z06.gz (可見光)、IR890715.Z06.gz (紅外)
  GMS-5    VISSR_GMS5_199806210400.tar 內的 VISSR_19980621_0331_IR1.A.IMG.gz (VIS、IR1、IR2、IR3)
  GOES-9   VISSR_GOES9_200406160000.tar (2003/05~2005/06 代替 GMS-5，格式同 GMS-5)
  衛星名稱與時間取自檔頭 (Mode Block 的觀測開始時刻，同 MTSAT 等以掃描開始時刻命名)，
  與 tar 檔名的整點不一定相同 (例如 GMS-3 06:00 的 tar 為 05:31 開始掃描)。
波段 (JMA 名稱；星下點解析度 可見光 1.25 km、紅外 5 km)：
  GMS-1~4  VIS 0.50~0.75 µm、IR 10.5~12.5 µm
  GMS-5    VIS 0.55~0.90、IR1 10.5~11.5、IR2 11.5~12.5、IR3 6.5~7.0 µm (水氣)
  GOES-9   VIS 0.55~0.75、IR1 10.2~11.2、IR2 11.5~12.5、IR3 6.47~7.02 µm (GOES Imager 第 1、4、5、3 頻道)
格式 (結構與欄位名稱同舊版 _VISSR_GMS_IR/_VIS.py、_VISSR_GMS5_IR/_VIS.py；參數段的欄位定義、校正表、導航參數與
地球邊緣修正在 sensors/vissr_common.py，與 GMS-TAPE 共用)：
  檔案由固定長度的紀錄組成：GMS-1~4 紅外 7008、可見光 13504 bytes (一個區塊 = 2 筆紀錄)；
  GMS-5 紅外 3664、可見光 13504 bytes (一個區塊 = 1 筆紀錄)。控制區塊 (Control Block) 記錄參數區與影像區的起始區塊。
  參數區由 2688 bytes 的參數段組成 (每筆紀錄放 1、2 或 4 段，其餘補空)，依序為：
    GMS-1~4  Mode、S/DB、IR 校正、VIS 校正、座標轉換、姿態預報、軌道預報 1、軌道預報 2
    GMS-5    Mode、S/DB、座標轉換、姿態預報、軌道預報 1、軌道預報 2、DCD、VIS 校正、IR1 校正、IR2 校正、WV 校正、
             Split window 校正、(保留) x 2、簡易座標轉換表、β 角取樣
  影像區每筆紀錄 = 一條掃描線：Line Control Word (64 bytes) + 紅外 256 / 可見光 64 bytes + 計數值 (每像素 1 byte)。
  實際資料常缺行、重複或在尾端附加無效紀錄，所以依每筆紀錄 LCW 的行號 (Line number) 放到影像中，缺行為 NaN。
校正：紅外 = 檔頭的計數值->亮溫表 (256 項)，計數值 0 (缺資料；GOES-9 的太空) 為 NaN。
      可見光 = 4 個感測器各自的計數值->反照率表 (64 項)，每條線依 LCW 的 Data ID 選用 (同舊版)。
導航：JMA MGIVSR (sensors/vissr_nav.py)，參數取自檔頭的座標轉換參數與姿態、軌道預報。
  中心像素號減去「Pixel difference of VISSR center from the normal position」(GOES-9 紅外 1.5，其餘 0)。
  [navigation] earth_edge_correction = true 時，再以各掃描線 LCW 的實際地球東西邊緣修正偏移
  (MSC Technical Note No.26；由紅外波段求出，換算後套用到同一時次的所有波段)。
完整檔頭：reader.header(scene, 'IR1') 回傳 {紀錄名稱: numpy 結構}，名稱與欄位同舊版，例如
    h['Coordinate_Conversion']['Longitude of SSP']、h['Orbit_prediction_data-1']['EFCS-X-component of position']、
    h['VIS2_Calibration_Table']、h['Line_Control_Word']['MJD of VISSR scan time'] (每條掃描線)。
  尚未定義格式的參數段 (S/DB、DCD、Split window 等) 以原始 bytes 保留。
"""
import gzip
import os
import re
import tarfile

import numpy as np

from satima.core.i18n import tr
from satima.geo.reader import FileCache, Scene
from vissr_common import (CONTROL_BLOCK, LINE_CONTROL_WORD, MODE, SEGMENT, SEGMENTS, Bytes, VissrBase, observation_time,
                          parse_segment, satellite_name)

# (世代, 種類) -> (紀錄長度, 每區塊紀錄數, 每紀錄參數段數, 計數值在紀錄中的起點)
_LAYOUT = {('GMS', 'VIS'): (13504, 2, 4, 128), ('GMS', 'IR'): (7008, 2, 2, 320),
           ('GMS5', 'VIS'): (13504, 1, 4, 128), ('GMS5', 'IR'): (3664, 1, 1, 320)}

_NAME_GMS = re.compile(r'^(VS|IR)(\d{6})\.Z(\d{2})(\.gz)?$', re.I)
_NAME_GMS5 = re.compile(r'^VISSR_(\d{8})_(\d{4})_(VIS|IR[123])\.A\.IMG(\.gz)?$', re.I)


# ==================================================================================================
#  檔頭解析 (參數段的格式定義與 parse_segment 在 vissr_common.py)
# ==================================================================================================
def read_header(stream, gen, kind):
    """VISSR 串流 (已在檔案開頭) -> (檔頭 dict, 影像區起點 bytes)。讀完參數區為止。"""
    reclen, per_block, per_record, _ = _LAYOUT[gen, kind]
    block = reclen * per_block
    raw = stream.read(CONTROL_BLOCK.itemsize)
    cb = Bytes(raw).record(CONTROL_BLOCK)
    h = {'Control_Block': cb}
    par0 = (int(cb['Head block number of parameter block']) - 1) * block
    image0 = (int(cb['I Head block number of image data']) - 1) * block
    if not 0 < par0 < image0:
        raise ValueError(tr('控制區塊不正確 (參數區 %d、影像區 %d)',
                            'invalid control block (parameter area %d, image area %d)') % (par0, image0))
    raw += stream.read(image0 - len(raw))
    h['Address_table'] = np.frombuffer(raw[32:int(cb['Control_block_size']) * block], '>i2')
    for j, name in enumerate(SEGMENTS[gen]):
        o = par0 + (j // per_record) * reclen + (j % per_record) * SEGMENT
        if o + SEGMENT <= image0:
            parse_segment(name, Bytes(raw[o:o + SEGMENT]), gen, h)
    return h, image0


def read_mode(stream, gen, kind):
    """只讀到 Mode Block (資料列表用，不必讀整個參數區)。"""
    reclen, per_block, _, _ = _LAYOUT[gen, kind]
    raw = stream.read(CONTROL_BLOCK.itemsize)
    cb = Bytes(raw).record(CONTROL_BLOCK)
    par0 = (int(cb['Head block number of parameter block']) - 1) * reclen * per_block
    raw += stream.read(par0 + MODE[gen].itemsize - len(raw))
    return Bytes(raw[par0:]).record(MODE[gen])


class VissrReader(VissrBase):
    NAME = 'GMS-VISSR'
    DESCRIPTION = 'GMS-1~5 / GOES-9 VISSR (JMA .Z?? / .A.IMG，gzip / tar)'
    DESCRIPTION_EN = 'GMS-1~5 / GOES-9 VISSR (JMA .Z?? / .A.IMG, gzip / tar)'
    PATTERNS = ('VISSR_GMS*.tar', 'VISSR_GOES9*.tar', 'VISSR_*.A.IMG', 'VISSR_*.A.IMG.gz',
                'VS??????.Z??', 'VS??????.Z??.gz', 'IR??????.Z??', 'IR??????.Z??.gz')

    def __init__(self, cfg):
        super().__init__(cfg)
        self._tars = FileCache(tarfile.open, maxsize=8)
        self._info = {}

    # ---- 檔案 --------------------------------------------------------------
    @staticmethod
    def _identify(base):
        """檔名 -> (世代, 波段)；不認得時 None。"""
        m = _NAME_GMS.match(base)
        if m:
            return 'GMS', 'VIS' if m.group(1).upper() == 'VS' else 'IR'
        m = _NAME_GMS5.match(base)
        if m:
            return 'GMS5', m.group(3).upper()
        return None

    def _open(self, source):
        """回傳可依序讀取 (可往後 seek) 的串流。"""
        path, member = source
        f = open(path, 'rb') if member is None else self._tars.get(path).extractfile(member)
        name = member or path
        return gzip.GzipFile(fileobj=f) if name.lower().endswith('.gz') else f

    def scenes(self, folder='.'):
        groups = {}

        def add(base, source):
            ident = self._identify(base)
            if ident is None:
                return
            gen, band = ident
            try:
                with self._open(source) as f:
                    mode = read_mode(f, gen, 'VIS' if band == 'VIS' else 'IR')
            except (OSError, ValueError, EOFError, tarfile.TarError) as err:
                print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (base, err))
                return
            sat = satellite_name(mode['Satellite name'], 'GMS-5' if gen == 'GMS5' else 'GMS')
            t = observation_time(mode)
            key = (sat, t)
            if key not in groups:
                groups[key] = Scene(sat, t, 'FD', extra={'gen': gen})
            # 同一個檔案同時有單獨的檔案與 tar 內的成員時，用單獨的檔案
            if band not in groups[key].files or source[1] is None:
                groups[key].files[band] = source

        for name in self.find(folder):
            path = os.path.join(folder, name)
            if name.lower().endswith('.tar'):
                try:
                    members = [m.name for m in self._tars.get(path).getmembers() if m.isfile()]
                except (tarfile.TarError, OSError) as err:
                    print(tr('  [略過] %s：%s', '  [skip] %s: %s') % (name, err))
                    continue
                for member in members:
                    add(os.path.basename(member), (path, member))
            else:
                add(name, (path, None))
        return sorted(groups.values(), key=lambda s: (s.time, s.satellite))

    # ---- 檔頭與影像結構 ----------------------------------------------------
    def _load(self, scene, channel):
        source = scene.files[channel]
        if source in self._info:
            return self._info[source]
        gen = scene.extra.get('gen', 'GMS5')
        kind = 'VIS' if channel == 'VIS' else 'IR'
        reclen, per_block, _, count0 = _LAYOUT[gen, kind]
        with self._open(source) as f:
            h, image0 = read_header(f, gen, kind)
            # 影像區各紀錄的 LCW (需解壓整個檔案一次)
            nmax = int(h['Control_Block']['Total block size of image data']) * per_block
            parts, n = [], 0
            while n < nmax:
                buf = f.read(reclen * min(1024, nmax - n))
                k = len(buf) // reclen
                if k == 0:
                    break
                parts.append(np.frombuffer(buf[:k * reclen], np.uint8).reshape(k, reclen)[:, :64].copy())
                n += k
        lcw = np.concatenate(parts).tobytes() if parts else b''
        h['Line_Control_Word'] = np.frombuffer(lcw, LINE_CONTROL_WORD)
        info = self._index(h, gen, channel)
        info.update(reclen=reclen, image0=image0, count0=count0, npix=reclen - count0)
        self._info[source] = info
        return info

    # ---- 讀取 --------------------------------------------------------------
    def read(self, scene, channel, window):
        r0, r1, c0, c1 = window
        info = self._load(scene, channel)
        out = np.full((r1 - r0, c1 - c0), np.nan, dtype=np.float32)
        recs = info['rec_of_row'][r0:r1]
        have = recs >= 0
        if not have.any():
            return out
        k0, k1 = int(recs[have].min()), int(recs[have].max()) + 1
        reclen = info['reclen']
        with self._open(scene.files[channel]) as f:
            f.seek(info['image0'] + k0 * reclen)              # gzip：只解壓到需要的最後一條線
            buf = f.read((k1 - k0) * reclen)
        n = len(buf) // reclen
        block = np.frombuffer(buf[:n * reclen], np.uint8).reshape(n, reclen)
        idx = recs - k0
        ok = have & (idx < n)
        counts = block[idx[ok], info['count0'] + c0:info['count0'] + c1]
        if info['kind'] == 'VIS':
            out[ok] = info['lut'][info['sensor'][recs[ok]][:, None], counts]
        else:
            out[ok] = info['lut'][counts]
        return out


READER = VissrReader
