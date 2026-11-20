"""MetOp AMSU-A Level 1B (EUMETSAT 的 EPS 原生格式；Data Store / EUMETCast 的 AMSA_xxx_1B_M0x_*)

檔名  AMSA_xxx_1B_<M01|M02|M03>_<起>Z_<訖>Z_<模式>_<處置>_<處理時間>Z   (.nat，或 Data Store 下載的 .zip；
      zip 裡的 .nat 直接讀，不必解壓縮)      M02 = MetOp-A、M01 = MetOp-B、M03 = MetOp-C
結構  (依 EUMETSAT「ATOVS Level 1b Product Guide」EUM/OPS-EPS/MAN/04/0030 第 11 節；大端)
      一連串的記錄，每筆開頭是 20 bytes 的通用檔頭 GRH：
        類別 u1 (1 MPHR、3 IPR、4 GEADR、5 GIADR、8 MDR)、儀器群 u1、子類別 u1、版本 u1、記錄長度 u4、
        起始時間 (2000-01-01 起算的日數 u2 + 當日毫秒 u4)、結束時間 (同)
      MPHR (3307 bytes) 是「名稱 = 值」的文字。MDR (子類別 2，版本 3 / 4，3464 bytes) 一筆 = 一條掃描 (8 秒)，30 個視場：
        offset   22  SCENE_RADIANCE   int32 x 30 視場 x 15 頻道，1e-7 mW/(m² sr cm⁻¹)
               1842  ANGULAR_RELATION int16 x 30 x 4 (太陽天頂、衛星天頂、太陽方位、衛星方位；0.01 度)
               2082  EARTH_LOCATION   int32 x 30 x 2 (緯度、經度；1e-4 度)
      (陣列的排列順序 —— 先視場後頻道 —— 以樣本檔確認：2026/11。)
亮溫  輻射量 R 以 Planck 函數在頻道中心波數 ν 換成亮溫：T = C2 ν / ln(1 + C1 ν³ / R)，
      C1 = 1.191044e-5 mW/(m² sr cm⁻⁴)、C2 = 1.438769 K/cm⁻¹；頻帶修正係數 A = 0、B = 1 (產品指南附錄 A)。
      中心波數用附錄 A 的 MetOp-B 值 (各顆衛星只差在第 5~6 位，影響 < 0.1 K)。
頻道  15 個 (橫向掃描，星下點約 48 km，極化隨掃描角旋轉，標的是星下點的極化)：
        1 23.8V  2 31.4V  3 50.3V  4 52.8V  5 53.596±0.115H  6 54.4H  7 54.94V  8 55.5H
        9 57.29H (fo = 57.290344)  10 fo±0.217H  11~14 fo±0.3222±0.048 / 0.022 / 0.010 / 0.0045 H  15 89V
      3~14 是氧氣吸收帶的溫度探空頻道；看得到地表與降水的是 1、2、(3)、15。沒有雙極化，只有單頻道亮溫。
      取樣很粗 (30 點、48 km 以上)，畫颱風只看得出大致的結構。
故障的頻道  輻射量為填充值 (負值) 的取樣點是無資料。有效取樣不到一半的頻道整個視為無資料
      (樣本中 MetOp-B 的 89 GHz 已故障，只剩 3% 的掃描有值，而且是固定的 280.7 K)。
      雜訊變大但仍有值的頻道 (樣本中 MetOp-B 的 54.4 GHz：76~417 K) 照原樣顯示；
      MDR 的 CALIBRATION_QUALITY 旗標在正常的頻道上也常被設定 (MetOp-C 的 50.3、52.8、53.596、55.5 GHz)，所以不拿來遮資料。
"""
import datetime
import io
import os
import re
import zipfile

import numpy as np

from satima.core.i18n import tr
from satima.mw.reader import Channel, Geo, Scene, SwathReader, clean_latlon

_NAME = re.compile(r'^AMSA_\w{3}_1B_(M0\d)_(\d{14})Z_(\d{14})Z_\w_\w_\d{14}Z(?:\.nat|\.zip)?$', re.I)
_SATELLITES = {'M01': 'MetOp-B', 'M02': 'MetOp-A', 'M03': 'MetOp-C'}
_EPOCH2000 = (datetime.datetime(2000, 1, 1) - datetime.datetime(1970, 1, 1)).total_seconds()
_MDR_SIZE, _NFOV, _NCH = 3464, 30, 15
_RAD, _LOC = 22, 2082
C1, C2 = 1.191044e-5, 1.438769
DEAD = 0.5              # 有效取樣的比例低於此值的頻道視為故障
# (名稱, 中心頻率 GHz, 星下點極化, 中心波數 cm⁻¹)
CHANNELS = (
    ('23.8', 23.8, 'V', 0.793897), ('31.4', 31.4, 'V', 1.047421), ('50.3', 50.3, 'V', 1.677830),
    ('52.8', 52.8, 'V', 1.761235), ('53.596±0.115', 53.596, 'H', 1.787785), ('54.4', 54.4, 'H', 1.814590),
    ('54.94', 54.94, 'V', 1.832608), ('55.5', 55.5, 'H', 1.851295), ('57.29', 57.29, 'H', 1.911001),
    ('57.29±0.217', 57.29, 'H', 1.911001), ('57.29±0.322±0.048', 57.29, 'H', 1.911001),
    ('57.29±0.322±0.022', 57.29, 'H', 1.911001), ('57.29±0.322±0.01', 57.29, 'H', 1.911001),
    ('57.29±0.322±0.0045', 57.29, 'H', 1.911001), ('89', 89.0, 'V', 2.968887),
)


def load(path):
    """.nat 檔，或內含 .nat 的 .zip -> 整個檔的位元組。"""
    if zipfile.is_zipfile(path):
        with zipfile.ZipFile(path) as z:
            names = [n for n in z.namelist() if n.lower().endswith('.nat')] or \
                    [n for n in z.namelist() if os.path.basename(n).upper().startswith('AMSA_')]
            if not names:
                raise ValueError(tr('%s 裡沒有 AMSA_*.nat', 'no AMSA_*.nat inside %s') % os.path.basename(path))
            return z.read(names[0])
    with open(path, 'rb') as f:
        return f.read()


def parse(buf):
    """整個檔 -> (輻射量 (掃描, 30, 15) float64 [mW/(m² sr cm⁻¹)], 緯度, 經度 (掃描, 30), 掃描時間 (1970 起算的秒數))。"""
    mv = memoryview(buf)
    rad, loc, times = [], [], []
    pos, n = 0, len(buf)
    while pos + 20 <= n:
        cls, sub = mv[pos], mv[pos + 2]
        size = int.from_bytes(mv[pos + 4:pos + 8], 'big')
        if size < 20 or pos + size > n:
            break
        if cls == 8 and sub == 2 and size == _MDR_SIZE:
            day = int.from_bytes(mv[pos + 8:pos + 10], 'big')
            ms = int.from_bytes(mv[pos + 10:pos + 14], 'big')
            times.append(_EPOCH2000 + day * 86400. + ms / 1000.)
            rad.append(np.frombuffer(mv, '>i4', _NFOV * _NCH, pos + _RAD))
            loc.append(np.frombuffer(mv, '>i4', _NFOV * 2, pos + _LOC))
        pos += size
    if not rad:
        raise ValueError(tr('沒有 AMSU-A 1B 的量測記錄 (MDR，%d bytes)',
                            'no AMSU-A 1B measurement records (MDR, %d bytes)') % _MDR_SIZE)
    rad = np.stack(rad).reshape(-1, _NFOV, _NCH).astype(np.float64) * 1e-7
    loc = np.stack(loc).reshape(-1, _NFOV, 2).astype(np.float64) * 1e-4
    return rad, loc[..., 0], loc[..., 1], np.array(times)


def brightness(rad, wavenumber):
    """輻射量 -> 亮溫 (K)；輻射量 <= 0 (無效) 為 NaN。"""
    with np.errstate(invalid='ignore', divide='ignore'):
        t = C2 * wavenumber / np.log1p(C1 * wavenumber ** 3 / rad)
    return np.where(rad > 0., t, np.nan).astype(np.float32)


class AMSUAReader(SwathReader):
    NAME = 'AMSUA-L1B'
    SENSOR = 'AMSU-A'
    SOURCE = 'EUMETSAT'
    DESCRIPTION = 'MetOp AMSU-A Level 1B (EUMETSAT EPS 原生格式 AMSA_xxx_1B_M0x_*，.nat 或 .zip)'
    DESCRIPTION_EN = 'MetOp AMSU-A Level 1B (EUMETSAT EPS native format AMSA_xxx_1B_M0x_*, .nat or .zip)'
    PATTERNS = ('AMSA_*_1B_M0*',)

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            m = _NAME.match(os.path.basename(name))
            if not m:
                continue
            start, end = (datetime.datetime.strptime(m.group(k), '%Y%m%d%H%M%S') for k in (2, 3))
            out.append(Scene(_SATELLITES.get(m.group(1).upper(), m.group(1).upper()), 'AMSU-A', start, end,
                             file=os.path.join(folder, name), format='L1B'))
        return out

    def channels(self, scene):
        return [Channel(name + pol, freq, pol, 'S', key=k) for k, (name, freq, pol, _) in enumerate(CHANNELS)]

    def _data(self, scene):
        if 'data' not in scene.extra:
            scene.extra['data'] = parse(load(scene.file))
        return scene.extra['data']

    def geo(self, scene, group):
        _, lat, lon, time = self._data(scene)
        lat, lon = clean_latlon(lat, lon)
        return [Geo(lat, lon, time, (scene.file, group))]

    def tb(self, scene, channel):
        rad = self._data(scene)[0]
        data = brightness(rad[:, :, channel.key], CHANNELS[channel.key][3])
        if np.isfinite(data).mean() < DEAD:
            data[:] = np.nan                    # 故障的頻道只剩零星幾條未校正的掃描：整個頻道視為無資料
        return [data]


READER = AMSUAReader
