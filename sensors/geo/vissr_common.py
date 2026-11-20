"""JMA VISSR (GMS-1~5、GOES-9) 的參數區格式與共用讀取流程 — GMS-VISSR 與 GMS-TAPE 共用
(檔名沒有 '-'，不會被當成讀取器載入)

格式 (結構與欄位名稱同舊版 _VISSR_GMS_IR/_VIS.py、_VISSR_GMS5_IR/_VIS.py；JMA「FORMAT OF VISSR ARCHIVE DATA」)：
  參數區由 2688 bytes (672 個 4-byte word) 的參數段組成，依序為 (SEGMENTS)：
    GMS-1~4  Mode、S/DB、IR 校正、VIS 校正、座標轉換、姿態預報、軌道預報 1、軌道預報 2
    GMS-5    Mode、S/DB、座標轉換、姿態預報、軌道預報 1、軌道預報 2、DCD、VIS 校正、IR1 校正、IR2 校正、WV 校正、
             Split window 校正、(保留) x 2、簡易座標轉換表、β 角取樣
  每條掃描線前有 Line Control Word (16 個 word)：Data ID、行號、掃描時刻、β 角、地球東西邊緣…
parse_segment() 透過「來源」物件取值 (record / array / raw)，同一份欄位定義可用於兩種儲存方式：
  Bytes    JMA 檔案：大端序 I*4 / R*4 / R*8 / ASCII (GMS-VISSR)
  Words36  磁帶映像：每個 4-byte word 存成一個 36-bit word (GMS-TAPE.py)
VissrBase：兩種讀取器共用的部分 — 掃描線索引、校正表、MGIVSR 導航參數與地球邊緣修正 (sensors/vissr_nav.py)。
  子類別只要實作 scenes()、_load() (讀檔頭與各線 LCW 後呼叫 _index()) 與 read()。
"""
import datetime
import re

import numpy as np

from satima.core.i18n import tr
from satima.geo.reader import Channel, GeoReader
from vissr_nav import EA, VissrNav, earth_edge_deviation, mjd_to_datetime

# ==================================================================================================
#  檔頭格式 (大端序；欄位名稱同舊版)
# ==================================================================================================
CONTROL_BLOCK = np.dtype([('Control_block_size', '>i2'), ('Head block number of parameter block', '>i2'),
                          ('Parameter block size', '>i2'), ('I Head block number of image data', '>i2'),
                          ('Total block size of image data', '>i2'), ('Available block size of image data', '>i2'),
                          ('Head valid line number', '>i2'), ('Final valid line number', '>i2'),
                          ('Final data block number', '>i2')])

_MODE_COMMON = [('Satellite number', '>i4'), ('Satellite name', 'S12'), ('Observation time-AD', 'S16'),
                ('Observation time-MJD', '>f8'), ('GMS operation mode', '>i4'), ('DPC operation mode', '>i4'),
                ('VISSR observation mode', '>i4'), ('Scanner selection', '>i4'), ('Sensor selection', '>i4'),
                ('Sensor mode', '>i4'), ('Scan frame mode', '>i4'), ('Scan mode', '>i4'),
                ('Upper limit of scan number', '>i4'), ('Lower limit of scan number', '>i4'),
                ('Equatorial scan line number', '>i4'), ('Spin rate', '>f4'),
                ('VIS-Bit length', '>i4'), ('VIS-Number of lines', '>i4'), ('VIS-Number of pixels', '>i4'),
                ('VIS-Stepping angle', '>f4'), ('VIS-Sampling angle', '>f4'), ('VIS-LCW-pixel size', '>i4'),
                ('VIS-DOC-pixel size', '>i4'), ('VIS-Reserved', '>i4'),
                ('IR-Bit length', '>i4'), ('IR-Number of lines', '>i4'), ('IR-Number of pixels', '>i4'),
                ('IR-Stepping angle', '>f4'), ('IR-Sampling angle', '>f4'), ('IR-LCW-pixel size', '>i4'),
                ('IR-DOC-pixel size', '>i4'), ('IR-Reserved', '>i4'),
                ('Satellite height', '>f4'), ('Earth radius', '>f4'), ('SSP-longitude', '>f4'), ('Reserved', 'S36'),
                ('VIS primary 1ch', '>i4'), ('VIS primary 2ch', '>i4'), ('VIS primary 3ch', '>i4'),
                ('VIS primary 4ch', '>i4'), ('VIS redundant 1ch', '>i4'), ('VIS redundant 2ch', '>i4'),
                ('VIS redundant 3ch', '>i4'), ('VIS redundant 4ch', '>i4'), ('IR1 primary', '>i4'),
                ('IR1 redundant', '>i4')]
MODE = {'GMS': np.dtype(_MODE_COMMON + [('Reserved2', 'S160')]),
        'GMS5': np.dtype(_MODE_COMMON + [('IR2 primary', '>i4'), ('IR2 redundant', '>i4'), ('WV primary', '>i4'),
                                         ('WV redundant', '>i4'), ('Reserved2', 'S144')])}
_Table_Status = np.dtype([('Data segment', '>i4'), ('Data presence', '>i4'), ('Date', '>i4'), ('Time', '>i4')])


def _per_sensor(item, sensors):
    return [('%s-%s' % (s, item), '>f4') for s in sensors]


_CC_ITEMS = ('Stepping angle along line', 'Sampling angle along pixel', 'Central line number of VISSR frame',
             'Central pixel number of VISSR frame', 'Pixel difference of VISSR center from the normal position',
             'Number of sensor elements', 'Total number of VISSR frame lines', 'Total number of VISSR frame pixels')
_CC_HEAD = [('Data segment', '>i4'), ('Data validity', '>i4'), ('Date', '>i4'), ('Time', '>i4'), ('Scheduled', '>f8')]
_CC_MISALIGN = ([('x-component', '>f4'), ('y-component', '>f4'), ('z-component', '>f4')] +
                [('ELMIS(%d,%d)' % (k, l), '>f4') for l in (1, 2, 3) for k in (1, 2, 3)] +
                [(n, '>f4') for n in ('Judgement of observation convergence time', 'Judgement of line convergence',
                                      'E-W angle of Sun-light condense prism', 'N-S angle of Sun-light condense prism',
                                      'pi', 'pi/180', '180/pi', 'Equatorial radius', 'Oblateness of earth',
                                      'Eccentricity of the earth orbit', 'First angle of VISSR observation in S/DB',
                                      'Upper limited line-2nd', 'Lower limited line-1st', 'Upper limited line-3rd',
                                      'Lower limited line-2nd')])
_CC_ORBIT = ([(n, '>f8') for n in ('Epoch time', 'Semi-major axis', 'Eccentricity', 'Orbital inclination',
                                   'Longitude of the ascending node', 'Argument of perigee', 'Mean anomaly',
                                   'Longitude of SSP', 'Latitude of SSP')] + [('Reserved2', 'S8')] +
             [(n, '>f8') for n in ('Attitude-Epoch time', 'Angle between Z-axis and satellite spin axis',
                                   'Angle change rate between spin axis and Z-axis', 'Angle between spin axis and ZY-axis',
                                   'Angle change rate between spin axis and ZT-axis', 'Daily mean of spin rate')])
_CC_SOLAR5 = [(n, '>f4') for n in (
    'VIS-Stepping angle along line-solar', 'IR-Stepping angle along line-solar', 'VIS-Stepping angle along pixel-solar',
    'IR-Stepping angle along pixel-solar', 'VIS-Center line of VISSR frame-solar', 'IR-Center line of VISSR frame-solar',
    'VIS-Center line of VISSR pixel-solar', 'IR-Center line of VISSR pixel-solar',
    'VIS-pixel difference of VISSR center-solar', 'IR-pixel difference of VISSR center-solar',
    'VIS-Sensor elements number-solar', 'IR-Sensor elements number-solar',
    'VIS-Total number of VISSR frame lines-solar', 'IR-Total number of VISSR frame lines-solar',
    'VIS-Total number of Pixels/lines VISSR frame-solar', 'IR-Total number of Pixels/lines VISSR frame-solar')]
_CC_TAIL5 = [(n, '>f4') for n in (
    'Stepping angle along line of IR1', 'Stepping angle along line of IR2', 'Stepping angle along line of WV',
    'Stepping angle along line of VIS', 'Sampling angle along pixel of IR1', 'Sampling angle along pixel of IR2',
    'Sampling angle along pixel of WV', 'Sampling angle along pixel of VIS', 'X component VISSR misalignment',
    'Y component VISSR misalignment', 'Z component VISSR misalignment')]
_COORDINATE_CONVERSION = {
    # GMS-1~4：VIS、IR、VIS-solar、IR-solar 各一組
    'GMS': np.dtype(_CC_HEAD + [f for item in _CC_ITEMS for f in _per_sensor(item, ('VIS', 'IR', 'VIS-solar', 'IR-solar'))] +
                    _CC_MISALIGN + [('Reserved1', 'S140')] + _CC_ORBIT + [('Reserved3', 'S2160')]),
    # GMS-5：VIS、IR1、IR2、WV 各一組
    'GMS5': np.dtype(_CC_HEAD + [f for item in _CC_ITEMS for f in _per_sensor(item, ('VIS', 'IR1', 'IR2', 'WV'))] +
                     _CC_MISALIGN + _CC_SOLAR5 + [('Reserved1', 'S76')] + _CC_ORBIT + [('Reserved3', 'S2116')] + _CC_TAIL5),
}
_PRED_HEAD = [('Data segment', '>i4'), ('Data validity', '>i4'), ('Date', '>i4'), ('Time', '>i4')]
_Attitude_prediction = np.dtype(_PRED_HEAD + [
    ('Start time of attitude prediction', '>f8'), ('End time of attitude prediction', '>f8'),
    ('Interval time of attitude prediction', '>f8'), ('Number of attitude prediction', '>i4'),
    ('Number of attitude prediction data set', '>i4')])
_Attitude_prediction_data = np.dtype([
    ('MJD', '>f8'), ('Date', '>i4'), ('Time', '>i4'), ('Predicted right ascension of attitude', '>f8'),
    ('Predicted declination of attitude', '>f8'), ('Sun-earth angle at prediction time', '>f8'),
    ('Satellite spin rate at prediction time', '>f8'), ('Right ascension of orbital plane at prediction time', '>f8'),
    ('Declination of orbital plane at prediction time', '>f8'), ('Reserved', '>f8'), ('Eclipse flag', '>i4'),
    ('Spin axis flag', '>i4')])
_Orbit_prediction = np.dtype(_PRED_HEAD + [
    ('Start time of orbit prediction', '>f8'), ('End time of orbit prediction', '>f8'),
    ('Interval time of orbit prediction', '>f8'), ('Number of orbit prediction', '>i4'),
    ('Number of orbit prediction data set', '>i4')])
_Orbit_prediction_data = np.dtype(
    [('MJD', '>f8'), ('Date', '>i4'), ('Time', '>i4')] +
    [('%s-%s-component of %s' % (frame, axis, what), '>f8') for frame in ('MICS', 'EFCS')
     for what in ('position', 'velocity') for axis in 'XYZ'] +
    [(n, '>f8') for n in ('Greenwich sidereal time', 'MICS-Azimuth', 'MICS-Elevation', 'EFCS-Azimuth', 'EFCS-Elevation')] +
    [('A%d' % i, '>f8') for i in range(1, 10)] +
    [(n, '>f8') for n in ('X-component of vector', 'Y-component of vector', 'Z-component of vector',
                          'Latitude of SSP', 'Longitude of SSP', 'Satellite height')] +
    [('Eclipse period flag', '>i4'), ('Reserved', '>i4')])
_VIS_Calibration = np.dtype(_PRED_HEAD + [('Sensor group', '>i4')])
_VIS_Calibration_P1 = np.dtype([('Channel number', '>i4'), ('Data validity', '>i4'), ('Date', '>i4'), ('Time', '>i4'),
                                ('Table ID', '>i4')])
_VIS_Calibration_P2 = np.dtype(
    [('VIS channel staircase brightness data-%d' % i, '>f4') for i in range(1, 7)] +
    [('Coefficients table of VIS staircase regression curve-%d' % i, '>f4') for i in range(1, 11)] +
    [(n, '>f4') for n in ('BT-Universal space brightness', 'BT-Solar brightness', 'BC-Universal space brightness',
                          'BC-Solar brightness', 'G', 'V0')] + [('Reserved', 'S36')])
_IR_Calibration_P1 = np.dtype(_PRED_HEAD + [('Sensor group', '>i4'), ('Table ID', '>i4'), ('Reserved', 'S8')])
_IR_Calibration_P2 = np.dtype(
    [('Staircase brightness data-%d' % i, '>f4') for i in range(1, 7)] +
    [('Coefficients table of staircase regression curve-%d' % i, '>f4') for i in range(1, 11)] +
    [(n, '>f4') for n in ('Brightness of space', 'Brightness of black body shutter', 'Reserved1', 'Voltage of space',
                          'Voltage of black body shutter', 'Reserved2', 'G', 'V0', 'Valid shutter temperature',
                          'Valid shutter radiation', 'Shutter temp. 1', 'Redundant mirror temp.', 'Primary mirror temp.',
                          'Baffle FW temp.', 'Baffle AF temp.', '+15 volt auxiliary power supply',
                          'Radiative cooler temp.1', 'Radiative cooler temp.2', 'Electronics module temp.',
                          'Scan mirror temp.', 'Shutter cavity temp.', 'Primary mirror sealed temp.',
                          'Redundant mirror sealed temp.', 'Shutter temp.2')] +
    [('Reserved', 'S8'), ('Flag of valid shutter temperature calculation', '>f4')])
_Pixel_Number = np.dtype([('Line', '>i2'), ('Pixel', '>i2')])
# 簡易座標轉換表的 15 個數為 IBM 浮點數 (舊版以 >u4 讀入未轉換)，這裡轉成一般浮點數
_SIMPLE_NAMES = ('Earth equator Radius', 'Satellite height', 'Stepping angle', 'Sampling angle', 'SSP latitude',
                 'SSP longitude', 'SSP line number', 'SSP pixel number', 'pi', 'Line correction(X)IR1-VIS',
                 'Pixel correction(Y)IR1-VIS', 'Line correction(X)IR1-IR2', 'Pixel correction(Y)IR1-IR2',
                 'Line correction(X)IR1-WV', 'Pixel correction(Y)IR1-WV')
_Simple_coordinate = np.dtype([(n, '>f8') for n in _SIMPLE_NAMES])
LINE_CONTROL_WORD = np.dtype([
    ('Data ID', '>i4'), ('Line number', '>i4'), ('Line name', '>i4'), ('Error line flag', '>i4'),
    ('Error message', '>i4'), ('Mode error flag', '>i4'), ('MJD of VISSR scan time', '>f8'), ('β angle', '>f4'),
    ('West side earth edge', '>i4'), ('East side earth edge', '>i4'), ('Received time (1)', '>i4'),
    ('Received time (1-1)', '>i4'), ('Received time (2)', '>i4'), ('Reserved', 'S8')])

SEGMENT = 2688
# 參數段的順序 (名稱以 '*' 開頭者尚未定義格式，以原始資料保留)
SEGMENTS = {
    'GMS': ('Mode', '*Information_of_SDB_operation', 'IR_calibration', 'VIS_calibration', 'Coordinate_Conversion',
            'Attitude', 'Orbit-1', 'Orbit-2'),
    'GMS5': ('Mode', '*Information_of_SDB_operation', 'Coordinate_Conversion', 'Attitude', 'Orbit-1', 'Orbit-2',
             '*DCD_communication', 'VIS_calibration', 'IR1_calibration', 'IR2_calibration', 'WV_calibration',
             '*Split_window_calibration', '*Reserved_segment_15', '*Reserved_segment_16', 'Simple_coordinate',
             '*Beta_angle_sampling'),
}
# 可見光各線的感測器 (LCW Data ID 低 16 位元)
_VIS_SENSOR = {'GMS': {0x0002: 0, 0x0004: 1, 0x0008: 2, 0x0010: 3}, 'GMS5': {0x0008: 0, 0x0010: 1, 0x0020: 2, 0x0040: 3}}

# ==================================================================================================
#  波段
# ==================================================================================================
# 波段 -> (座標轉換參數的前綴, 校正表前綴)
_PREFIX = {('GMS', 'VIS'): ('VIS', 'VIS'), ('GMS', 'IR'): ('IR', 'IR'), ('GMS5', 'VIS'): ('VIS', 'VIS'),
           ('GMS5', 'IR1'): ('IR1', 'IR1'), ('GMS5', 'IR2'): ('IR2', 'IR2'), ('GMS5', 'IR3'): ('WV', 'WV')}
# (名稱, 中心波長 µm, 解析度 km, 種類, 功能代稱)
_BANDS = {
    'GMS': (('VIS', 0.63, 1.25, 'vis', 'Visible'), ('IR', 11.5, 5.0, 'ir', 'LW Window')),
    'GMS5': (('VIS', 0.73, 1.25, 'vis', 'Visible'), ('IR1', 11.0, 5.0, 'ir', 'LW Window'),
             ('IR2', 12.0, 5.0, 'ir', 'Dirty LW Window'), ('IR3', 6.75, 5.0, 'wv', 'Water Vapor')),
    'GOES9': (('VIS', 0.65, 1.25, 'vis', 'Visible'), ('IR1', 10.7, 5.0, 'ir', 'LW Window'),
              ('IR2', 12.0, 5.0, 'ir', 'Dirty LW Window'), ('IR3', 6.75, 5.0, 'wv', 'Water Vapor')),
}
CHANNELS = {key: {b[0]: Channel(*b[:4], alias=b[4]) for b in bands} for key, bands in _BANDS.items()}


def _ibm32(words):
    """IBM 32 位元浮點數 -> float64。"""
    u = np.asarray(words, dtype=np.uint32).astype(np.int64)
    sign = np.where(u >> 31, -1.0, 1.0)
    return sign * ((u & 0xFFFFFF) / float(1 << 24)) * 16.0 ** (((u >> 24) & 0x7F) - 64)


def satellite_name(raw, fallback):
    m = re.search(r'(GMS|GOES)\s*-?\s*(\d+)', raw.decode('ascii', 'replace').upper())
    return '%s-%s' % m.groups() if m else fallback


def observation_time(mode):
    text = mode['Observation time-AD'].decode('ascii', 'replace')
    try:
        return datetime.datetime.strptime(text[:12], '%Y%m%d%H%M')
    except ValueError:
        return mjd_to_datetime(mode['Observation time-MJD']).replace(second=0, microsecond=0)


# ==================================================================================================
#  檔頭解析
# ==================================================================================================
class Bytes:
    """參數段的來源：JMA 檔案的 bytes (大端序)。offset 為 bytes；磁帶的 36-bit word 見 GMS-TAPE.py 的 Words36。"""

    def __init__(self, buf):
        self.buf = bytes(buf)

    def record(self, dtype, offset=0):
        return np.frombuffer(self.buf, dtype, 1, offset)[0]

    def array(self, dtype, count, offset=0):
        return np.frombuffer(self.buf, dtype, count, offset)

    def raw(self):
        return self.buf


def parse_segment(name, src, gen, h):
    """一個參數段 (2688 bytes = 672 個 word) -> 放進 h (名稱同舊版)。src：Bytes 或 Words36。"""
    if name.startswith('*'):
        h[name[1:]] = src.raw()
    elif name == 'Mode':
        h['Mode_Block'] = src.record(MODE[gen])
        h['Table_Status'] = src.array(_Table_Status, 10 if gen == 'GMS' else 15, 400)
    elif name == 'Coordinate_Conversion':
        h['Coordinate_Conversion'] = src.record(_COORDINATE_CONVERSION[gen])
    elif name == 'Attitude':
        h['Attitude_prediction'] = src.record(_Attitude_prediction)
        h['Attitude_prediction_data'] = src.array(_Attitude_prediction_data, 33, 48)
    elif name.startswith('Orbit-'):
        k = name[-1]
        h['Orbit_prediction-' + k] = src.record(_Orbit_prediction)
        h['Orbit_prediction_data-' + k] = src.array(_Orbit_prediction_data, 9, 48)
    elif name == 'VIS_calibration':
        h['VIS_Calibration'] = src.record(_VIS_Calibration)
        for i in range(4):
            o = 20 + 400 * i
            h['VIS%d_Calibration_P1' % (i + 1)] = src.record(_VIS_Calibration_P1, o)
            h['VIS%d_Calibration_Table' % (i + 1)] = src.array('>f4', 64, o + 20)
            h['VIS%d_Calibration_P2' % (i + 1)] = src.record(_VIS_Calibration_P2, o + 276)
    elif name.endswith('_calibration'):              # IR / IR1 / IR2 / WV
        p = name.split('_')[0]
        h[p + '_Calibration_P1'] = src.record(_IR_Calibration_P1)
        h[p + '_Radiation_Table'] = src.array('>f4', 256, 32)
        h[p + '_BT_Table'] = src.array('>f4', 256, 1056)
        h[p + '_Calibration_P2'] = src.record(_IR_Calibration_P2, 2080)
    elif name == 'Simple_coordinate':
        h['Pixel Number'] = src.array(_Pixel_Number, 625, 0)
        values = _ibm32(src.array('>u4', 15, 2500))
        h['Simple_coordinate'] = np.array(tuple(values), dtype=_Simple_coordinate)


# ==================================================================================================
#  共用讀取流程
# ==================================================================================================
class VissrBase(GeoReader):
    SENSOR = 'VISSR'
    CHANNELS = {**CHANNELS['GMS'], **CHANNELS['GMS5']}
    TRUECOLOR = {}                  # 只有一個寬頻可見光

    def channels(self, scene):
        table = CHANNELS['GOES9' if scene.satellite.startswith('GOES') else scene.extra.get('gen', 'GMS5')]
        return [ch for name, ch in table.items() if name in scene.files]

    # ---- 檔頭與影像結構 ----------------------------------------------------
    def header(self, scene, channel):
        """完整檔頭 {紀錄名稱: numpy 結構} (名稱同舊版)，含每條掃描線的 Line_Control_Word。"""
        return self._load(scene, channel)['header']

    def _load(self, scene, channel):
        """該波段的影像結構 (快取)：讀檔頭與各線 LCW 放進 h 後呼叫 _index()，再加上讀取影像所需的位置資訊。"""
        raise NotImplementedError

    def _index(self, h, gen, channel, usable=None):
        """檔頭 h (含各紀錄的 Line_Control_Word) -> 掃描線索引、校正表與導航的掃描參數。
        實際資料常缺行、重複或在尾端附加無效紀錄，所以依每筆紀錄 LCW 的行號放到影像中 (重複的行取第一筆)。
        usable：額外的有效紀錄條件 (布林陣列)。
        回傳 dict：header、gen、kind、line0 (第 0 列的行號)、rec_of_row (各列的紀錄索引，缺行 -1)、cc_prefix、
                   lut (計數值 -> 亮溫 °C / 反照率)、sensor (可見光各紀錄的感測器 0~3)、scan。"""
        kind = 'VIS' if channel == 'VIS' else 'IR'
        lcw = h['Line_Control_Word']
        cc, mode = h['Coordinate_Conversion'], h['Mode_Block']
        cc_p, cal_p = _PREFIX[gen, channel]
        frame = max(int(cc[cc_p + '-Total number of VISSR frame lines']), int(mode[kind + '-Number of lines']))
        lines = lcw['Line number'].astype(np.int64)
        seg = lcw['Data ID'] & 0xFFFF
        valid = (lines >= 1) & (lines <= frame) & (seg != 0) & (lcw['Error message'] == 0)
        if usable is not None:
            valid &= usable
        if not valid.any():
            raise ValueError(tr('%s 沒有有效的掃描線', '%s has no valid scan lines') % channel)
        uniq, first = np.unique(lines[valid], return_index=True)          # 重複的行取第一筆
        line0 = int(uniq[0])
        rec_of_row = np.full(int(uniq[-1]) - line0 + 1, -1, dtype=np.int64)
        rec_of_row[uniq - line0] = np.flatnonzero(valid)[first]
        info = dict(header=h, gen=gen, kind=kind, line0=line0, rec_of_row=rec_of_row, cc_prefix=cc_p)
        if kind == 'VIS':
            sensor_map = _VIS_SENSOR[gen]
            sensor = np.array([sensor_map.get(int(s), (int(l) - 1) % 4) for s, l in zip(seg, lines)], dtype=np.intp)
            info['sensor'] = sensor
            lut = np.full((4, 256), np.nan, dtype=np.float32)
            for i in range(4):
                lut[i, :64] = h['VIS%d_Calibration_Table' % (i + 1)]
            info['lut'] = lut
        else:
            lut = (h[cal_p + '_BT_Table'] - 273.15).astype(np.float32)
            lut[0] = np.nan                                  # 計數值 0：缺資料 (GOES-9 的太空)
            info['lut'] = lut
        info['scan'] = self._scan_params(h, cc_p)
        return info

    @staticmethod
    def _scan_params(h, p):
        """MGIVSR 的掃描參數 (同 JMA 程式與舊版：觀測時刻以 Scheduled 起算、自旋速率取每日平均)。
        LCW 記錄的各線掃描時刻在 GMS-5 與此一致；GMS-1~4 差約 30~45 秒、GOES-9 差約 3.5 分鐘，
        但導航參數是配合 MGIVSR 產生的，地球邊緣比對兩者結果相同，所以依 JMA 的公式。"""
        cc = h['Coordinate_Conversion']
        elmis = np.array([[cc['ELMIS(%d,%d)' % (k, l)] for l in (1, 2, 3)] for k in (1, 2, 3)], dtype=np.float64)
        # 像素號中心要減去「VISSR 中心與標準位置的像素差」(GOES-9 紅外為 1.5；GMS 為 0)：
        # 以每條線實際的地球東西邊緣 (LCW) 比對，GOES-9 不減時東西偏 -1.54 像素，減去後為 0
        pixel0 = float(cc[p + '-Central pixel number of VISSR frame']) - \
            float(cc[p + '-Pixel difference of VISSR center from the normal position'])
        return dict(step=float(cc[p + '-Stepping angle along line']), samp=float(cc[p + '-Sampling angle along pixel']),
                    line0=float(cc[p + '-Central line number of VISSR frame']), pixel0=pixel0,
                    sensors=float(cc[p + '-Number of sensor elements']) or 1., elmis=elmis,
                    t0=float(cc['Scheduled']), spin=float(cc['Daily mean of spin rate']))

    # ---- 地球邊緣修正 --------------------------------------------------------
    def edge_correction(self, scene):
        """MSC Technical Note No.26 的地球邊緣修正：由場景紅外波段各掃描線的實際地球東西邊緣 (LCW) 與導航預測的
        邊緣比較，求影像的南北 (行)、東西 (像素) 偏移 (紅外單位)。結果存在 scene.extra，各波段共用。"""
        if 'edge_correction' in scene.extra:
            return scene.extra['edge_correction']
        res = None
        ir = next((c for c in ('IR', 'IR1', 'IR2', 'IR3') if c in scene.files), None)
        if ir is not None:
            info = self._load(scene, ir)
            rows = np.flatnonzero(info['rec_of_row'] >= 0)[::2]
            lcw = info['header']['Line_Control_Word'][info['rec_of_row'][rows]]
            west, east = lcw['West side earth edge'], lcw['East side earth edge']
            ok = (west > 0) & (east > west)
            res = earth_edge_deviation(self._nav(scene, ir, None), rows[ok], west[ok], east[ok])
            if res is not None and (abs(res['line']) > 10. or abs(res['pixel']) > 10.):
                print(tr('  [略過地球邊緣修正] 估計的偏移過大 (南北 %+.1f 行、東西 %+.1f 像素)，邊緣資料可能有問題',
                         '  [skip earth edge correction] the estimated offset is too large (N-S %+.1f lines, '
                         'E-W %+.1f pixels); the edge data may be bad') % (
                    res['line'], res['pixel']))
                res = None
            if res is not None:
                res.update(channel=ir, step=info['scan']['step'], samp=info['scan']['samp'])
                print(tr('  地球邊緣修正 (MSC Tech. Note No.26，%s)：南北 %+.2f 行、東西 %+.2f 像素 (%d 條掃描線)',
                         '  Earth edge correction (MSC Tech. Note No.26, %s): N-S %+.2f lines, E-W %+.2f pixels '
                         '(%d scan lines)') % (
                    ir, res['line'], res['pixel'], res['n']))
        scene.extra['edge_correction'] = res
        return res

    # ---- 導航 --------------------------------------------------------------
    def nav(self, scene, channel):
        corr = self.edge_correction(scene) if self.cfg.navigation.earth_edge_correction else None
        return self._nav(scene, channel, corr)

    def _nav(self, scene, channel, corr):
        info = self._load(scene, channel)
        scan = dict(info['scan'])
        if corr is not None:                    # 紅外的偏移換算成該波段的行、像素
            scan['line0'] += corr['line'] * corr['step'] / scan['step']
            scan['pixel0'] += corr['pixel'] * corr['samp'] / scan['samp']
        h = info['header']
        att = h['Attitude_prediction_data']
        orb = np.concatenate([h['Orbit_prediction_data-1'], h['Orbit_prediction_data-2']])
        npa = np.stack([orb['A%d' % i] for i in range(1, 10)], axis=-1).reshape(-1, 3, 3).transpose(0, 2, 1)
        return VissrNav(
            (len(info['rec_of_row']), info['npix']), scan,
            attitude=dict(mjd=att['MJD'], alpha=att['Predicted right ascension of attitude'],
                          delta=att['Predicted declination of attitude'], beta=att['Sun-earth angle at prediction time']),
            orbit=dict(mjd=orb['MJD'], sat=np.stack([orb['EFCS-%s-component of position' % a] for a in 'XYZ'], axis=-1),
                       gst=orb['Greenwich sidereal time'], sun_az=orb['EFCS-Azimuth'], sun_el=orb['EFCS-Elevation'], npa=npa),
            first_line=info['line0'], key=(scene.files[channel], channel, corr is not None),
            engine=str(self.cfg.navigation.engine).lower())

    def satellite_position(self, scene):
        cc = self.header(scene, next(iter(scene.files)))['Coordinate_Conversion']
        return float(cc['Longitude of SSP']), float(cc['Semi-major axis']) - EA / 1000.
