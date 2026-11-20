"""TRMM TMI Level 1B 亮溫 (NASA PPS / GES DISC 的 1B.TRMM.TMI.*，GPM 格式的 HDF5 / netCDF-4)

檔名  1B.TRMM.TMI.<校正版本>.<YYYYMMDD>-S<HHMMSS>-E<HHMMSS>.<軌道>.<版本>.nc (或 .HDF5)
來源  https://arthurhouhttps.pps.eosdis.nasa.gov/gpmdata/年/月/日/1B/ (與 L1C 同一個帳號)，或 GES DISC 的 GPM_1BTMI
結構  同 PPS 的 swath 檔 (pps_common.py)，亮溫為 Tb (未經 XCAL 交互校正)，沒有 Quality：
        S1  10.65 GHz V / H              每條掃描 104 點
        S2  19.35 V / H、21.3 V、37.0 V / H   104 點
        S3  85.5 GHz V / H               208 點
      頻道順序固定，寫在下面的 TABLE (2026/11 與同一軌道 L1C 的 long_name 核對)。
舊的 HDF4 格式 (1B11.*) 不支援。
"""
import re

from pps_common import PPSReader


class TMIReader(PPSReader):
    NAME = 'TMI-L1B'
    SENSOR = 'TMI'
    DESCRIPTION = 'TRMM TMI Level 1B 亮溫 (NASA PPS / GES DISC)'
    DESCRIPTION_EN = 'TRMM TMI Level 1B brightness temperature (NASA PPS / GES DISC)'
    PATTERNS = ('1B.TRMM.TMI.*.nc', '1B.TRMM.TMI.*.HDF5')
    VAR = 'Tb'
    FORMAT = 'L1B'
    FILENAME = re.compile(r'^1B\.(?P<sat>TRMM)\.(?P<sensor>TMI)\.[^.]+\.(?P<date>\d{8})-S(?P<s>\d{6})'
                          r'-E(?P<e>\d{6})(?:\.(?P<orbit>\d+))?\.', re.I)
    TABLE = {
        'S1': [(10.65, None, 'V', ''), (10.65, None, 'H', '')],
        'S2': [(19.35, None, 'V', ''), (19.35, None, 'H', ''), (21.3, None, 'V', ''), (37.0, None, 'V', ''),
               (37.0, None, 'H', '')],
        'S3': [(85.5, None, 'V', ''), (85.5, None, 'H', '')],
    }


READER = TMIReader
