"""GCOM-W1 AMSR2 Level 1B 亮溫 (JAXA G-Portal，HDF5)

檔名  GW1AM2_<YYYYMMDDhhmm>_<路徑><A/D>_L1SGBTBR_<版本>.h5   (A = 升軌、D = 降軌，一個檔約半圈)
頻道  6.925、7.3、10.65、18.7、23.8、36.5、89.0 GHz 的 V / H
結構與低頻經緯度的推算見 amsr_common.py。
(2026/11 以同一軌道的 NASA L1C 驗證：10.65~36.5 GHz 四個頻段的位置相差不到 1 公尺。)
"""
import re

from amsr_common import AMSRReader


class AMSR2Reader(AMSRReader):
    NAME = 'AMSR2-L1B'
    SENSOR = 'AMSR2'
    SATELLITE = 'GCOM-W1'
    DESCRIPTION = 'GCOM-W1 AMSR2 Level 1B 亮溫 (JAXA G-Portal)'
    DESCRIPTION_EN = 'GCOM-W1 AMSR2 Level 1B brightness temperature (JAXA G-Portal)'
    PATTERNS = ('GW1AM2_*_L1SGBTBR_*.h5',)
    FILENAME = re.compile(r'^GW1AM2_(\d{12})_(\d{3})([AD])_L1SGBTBR_', re.I)


READER = AMSR2Reader
