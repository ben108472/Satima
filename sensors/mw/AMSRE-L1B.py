"""Aqua AMSR-E Level 1B 亮溫 —— JAXA 重新處理成「AMSR2 格式」的版本 (G-Portal，HDF5)

檔名  PM1AME_<YYYYMMDDhhmm>_<路徑><A/D>_L1SGBTBR_<版本>.h5   (A = 升軌、D = 降軌，一個檔約半圈)
頻道  6.925、10.65、18.7、23.8、36.5、89.0 GHz 的 V / H (沒有 7.3 GHz)
      89 GHz 的 A 掃描在 2004/11 故障，之後的檔案 89A 全部是無效值，只有 B 掃描 (沿軌 10 km)；
      合併的 89V / 89H 自動只用有資料的那一組。
結構與低頻經緯度的推算同 AMSR2，見 amsr_common.py。舊的 HDF4 格式 (P1AME*.00) 不支援。
"""
import re

from amsr_common import AMSRReader


class AMSREReader(AMSRReader):
    NAME = 'AMSRE-L1B'
    SENSOR = 'AMSR-E'
    SATELLITE = 'Aqua'
    DESCRIPTION = 'Aqua AMSR-E Level 1B 亮溫，AMSR2 格式 (JAXA G-Portal)'
    DESCRIPTION_EN = 'Aqua AMSR-E Level 1B brightness temperature, AMSR2 format (JAXA G-Portal)'
    PATTERNS = ('PM1AME_*_L1SGBTBR_*.h5',)
    FILENAME = re.compile(r'^PM1AME_(\d{12})_(\d{3})([AD])_L1SGBTBR_', re.I)


READER = AMSREReader
