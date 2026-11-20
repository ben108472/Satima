"""NASA PPS 的 GPM 星座 Level 1C (XCAL 交互校正亮溫 Tc)

GPM 星座所有輻射計的 L1C 都是同一種格式 (HDF5 / netCDF-4)，所以這一支讀取器涵蓋 PPS 上所有的 1C 產品，
不需逐儀器寫頻道表：頻道由各 swath 的 Tc long_name 解析 (見 pps_common.py)。

PPS 封存 (arthurhouhttps.pps.eosdis.nasa.gov/gpmdata/年/月/日/1C/) 上的產品 (2026/11 盤點 1987~2026)：
  圓錐掃描  SSM/I (F08、F10、F11、F13、F14、F15)、SSMIS (F16~F19)、TMI (TRMM)、AMSR-E (Aqua)、
            AMSR2 (GCOM-W1)、AMSR3 (GOSAT-GW)、GMI (GPM)、MWI (WSF-M1)
  橫向掃描  AMSU-B (NOAA-15/16/17)、MHS (NOAA-18/19、MetOp-A/B/C)、ATMS (S-NPP、NOAA-20/21)、
            SAPHIR (Megha-Tropiques)；沒有雙極化，只有單頻道亮溫
  1C-R      GMI 的另一個版本：高頻 (S2：166、183 GHz) 重新對位到低頻 (S1) 的取樣位置，S1 與 1C 相同。
            列表的格式欄顯示 L1C-R

檔名  1C.<衛星>.<儀器>.<校正版本>.<YYYYMMDD>-S<HHMMSS>-E<HHMMSS>.<軌道>.<版本>.nc (或 .HDF5)
      即時產品沒有軌道編號，副檔名為 .RT-NC / .RT-H5
品質  依 L1C ATBD (V08)：Quality 0 = 正常；正值 = 警告 (太陽耀光、RFI…)，亮溫保留；負值 = 嚴重錯誤，
      出錯頻道的 Tc 已設成無資料 (-9999.9)，所以讀取時只要把無資料的 Tc 去掉即可。
      Quality 是同一個 swath 的各頻道共用的，許多負值其實只指其中一個頻道 (-4 一個頻道缺資料、
      SSMIS -125 150H 故障、AMSU-B -100~-104 某一頻道不可用…)，其他頻道的 Tc 仍然有效：
      [swath] mask_bad_quality = true 會把這些取樣點整點丟掉 (F18 SSMIS 的 183 GHz、NOAA-15 AMSU-B 2008 年的
      89 / 150 GHz 會整個不見)，所以預設為 false。
"""
import re

from pps_common import PPSReader


class L1CReader(PPSReader):
    NAME = 'PPS-L1C'
    SENSOR = 'GPM constellation'
    DESCRIPTION = 'NASA PPS Level 1C 交互校正亮溫 (SSM/I、TMI、AMSR-E、SSMIS、AMSR2/3、GMI、MWI、AMSU-B、MHS、ATMS、SAPHIR)'
    DESCRIPTION_EN = ('NASA PPS Level 1C intercalibrated brightness temperature (SSM/I, TMI, AMSR-E, SSMIS, AMSR2/3, '
                      'GMI, MWI, AMSU-B, MHS, ATMS, SAPHIR)')
    PATTERNS = ('1C*.nc', '1C*.nc4', '1C*.HDF5', '1C*.RT-NC', '1C*.RT-H5')
    VAR = 'Tc'
    FORMAT = 'L1C'
    FILENAME = re.compile(r'^1C(?P<variant>-R)?\.(?P<sat>[^.]+)\.(?P<sensor>[^.]+)\.[^.]+\.(?P<date>\d{8})-S(?P<s>\d{6})'
                          r'-E(?P<e>\d{6})(?:\.(?P<orbit>\d+))?\.', re.I)


READER = L1CReader
