"""GPM GMI Level 1B 亮溫 (JAXA G-Portal 的 GPMCOR_GMI_*_1BS_*；NASA PPS 的 1B.GPM.GMI.* 為同一份檔案)

檔名  JAXA：GPMCOR_GMI_<YYMMDDhhmm>_<hhmm>_<軌道>_1BS_G1B_<版本>.h5 (或 .nc)
      NASA：1B.GPM.GMI.TB2021.<YYYYMMDD>-S<HHMMSS>-E<HHMMSS>.<軌道>.<版本>.HDF5
結構  同 PPS 的 swath 檔 (pps_common.py)：S1 (10.65~89 GHz，9 個頻道)、S2 (166、183 GHz，4 個頻道)，
      亮溫為 Tb (未經 XCAL 交互校正)。頻道順序固定，寫在下面的 TABLE。
"""
import datetime
import os
import re

from pps_common import PPSReader, granule_times
from satima.mw.reader import Scene

_JAXA = re.compile(r'^GPMCOR_GMI_(?P<start>\d{10})_(?P<end>\d{4})_(?P<orbit>\d+)_1B', re.I)
_NASA = re.compile(r'^1B\.GPM\.GMI\.[^.]+\.(?P<date>\d{8})-S(?P<s>\d{6})-E(?P<e>\d{6})(?:\.(?P<orbit>\d+))?\.', re.I)


class GMIReader(PPSReader):
    NAME = 'GMI-L1B'
    SENSOR = 'GMI'
    DESCRIPTION = 'GPM GMI Level 1B 亮溫 (JAXA G-Portal / NASA PPS)'
    DESCRIPTION_EN = 'GPM GMI Level 1B brightness temperature (JAXA G-Portal / NASA PPS)'
    PATTERNS = ('GPMCOR_GMI_*_1B*.h5', 'GPMCOR_GMI_*_1B*.nc', '1B.GPM.GMI.*.HDF5', '1B.GPM.GMI.*.nc')
    VAR = 'Tb'
    FORMAT = 'L1B'
    TABLE = {
        'S1': [(10.65, None, 'V', ''), (10.65, None, 'H', ''), (18.7, None, 'V', ''), (18.7, None, 'H', ''),
               (23.8, None, 'V', ''), (36.64, None, 'V', ''), (36.64, None, 'H', ''), (89.0, None, 'V', ''),
               (89.0, None, 'H', '')],
        'S2': [(166.0, None, 'V', ''), (166.0, None, 'H', ''), (183.31, 3., 'V', ''), (183.31, 7., 'V', '')],
    }

    def scenes(self, folder='.'):
        out = []
        for name in self.find(folder):
            base = os.path.basename(name)
            m = _JAXA.match(base)
            source = 'JAXA' if m else 'PPS'
            if m:
                start = datetime.datetime.strptime(m['start'], '%y%m%d%H%M')
                end = datetime.datetime.strptime(m['start'][:6] + m['end'], '%y%m%d%H%M')
                if end < start:
                    end += datetime.timedelta(days=1)
            else:
                m = _NASA.match(base)
                if not m:
                    continue
                start, end = granule_times(m['date'], m['s'], m['e'])
            out.append(Scene('GPM', 'GMI', start, end, file=os.path.join(folder, name), orbit=m['orbit'] or '',
                             format=self.FORMAT, extra={'source': source}))
        return out


READER = GMIReader
