"""列出資料夾裡 VIIRS 定位檔的涵蓋範圍與所含的資料集 (挑測試用的地圖中心、確認檔案內容用)。

  python tests/viirs_peek.py <資料夾>
"""
import glob
import os
import sys

import h5py
import numpy as np

folder = sys.argv[1] if len(sys.argv) > 1 else '.'
for path in sorted(glob.glob(os.path.join(folder, 'G*_ops.h5'))):
    name = os.path.basename(path)
    with h5py.File(path, 'r') as f:
        for group in f['All_Data']:
            g = f['All_Data'][group]
            lat = np.asarray(g['Latitude'][::32, ::32], dtype=np.float64)
            lon = np.asarray(g['Longitude'][::32, ::32], dtype=np.float64)
            ok = (np.abs(lat) <= 90) & (np.abs(lon) <= 180)
            mid = lat.shape[0] // 2, lat.shape[1] // 2
            extra = [k for k in ('Latitude_TC', 'Longitude_TC', 'Height', 'MoonIllumFraction') if k in g]
            print('%-48s %-22s %s  lat %6.2f~%6.2f lon %7.2f~%7.2f  中心 %.2f %.2f  %s' % (
                name[:48], group, g['Latitude'].shape, lat[ok].min(), lat[ok].max(), lon[ok].min(), lon[ok].max(),
                lat[mid], lon[mid], ' '.join(extra)))
