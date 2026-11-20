"""把分開的 VIIRS 定位檔與 SDR 檔併成 NOAA CLASS 那種「合併檔」(檔名前段列出所含的產品)，測試合併檔的讀取用。

  python tests/pack_viirs.py <來源資料夾> <輸出資料夾> <定位產品> <SDR 產品> [衛星_d日期 的篩選字串]
  例：python tests/pack_viirs.py D:/VIIRS D:/tmp/packed GITCO SVI05 j01_d20261027

每個片段產生一個 <定位>-<SDR>_<衛星>_d…_ops.h5，內容是兩個來源檔的 All_Data、Data_Products 群組原樣複製
(HDF5 內部的群組名稱與真正的合併檔相同)。這是用真實資料「組」出來的檔，不是 CLASS 下載的原檔。
"""
import glob
import os
import sys

import h5py

src, dst, geo, sdr = sys.argv[1:5]
match = sys.argv[5] if len(sys.argv) > 5 else ''
os.makedirs(dst, exist_ok=True)
n = 0
for gpath in sorted(glob.glob(os.path.join(src, geo + '_*%s*_ops.h5' % match))):
    key = os.path.basename(gpath).split('_c')[0][len(geo):]         # _j01_d…_t…_e…_b…
    spaths = glob.glob(os.path.join(src, sdr + key + '_c*_ops.h5'))
    if not spaths:
        continue
    tail = os.path.basename(spaths[0])[len(sdr):]
    out = os.path.join(dst, '%s-%s%s' % (geo, sdr, tail))
    with h5py.File(out, 'w') as fo:
        for path in (gpath, spaths[0]):
            with h5py.File(path, 'r') as fi:
                for top in fi:
                    g = fo.require_group(top)
                    for name in fi[top]:
                        if name not in g:
                            fi.copy('%s/%s' % (top, name), g)
    n += 1
    print(os.path.basename(out))
print('%d packaged files' % n)
