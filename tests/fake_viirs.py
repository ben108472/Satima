"""產生「合成的」VIIRS SDR 測試檔 (不是真的觀測資料)，用來驗證 VIIRS 讀取器的流程：

  python tests/fake_viirs.py <輸出資料夾>

手邊沒有真正的 VIIRS 檔時，至少能確認檔名解析、片段接續、合併檔、定位檔的退路這幾條路徑走得通。
結構只做到讀取器會讀的資料集 (All_Data/VIIRS-<頻道>-SDR_All、VIIRS-MOD-GEO[-TC]_All、VIIRS-DNB-GEO_All)，
數值是人造的平滑場；不能用來驗證校正或定位的正確性。

產生的片段 (S-NPP 2023/05/26，每個 85.4 秒，共 5 個時間，第 4 個故意缺)：
  1  SVM15 + GMTCO 分開的檔；SVDNB + GDNBO (含 Latitude_TC)
  2  GMTCO-SVM15 合併檔 (定位與資料在同一個檔)；GDNBO-SVDNB 合併檔
  3  SVM15 + GMODO (只有未做地形校正的定位)；另有一個較舊的重複檔 (產製時間較早，應被忽略)
  4  (缺)
  5  SVM15 + GMTCO
"""
import datetime
import os
import sys

import h5py
import numpy as np

NSCAN, PER_SCAN, NCOL = 12, 16, 400         # 每個片段 12 條掃描 x 16 列，每列 400 點 (真的是 48 x 16 x 3200)
GRAN = 85.4
T0 = datetime.datetime(2023, 5, 26, 4, 3, 20, 400000)
IET0 = datetime.datetime(1958, 1, 1)
LEAP = 37                                   # 2017 年以後 TAI - UTC


def _stamp(t):
    return '%s%d' % (t.strftime('%H%M%S'), t.microsecond // 100000)


def name(prods, k, created='20230526050000000000'):
    start = T0 + datetime.timedelta(seconds=GRAN * k)
    end = start + datetime.timedelta(seconds=GRAN - 0.2)
    return '%s_npp_d%s_t%s_e%s_b59876_c%s_noac_ops.h5' % ('-'.join(prods), start.strftime('%Y%m%d'), _stamp(start),
                                                          _stamp(end), created)


def geometry(k, shift=0.):
    """第 k 個片段的 (緯度, 經度)：由南往北，掃描帶中心在 122E；shift 模擬地形校正造成的位移 (度)。"""
    ny = NSCAN * PER_SCAN
    rows = (k * ny + np.arange(ny))[:, None]
    cols = np.arange(NCOL)[None, :]
    lat = 15. + rows * 0.03 + (cols - NCOL / 2.) * 0.0015 + shift
    lon = 122. + (cols - NCOL / 2.) * 0.035 - rows * 0.003 + shift
    return lat.astype(np.float32), lon.astype(np.float32)


def write_geo(f, group, k, dnb=False, tc=False):
    g = f.require_group('All_Data/' + group)
    lat, lon = geometry(k)
    g.create_dataset('Latitude', data=lat, fillvalue=-999.3)
    g.create_dataset('Longitude', data=lon, fillvalue=-999.3)
    if tc:
        lat2, lon2 = geometry(k, shift=0.02)
        g.create_dataset('Latitude_TC', data=lat2, fillvalue=-999.3)
        g.create_dataset('Longitude_TC', data=lon2, fillvalue=-999.3)
    for key, value in (('SolarZenithAngle', 120. if dnb else 40.), ('SolarAzimuthAngle', 100.),
                       ('SatelliteZenithAngle', 30.), ('SatelliteAzimuthAngle', 80.), ('Height', 0.)):
        g.create_dataset(key, data=np.full(lat.shape, value, np.float32), fillvalue=-999.3)
    start = T0 + datetime.timedelta(seconds=GRAN * k)
    sec = (start - IET0).total_seconds() + LEAP + (np.arange(NSCAN) + 0.5) * GRAN / NSCAN
    g.create_dataset('MidTime', data=(sec * 1e6).astype(np.int64), fillvalue=-993)
    if dnb:
        g.create_dataset('LunarZenithAngle', data=np.full(lat.shape, 50., np.float32), fillvalue=-999.3)
        g.create_dataset('MoonIllumFraction', data=np.array([80.], np.float32))
        g.create_dataset('MoonPhaseAngle', data=np.array([50.], np.float32))


def write_m15(f, k):
    g = f.require_group('All_Data/VIIRS-M15-SDR_All')
    lat, lon = geometry(k)
    bt = 290. - 70. * np.exp(-((lat - 17.5) ** 2 + (lon - 122.) ** 2) / 1.5) + 3. * np.sin(lon * 4.) * np.cos(lat * 5.)
    scale, offset = 0.0025, 180.
    g.create_dataset('BrightnessTemperature', data=np.round((bt - offset) / scale).astype(np.uint16), fillvalue=65535)
    g.create_dataset('BrightnessTemperatureFactors', data=np.array([scale, offset], np.float32))


def write_dnb(f, k):
    g = f.require_group('All_Data/VIIRS-DNB-SDR_All')
    lat, lon = geometry(k)
    rad = 2e-9 + 3e-8 * np.exp(-((lat - 17.5) ** 2 + (lon - 122.) ** 2) / 1.5)
    g.create_dataset('Radiance', data=rad.astype(np.float32), fillvalue=-999.3)


def main(out):
    os.makedirs(out, exist_ok=True)

    def make(prods, k, writers, **kw):
        path = os.path.join(out, name(prods, k, **kw))
        with h5py.File(path, 'w') as f:
            for w in writers:
                w(f)

    geo = lambda group, k, **kw: (lambda f: write_geo(f, group, k, **kw))
    m15 = lambda k: (lambda f: write_m15(f, k))
    dnb = lambda k: (lambda f: write_dnb(f, k))
    # 1：分開的檔
    make(['SVM15'], 0, [m15(0)])
    make(['GMTCO'], 0, [geo('VIIRS-MOD-GEO-TC_All', 0)])
    make(['SVDNB'], 0, [dnb(0)])
    make(['GDNBO'], 0, [geo('VIIRS-DNB-GEO_All', 0, dnb=True, tc=True)])
    # 2：合併檔
    make(['GMTCO', 'SVM15'], 1, [geo('VIIRS-MOD-GEO-TC_All', 1), m15(1)])
    make(['GDNBO', 'SVDNB'], 1, [geo('VIIRS-DNB-GEO_All', 1, dnb=True, tc=True), dnb(1)])
    # 3：只有未做地形校正的定位；外加一個較舊的重複檔
    make(['SVM15'], 2, [m15(2)])
    make(['SVM15'], 2, [m15(2)], created='20230526040000000000')
    make(['GMODO'], 2, [geo('VIIRS-MOD-GEO_All', 2)])
    # 4 缺；5：分開的檔
    make(['SVM15'], 4, [m15(4)])
    make(['GMTCO'], 4, [geo('VIIRS-MOD-GEO-TC_All', 4)])
    print('%d files in %s' % (len(os.listdir(out)), out))


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'fake_viirs')
