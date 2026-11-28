"""由 FY-3D MERSI-II 的真實檔組出「合成的」FY-3A/B/C 檔 (VIRR、MERSI-1)，用來驗證沒有樣本的讀取器流程：

  python tests/fake_fy3.py <FY-3D 1000M 檔> <FY-3D GEO1K 檔> <輸出資料夾>

結構照 NSMC 的格式說明 (FY-3A/B L1 Product Format、FY-3C 扫描辐射计 / 中分辨率光谱成像仪 L1 数据格式说明)，
只做到讀取器會讀的資料集與屬性；影像內容是 MERSI-II 換算過去的，所以畫出來應該和 FY-3D 的圖一樣 (同時間、同位置)。
只能確認檔名解析、資料集位置、定標公式的串接走得通，不能驗證 FY-3A/B/C 真實檔案的格式細節。

產生的檔 (時間沿用 FY-3D 檔的時間)：
  FY3C_VIRRX_GBAL_L1_<時間>_1000M_MS.HDF + _GEOXX_MS.HDF     FY-3C VIRR (Data/、Geolocation/、Timedata/ 群組)
  FY3A_VIRRX_GBAL_L1_<時間>_1000M_MS.HDF                     FY-3A VIRR (全部在最上層，屬性拼成 Emmisive_…)
  FY3B_MERSI_GBAL_L1_<時間>_1000M_MS.HDF                     FY-3B MERSI-1 (最上層，經緯度與角度在檔內)
"""
import os
import re
import sys

import h5py
import numpy as np

C1, C2 = 1.191042e-5, 1.4387752


def _attrs(f, src, sat, name):
    for k in ('Observing Beginning Date', 'Observing Beginning Time', 'Observing Ending Date', 'Observing Ending Time',
              'EarthSun Distance Ratio', 'Day Or Night Flag'):
        if k in src.attrs:
            f.attrs[k] = src.attrs[k]
    f.attrs['Satellite Name'] = np.bytes_(sat)
    f.attrs['File Name'] = np.bytes_(name)


def _ds(f, path, data, **attrs):
    d = f.create_dataset(path, data=data, compression='gzip', compression_opts=1)
    for k, v in attrs.items():
        d.attrs[k] = v
    return d


def main(l1, geo, out):
    os.makedirs(out, exist_ok=True)
    stamp = re.search(r'_(\d{8}_\d{4})_', os.path.basename(l1)).group(1)
    with h5py.File(l1, 'r') as f, h5py.File(geo, 'r') as g:
        vis = f['Calibration/VIS_Cal_Coeff'][:].astype(np.float64)
        rows = slice(100, 1900)                                     # VIRR：1800 條掃描

        def ref_pct(b, rr=rows):                                    # MERSI-II 第 b 頻道的反射率 %
            if b <= 4:
                dn = f['Data/EV_250_Aggr.1KM_RefSB'][b - 1, rr].astype(np.float64)
            else:
                dn = f['Data/EV_1KM_RefSB'][b - 5, rr].astype(np.float64)
            c = vis[b - 1]
            p = c[0] + c[1] * dn + c[2] * dn * dn
            p[dn >= 65533] = np.nan
            return p

        def rad(b, rr=rows):                                        # MERSI-II 第 b 頻道的輻亮度
            if b <= 23:
                ds = f['Data/EV_1KM_Emissive']
                k = b - 20
            else:
                ds = f['Data/EV_250_Aggr.1KM_Emissive']
                k = b - 24
            raw = ds[k, rr].astype(np.float64)
            r = raw * float(np.ravel(ds.attrs['Slope'])[k])
            r[raw >= 65533] = np.nan
            return r

        ecw = np.asarray(f.attrs['Effect_Center_WaveLength'], dtype=np.float64)
        wn = 1e4 / ecw[[19, 23, 24]]                                # CH20、24、25 的等效中心波數
        lat, lon = g['Geolocation/Latitude'][rows], g['Geolocation/Longitude'][rows]
        ang = {k: g['Geolocation/' + k][rows] for k in ('SolarZenith', 'SolarAzimuth', 'SensorZenith', 'SensorAzimuth')}
        dem = g['Geolocation/DEM'][rows]
        day = g['Timedata/Day_Count'][:].astype(np.float64).ravel()
        ms = g['Timedata/Millisecond_Count'][:].astype(np.float64).ravel()
        t = day * 86400e3 + ms                                      # 毫秒；逐掃描帶 (10 列) -> 逐列
        t_rows = np.interp(np.arange(2000), np.arange(200) * 10 + 4.5, t)[rows]
        day_r = np.floor(t_rows / 86400e3)
        ms_r = t_rows - day_r * 86400e3

        # ---- VIRR ----
        virr_bands = (3, 4, 6, 1, 10, 2, 5)                         # VIRR 1、2、6、7、8、9、10 <- MERSI-II
        k_ref, b_ref = 0.1, -1.0
        refsb = np.stack([np.clip(np.nan_to_num((ref_pct(b) - b_ref) / k_ref, nan=0.), 1, 32000) for b in virr_bands])
        refsb = refsb.astype(np.uint16)
        emis = np.stack([np.clip(np.nan_to_num(rad(b) / 0.01, nan=0.), 0, 32000) for b in (20, 24, 25)]).astype(np.uint16)
        scales = np.full((1800, 3), 0.01, np.float32)
        offsets = np.zeros((1800, 3), np.float32)
        coef = np.array([k_ref, b_ref] * 7, np.float32)
        btc = np.array([1., 0.] * 3, np.float32)
        for sat, grouped in (('FY3C', True), ('FY3A', False)):
            name = '%s_VIRRX_GBAL_L1_%s_1000M_MS.HDF' % (sat, stamp)
            with h5py.File(os.path.join(out, name), 'w') as o:
                _attrs(o, f, 'FY-' + sat[2:], name)
                pre = 'Data/' if grouped else ''
                _ds(o, pre + 'EV_RefSB', refsb, Slope=np.ones(7, np.float32), Intercept=np.zeros(7, np.float32),
                    valid_range=np.array([0, 32767], np.int32))
                _ds(o, pre + 'EV_Emissive', emis, Slope=np.ones(3, np.float32), Intercept=np.zeros(3, np.float32),
                    valid_range=np.array([0, 32767], np.int32))
                _ds(o, pre + 'Emissive_Radiance_Scales', scales)
                _ds(o, pre + 'Emissive_Radiance_Offsets', offsets)
                o.attrs['RefSB_Cal_Coefficients'] = coef
                o.attrs['Emissive_Centroid_Wave_Number' if grouped else 'Emmisive_Centroid_Wave_Number'] = \
                    wn.astype(np.float32)
                o.attrs['Emissive_BT_Coefficients' if grouped else 'Emmisive_BT_Coefficients'] = btc
                geo_f = o
                if grouped:
                    gname = '%s_VIRRX_GBAL_L1_%s_GEOXX_MS.HDF' % (sat, stamp)
                    geo_f = h5py.File(os.path.join(out, gname), 'w')
                    _attrs(geo_f, f, 'FY-' + sat[2:], gname)
                gp = 'Geolocation/' if grouped else ''
                tp = 'Timedata/' if grouped else ''
                _ds(geo_f, gp + 'Latitude', lat)
                _ds(geo_f, gp + 'Longitude', lon)
                for k, v in ang.items():
                    _ds(geo_f, gp + k, v, Slope=np.array([0.01], np.float32), Intercept=np.array([0.], np.float32))
                _ds(geo_f, gp + 'Height', dem)
                _ds(geo_f, tp + 'Day_Count', day_r.astype(np.int32)[:, None])
                _ds(geo_f, tp + ('Millisecond_Count' if grouped else 'Msec_Count'), ms_r.astype(np.int32)[:, None])
                if grouped:
                    geo_f.close()

        # ---- MERSI-1 (FY-3B，最上層) ----
        allr = slice(0, 2000)
        name = 'FY3B_MERSI_GBAL_L1_%s_1000M_MS.HDF' % stamp
        m1 = (1, 2, 3, 4)                                            # MERSI-1 的 1~4 = MERSI-II 的 1~4
        m1_1km = (6, 7, 8, 9, 10, 11, 11, 3, 12, 14, 15, 16, 18, 17, 19)   # MERSI-1 的 6~20 (近似的 MERSI-II 頻道)
        k1 = 0.025
        with h5py.File(os.path.join(out, name), 'w') as o:
            _attrs(o, f, 'FY-3B', name)
            dn = lambda b: np.clip(np.nan_to_num(ref_pct(b, allr) / k1, nan=0.), 1, 4095).astype(np.uint16)
            _ds(o, 'EV_250_Aggr.1KM_RefSB', np.stack([dn(b) for b in m1]), Slope=np.ones(4, np.float32),
                Intercept=np.zeros(4, np.float32), valid_range=np.array([0, 4095], np.int32))
            _ds(o, 'EV_1KM_RefSB', np.stack([dn(b) for b in m1_1km]), Slope=np.ones(15, np.float32),
                Intercept=np.zeros(15, np.float32), valid_range=np.array([0, 4095], np.int32))
            # 第 5 頻道 (11.25 µm)：輻亮度 x 100 存成整數，Slope 誤植為 1 (同真實檔案)
            r5 = rad(24, allr)
            _ds(o, 'EV_250_Aggr.1KM_Emissive', np.clip(np.nan_to_num(r5 * 100, nan=0.), 0, 30000).astype(np.uint16),
                Slope=np.array([1.], np.float32), Intercept=np.array([0.], np.float32),
                valid_range=np.array([0, 4095], np.int32))
            o.attrs['VIR_Cal_Coeff'] = np.array([0., k1, 0.] * 19, np.float32)
            _ds(o, 'Latitude', g['Geolocation/Latitude'][:])
            _ds(o, 'Longitude', g['Geolocation/Longitude'][:])
            for k in ang:
                _ds(o, k, g['Geolocation/' + k][:], Slope=np.array([0.01], np.float32))
            _ds(o, 'EV_start_time', f['Calibration/EV_start_time'][:])
    print('written to', out)


if __name__ == '__main__':
    if len(sys.argv) != 4:
        sys.exit(__doc__)
    main(*sys.argv[1:])
