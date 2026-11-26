"""
SLSTR-SAFE — Sentinel-3A/B SLSTR L1B RBT，SAFE 格式
原始腳本：SatimaV3-S3-loop.py

所需檔案
  SAFE 產品 S3A_SL_1_RBT____<開始>_<結束>_..._.SEN3 (資料夾內為 xfdumanifest.xml + 各頻道 NetCDF)
  可為 .SEN3.tar、.SEN3.zip 或已解壓的 .SEN3 資料夾

頻道
  1.S1 0.555  2.S2 0.659  3.S3 0.865  4.S4 1.375  5.S5 1.61  6.S6 2.25   (a 視角 500 m，天底，反射率)
  7.S7 3.74   8.S8 10.85  9.S9 12.0                                       (i 視角 1 km，天底，亮溫)

可見光 / 短波紅外反射率 (參考 satpy slstr_l1b 與 ESA SLSTR L1 產品說明)
  ρ = π · L · adj / (F0 · cos θs)
  L    Sn_radiance_an (mW m-2 sr-1 nm-1)
  F0   Sn_quality_an.nc 的 Sn_solar_irradiance_an，依 indices_an.nc 的 detector_an 逐像元選偵測器
       (已含當日的日地距離修正：同一年 9/22 與 10/29 的 F0 差 2.1%，與日地距離的理論值 2.0% 相符)
  θs   geometry_tn.nc 的 solar_zenith_tn (16 km × 1 km tie-point 網格)，
       依 cartesian_tx / cartesian_an 的影像平面座標 (x, y) 雙線性內插到 500 m 網格 (同 satpy 作法)；
       1/cos θs 使用 satima.truecolor.sunz_corr (88° 起漸減，晨昏線附近不會爆亮)
  adj  ESA 建議的替代校正係數 (S3.PN-SLSTR-L1.08；satpy CHANCALIB_FACTORS 天底值)。
       新版產品 (quality 檔有 Sn_vic_cal_factor_an) 已在輻射值中套用，不再重複乘；
       舊版產品 (例如 processing baseline 003) 沒有此變數時才乘上

與原腳本的差異
  - 反射率由 πL/1500 的近似改為上式 (原本三個頻道都除 1500 且未除 cos θs；
    例如 θs = 40° 時 S1 約為正確值的 0.97 倍、S2 0.79 倍、S3 只有 0.49 倍)
  - 移除可見光繪圖前「反射率 < 0.01 改成最大值」的處理 (正確反射率下會把暗的海面畫成白色；
    填充值本來就已由 NetCDF 自動遮蔽)
  - 選單加入 S4–S6 (短波紅外反射率) 與 S7 (亮溫)
  - 兩種 dpi 的檔名頻道格式不一致 (B8 / B08)，統一為 B08
  - 不再解壓縮：直接從 zip / tar 讀取 NetCDF (原本每次整包解壓到目前目錄)，也接受已解壓的 .SEN3 資料夾
  - 中心時間：取最接近輸入中心點的像元，用該像元的實際觀測時間 (由掃描編號換算，GPS 時間已換成 UTC)
    (原本以檔名起訖時間依裁切範圍的中間列線性推算；錐狀掃描使掃描邊緣比列時間晚可達 1 分鐘)
"""
from datetime import datetime

import numpy as np
from scipy.interpolate import RectBivariateSpline

from satima.po import SatReader, ui, archive, swath, obstime, truecolor as tc
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr

CHANNEL_MENU = ('1.VIS(0.555) 2.VIS(0.659) 3.VIS(0.865) 4.SWIR(1.375) 5.SWIR(1.61) 6.SWIR(2.25) '
                '7.IR(3.74) 8.IR(10.85) 9.IR(12.0)')
REFL_CHANNELS = (1, 2, 3, 4, 5, 6)
BT_CHANNELS = (7, 8, 9)
# ESA S3.PN-SLSTR-L1.08 建議的天底視角輻射值調整 (僅用於沒有 vic_cal_factor 的舊版產品)
RAD_ADJUST_NADIR = {1: 0.97, 2: 0.98, 3: 0.98, 4: 1.0, 5: 1.11, 6: 1.13}


# ===== 單一檔案 =====
class SlstrFile:

    def __init__(self, fullname, channel=None, time=None):
        self.fullname = fullname
        self.channel = channel
        self.satellite = 'Sentinel-%s' % (self.fullname[1:3])
        if time:
            self.time = datetime.strptime(time, '%Y%m%dT%H%M%S')
        else:
            self.time = datetime.strptime(self.fullname[16:31], '%Y%m%dT%H%M%S')
        try:                    # 檔名第 2 個時間是結束時間：S3A_xx_1_xxx____<起>_<訖>_…
            self.end = datetime.strptime(self.fullname[32:47], '%Y%m%dT%H%M%S')
        except ValueError:
            self.end = None

    def __lt__(self, other):
        return self.time < other.time

    def _dataset(self, key):
        '''開啟 SAFE 產品中名稱包含 key 的 NetCDF (zip / tar 直接讀進記憶體，不解壓縮)。'''
        if getattr(self, 'safe', None) is None:
            self.safe = archive.SafeProduct(self.fullname)
        return self.safe.dataset(key)

    def extract(self):
        print(self.fullname)
        if self.channel in REFL_CHANNELS:
            name = 'S%s_radiance_an' % (int(self.channel))
            nc_obj = self._dataset(name)
            radiance = nc_obj.variables[name][:]
            nc_obj.close()
            data = self.rad_to_refl(radiance, int(self.channel))
        else:
            name = 'S%s_BT_in' % (int(self.channel))
            nc_obj = self._dataset(name)
            data = nc_obj.variables[name][:] - 273.15
            nc_obj.close()
        print(name, 'max', np.nanmax(data))
        return data

    def rad_to_refl(self, radiance, channel):
        """輻射 (a 視角天底) -> 反射率 πL·adj / (F0·cosθs)。"""
        q = self._dataset('S%d_quality_an' % channel)
        f0 = np.ma.filled(q.variables['S%d_solar_irradiance_an' % channel][:].astype(np.float64), np.nan)
        has_vic = ('S%d_vic_cal_factor_an' % channel) in q.variables
        vic = float(q.variables['S%d_vic_cal_factor_an' % channel][0]) if has_vic else None
        q.close()
        ind = self._dataset('indices_an')
        det = np.ma.filled(ind.variables['detector_an'][:], -1).astype(np.int32)
        ind.close()

        adj = 1.0 if has_vic else RAD_ADJUST_NADIR[channel]
        print('S%d F0 %s mW m-2 nm-1, %s' % (channel, np.round(f0, 2),
              tr('已套用替代校正係數 %.4f', 'vicarious factor %.4f already applied') % vic if has_vic
              else tr('套用 ESA 建議的調整係數 %.2f', 'apply ESA adjustment %.2f') % adj))
        f0_px = np.where(det >= 0, f0[np.clip(det, 0, f0.size - 1)], np.nan)
        sza = self.solar_zenith(radiance.shape)
        rad = swath.as_float(radiance)
        with np.errstate(invalid='ignore', divide='ignore'):
            refl = np.pi * rad * adj / f0_px * tc.sunz_corr(sza)
        return np.ma.masked_invalid(refl)

    def solar_zenith(self, shape):
        """太陽天頂角：tie-point 網格 (tx) 依影像平面座標雙線性內插到 a 視角 500 m 網格。"""
        geom = self._dataset('geometry_tn')
        sza = swath.fill_invalid(swath.as_float(geom.variables['solar_zenith_tn'][:]))
        geom.close()
        tx = self._dataset('cartesian_tx')
        x_tx = swath.as_float(tx.variables['x_tx'][0, :])
        y_tx = swath.as_float(tx.variables['y_tx'][:, 0])
        tx.close()
        an = self._dataset('cartesian_an')
        x_an = swath.as_float(an.variables['x_an'][:])
        y_an = swath.as_float(an.variables['y_an'][:])
        an.close()
        if x_tx[0] > x_tx[-1]:          # x 沿行遞減 -> 反轉成遞增
            x_tx, sza = x_tx[::-1], sza[:, ::-1]
        spline = RectBivariateSpline(y_tx, x_tx, sza, kx=1, ky=1)
        out = np.full(shape, np.nan)
        ok = np.isfinite(x_an) & np.isfinite(y_an)
        out[ok] = spline.ev(y_an[ok], x_an[ok])
        return out

    def geocoord(self):
        if self.channel in REFL_CHANNELS:
            nc_obj = self._dataset('geodetic_an.nc')
            lats = nc_obj.variables['latitude_an'][:]
            lons = nc_obj.variables['longitude_an'][:]
        else:
            nc_obj = self._dataset('geodetic_in.nc')
            lats = nc_obj.variables['latitude_in'][:]
            lons = nc_obj.variables['longitude_in'][:]
        nc_obj.close()
        print(np.amax(lats), np.amin(lats), np.amax(lons), np.amin(lons))
        return lats, lons

    def plot(self, canvas, lats, lons, data, variant):
        if data.size == 0:
            return ''
        if self.channel in REFL_CHANNELS:
            canvas.mesh(lons, lats, data, cmap='gray', vmin=0.0, vmax=1.0)
            return ''
        return canvas.bt_mesh(lons, lats, data, variant)

    def row_times(self):
        """每個像元的觀測時間 (UTC，與 geocoord 同大小的 2-D 陣列)，依 indices 檔的掃描編號換算：
            t = Nadir_Minimal_ts + (scan − Nadir_First_scan) × SCANSYNC + (NADIR_FIRST_PIXEL + pixel) × PIXSYNC
        SLSTR 為錐狀掃描，同一列的觀測時間並不相同：與 time_stamp (星下點通過該列的時間) 相差
        約 −12 s (軌道中線) 到 +60 s (掃描邊緣)。Sentinel-3 的時間欄位為 GPS 時間 (比 UTC 快 18 s)。
        產品沒有掃描編號時退回逐列的 time_stamp。"""
        g = 'a' if self.channel in REFL_CHANNELS else 'i'
        grid = g + 'n'
        tm = self._dataset('time_%s' % grid)
        try:
            ind = self._dataset('indices_%s' % grid)
            scan = np.ma.filled(ind.variables['scan_%s' % grid][:].astype(np.float64), np.nan)
            pix = np.ma.filled(ind.variables['pixel_%s' % grid][:].astype(np.float64), np.nan)
            ind.close()
            first = np.ma.filled(tm.variables['Nadir_First_scan_%s' % g][:].astype(np.float64), np.nan)
            ts0 = np.ma.filled(tm.variables['Nadir_Minimal_ts_%s' % g][:].astype(np.float64), np.nan)
            scansync = float(tm.variables['SCANSYNC'][0])
            pixsync = float(tm.variables['PIXSYNC_%s' % g][0])
            first_pix = float(tm.variables['NADIR_FIRST_PIXEL_%s' % g][0])
        except (KeyError, FileNotFoundError):
            t = obstime.gps2000(tm.variables['time_stamp_%s' % g][:])
            tm.close()
            return t
        tm.close()
        ref = np.nanmedian(ts0 - first * scansync)          # 掃描編號 0 的時間 (掃描週期固定)
        return obstime.gps2000(ref + scan * scansync + (first_pix + pix) * pixsync)


# ===== 讀取器 =====
class SlstrReader(SatReader):
    NAME = 'SLSTR-SAFE'
    SENSOR = 'SLSTR'
    DESCRIPTION = 'Sentinel-3 SLSTR L1B RBT (SAFE .SEN3 tar/zip/資料夾)'
    DESCRIPTION_EN = 'Sentinel-3 SLSTR L1B RBT (SAFE .SEN3 tar/zip/folder)'
    SOURCE = 'Copernicus'
    PATTERNS = ('S*_SL*tar', 'S*_SL*zip', 'S*_SL*.SEN3')
    MULTI_FILE = True

    @classmethod
    def find(cls, folder=None):
        return archive.prefer_archives(super().find(folder))

    def make_file(self, name):
        # S3A_SL_1_RBT____20181001T211249_20181001T211549_20181003T063018_0179_036_214_2700_MAR_O_NT_003.SEN3.tar
        return SlstrFile(name, channel=None, time=name[16:31])

    def listing(self, f):
        return (f.satellite, 'SLSTR', f.time)

    def select_band(self):
        self.channel = ui.ask_int(tr('輸入頻道編號 (%s)：', '%s : ') % CHANNEL_MENU, REFL_CHANNELS + BT_CHANNELS)
        for f in self.files:
            f.channel = self.channel

    def select_variants(self):
        if self.channel in BT_CHANNELS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare(self):
        datas, lats, lons = [], [], []
        for f in self.files:
            datas.append(f.extract())
            la, lo = f.geocoord()
            lats.append(la)
            lons.append(lo)
        data = self.concat(datas)
        lats = self.concat(lats)
        lons = self.concat(lons)
        rowtime = np.concatenate([f.row_times() for f in self.files], axis=0)
        if len(self.files) > 1:
            print(tr('已合併 %d 個片段：', 'Merged %d granules:') % len(self.files), np.shape(lats))
        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        self.data = self.crop(data)
        self.time = self.center_time(rowtime)
        if self.channel in BT_CHANNELS:
            self.report_eyetemp(self.data)
        return True

    def slice(self, georange, lats, lons):
        return self.slice_swath(lats, lons)

    def plot(self, canvas, variant):
        cm = self.files[0].plot(canvas, self.lats, self.lons, self.data, variant)
        return dict(band='BAND-%02d' % self.channel, colorbar=cm, label='BAND-%02d' % self.channel)


READER = SlstrReader
