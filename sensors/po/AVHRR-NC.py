"""
AVHRR-NC — Metop-A/B/C AVHRR/3 L1B (EUMETSAT NetCDF)
原始腳本：SatimaV3METOP-loop.py

所需檔案
  W_XX-EUMETSAT-Darmstadt,HIRES+RADIOMETER,METOPA+AVHR_C_EUMP_<時間>_<軌道>_eps_o_l1b.nc

頻道
  1.CH1 0.630  2.CH2 0.862  3.CH3A  (反射率 πL/F，同原程式)
  4.CH3B  5.CH4 10.8  6.CH5 12.0    (亮溫)

與原腳本的差異
  - NetCDF 只開一次並重複使用 (原本 extract / geocoord / timeprocess 各開一次且未關閉)
  - 中心時間：取最接近輸入中心點的像元，用該掃描線的 record_start_time
    (原本以檔案起訖時間依裁切範圍的中間列線性推算)
  - 亮溫的頻帶修正為 T = beta + alpha·T* (Product Guide §6.3.4 的 T = A + B·T*)，原本寫成
    (T* - beta) / alpha，亮溫偏高 (CH3B 約 1.6-2.4 K，CH4 約 0.1-0.6 K，CH5 約 0-0.3 K)
"""
from datetime import datetime

import netCDF4 as nc
import numpy as np

from satima.po import SatReader, ui, obstime
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr

BAND_MENU = '1.BAND-1(0.630) 2.BAND-2(0.862) 3.BAND-3A 4.BAND-3B 5.BAND-4(10.8) 6.BAND5(12)'
BAND_LABEL = {1: '01', 2: '02', 3: '03A', 4: '03B', 5: '04', 6: '05'}
BT_BANDS = (4, 5, 6)


# ===== 單一檔案 =====
class AvhrrNcFile:

    def __init__(self, fullname, channel=None, time=None):
        self.fullname = fullname
        self.satellite = 'METOP-%s' % (fullname[46:47])
        self.channel = 'IR/VIS'
        if time:
            self.time = datetime.strptime(time, '%Y%m%d%H%M%S')
        else:
            self.time = datetime.strptime(self.fullname[60:74], '%Y%m%d%H%M%S')

    def __lt__(self, other):
        return self.time < other.time

    def dataset(self):
        if not hasattr(self, '_nc'):
            self._nc = nc.Dataset(self.fullname)
        return self._nc

    def extract(self, band):
        nc_obj = self.dataset()
        if band == 1 or band == 2:
            raw = nc_obj.variables['scene_radiances%s' % (band)][:]
            out = self.VISConvert(nc_obj, raw, band)
        elif band == 3:
            raw = nc_obj.variables['scene_radiances3a'][:]
            out = self.VISConvert(nc_obj, raw, band)
        elif band == 4:
            raw = nc_obj.variables['scene_radiances3b'][:]
            out = self.BTConvert(nc_obj, raw, band)
        elif band == 5 or band == 6:
            raw = nc_obj.variables['scene_radiances%s' % (band - 1)][:]
            out = self.BTConvert(nc_obj, raw, band)
        return out

    def BTConvert(self, nc_obj, R, band):
        C1 = nc_obj.c1
        C2 = nc_obj.c2
        if band == 4:
            wnc, beta, alpha = nc_obj.channel_3b_wnc, nc_obj.channel_3b_beta, nc_obj.channel_3b_alpha
        if band == 5:
            wnc, beta, alpha = nc_obj.channel_4_wnc, nc_obj.channel_4_beta, nc_obj.channel_4_alpha
        if band == 6:
            wnc, beta, alpha = nc_obj.channel_5_wnc, nc_obj.channel_5_beta, nc_obj.channel_5_alpha
        print('C1 %s  C2 %s' % (C1, C2))
        print('wnc %s  beta %s  alpha %s' % (wnc, beta, alpha))
        # beta (K) = EPS 的 CONSTANT1 (A)，alpha (K/K) = CONSTANT2_SLOPE (B)
        # AVHRR L1B Product Guide §6.3.4：T* = C2·γ / ln(1 + C1·γ³/R)，T = A + B·T*
        # (原本寫成 (T* - beta) / alpha，CH3B 高估約 1.6-2.4 K、CH4 0.1-0.6 K、CH5 0-0.3 K)
        with np.errstate(invalid='ignore', divide='ignore'):
            bt = beta + alpha * (C2 * wnc) / np.log(1. + (C1 * wnc ** 3) / R) - 273.15
        return bt

    def VISConvert(self, nc_obj, R, band):
        if band == 1:
            f_sol = nc_obj.channel_1_f_sol
        if band == 2:
            f_sol = nc_obj.channel_2_f_sol
        if band == 3:
            f_sol = nc_obj.channel_3a_f_sol
        refl = (R * np.pi) / f_sol
        return refl

    def geocoord(self, band):
        nc_obj = self.dataset()
        lats = nc_obj.variables['lat'][:]
        lons = nc_obj.variables['lon'][:]
        return lats, lons

    def plot(self, canvas, lats, lons, data, variant, band):
        if data.size == 0:
            return ''
        if band in (1, 2, 3):
            canvas.mesh(lons, lats, data, cmap='gray', vmin=0., vmax=1.)
            return ''
        return canvas.bt_mesh(lons, lats, data, variant)

    def row_times(self):
        """每條掃描線的觀測時間 (UTC)：record_start_time (EPS MDR 的 GRH 紀錄時間，2000-01-01 起的秒數)。
        檔案只存到整秒 (每秒 6 條掃描線)，誤差 < 1 s。"""
        var = self.dataset().variables['record_start_time']
        var.set_auto_scale(False)
        return obstime.from_utc('2000-01-01', var[:], 's')


# ===== 讀取器 =====
class AvhrrNcReader(SatReader):
    NAME = 'AVHRR-NC'
    SENSOR = 'AVHRR3'
    DESCRIPTION = 'Metop AVHRR/3 L1B (EUMETSAT NetCDF)'
    SOURCE = 'EUMETSAT'
    PATTERNS = ('W_XX-EUMETSAT-Darmstadt,HIRES+RADIOMETER,METOP*.nc',)
    MULTI_FILE = False

    def make_file(self, name):
        # W_XX-EUMETSAT-Darmstadt,HIRES+RADIOMETER,METOPA+AVHR_C_EUMP_20190901140403_66771_eps_o_l1b.nc
        return AvhrrNcFile(name, time=name[60:74])

    def listing(self, f):
        return (f.satellite, f.channel, f.time)

    def select_band(self):
        self.file = self.files[0]
        print(self.file.fullname)
        print(BAND_MENU)
        self.band = ui.ask_int(tr('輸入頻道編號：', 'Input the band: '), range(1, 7))

    def select_variants(self):
        if self.band in BT_BANDS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    def prepare(self):
        f = self.file
        data = f.extract(self.band)
        lats, lons = f.geocoord(self.band)
        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        self.data = self.crop(data)
        self.time = self.center_time(f.row_times())
        if self.band in BT_BANDS:
            self.report_eyetemp(self.data)
        return True

    def slice(self, georange, lats, lons):
        return self.slice_swath(lats, lons)

    def plot(self, canvas, variant):
        cm = self.file.plot(canvas, self.lats, self.lons, self.data, variant, self.band)
        channel = BAND_LABEL[self.band]
        return dict(band='B%s' % channel, colorbar=cm, label='BAND%s' % channel)


READER = AvhrrNcReader
