"""
MODIS-HDF4 — Terra / Aqua MODIS L1B (HDF4)
原始腳本：SatimaV3MODIS-loop.py

所需檔案 (同一目錄)
  MOD021KM / MYD021KM   1 km (紅外線)
  MOD02HKM / MYD02HKM   500 m (可見光 B1–B7)
  MOD02QKM / MYD02QKM   250 m (可見光 B1–B2)
  MOD03 / MYD03         地理定位 (經緯度、角度、地形高度、掃描時間)

產品
  紅外線 (任一熱紅外波段，B27/B28 用水氣色表，其餘用紅外線色表) / 可見光單波段 /
  真實色 1.RGB (原始合成) 2.CORR (sunz + CREFL + 銳化 + CIRA) 3.WV (NASA Worldview)

與原腳本的差異
  - 紅外線不限 B27/B31 (亮溫係數本來就涵蓋所有熱紅外波段)
  - CREFL 中與其他衛星相同的部分 (chand、csalbr、sunz、增強曲線) 改呼叫 satima/truecolor.py
"""
import datetime
from glob import glob

import numpy as np
from pyhdf import SD
from scipy.interpolate import RectBivariateSpline

from satima.po import SatReader, ui, obstime, truecolor as tc
from satima.core.colormaps import ask_colormaps
from satima.core.i18n import tr

EMISSIVE_BANDS = [20, 21, 22, 23, 24, 25, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36]
WV_BANDS = (27, 28)
PRODUCT_RES = {'IR 1KM': '1km', 'VIS 500m': '500m', 'VIS 250m': '250m'}


class ModisGeo:
    """MODIS 地理定位資料與 1km -> 250m / 500m 經緯度插值 (純 numpy)。

    演算法對應 python-geotiepoints 1.9 (satpy modis_l1b reader) 的 modisinterpolator：
      1. 每個掃描帶 (10 條 1km 線) 各自插值，不跨掃描帶 (處理 bow-tie 重疊)
      2. 於地心直角座標 (x, y, z) 內插，避免換日線與極區問題
      3. 依衛星天頂角計算沿掃描方向的 expansion 修正 (掃描角非線性)
      4. 子像素偏移：細網格第 r 列對應 1km 列 (r - (f-1)/2) / f，第 c 行對應 c / f
    舊版 geotiepoints 另有 alignment 修正項；以 MYD021KM 5km -> 1km 對 MYD03 真值驗證，
    加入該項誤差反而變大 (平均 7.1 m vs 5.8 m)，且 1.9 版輸出亦不含此項，故不採用。

    用法：
        geo = ModisGeo.from_file('MYD03.A2006302.0500...hdf')
        lats250, lons250 = geo.to_fine(4)     # 250m
        lats500, lons500 = geo.to_fine(2)     # 500m
    """

    EARTH_R = 6371.0      # km
    SAT_H = 705.0         # km, Terra/Aqua 軌道高度
    SCAN_LEN_1KM = 10     # 每個掃描帶的 1km 線數

    def __init__(self, lats, lons, satz, scan_len=SCAN_LEN_1KM):
        """lats, lons, satz: (nlines, ncols) 經度、緯度、衛星天頂角 (度)。
        scan_len: 每個掃描帶的線數 (1km 為 10；5km 為 2)。"""
        if lats.shape[0] % scan_len:
            raise ValueError(tr('線數 %d 不是 %d 的整數倍', 'Number of lines %d is not a multiple of %d')
                             % (lats.shape[0], scan_len))
        self.lats = np.asarray(lats)
        self.lons = np.asarray(lons)
        self.satz = np.asarray(satz)
        self.scan_len = scan_len

    @classmethod
    def from_file(cls, geofile):
        """由 MOD03/MYD03 讀入 1km Latitude、Longitude、SensorZenith。"""
        hdf = SD.SD(geofile, SD.SDC.READ)
        try:
            lats = hdf.select('Latitude')[:].astype(np.float32)
            lons = hdf.select('Longitude')[:].astype(np.float32)
            sds = hdf.select('SensorZenith')
            satz = (sds[:] * sds.attributes().get('scale_factor', 0.01)).astype(np.float32)
        finally:
            hdf.end()
        return cls(lats, lons, satz)

    @property
    def nscan(self):
        return self.lats.shape[0] // self.scan_len

    # ---- 公開介面 ----------------------------------------------------------
    def to_fine(self, factor, scans_per_chunk=20):
        """1km -> 細網格。factor: 4 -> 250m, 2 -> 500m。
        回傳 (lats, lons)，float32，shape (nlines*factor, ncols*factor)。"""
        # 沿軌：1km 像元中心位於 factor 條細線的中央 -> 偏移 (f-1)/2；沿掃描：對齊第一個細像元
        return self.interpolate(factor, (factor - 1) / 2., 0.,
                                self.lats.shape[1] * factor, scans_per_chunk)

    def iter_fine(self, fields, factor, scans_per_chunk=20):
        """把與 1km 經緯度同網格的任意欄位 (如角度、校正係數) 以相同幾何插值到細網格。
        逐掃描帶區塊產生 (細網格列 slice, [各欄位的細網格區塊 float64])，避免整張大陣列佔用記憶體。"""
        return self.iter_chunks(fields, factor, (factor - 1) / 2., 0.,
                                self.lats.shape[1] * factor, scans_per_chunk)

    def interpolate(self, f, row_off, col_off, ncol_fine, scans_per_chunk=20):
        """通用經緯度插值 (於地心直角座標內插)，參數意義見 iter_chunks。"""
        out_lat = np.empty((self.nscan * self.scan_len * f, ncol_fine), dtype=np.float32)
        out_lon = np.empty_like(out_lat)
        xyz = self._lonlat_to_xyz(self.lons, self.lats)
        for rows, (x, y, z) in self.iter_chunks(xyz, f, row_off, col_off, ncol_fine, scans_per_chunk):
            out_lat[rows], out_lon[rows] = self._xyz_to_latlon(x, y, z)
        return out_lat, out_lon

    def iter_chunks(self, fields, f, row_off, col_off, ncol_fine, scans_per_chunk=20):
        """通用掃描帶插值。f: 倍率；細網格第 r 列 (掃描帶內) 對應粗網格列 (r-row_off)/f，
        第 c 行對應 (c-col_off)/f；超出範圍的列/行以兩端格點線性外插。"""
        L = self.scan_len
        nscan, ncol = self.nscan, self.lats.shape[1]

        # 細網格位置 -> 使用的格點對 (i0, i0+1) 與分數權重
        rpos = (np.arange(L * f) - row_off) / f
        ri = np.clip(np.floor(rpos).astype(int), 0, L - 2)
        st = (rpos - ri)[:, None]                              # (L*f, 1)
        cpos = (np.arange(ncol_fine) - col_off) / f
        ci = np.clip(np.floor(cpos).astype(int), 0, ncol - 2)
        ss = (cpos - ci)[None, :]                              # (1, ncol_fine)

        fields = [np.asarray(v, dtype=np.float64).reshape(nscan, L, ncol) for v in fields]
        sz = np.deg2rad(self.satz.astype(np.float64)).reshape(nscan, L, ncol)

        step = L * f
        for s0 in range(0, nscan, scans_per_chunk):
            s1 = min(s0 + scans_per_chunk, nscan)
            c_exp = self._expansion(sz[s0:s1, :, :-1], sz[s0:s1, :, 1:])
            a_scan = ss + ss * (1 - ss) * c_exp[:, ri][:, :, ci]   # 修正後的沿掃描權重
            yield (slice(s0 * step, s1 * step),
                   [self._bilinear(v[s0:s1], ri, ci, st, a_scan).reshape(-1, ncol_fine) for v in fields])

    # ---- 內部計算 ----------------------------------------------------------
    @staticmethod
    def _lonlat_to_xyz(lons, lats):
        lat_r = np.deg2rad(lats.astype(np.float64))
        lon_r = np.deg2rad(lons.astype(np.float64))
        return np.cos(lat_r) * np.cos(lon_r), np.cos(lat_r) * np.sin(lon_r), np.sin(lat_r)

    @staticmethod
    def _xyz_to_latlon(x, y, z):
        return np.rad2deg(np.arctan2(z, np.hypot(x, y))), np.rad2deg(np.arctan2(y, x))

    @staticmethod
    def _bilinear(v, ri, ci, st, a_scan):
        """v: (nscan, L, ncol) 單一座標分量；沿軌權重 st，沿掃描權重 a_scan。"""
        top, bot = v[:, ri], v[:, ri + 1]
        v1 = (1 - a_scan) * top[:, :, ci] + a_scan * top[:, :, ci + 1]
        v2 = (1 - a_scan) * bot[:, :, ci] + a_scan * bot[:, :, ci + 1]
        return (1 - st) * v1 + st * v2

    @classmethod
    def _view_angle(cls, zeta):
        """衛星天頂角 zeta -> 地心角對應的 phi (球面幾何，rad)。"""
        return np.arcsin(cls.EARTH_R * np.sin(zeta) / (cls.EARTH_R + cls.SAT_H))

    @classmethod
    def _zenith_angle(cls, phi):
        """_view_angle 的反函數。"""
        return np.arcsin((cls.EARTH_R + cls.SAT_H) * np.sin(phi) / cls.EARTH_R)

    @classmethod
    def _expansion(cls, satz_a, satz_b):
        """由同一列相鄰兩個粗網格點的衛星天頂角 (rad) 計算 expansion 係數。"""
        phi_a, phi_b = cls._view_angle(satz_a), cls._view_angle(satz_b)
        theta_a, theta_b = satz_a - phi_a, satz_b - phi_b
        phi = (phi_a + phi_b) / 2
        theta = cls._zenith_angle(phi) - phi
        # 星下點兩側對稱時分母為 0 的處理
        denom = np.where(theta_a == theta_b, theta_a * 2, theta_a - theta_b)
        return 4 * (((theta_a + theta_b) / 2 - theta) / denom)


def read_refsb(filename, sds_name, band):
    """讀取 MODIS L1B 反射率 SDS 中的單一波段，回傳 float32 反射率。

    注意 L1B 的 reflectance_scales 給的是 ρ·cos(θs) (未做太陽天頂角校正)。
    DN > 32767 為品質旗標：65533 (飽和) 設為該波段最大可表示反射率，其餘設 NaN。
    """
    hdf = SD.SD(filename, SD.SDC.READ)
    try:
        sds = hdf.select(sds_name)
        attrs = sds.attributes()
        names = [b.strip() for b in attrs['band_names'].split(',')]
        i = names.index(str(band))
        raw = sds[i, :, :]
    finally:
        hdf.end()
    scale = attrs['reflectance_scales'][i]
    offset = attrs['reflectance_offsets'][i]
    refl = (raw.astype(np.float32) - offset) * scale
    refl[raw > 32767] = np.nan
    refl[raw == 65533] = (32767 - offset) * scale
    return refl


class ModisTrueColor:
    """MODIS 真實色 (R = B1 0.645µm, G = B4 0.555µm, B = B3 0.469µm)。

    mode
      1 'RGB'  純 RGB 合成：L1B 反射率 (ρ·cosθs) 直接線性 0–1，不做任何校正與增強
      2 'CORR' 全校正：太陽天頂角校正 + CREFL (Rayleigh 散射、O3/H2O 吸收、地形高度)
               + B1 比值銳化 (250m) + CIRA 對數增強 (satpy true_color 預設增強)
      3 'WV'   NASA Worldview 風格：校正與銳化同 CORR，增強改用 MODIS Rapid Response
               非線性增強 (反射率 0–1.1 線性對應 0–255 後查表)
    解析度跟隨所選檔案：QKM -> 250m (B3/B4 由 500m 升採樣)，HKM -> 500m (B1 用 EV_250_Aggr500)。
    需要同時間的 HKM (B3/B4)、MOD03/MYD03 (經緯度、角度、地形高度)，250m 另需 QKM。
    """

    MODES = {1: 'RGB', 2: 'CORR', 3: 'WV'}
    MODE_DESC = {1: 'RGB (raw composite)', 2: 'CORR (sunz + CREFL + sharpening + CIRA stretch)',
                 3: 'WV (NASA Worldview style)'}

    # CREFL 1.7.1 MODIS 係數 (NASA CREFL SPA；與 satpy _crefl_utils 相同)
    #   band: (aH2O, bH2O, aO3, taur0)
    CREFL_COEF = {1: (-5.60723, 0.820175, 0.0715289, 0.05100),
                  3: (0., 0., 0.00743232, 0.19325),
                  4: (0., 0., 0.089691, 0.09536)}
    UO3, UH2O = 0.319, 2.93          # 氣候值臭氧 (cm-atm)、水氣 (g/cm2)
    MAXAIRMASS = tc.MAXAIRMASS
    SCALEHEIGHT = tc.SCALE_HEIGHT    # m
    REFLMIN, REFLMAX = tc.REFLMIN, tc.REFLMAX
    SUNZ_LIMIT, SUNZ_MAX = 88., 95.  # 太陽天頂角校正上限與漸減終點 (同 satpy sunz_corrected)

    # MODIS Rapid Response true color 增強表 (Worldview Corrected Reflectance 外觀)
    WV_REFL_MAX = tc.WV_REFL_MAX
    WV_XP = list(tc.WV_XP)
    WV_FP = list(tc.WV_FP)

    def __init__(self, l1b_file):
        prefix, tstr = l1b_file[0:3], l1b_file[9:22]
        self.factor = 4 if l1b_file[3:8] == '02QKM' else 2
        self.hkm = self._find(prefix + '02HKM', tstr)
        self.qkm = self._find(prefix + '02QKM', tstr) if self.factor == 4 else None
        self.geofile = self._find(prefix + '03', tstr)
        self._linear = {}      # 快取線性 RGB：False -> 未校正, True -> CREFL 校正 (+銳化)

    @staticmethod
    def _find(product, tstr):
        files = glob('%s.%s*.hdf' % (product, tstr))
        if not files:
            raise FileNotFoundError(tr('真實色需要 %s.%s*.hdf', 'True color requires %s.%s*.hdf') % (product, tstr))
        return files[0]

    # ---- 公開介面 ----------------------------------------------------------
    def compute(self, mode):
        """回傳 (ny, nx, 3) float32 RGB，數值 0–1；無效像元為 NaN。
        CORR 與 WV 共用同一份 CREFL 校正結果 (只差增強)，同一物件連續呼叫時不重算。"""
        if mode not in self.MODES:
            raise ValueError(tr('真實色模式須為 %s 其中之一', 'True color mode must be one of %s') % list(self.MODES))
        rgb = self.linear(corrected=(mode != 1))
        if mode == 1:
            return np.clip(rgb, 0., 1.)
        out = self.cira_stretch(rgb) if mode == 2 else self.worldview_stretch(rgb)
        return out.astype(np.float32, copy=False)

    def linear(self, corrected):
        """線性反射率 RGB (未增強)。corrected=False：L1B 原值；True：太陽天頂角 + CREFL + 銳化。"""
        if corrected in self._linear:
            return self._linear[corrected]
        if self.factor == 4:
            red = read_refsb(self.qkm, 'EV_250_RefSB', 1)
        else:
            red = read_refsb(self.hkm, 'EV_250_Aggr500_RefSB', 1)
        green = read_refsb(self.hkm, 'EV_500_RefSB', 4)
        blue = read_refsb(self.hkm, 'EV_500_RefSB', 3)

        if corrected:
            geo = ModisGeo.from_file(self.geofile)
            ang = self._read_angles(self.geofile)
            red = self.crefl(red, 1, geo, ang, self.factor)
            green = self.crefl(green, 4, geo, ang, 2)
            blue = self.crefl(blue, 3, geo, ang, 2)
            if self.factor == 4:
                green, blue = self.ratio_sharpen(red, green, blue)
        elif self.factor == 4:
            green, blue = self.up2(green), self.up2(blue)
        rgb = np.dstack((red, green, blue))
        self._linear[corrected] = rgb
        return rgb

    # ---- 輻射校正 ----------------------------------------------------------
    @staticmethod
    def _read_angles(geofile):
        """讀 MOD03/MYD03 1km 角度 (度) 與地形高度 (m)。"""
        hdf = SD.SD(geofile, SD.SDC.READ)
        out = {}
        try:
            for key, name in (('sza', 'SolarZenith'), ('saz', 'SolarAzimuth'),
                              ('vza', 'SensorZenith'), ('vaz', 'SensorAzimuth')):
                sds = hdf.select(name)
                raw = sds[:]
                val = raw * sds.attributes().get('scale_factor', 0.01)
                val[raw == sds.attributes().get('_FillValue', -32767)] = np.nan
                out[key] = val.astype(np.float64)
            sds = hdf.select('Height')
            raw = sds[:]
            h = raw.astype(np.float64)
            h[(raw == sds.attributes().get('_FillValue', -32767)) | (h < 0)] = 0.
            out['height'] = h
        finally:
            hdf.end()
        return out

    def atm_fields(self, band, ang):
        """於 1km 網格計算 CREFL 所需欄位，校正式整理為
            ρ_toa = refl_L1B × sunz_corr
            ρ_c   = (ρ_toa / tOG − ρ_ray) / T  ;  ρ_c /= 1 + ρ_c × S
        回傳 (gain, offset, sphalb)，使 ρ_c' = refl_L1B × gain − offset，ρ_c = ρ_c' / (1 + ρ_c' × sphalb)。"""
        ah2o, bh2o, ao3, taur0 = self.CREFL_COEF[band]
        sza = ang['sza']
        with np.errstate(invalid='ignore', divide='ignore'):
            mus = np.cos(np.deg2rad(sza))
            mus = np.where(mus >= 0, mus, np.nan)          # 夜間 -> NaN
            muv = np.cos(np.deg2rad(ang['vza']))
            phi = ang['saz'] - ang['vaz']
            # satpy 對 air mass > 18 設為 -1 (會使 B1 水氣項變 NaN)；此處改為截斷於 18，避免晨昏線附近缺值
            air_mass = np.minimum(1. / mus + 1. / muv, self.MAXAIRMASS)
            t_o3 = tc.o3_transmittance(air_mass, self.UO3, ao3)
            t_h2o = tc.h2o_transmittance_modis(air_mass, self.UH2O, ah2o, bh2o)
            taur = taur0 * np.exp(-ang['height'] / self.SCALEHEIGHT)
            rhoray, trdown, trup = self._chand(phi, muv, mus, taur)
            t_up = ((2 / 3. + muv) + (2 / 3. - muv) * trup) / (4 / 3. + taur)
            t_down = ((2 / 3. + mus) + (2 / 3. - mus) * trdown) / (4 / 3. + taur)
            t_total = t_up * t_down * t_h2o
            gain = self._sunz_corr(sza) / (t_o3 * t_total)
            offset = rhoray / t_total
        return gain, offset, self._csalbr(taur)

    def crefl(self, refl, band, geo, ang, factor):
        """對細網格反射率做太陽天頂角校正 + CREFL；係數在 1km 計算後依掃描幾何插值。"""
        expect = (geo.lats.shape[0] * factor, geo.lats.shape[1] * factor)
        if refl.shape != expect:
            raise ValueError(tr('B%d 的形狀 %s 與定位資料不符 (應為 %s)',
                                'B%d shape %s does not match geolocation (expected %s)') % (band, refl.shape, expect))
        out = np.empty_like(refl)
        for rows, (gain, offset, sphalb) in geo.iter_fine(self.atm_fields(band, ang), factor):
            r = refl[rows] * gain - offset
            r /= 1. + r * sphalb
            out[rows] = np.clip(r, self.REFLMIN, self.REFLMAX)
        return out

    @classmethod
    def _sunz_corr(cls, sza):
        """太陽天頂角校正係數 1/cos(θs)；88° 以上逐漸減弱，95° 以上為 0 (同 satpy sunzen_corr_cos)。"""
        return tc.sunz_corr(sza, cls.SUNZ_LIMIT, cls.SUNZ_MAX)

    @staticmethod
    def _csalbr(tau):
        """Rayleigh 大氣球面反照率 (CREFL csalbr，共用版本見 satima/truecolor.py)。"""
        return tc.csalbr(tau)

    @staticmethod
    def _chand(phi, muv, mus, taur):
        """Rayleigh 路徑反射率與上下行直射穿透率 (CREFL chand，共用版本見 satima/truecolor.py)。
        phi: 太陽與衛星方位角差 (度)；mus/muv: 太陽/衛星天頂角餘弦；taur: 分子光學厚度。"""
        return tc.chand(phi, muv, mus, taur)

    # ---- 500m <-> 250m 與銳化 ----------------------------------------------
    # 幾何關係 (同 ModisGeo 的子像素偏移)：沿軌 500m 第 j 列 = 250m 第 2j、2j+1 列的中心；
    # 沿掃描 500m 第 k 行中心 = 250m 第 2k 行，250m 第 2k+1 行位於 500m 第 k、k+1 行之間。
    @staticmethod
    def up2(a):
        """500m -> 250m：沿軌複製，沿掃描線性內插 (右端複製)。"""
        a = np.repeat(a, 2, axis=0)
        out = np.empty((a.shape[0], a.shape[1] * 2), dtype=a.dtype)
        out[:, 0::2] = a
        out[:, 1:-1:2] = 0.5 * (a[:, :-1] + a[:, 1:])
        out[:, -1] = a[:, -1]
        return out

    @staticmethod
    def down2(a):
        """250m -> 500m：沿軌 2 列平均，沿掃描 [1/4, 1/2, 1/4] 三角權重 (與 up2 對應)。"""
        a = 0.5 * (a[0::2] + a[1::2])
        p = np.pad(a, ((0, 0), (1, 0)), mode='edge')
        return 0.25 * p[:, 0:-1:2] + 0.5 * p[:, 1::2] + 0.25 * p[:, 2::2]

    @classmethod
    def ratio_sharpen(cls, red, green, blue):
        """以 250m B1 與其降解析版本的比值，銳化升採樣後的 B4/B3 (同 satpy RatioSharpenedRGB)。"""
        low = cls.up2(cls.down2(red))
        with np.errstate(invalid='ignore', divide='ignore'):
            ratio = red / low
        ratio[~np.isfinite(ratio) | (ratio < 0)] = 1.
        np.clip(ratio, 0., 1.5, out=ratio)
        return cls.up2(green) * ratio, cls.up2(blue) * ratio

    # ---- 增強 ----------------------------------------------------------------
    @staticmethod
    def cira_stretch(rgb):
        """CIRA 對數增強 (satpy true_color 預設)，輸入反射率 0–1。"""
        return tc.cira_stretch(rgb)

    @classmethod
    def worldview_stretch(cls, rgb):
        """MODIS Rapid Response 非線性增強：反射率 0–1.1 -> 0–255 後分段線性查表。"""
        return tc.worldview_stretch(rgb)


# ===== 單一檔案 =====
class ModisFile:

    SCAN_PERIOD = 1.4771   # s，MODIS 每個掃描帶的週期

    def __init__(self, fullname, channel=None, time=None):
        self.fullname = fullname
        channel = channel or fullname[0:8]
        if channel[3:8] == '021KM':
            self.channel = 'IR 1KM'
        elif channel[3:8] == '02QKM':
            self.channel = 'VIS 250m'
        elif channel[3:8] == '02HKM':
            self.channel = 'VIS 500m'
        else:
            raise ValueError(tr('不認得的 MODIS 產品 %s', 'Unknown MODIS product %s') % channel)
        self.satellite = 'TERRA' if channel[0:3] == 'MOD' else 'AQUA'
        self.tcmode = None   # 真實色模式 (band 0 時由讀取器設定)
        self.center = None   # (clat, clon)，由讀取器設定

        if time:
            self.time = datetime.datetime.strptime(time, 'A%Y%j.%H%M')
        else:
            self.time = datetime.datetime.strptime(self.fullname[9:22], 'A%Y%j.%H%M')

    def __lt__(self, other):
        return self.time < other.time

    def extract(self, band):
        if band == 0:
            # 真實色：回傳 (ny, nx, 3) RGB；物件保留在 self，畫多種模式時共用快取
            if getattr(self, '_truecolor', None) is None:
                self._truecolor = ModisTrueColor(self.fullname)
            return self._truecolor.compute(self.tcmode)
        if self.channel in ('VIS 250m', 'VIS 500m'):
            # 500m 檔的 B1/B2 為 250m 聚合到 500m 的 EV_250_Aggr500_RefSB
            if self.channel == 'VIS 250m':
                sds_name = 'EV_250_RefSB'
            else:
                sds_name = 'EV_500_RefSB' if band >= 3 else 'EV_250_Aggr500_RefSB'
            # 單波段沿用原本作法：無效值 (fill 等旗標) 設 0；飽和值改為最大反射率 (原本會被設 0 變黑)
            return np.nan_to_num(read_refsb(self.fullname, sds_name, band), nan=0.)
        if self.channel == 'IR 1KM':
            # 輻射 (W/m2/um/sr) = (DN - offset) * scale；DN > 32767 為品質旗標 (fill、飽和等) -> NaN
            hdf = SD.SD(self.fullname, SD.SDC.READ)
            try:
                sds = hdf.select('EV_1KM_Emissive')
                attrs = sds.attributes()
                i = [int(b) for b in attrs['band_names'].split(',')].index(band)
                raw = sds[i, :, :]
            finally:
                hdf.end()
            radiance = (raw.astype(np.float64) - attrs['radiance_offsets'][i]) * attrs['radiance_scales'][i]
            radiance[raw > 32767] = np.nan
            with np.errstate(invalid='ignore', divide='ignore'):
                bt = self.brightness_temperature(band, radiance) - 273.15
            return bt

    def brightness_temperature(self, band, radiance):
        # 熱紅外亮溫係數 (MODIS 官方 modis_bright 作法)：
        #   T = (T_mono(cwn) - tci) / tcs，T_mono 為有效中心波數 cwn (cm-1) 下的 Planck 反演
        # 係數取自 NASA OceanColor SeaDAS/OCSSW 的 planck_functions.f90 (cwn/tcs/tci_terra、_aqua)：
        #   https://oceancolor.gsfc.nasa.gov/docs/ocssw/planck__functions_8f90_source.html
        # 原程式以單一波長 (B27 6.715 µm、B31 11.03 µm) 反演，B27 在 Aqua 偏暖約 1.1 K、Terra 偏暖 0.7–0.9 K。
        BT_COEF = {
            'TERRA': {
                'cwn': [2.641775e+03, 2.505277e+03, 2.518028e+03, 2.465428e+03, 2.235815e+03, 2.200346e+03,
                        1.477967e+03, 1.362737e+03, 1.173190e+03, 1.027715e+03, 9.080884e+02, 8.315399e+02,
                        7.483394e+02, 7.308963e+02, 7.188681e+02, 7.045367e+02],
                'tcs': [9.993411e-01, 9.998646e-01, 9.998585e-01, 9.998682e-01, 9.998820e-01, 9.998845e-01,
                        9.994878e-01, 9.994918e-01, 9.995496e-01, 9.997399e-01, 9.995607e-01, 9.997256e-01,
                        9.999161e-01, 9.999167e-01, 9.999192e-01, 9.999282e-01],
                'tci': [4.770522e-01, 9.264053e-02, 9.756834e-02, 8.928794e-02, 7.309468e-02, 7.061112e-02,
                        2.204659e-01, 2.045902e-01, 1.599076e-01, 8.249532e-02, 1.302885e-01, 7.181662e-02,
                        1.970605e-02, 1.912743e-02, 1.816222e-02, 1.579983e-02]},
            'AQUA': {
                'cwn': [2.647409e+03, 2.511760e+03, 2.517908e+03, 2.462442e+03, 2.248296e+03, 2.209547e+03,
                        1.474262e+03, 1.361626e+03, 1.169626e+03, 1.028740e+03, 9.076813e+02, 8.308411e+02,
                        7.482978e+02, 7.307766e+02, 7.182094e+02, 7.035007e+02],
                'tcs': [9.993363e-01, 9.998626e-01, 9.998627e-01, 9.998707e-01, 9.998737e-01, 9.998770e-01,
                        9.995694e-01, 9.994867e-01, 9.995270e-01, 9.997382e-01, 9.995270e-01, 9.997271e-01,
                        9.999173e-01, 9.999070e-01, 9.999198e-01, 9.999233e-01],
                'tci': [4.818401e-01, 9.426663e-02, 9.458604e-02, 8.736613e-02, 7.873285e-02, 7.550804e-02,
                        1.848769e-01, 2.064384e-01, 1.674982e-01, 8.304364e-02, 1.343433e-01, 7.135051e-02,
                        1.948513e-02, 2.131043e-02, 1.804156e-02, 1.683156e-02]},
        }

        """MODIS 熱紅外波段輻射 (W/m2/um/sr) -> 亮溫 (K)，依衛星 (self.satellite) 選係數。"""
        i = EMISSIVE_BANDS.index(band)
        coef = BT_COEF[self.satellite]
        wavel = 1.e4 / coef['cwn'][i]                    # cm-1 -> µm
        return (self.planckInvert(wavel, radiance) - coef['tci'][i]) / coef['tcs'][i]

    def planckInvert(self, wavel, Llambda):
        """input wavelength in microns and Llambda in W/m^2/micron/sr, output
        output brightness temperature in K  (note that we've remove the factor
        of pi because we are working with radiances, not fluxes)
        """
        c = 2.99792458e+08  # m/s -- speed of light in vacumn
        h = 6.62606876e-34  # J s  -- Planck's constant
        kb = 1.3806503e-23  # J/K  -- Boltzman's constant
        c1 = 2. * h * c ** 2.
        c2 = h * c / kb

        Llambda = Llambda * 1.e6  # convert to W/m^2/m/sr
        wavel = wavel * 1.e-6  # convert wavelength to m
        Tbright = c2 / (wavel * np.log(c1 / (wavel ** 5. * Llambda) + 1.))

        return Tbright

    def bilin_2d(self, data, interval, lon_mode=False):
        # (目前未使用，保留)
        data = data.copy()

        if interval == 1:
            return data

        if lon_mode is True:
            max_diff = np.nanmax(np.abs(data[:, :-1] - data[:, 1:]))
            if max_diff > 180.:
                data[data < 0] = 360. + data[data < 0]

        data = np.concatenate((data, data[-1].reshape(1, -1)), axis=0)
        data = np.concatenate((data, data[:, -1].reshape(-1, 1)), axis=1)

        ratio_horizontal = np.tile(np.linspace(0, (interval - 1) / interval, interval, dtype=np.float32),
                                   (data.shape[0] * interval, data.shape[1] - 1))
        ratio_vertical = np.tile(np.linspace(0, (interval - 1) / interval, interval, dtype=np.float32).reshape(-1, 1),
                                 (data.shape[0] - 1, (data.shape[1] - 1) * interval))
        repeat_data = np.repeat(data, interval, axis=0)
        repeat_data = np.repeat(repeat_data, interval, axis=1)

        horizontal_interp = (1. - ratio_horizontal) * repeat_data[:, :-interval] + ratio_horizontal * repeat_data[:, interval:]
        ret = (1. - ratio_vertical) * horizontal_interp[:-interval, :] + ratio_vertical * horizontal_interp[interval:, :]

        return ret

    def _geofile(self):
        geofile = glob('%s03.%s*.hdf' % (self.fullname[0:3], self.fullname[9:22]))
        if not geofile:
            raise FileNotFoundError(tr('找不到定位檔 %s03.%s*.hdf', 'Geolocation file %s03.%s*.hdf not found')
                                    % (self.fullname[0:3], self.fullname[9:22]))
        return geofile[0]

    def geocoord(self, band):
        if self.channel == 'IR 1KM':
            # 1km 經緯度直接取自 MOD03/MYD03 (原本只處理 B27/B31)
            geo = SD.SD(self._geofile(), SD.SDC.READ)
            try:
                lats = geo.select('Latitude')[:]
                lons = geo.select('Longitude')[:]
            finally:
                geo.end()
        else:
            # 1km 經緯度與衛星天頂角取自 MOD03/MYD03，以 satz 掃描幾何插值至 250m/500m
            # (見 ModisGeo，與 satpy/geotiepoints 相同演算法，可處理 bow-tie 與掃描帶邊界)
            # 經緯度只跟檔案解析度有關，快取起來；畫多種真實色模式時不必重算
            if getattr(self, '_fine_geo', None) is None:
                factor = 4 if self.channel == 'VIS 250m' else 2
                self._fine_geo = ModisGeo.from_file(self._geofile()).to_fine(factor)
            lats, lons = self._fine_geo
        return lats, lons

    def row_times(self, nrows):
        """每一列資料的觀測時間 (UTC datetime64)，shape (nrows,)。
        取自 MOD03/MYD03 逐掃描帶的 'EV center time' (TAI93，換 UTC 時扣除 1993 年後的閏秒)；
        同一掃描帶 (1km 10 列、500m 20 列、250m 40 列) 同一時間。缺掃描以前後掃描內插；
        沒有 MOD03 時以檔名時間加掃描週期推算。"""
        geofile = glob('%s03.%s*.hdf' % (self.fullname[0:3], self.fullname[9:22]))
        t = None
        if geofile:
            hdf = SD.SD(geofile[0], SD.SDC.READ)
            try:
                t = hdf.select('EV center time')[:].astype(np.float64)
            finally:
                hdf.end()
            t[t <= 0] = np.nan                                 # 缺掃描 (fill value)
            ok = np.isfinite(t)
            if ok.any():
                k = np.arange(t.size)
                t = obstime.tai93(np.interp(k, k[ok], t[ok]))
            else:
                t = None
        if t is None:
            nscan = nrows // 10 if self.channel == 'IR 1KM' else nrows // (40 if self.channel == 'VIS 250m' else 20)
            t = obstime.from_utc(self.time, (np.arange(nscan) + 0.5) * self.SCAN_PERIOD)
        return obstime.per_scan(t, nrows)

    @staticmethod
    def plot_rgb(canvas, lats, lons, rgb):
        """畫 RGB 影像 (NaN 像元為透明)，實作見 satima/canvas.py。"""
        return canvas.rgb(lons, lats, rgb)

    def plot(self, canvas, lats, lons, data, variant, band):
        """variant：真實色為模式編號，紅外線為色表。回傳色表/模式名稱。"""
        if data.size == 0:
            return ''
        print(data.shape)
        if band == 0:
            self.plot_rgb(canvas, lats, lons, data)
            return ModisTrueColor.MODES[self.tcmode]
        if self.channel.startswith('V'):
            canvas.mesh(lons, lats, data, cmap='gray', vmin=0., vmax=1.)
            return ''
        return canvas.bt_mesh(lons, lats, data, variant)


# ===== 讀取器 =====
class ModisReader(SatReader):
    NAME = 'MODIS-HDF4'
    SENSOR = 'MODIS'
    DESCRIPTION = 'Terra / Aqua MODIS L1B (HDF4)'
    SOURCE = 'NASA'
    GRANULE_SECONDS = 300       # 檔名只有起始時間：一個 granule 標稱 5 分鐘
    PATTERNS = ('M*D0*KM.A*.hdf',)
    MULTI_FILE = True

    def make_file(self, name):
        # MYD021KM.A2023160.0520.061.2023160192407.hdf
        return ModisFile(name, channel=name[0:8], time=name[9:22])

    def listing(self, f):
        return (f.satellite, f.channel, f.time)

    # ---- 選單 (原 getband) ---------------------------------------------------
    def select_band(self):
        kinds = set(f.channel for f in self.files)
        if len(kinds) > 1:
            print(tr('注意：選到的檔案混有不同產品 (%s)，只使用 %s 的檔案。',
                     'WARNING: mixed products selected (%s); only %s files are used.')
                  % (', '.join(sorted(kinds)), self.files[0].channel))
            self.files = [f for f in self.files if f.channel == self.files[0].channel]
        self.multifilesflag = len(self.files) >= 2
        self.tcmode = None
        channel = self.files[0].channel
        self.res = PRODUCT_RES[channel]
        if channel == 'IR 1KM':
            print('27.VAPOR  31.IR  (any emissive band %d-%d except 26)' % (EMISSIVE_BANDS[0], EMISSIVE_BANDS[-1]))
            valid = EMISSIVE_BANDS
        elif channel == 'VIS 250m':
            print('1.VIS 0.645  2.VIS 0.858  0.TRUE COLOR (250m)')
            valid = (0, 1, 2)
        else:
            print('1.VIS 0.645  2.VIS 0.858  3.VIS 0.469  4.VIS 0.555  5.NIR 1.24  6.NIR 1.64  7.NIR 2.13  0.TRUE COLOR (500m)')
            valid = range(0, 8)
        self.band = ui.ask_int(tr('輸入頻道編號：', 'Input the band: '), valid)

    def select_variants(self):
        if self.band == 0:
            return ui.ask_menu('', ModisTrueColor.MODE_DESC, tr('選擇真實色模式 (可複選，0 = 全部)：',
                                                                'Select true color mode: '))
        if self.band in WV_BANDS:
            return ask_colormaps(self.cfg, 'wv')
        if self.band in EMISSIVE_BANDS:
            return ask_colormaps(self.cfg, 'ir')
        return [None]

    # ---- 資料 --------------------------------------------------------------
    def prepare(self):
        for f in self.files:
            f.center = self.center

        # 多檔時沿軌依序串接；rowtime 為每一列的觀測時間 (供中心時間計算)
        data, lats, lons, rowtime = [], [], [], []
        for f in self.files:
            if self.band != 0:          # 真實色的 RGB 依模式在 plot 時才計算
                data.append(f.extract(self.band))
            la, lo = f.geocoord(self.band)
            lats.append(la)
            lons.append(lo)
            rowtime.append(f.row_times(la.shape[0]))
        self.data = np.concatenate(data, axis=0) if data else None
        lats = np.concatenate(lats, axis=0)
        lons = np.concatenate(lons, axis=0)
        rowtime = np.concatenate(rowtime)
        if self.multifilesflag:
            print(tr('已合併 %d 個片段：', 'Merged %d granules:') % len(self.files), lats.shape)

        timeindex = self.slice(self.georange, lats, lons)
        if timeindex is None:
            return False
        if self.data is not None:
            self.data = self.crop(self.data)
        self.time = self.center_time(rowtime)
        if self.band in EMISSIVE_BANDS:
            self.report_eyetemp(self.data)
        return True

    def slice(self, georange, lats, lons):
        """裁出包住 georange 的最小矩形 (含最後一列/行)；範圍內沒有資料時回傳 None。"""
        return self.slice_swath(lats, lons)

    def interpolation(self, lons, lats, data, ip=2):
        # (目前未使用，保留)
        ny, nx = data.shape
        x = np.arange(nx)
        y = np.arange(ny)
        spline = RectBivariateSpline(y, x, data)
        newx = np.arange(0, nx, 1 / ip)
        newy = np.arange(0, ny, 1 / ip)
        return spline(newy, newx)

    def plot(self, canvas, variant):
        f0 = self.files[0]
        if self.band == 0:
            # 每種模式資料不同；經緯度與 CREFL 結果已快取在 ModisFile 中，不會重算
            for f in self.files:
                f.tcmode = variant
            data = self.crop(np.concatenate([f.extract(0) for f in self.files], axis=0))
        else:
            data = self.data
        cm = f0.plot(canvas, self.lats, self.lons, data, variant, self.band)
        if self.band == 0:
            return dict(band='TrueColor', colorbar=cm, label='TRUE COLOR %s' % cm)
        return dict(band='BAND-%02d' % self.band, colorbar=cm, label='BAND-%02d' % self.band)


READER = ModisReader
