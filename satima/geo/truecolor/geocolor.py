"""GeoColor：白天真實色 + 夜間紅外高雲、低雲與城市燈光 (CIRA GeoColor V2.0；Miller et al. 2020 JTECH 第 5 節)

白天圖層 DL：本套件的 SHAC 真實色 (composite.py)，含 SZA 80°~90° 的晨昏線淡出：
  CIRA 的白天圖層在晨昏帶同樣變暗 (反推雲頂色階 81°/83°/85°：CIRA 201/172/151，本套件 200/177/157)。
夜間圖層 NL (由上而下三層，Eq. 12)：
  高雲   N_IR = 1 - N(BT10.4)[IRmin, IRmax]，IRmin 在 |緯度| 30°~60° 之間增加 20 K (Eq. 9)，白色
  低雲   LC = N(BT10.4 - BT3.9)[1.0, 4.5] (陸地) / [0.0, 4.0] (水面)；BT10.4 < 230 K 時為 0，顏色 (0.55, 0.75, 0.98)
  地表   城市燈光 N_L = N(log10 輻射量)[-0.5, 2.0]，N_L > 0.2 處為 (0.8 x N_L)^(0.75, 1.25, 2.0) (鈉燈的橘黃色)；
         其餘陸地為深紫 (0.06, 0.03, 0.13)，水面為黑色
  NL = N_IR + (1 - N_IR) x [A x LC + (1 - LC) x S]
日夜混合 (Eq. 13)：C = Nmo x DL + (1 - Nmo) x NL，Nmo = N(cos SZA)[c0, c1]^p
  3.9 µm 在晨昏時含太陽反射，低雲會在日出日落前後短暫變淡；這段時間白天圖層佔多數，CIRA 的說明書也提到同樣現象。
地表海拔的浮雕 (論文以 0~50 km 正規化，幾乎看不出來) 省略。

預設參數取自 CIRA 目前 (2025 年更新後) 的 SLIDER 圖片，與論文不同的部分
(由 GOES-19 2026-11-12 23:01 UTC 與 MTG-I1 2026-11-09 13:20 UTC 逐像素反推，兩顆衛星結果一致)：
  紅外範圍  [178, 301] K (論文 [200, 280])；冷端同樣隨緯度增加 20 K (反推：熱帶約 178、40° 約 186、50° 約 195 K)
  日夜混合  cos SZA [0.00, 0.18]、次方 1.75，SZA 約 79.6°~90° 漸變 (論文 [0.1, 0.3]、1.5，即 72.5°~84.3°)；
            CIRA 2025 年的說明：縮窄晨昏帶、白天延伸最多約一小時
  水面背景  黑色 (CIRA 夜間晴空海面的灰色來自紅外圖層)

需要的波段：真實色的 red / green / blue / nir，加上 ir (10.4 µm) 與 swir (3.9 µm)。
"""
import time

import numpy as np

from .. import resample
from ...core.i18n import tr
from . import astronomy, background
from .composite import CHUNK, TrueColor

DEFAULTS = {
    'nightlights': 'nightlights',
    'twilight_cos': [0.0, 0.18],
    'twilight_power': 1.75,
    'ir_range': [178., 301.],
    'low_cloud_color': [0.55, 0.75, 0.98],
    'land_color': [0.06, 0.03, 0.13],
    'water_color': [0.0, 0.0, 0.0],
}
IR_LAT_RISE = 20.           # K：IRmin 在 |緯度| 30°~60° 之間增加的量 (論文 Eq. 9)
LIGHT_EXPONENTS = np.array([0.75, 1.25, 2.0])
LIGHT_DIM = 0.8
LIGHT_THRESHOLD = 0.2
LOW_CLOUD_MIN_BT = 230.     # K
BTD_LAND = (1.0, 4.5)       # K
BTD_WATER = (0.0, 4.0)


def _norm(x, lo, hi):
    """CIRA Eq. 1：線性正規化到 0~1"""
    with np.errstate(invalid='ignore'):
        return np.clip((x - lo) / (hi - lo), 0., 1.)


def day_weight(cos_sza, opts):
    lo, hi = (float(v) for v in opts['twilight_cos'])
    return _norm(cos_sza, lo, hi) ** float(opts['twilight_power'])


def high_cloud(bt_k, lat, opts):
    """高雲圖層的不透明度 N_IR = 1 - N(BT)[IRmin, IRmax] (Eq. 9)，IRmin 在 |緯度| 30°~60° 之間增加 20 K。
    bt_k 為 K (與 lat 同形狀)；無資料為 0。夜光底圖 (satima/landsea.py) 也使用。"""
    ir_lo, ir_hi = (float(v) for v in opts['ir_range'])
    irmin = ir_lo + IR_LAT_RISE * np.clip((np.abs(lat) - 30.) / 30., 0., 1.)
    n_ir = np.nan_to_num(1. - _norm(bt_k, irmin, ir_hi))
    return np.where(np.isfinite(bt_k), n_ir, 0.)


def add_lights(surf, rad):
    """城市燈光：N_L = N(log10 輻射量)[-0.5, 2.0]，N_L > 0.2 處的地表改為 (0.8 x N_L)^(0.75, 1.25, 2.0)
    (鈉燈的橘黃色)。surf (..., 3) 就地修改並回傳；rad 與 surf[..., 0] 同形狀。夜光底圖也使用。"""
    nl = _norm(np.log10(np.maximum(rad, 1e-10)), -0.5, 2.0)
    lit = nl > LIGHT_THRESHOLD
    surf[lit] = (LIGHT_DIM * nl[lit, None]) ** LIGHT_EXPONENTS
    return surf


def night_layer(bt, swir, lon, lat, lights, opts):
    """夜間圖層 (n, 3)，數值 0~1。bt、swir 為 °C。"""
    bt_k, sw_k = bt.astype(np.float64) + 273.15, swir.astype(np.float64) + 273.15
    n_ir = high_cloud(bt_k, lat, opts)

    land = background.land_mask(lon, lat)
    if land is None:
        land = np.ones(lon.shape, dtype=bool)
    btd = np.where(bt_k < LOW_CLOUD_MIN_BT, 0., bt_k - sw_k)
    lo = np.where(land, BTD_LAND[0], BTD_WATER[0])
    hi = np.where(land, BTD_LAND[1], BTD_WATER[1])
    lc = np.nan_to_num(_norm(btd, lo, hi))

    surf = np.where(land[:, None], np.asarray(opts['land_color'], float), np.asarray(opts['water_color'], float))
    if lights is not None:
        add_lights(surf, lights.sample(lon, lat))
    low =np.asarray(opts['low_cloud_color'], float)
    return n_ir[:, None] + (1. - n_ir[:, None]) * (low * lc[:, None] + (1. - lc[:, None]) * surf)


class GeoColor(TrueColor):
    """GeoColor：TrueColor 加上夜間圖層與日夜混合 (分條結果與整張相同)。"""

    def __init__(self, reader, scene, recipe, grid, options=None, geo_options=None, low_memory=False, home='.',
                 interp='nearest'):
        super().__init__(reader, scene, recipe, grid, options, low_memory)
        self.geo = dict(DEFAULTS)
        self.geo.update(geo_options or {})
        self.interp = interp            # 夜間紅外 (解析度通常比輸出粗) 的取值方式；白天波段維持最近鄰 + 比值銳化
        self.night = {}
        for role in ('ir', 'swir'):
            name = recipe[role]
            if role == 'ir' and 'bt' in self.sources:
                self.night[role] = self.sources['bt']
            else:
                k = self.factor(name)
                self.night[role] = (resample.Source(reader, scene, name, grid, k, coarse=low_memory, interp=interp), k)
        self.lights = background.night_lights(self.geo['nightlights'], home, grid.pixel_km)
        background.land_mask(np.zeros(1), np.zeros(1))     # 海陸遮罩在平行計算前先載入
        self.n_day = self.n_night = 0

    def compute(self, grid):
        """回傳 (列, 行, 3) 的 float32 陣列，數值 0~1，第一列為北邊。"""
        bands = self.sample_bands(grid)
        night = {role: src.sample(grid, k, interp=self.interp) for role, (src, k) in self.night.items()}

        vis_ok = np.ones(grid.shape, dtype=bool)
        for role in ('blue', 'green', 'nir', 'red_lo', 'red_hi'):
            if role in bands:
                vis_ok &= np.isfinite(bands[role])
        ir_ok = np.isfinite(night['ir'])
        index = np.flatnonzero(vis_ok | ir_ok)
        vis_ok, ir_ok = vis_ok.ravel(), ir_ok.ravel()
        n_day = n_night = 0

        flat = {role: band.ravel() for role, band in bands.items()}
        nflat = {role: band.ravel() for role, band in night.items()}
        rgb = np.zeros((grid.size, 3), dtype=np.float32)
        for start in range(0, index.size, CHUNK):
            sel = index[start:start + CHUNK]
            lon_t, lat_t = grid.lonlat_at(sel)
            w = day_weight(astronomy.cos_sun_zenith(self.scene.time, lon_t, lat_t), self.geo)
            out = np.zeros((sel.size, 3))
            day = (w > 0) & vis_ok[sel]
            if day.any():
                ds = sel[day]
                values = {role: band[ds] for role, band in flat.items()}
                out[day] = w[day, None] * self.day(values, lon_t[day], lat_t[day])
                n_day += int(day.sum())
            nt = (w < 1) & ir_ok[sel]
            if nt.any():
                ns = sel[nt]
                out[nt] += (1. - w[nt, None]) * night_layer(nflat['ir'][ns], nflat['swir'][ns], lon_t[nt], lat_t[nt],
                                                             self.lights, self.geo)
                n_night += int(nt.sum())
            rgb[sel] = out
        with self._lock:
            self.valid += index.size
            self.n_day += n_day
            self.n_night += n_night
        return rgb.reshape(grid.shape + (3,))

    def finish(self):
        if self.valid == 0:
            raise ValueError(tr('指定範圍內沒有有效資料。', 'No valid data in the requested range.'))
        print(tr('  GeoColor 計算完成 (%.1f 秒；白天/晨昏 %d、夜間/晨昏 %d 像素)',
                 '  GeoColor done (%.1f s; day/twilight %d, night/twilight %d pixels)') % (
            time.time() - self.t0, self.n_day, self.n_night))
