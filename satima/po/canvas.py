"""PO 的畫布：掃描帶取樣到共用地圖的格點 (最近鄰)，再交給共用畫布 (satima/core/canvas.py) 畫線、標註、存檔

讀取器 (sensors/po/*.py) 的 plot() 仍然呼叫 mesh() / bt_mesh() / rgb()，介面與以前相同；這裡只把這次要畫的
資料記下來，等主流程決定圖的大小 (原始解析度或 dpi) 之後才取樣 (render)。dpi 模式輸出多張時每張各取樣一次，
讀取器只需要算一次。影像不再經過 matplotlib 的 pcolormesh。

取樣 (同 Satima-PO 的原始解析度模式)：
  1. 輸出格點：共用地圖 (merc / cyl / lcc，WGS84) 的規則格點。原始解析度 = 這次所畫資料的取樣間距，由經緯度量出
     (native_res)：裁切範圍內沿掃描、沿軌兩個方向的間距中較細者，所以不需要各格式另外宣告解析度，
     真實色銳化後的 250 m / 375 m 格點也自動適用。
  2. 每個輸出像素取最近的有效取樣點 (最近鄰)，數值不經內插。於地圖平面 (像素座標) 上以 KD-tree 尋找；
     麥卡托與 LCC 是正形投影，平面上最近即實際上最近 (等距圓柱 cyl 先把東西向乘上中央緯度的 cos，近似等向)。
     離最近取樣點超過該處取樣間距約半格 (FAR 倍) 的像素在掃描帶外，為透明 (顯示 [image] background)。
     資料無效的取樣點不參與，所以 VIIRS bow-tie deletion、掃描帶重疊處會由相鄰掃描的有效像元補上；
     缺掃描 (缺檔) 處超過距離上限，留空不補。輸出像素比取樣粗時等於抽樣。
同一組經緯度的取樣結果 (各輸出像素對應哪個取樣點) 會留著，選「0 = 全部色表」時只算一次。
"""
import warnings

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import Normalize
from scipy.spatial import cKDTree

from . import swath
from ..core.canvas import Canvas as _Canvas
from ..core.i18n import tr

EARTH_KM = 6371.0
FAR = 1.25                  # 離最近取樣點超過「該處取樣格子半對角線」的幾倍視為掃描帶外
STRIP = 1 << 21             # 分條處理：每一條的輸出像素數
SPACING_ROWS = 2000         # 量取樣間距時最多取幾列
GAP_PCT = 99.               # 掃描與掃描之間的間距：沿軌位移的第幾百分位 (一次掃描最多 40 列 = 2.5% 的列是掃描交界)
GAP_MAX = 6.                # 掃描之間的間距最多算到一般間距的幾倍 (更大的是缺掃描，不補)
NORTH_STEP = 0.01           # 量各取樣點處一個像素的實際長度時，往北移動的緯度 (度)

_CACHE = {}                 # 最近的取樣結果 (見 Canvas._locate)


def _distance_km(lat1, lon1, lat2, lon2):
    """大圓距離 (km，球面；輸入為弧度)。"""
    h = np.sin((lat2 - lat1) / 2.) ** 2 + np.cos(lat1) * np.cos(lat2) * np.sin((lon2 - lon1) / 2.) ** 2
    return 2. * EARTH_KM * np.arcsin(np.sqrt(np.clip(h, 0., 1.)))


def spacing_km(lats, lons):
    """各行 (沿掃描位置) 的取樣間距 (km)：回傳 (沿掃描 sx, 沿軌 sy, 沿軌最大間距 sgap)，
    各為長度 = 行數的陣列；無法計算時 None。

    sx、sy 取跨列的中位數：極軌成像儀的像元大小主要隨掃描位置變化；MODIS / VIIRS 一次掃描有 10~40 列，
    掃描與掃描之間的跳動每 10~40 列才出現一次，不影響中位數。
    sgap (>= sy) 是掃描與掃描之間的間距：相鄰兩列的位移投影到沿軌方向 (各行位移的中位數方向) 後取第 GAP_PCT 百分位。
    掃描帶邊緣相鄰掃描互相重疊 (bow-tie)，投影為負，sgap = sy；星下點附近兩次掃描之間可能有空隙
    (軌道降低後的 Aqua MODIS：掃描內 0.98 km、掃描間 1.57 km)，sgap 即這個較大的間距，取樣時才不會留下細縫。
    列數很多時只取其中 SPACING_ROWS 列。"""
    ny, nx = lats.shape
    if ny < 2 or nx < 2:
        return None
    rows = np.arange(ny - 1)
    if rows.size > SPACING_ROWS:
        rows = np.unique(np.random.default_rng(0).integers(0, ny - 1, SPACING_ROWS))
    la0, lo0 = np.deg2rad(lats[rows]), np.deg2rad(lons[rows])
    la1, lo1 = np.deg2rad(lats[rows + 1]), np.deg2rad(lons[rows + 1])
    with warnings.catch_warnings(), np.errstate(invalid='ignore'):
        warnings.simplefilter('ignore', RuntimeWarning)             # 整行都是 NaN
        # 相鄰兩列的位移 (km；北、東分量)
        north = (la1 - la0) * EARTH_KM
        east = (lo1 - lo0) * EARTH_KM * np.cos(0.5 * (la0 + la1))
        mn, me = np.nanmedian(north, axis=0), np.nanmedian(east, axis=0)
        sy = np.hypot(mn, me)
        sgap = np.nanpercentile((north * mn + east * me) / sy, GAP_PCT, axis=0)
        dx = np.nanmedian(_distance_km(la0[:, :-1], lo0[:, :-1], la0[:, 1:], lo0[:, 1:]), axis=0)
        # 行與行之間的間距 (少一個) -> 各行：左右兩個的平均，邊緣用單邊
        sx = np.nanmean(np.stack([np.append(dx, np.nan), np.append(np.nan, dx)]), axis=0)
    out = []
    for s in (sx, sy):
        good = np.isfinite(s) & (s > 1e-3)
        if not good.any():
            return None
        out.append(np.where(good, s, np.median(s[good])))
    sx, sy = out
    with np.errstate(invalid='ignore'):
        sgap = np.where(np.isfinite(sgap), np.clip(sgap, sy, GAP_MAX * sy), sy)
    return sx, sy, sgap


class Canvas(_Canvas):
    """讀取器使用的介面：mesh / bt_mesh / rgb；主流程使用：clear / empty / native_res / render，其餘同共用畫布。"""

    def __init__(self, view, cfg):
        super().__init__(view, cfg)
        latmin, latmax, lonmin, lonmax = (float(v) for v in view.georange)
        self.clon = 0.5 * (lonmin + lonmax)         # 經度一律展開到以地圖中央經度為中心的 ±180 度
        self.layer = None

    # ---- 讀取器呼叫的介面 (介面同以前的 MapCanvas) ---------------------------
    def clear(self):
        self.layer = None

    @property
    def empty(self):
        """讀取器沒有交出可畫的資料。"""
        return self.layer is None

    def mesh(self, lons, lats, data, cmap=None, vmin=None, vmax=None, norm=None, **kwargs):
        """純量場 (反射率、亮溫、DNB)；cmap / vmin / vmax / norm 的意義同 matplotlib。"""
        values = swath.as_float(np.ma.masked_invalid(np.ma.asarray(data)))
        valid = np.isfinite(values)
        if not self._check(lons, lats, valid):
            return None
        if norm is None:
            norm = Normalize(vmin, vmax)
        if norm.vmin is None or norm.vmax is None:
            norm.autoscale_None(values[valid])
        cmap = plt.get_cmap(cmap)
        flat = values.ravel()
        self.layer = dict(lons=lons, lats=lats, valid=valid, colorize=lambda k: cmap(norm(flat[k]), bytes=True))
        return self.layer

    def bt_mesh(self, lons, lats, data, entry):
        """亮溫 (°C) 配色表；回傳色表名稱。"""
        data = np.ma.masked_outside(np.ma.masked_invalid(np.ma.asarray(data)), -150., 100.)
        cmap, vmin, vmax = entry.load()
        self.mesh(lons, lats, data, cmap=cmap, vmin=vmin, vmax=vmax)
        return entry.name

    def rgb(self, lons, lats, rgb):
        """RGB 影像 (數值 0–1，NaN 為透明)。"""
        rgb = np.asarray(np.ma.filled(rgb, np.nan) if np.ma.isMaskedArray(rgb) else rgb, dtype=np.float32)
        valid = np.isfinite(rgb).all(axis=-1)
        if not self._check(lons, lats, valid):
            return None
        flat = rgb.reshape(-1, 3)

        def colorize(k):
            out = np.full((k.size, 4), 255, dtype=np.uint8)
            out[:, :3] = np.rint(np.clip(flat[k], 0., 1.) * 255.)
            return out

        self.layer = dict(lons=lons, lats=lats, valid=valid, colorize=colorize)
        return self.layer

    def _check(self, lons, lats, valid):
        if np.shape(lats) != valid.shape or np.shape(lons) != valid.shape:
            raise ValueError(tr('資料 %s 與經緯度 %s 的形狀不同', 'data %s and lat/lon %s have different shapes')
                             % (valid.shape, np.shape(lats)))
        if not valid.any():
            print(tr('注意：沒有有效的資料可畫。', 'WARNING: no valid data to plot.'))
            return False
        return True

    # ---- 主流程呼叫的介面 ---------------------------------------------------
    def _geometry(self):
        """這次資料的經緯度 (經度已展開) 與取樣間距；結果留在 layer 裡。"""
        L = self.layer
        if 'sp' not in L:
            L['lat'] = swath.as_float(L['lats'])
            L['lon'] = swath.wrap_lon(L['lons'], self.clon)
            L['sp'] = spacing_km(L['lat'], L['lon'])
        return L['lat'], L['lon'], L['sp']

    def native_res(self):
        """這次資料的原始取樣間距 (km)：沿掃描、沿軌兩個方向中較細者；量不出來時 None。"""
        _, _, sp = self._geometry()
        if sp is None:
            print(tr('注意：無法量出取樣間距 (有效的經緯度太少)。',
                     'WARNING: cannot measure the sampling distance (too few valid coordinates).'))
            return None
        return float(min(sp[0].min(), sp[1].min()))

    def render(self, width):
        """取樣成 (列, width, 4) uint8 的影像 (第一列為北邊；無資料處 alpha = 0)。掃描帶沒有涵蓋地圖時回傳 None。"""
        src = self._locate(width)
        if src is None:
            return None
        ny, nx = src.shape
        image = np.zeros((ny, nx, 4), dtype=np.uint8)
        colorize = self.layer['colorize']
        step = max(1, STRIP // nx)
        for i0 in range(0, ny, step):
            k = src[i0:i0 + step]
            hit = k >= 0
            if hit.any():
                image[i0:i0 + step][hit] = colorize(k[hit])
        return image

    def _locate(self, width):
        """各輸出像素最近的有效取樣點：(列, 行) int32，值為取樣點在攤平陣列中的位置，沒有為 -1。
        範圍內沒有可用的取樣點時回傳 None。"""
        L = self.layer
        valid = L['valid']
        key = (id(L['lons']), id(L['lats']), self.view, int(width))
        hit = _CACHE.get(key)
        if hit is not None and np.array_equal(hit['valid'], valid):
            return hit['src']
        lats, lons, sp = self._geometry()
        if sp is None:
            return None
        sx, sy, sgap = sp
        ny = max(1, int(round(width * self.aspect)))

        # 可用的取樣點：經緯度與資料都有效，且在地圖範圍附近 (裁切框是索引空間的矩形，斜的掃描帶會帶進很多範圍外的點)
        latmin, latmax, lonmin, lonmax = self.frame_range()
        # 各行的距離上限：取樣格子的半對角線 x FAR (沿軌取掃描之間的間距，星下點的掃描間空隙才補得起來)
        reach_km = FAR * 0.5 * np.hypot(sx, sgap)
        mlat = 1.5 * float(reach_km.max()) / 111.
        mlon = mlat / max(0.05, np.cos(np.deg2rad(min(89., max(abs(latmin), abs(latmax))))))
        with np.errstate(invalid='ignore'):
            near = (valid & (np.abs(lats) < 89.95) & (lats > latmin - mlat) & (lats < latmax + mlat)
                    & (lons > lonmin - mlon) & (lons < lonmax + mlon))
        idx = np.flatnonzero(near)
        del near
        if idx.size == 0:
            print(tr('注意：地圖範圍內沒有有效的資料。', 'WARNING: no valid data inside the map range.'))
            return None
        lat_s = lats.ravel()[idx].astype(np.float64)
        lon_s = lons.ravel()[idx].astype(np.float64)
        x0, x1, y0, y1 = self.map.extent
        x, y = self.map.forward(lon_s, lat_s)
        # 各取樣點處一個像素的實際長度 (km)：往北移一小段，看在圖上移動了幾個像素 (適用各種投影)
        x2, y2 = self.map.forward(lon_s, np.minimum(lat_s + NORTH_STEP, 89.99))
        with np.errstate(invalid='ignore', divide='ignore'):
            moved = np.hypot((x2 - x) / (x1 - x0) * width, (y2 - y) / (y1 - y0) * ny)
            thr = reach_km[idx % lats.shape[1]] / (np.deg2rad(NORTH_STEP) * EARTH_KM / moved)
        del x2, y2, moved, lat_s, lon_s
        # 等距圓柱不是正形投影 (東西向的像素比較短)：平面座標的東西向乘上中央緯度的 cos，近似等向
        kx = np.cos(np.deg2rad(0.5 * (latmin + latmax))) if self.view.projection == 'cyl' else 1.
        pts = np.empty((idx.size, 2))
        pts[:, 0] = ((x - x0) / (x1 - x0) * width - 0.5) * kx      # 像素座標：整數 = 像素中心
        pts[:, 1] = (y1 - y) / (y1 - y0) * ny - 0.5                # 第 0 列為北邊
        del x, y
        good = np.isfinite(pts).all(axis=1) & np.isfinite(thr)
        if not good.all():
            idx, pts, thr = idx[good], pts[good], thr[good]
        if idx.size == 0:
            return None
        thr = thr.astype(np.float32)
        tree = cKDTree(pts, leafsize=16, compact_nodes=False, balanced_tree=False)
        bound = float(thr.max())

        src = np.full((ny, width), -1, dtype=np.int32)
        cols = np.arange(width, dtype=np.float64) * kx
        step = max(1, STRIP // width)
        for i0 in range(0, ny, step):
            i1 = min(ny, i0 + step)
            q = np.empty((i1 - i0, width, 2))
            q[..., 0] = cols
            q[..., 1] = np.arange(i0, i1, dtype=np.float64)[:, None]
            dist, k = tree.query(q.reshape(-1, 2), distance_upper_bound=bound, workers=-1)
            hit = np.flatnonzero(k < idx.size)                      # 超過距離上限的像素沒有結果
            k = k[hit]
            ok = dist[hit] <= thr[k]
            src[i0:i1].reshape(-1)[hit[ok]] = idx[k[ok]]
        if not (src >= 0).any():
            print(tr('注意：掃描帶沒有涵蓋這個地圖範圍。', 'WARNING: the swath does not cover the selected range.'))
            return None
        _CACHE.clear()                  # 只留最近一組幾何 (含不同寬度會互相清掉；dpi 模式的張數很少，重算即可)
        _CACHE[key] = dict(valid=valid.copy(), src=src, refs=(L['lons'], L['lats']))    # refs：留住物件，id 才不會重複
        return src

    def frame_range(self):
        """圖框涵蓋的經緯度範圍 (南, 北, 西, 東；經度以地圖中央經度為中心展開)。
        merc / cyl 就是地圖範圍；LCC 的圖框是涵蓋整個範圍的最小矩形，四角會超出輸入的經緯度範圍。"""
        if not hasattr(self, '_frame'):
            latmin, latmax, lonmin, lonmax = (float(v) for v in self.view.georange)
            if self.view.projection == 'lcc':
                x0, x1, y0, y1 = self.map.extent
                t = np.linspace(0., 1., 201)
                bx = np.concatenate([x0 + (x1 - x0) * t, np.full(t.size, x1), x0 + (x1 - x0) * t, np.full(t.size, x0)])
                by = np.concatenate([np.full(t.size, y0), y0 + (y1 - y0) * t, np.full(t.size, y1), y0 + (y1 - y0) * t])
                lon, lat = self.map.inverse(bx, by)
                lon = swath.wrap_lon(lon, self.clon)
                ok = np.isfinite(lon) & np.isfinite(lat)
                if ok.any():
                    latmin, latmax = min(latmin, float(lat[ok].min())), max(latmax, float(lat[ok].max()))
                    lonmin, lonmax = min(lonmin, float(lon[ok].min())), max(lonmax, float(lon[ok].max()))
            self._frame = (latmin, latmax, lonmin, lonmax)
        return self._frame
