"""海陸底圖 (GRAY 的第 2 層)：灰階雲疊在單色的陸地與海洋上 (CWA，仿 CWA 官網的「紅外線彩色雲圖」)，
或疊在城市燈光上 (LIGHTS 夜光底圖，雲的透明度與顏色依 GeoColor 夜間圖層，見 compose 的 night)

以下為 CWA 底圖：

合成：輸出 = a x 灰階 + (1 - a) x 底色
  底色    陸地與海洋兩色，依陸地比例 (海岸處的反鋸齒比例) 混合；地形紋理來自紅外本身 (地表溫度差)
  紅外    灰階 = 亮溫的分段線性對照 (ir_gray)；a = clip((T0 - T) / (T0 - T1), 0, 1) ^ alpha_power
  可見光  灰階 = 反射率 (同單波段 GRAY，[visible] 的倍數與 gamma 也適用)；a 依反射率 vis_alpha 範圍同樣漸變
預設值是以 CWA 官網 2026-11-13 11:30 UTC 的 Himawari-9 全圓盤「紅外線彩色」與「紅外線灰階」圖，
對照同時次的 B13 亮溫逐像素反推 (透明度擬合誤差 0.009)。

海陸比例：basemap 內建的 GSHHS 全解析度 (f) 海岸線多邊形 (湖泊算水、湖中島算陸地)，
以 matplotlib Agg 反鋸齒柵格化成 30 角秒 (約 0.93 km) 的陸地比例 (0~255)，
存成分塊 DEFLATE GeoTIFF (含 2、4、8… 倍的縮圖層)。第一次使用時自動產生，之後直接讀檔。
取樣時選不比輸出像素粗的最粗一層，雙線性內插。
"""
import os
import threading
import time

import numpy as np

from ..core.i18n import tr

ARCSEC = 30                    # 遮罩解析度 (角秒)
TILE = 512
STRIP = 1024                    # 產生時每次柵格化的列數
_LOCK = threading.Lock()
_OPEN = {}


# ---- 產生遮罩 ------------------------------------------------------------------
def _gshhs_polygons():
    """GSHHS 全解析度多邊形 -> {類別: [(南, 北, 經緯度陣列)]}；1 陸地、2 湖泊、3 湖中島、5 南極。"""
    import mpl_toolkits.basemap_data as bd
    folder = list(bd.__path__)[0]
    polys = {}
    with open(os.path.join(folder, 'gshhs_f.dat'), 'rb') as f:
        for line in open(os.path.join(folder, 'gshhsmeta_f.dat')):
            s = line.split()
            typ, south, north = int(s[0]), float(s[3]), float(s[4])
            f.seek(int(s[5]))
            xy = np.frombuffer(f.read(int(s[6])), '<f4').reshape(-1, 2).astype(np.float64)
            polys.setdefault(typ, []).append((south, north, xy))
    return polys


def _compound(items, south, north):
    """與緯度範圍相交的多邊形合成一條複合路徑 (相鄰多邊形的共用邊不會留下反鋸齒的接縫)。"""
    from matplotlib.path import Path
    sel = [xy for s, n, xy in items if n >= south and s <= north]
    if not sel:
        return None
    verts = np.concatenate(sel)
    codes = np.full(len(verts), Path.LINETO, dtype=Path.code_type)
    start = 0
    for xy in sel:
        codes[start] = Path.MOVETO
        start += len(xy)
    return Path(verts, codes)


SUB = 4                         # 柵格化時每列的垂直取樣數 (水平方向為精確的覆蓋比例)

try:
    import numba as _nb
except ImportError:
    _nb = None

if _nb is not None:
    @_nb.njit(cache=True)
    def _crossings(x0, y0, x1, y1, sub, nsub):
        """各邊與子掃描線 (中心 y = (k + 0.5) / sub) 的交點 -> (子列, x)。端點採半開區間，頂點不會重複計算。"""
        n = 0
        for i in range(x0.size):
            ya, yb = min(y0[i], y1[i]), max(y0[i], y1[i])
            k0 = max(0, int(np.ceil(ya * sub - 0.5)))
            k1 = min(nsub, int(np.ceil(yb * sub - 0.5)))
            if k1 > k0:
                n += k1 - k0
        ks = np.empty(n, np.int64)
        xs = np.empty(n, np.float64)
        n = 0
        for i in range(x0.size):
            if y0[i] == y1[i]:
                continue
            ya, yb = min(y0[i], y1[i]), max(y0[i], y1[i])
            k0 = max(0, int(np.ceil(ya * sub - 0.5)))
            k1 = min(nsub, int(np.ceil(yb * sub - 0.5)))
            slope = (x1[i] - x0[i]) / (y1[i] - y0[i])
            for k in range(k0, k1):
                ks[n] = k
                xs[n] = x0[i] + ((k + 0.5) / sub - y0[i]) * slope
                n += 1
        return ks[:n], xs[:n]

    @_nb.njit(cache=True)
    def _fill(ks, xs, k_first, row0, sub, sign, cov):
        """依交點 (已依子列、x 排序) 以奇偶規則填入覆蓋比例 cov[列 - row0, 行] (乘上 sign 累加)。"""
        nrow, nx = cov.shape
        diff = np.zeros(nx + 1, np.float32)
        w = sign / sub
        i = k_first
        n = ks.size
        for r in range(nrow):
            kend = (row0 + r + 1) * sub
            any_ = False
            while i < n and ks[i] < kend:
                j = i
                while j < n and ks[j] == ks[i]:
                    j += 1
                for p in range(i, j - 1, 2):            # 成對的交點為一段 (奇數時略過最後一個)
                    xa = min(max(xs[p], 0.), float(nx))
                    xb = min(max(xs[p + 1], 0.), float(nx))
                    if xb <= xa:
                        continue
                    ia, ib = int(xa), int(xb)
                    if ia == ib:
                        cov[r, min(ia, nx - 1)] += (xb - xa) * w
                    else:
                        cov[r, ia] += (ia + 1 - xa) * w
                        if ib < nx:
                            cov[r, ib] += (xb - ib) * w
                        diff[ia + 1] += w
                        diff[ib] -= w
                    any_ = True
                i = j
            if any_:
                acc = 0.
                for c in range(nx):
                    acc += diff[c]
                    cov[r, c] += acc
                    diff[c] = 0.
                diff[nx] = 0.
        return i


def _edges(items, d):
    """多邊形 -> 各邊端點的像素座標 (x 向東、y 向南，(0, 0) 為 180°W、90°N 的角落)。"""
    x0, y0, x1, y1 = [], [], [], []
    for _, _, xy in items:
        px = (xy[:, 0] + 180.) / d
        py = (90. - xy[:, 1]) / d
        if xy[0, 0] != xy[-1, 0] or xy[0, 1] != xy[-1, 1]:      # 確保封閉
            px, py = np.r_[px, px[0]], np.r_[py, py[0]]
        x0.append(px[:-1]); y0.append(py[:-1]); x1.append(px[1:]); y1.append(py[1:])
    if not x0:
        return None
    return tuple(np.concatenate(a) for a in (x0, y0, x1, y1))


def _rasterize_numba(layers, nx, ny, d):
    frac = np.zeros((ny, nx), dtype=np.uint8)
    prepared = []
    for items, sign in layers:
        e = _edges(items, d)
        if e is None:
            continue
        ks, xs = _crossings(*e, SUB, ny * SUB)
        order = np.lexsort((xs, ks))
        prepared.append([ks[order], xs[order], 0, sign])
    for r0 in range(0, ny, STRIP):
        r1 = min(ny, r0 + STRIP)
        cov = np.zeros((r1 - r0, nx), dtype=np.float32)
        for item in prepared:
            item[2] = _fill(item[0], item[1], item[2], r0, SUB, item[3], cov)
        frac[r0:r1] = (np.clip(cov, 0., 1.) * 255. + 0.5).astype(np.uint8)
        print('\r    %3d%%' % (100 * r1 // ny), end='', flush=True)
    print()
    return frac


def build(path):
    """以 GSHHS 全解析度海岸線產生陸地比例 GeoTIFF (30 角秒；numba 約 1 分鐘，沒有 numba 時約 10 分鐘)。"""
    t0 = time.time()
    print(tr('  第一次使用海陸底圖：以 GSHHS 全解析度海岸線產生海陸遮罩 (只需一次) ...',
             '  First use of the land/sea basemap: building the land/sea mask from the GSHHS full-resolution '
             'coastlines (one time only) ...'))
    polys = _gshhs_polygons()
    d = ARCSEC / 3600.
    nx, ny = int(round(360. / d)), int(round(180. / d))
    if _nb is not None:
        # 陸地與南極 (+)、湖泊 (-)、湖中島 (+)
        frac = _rasterize_numba([(polys.get(1, []) + polys.get(5, []), 1.), (polys.get(2, []), -1.),
                                 (polys.get(3, []), 1.)], nx, ny, d)
    else:
        frac = _rasterize_agg(polys, nx, ny, d)
    _write(path, frac, d)
    print(tr('  海陸遮罩已存檔 %s (%dx%d，%.1f MB，%.0f 秒)', '  Land/sea mask saved %s (%dx%d, %.1f MB, %.0f s)')
          % (path, nx, ny, os.path.getsize(path) / 1e6, time.time() - t0))


def _rasterize_agg(polys, nx, ny, d):
    """沒有 numba 時：matplotlib Agg 反鋸齒繪製 (每個頂點約 20 µs，全球約 10 分鐘)。"""
    import matplotlib
    matplotlib.use('Agg')
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from matplotlib.figure import Figure
    from matplotlib.patches import PathPatch

    frac = np.zeros((ny, nx), dtype=np.uint8)
    # 繪製順序：陸地與南極 (白) -> 湖泊 (黑) -> 湖中島 (白)
    layers = [(polys.get(1, []) + polys.get(5, []), 'white'), (polys.get(2, []), 'black'),
              (polys.get(3, []), 'white')]
    dpi = 100.
    for r0 in range(0, ny, STRIP):
        r1 = min(ny, r0 + STRIP)
        north, south = 90. - r0 * d, 90. - r1 * d
        fig = Figure(figsize=((nx + 1e-3) / dpi, (r1 - r0 + 1e-3) / dpi), dpi=dpi)
        fig.patch.set_facecolor('black')
        ax = fig.add_axes([0, 0, 1, 1])
        ax.set_axis_off()
        ax.set_facecolor('black')
        ax.set_xlim(-180., 180.)
        ax.set_ylim(south, north)
        for z, (items, color) in enumerate(layers):
            p = _compound(items, south - d, north + d)
            if p is not None:
                ax.add_patch(PathPatch(p, facecolor=color, edgecolor='none', antialiased=True, zorder=z + 1))
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        frac[r0:r1] = np.asarray(canvas.buffer_rgba())[:r1 - r0, :nx, 0]
        print('\r    %3d%%' % (100 * r1 // ny), end='', flush=True)
    print()
    return frac


def _write(path, frac, d):
    """分塊 DEFLATE GeoTIFF，含 2、4、8… 倍的縮圖層 (區塊平均)。"""
    import tifffile
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    tmp = path + '.tmp'
    geokeys = (1, 1, 0, 3, 1024, 0, 1, 2, 1025, 0, 1, 1, 2048, 0, 1, 4326)    # 經緯度 (WGS 84)、PixelIsArea
    with tifffile.TiffWriter(tmp) as tw:
        level, k = frac, 0
        while True:
            extra = [] if k else [(33550, 'd', 3, (d, d, 0.)), (33922, 'd', 6, (0., 0., 0., -180., 90., 0.)),
                                  (34735, 'H', len(geokeys), geokeys)]
            tw.write(level, tile=(TILE, TILE), compression='zlib', photometric='minisblack',
                     subfiletype=1 if k else 0, extratags=extra, metadata=None)
            if level.shape[1] <= 2 * TILE or level.shape[0] % 2 or level.shape[1] % 2:
                break
            level = level.reshape(level.shape[0] // 2, 2, level.shape[1] // 2, 2).mean((1, 3)).round().astype(np.uint8)
            k += 1
    os.replace(tmp, path)


# ---- 取樣 ----------------------------------------------------------------------
class LandSea:
    """依經緯度取陸地比例 (0~1)，雙線性內插。
    經緯度格點的分塊 GeoTIFF (含縮圖層) 共用：多波段 (RGB) 的檔案見 satima/bluemarble.py。"""

    LABEL = '海陸遮罩'
    LABEL_EN = 'Land/sea mask'

    def __init__(self, path, pixel_km, cache_tiles=256):
        from .truecolor.background import TiledTiff
        self.tif = TiledTiff(path, cache_tiles)
        base = self.tif.levels[0]
        base_km = self.tif.dlon * 111.32
        pick = 0          # 不比輸出像素粗的最粗一層
        for i, lv in enumerate(self.tif.levels):
            if base_km * base['width'] / lv['width'] <= pixel_km * 1.001:
                pick = i
        self.level = pick
        lv = self.tif.levels[pick]
        self.nx, self.ny, self.tw, self.th = lv['width'], lv['height'], lv['tw'], lv['th']
        self.lon0, self.lat0 = self.tif.lon0, self.tif.lat0
        self.dlon = self.tif.dlon * base['width'] / self.nx
        self.dlat = self.tif.dlat * base['height'] / self.ny
        self.bands = lv['spp']
        print(tr('  %s：%s (%dx%d，約 %.2f km)', '  %s: %s (%dx%d, about %.2f km)')
              % (tr(self.LABEL, self.LABEL_EN), os.path.basename(path), self.nx, self.ny, self.dlon * 111.32))

    def _values(self, row, col):
        """各點的值 (依所在分塊分組讀取，跨換日線也只讀用到的分塊)；多波段時多一維。"""
        out = np.empty(row.shape + ((self.bands,) if self.bands > 1 else ()), dtype=np.float32)
        across = -(-self.nx // self.tw)
        key = (row // self.th) * across + col // self.tw
        order = np.argsort(key, kind='stable')
        ks = key[order]
        bounds = np.flatnonzero(np.diff(ks)) + 1
        for a, b in zip(np.r_[0, bounds], np.r_[bounds, ks.size]):
            idx = order[a:b]
            ty, tx = divmod(int(ks[a]), across)
            tile = self.tif._tile(self.level, ty, tx)
            out[idx] = tile[row[idx] - ty * self.th, col[idx] - tx * self.tw]
        return out

    def bilinear(self, lon, lat):
        """雙線性內插的檔案值 (float32；經緯度為 NaN 處為 0)；多波段時最後多一維。"""
        lon, lat = np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)
        extra = (self.bands,) if self.bands > 1 else ()
        out = np.zeros(lon.shape + extra, dtype=np.float32)
        ok = np.isfinite(lon) & np.isfinite(lat)
        if not ok.any():
            return out
        x = np.mod(lon[ok] - self.lon0, 360.) / self.dlon - 0.5           # 像素中心
        y = np.clip((self.lat0 - lat[ok]) / self.dlat - 0.5, 0., self.ny - 1.)
        c0, r0 = np.floor(x).astype(np.int64), np.floor(y).astype(np.int64)
        fx, fy = (x - c0).astype(np.float32), (y - r0).astype(np.float32)
        if extra:
            fx, fy = fx[:, None], fy[:, None]
        r1 = np.minimum(r0 + 1, self.ny - 1)
        c0, c1 = c0 % self.nx, (c0 + 1) % self.nx
        out[ok] = ((1 - fy) * ((1 - fx) * self._values(r0, c0) + fx * self._values(r0, c1)) +
                   fy * ((1 - fx) * self._values(r1, c0) + fx * self._values(r1, c1)))
        return out

    def fraction(self, lon, lat):
        """陸地比例 0~1 (經緯度為 NaN 處為 0)。"""
        return self.bilinear(lon, lat) / 255.


def open_mask(path, pixel_km):
    """取得遮罩 (檔案不存在時先產生)；同一檔案與解析度共用。"""
    with _LOCK:
        if not os.path.isfile(path):
            build(path)
        key = (path, round(float(pixel_km), 6))
        if key not in _OPEN:
            _OPEN[key] = LandSea(path, pixel_km)
        return _OPEN[key]


# ---- 合成 ----------------------------------------------------------------------
def compose(data, lon, lat, mask, reflective, c, out, night=None):
    """data：亮溫 (°C) 或反射率 (已依 [visible] 處理)；寫入 out (列, 行, 4) uint8 RGBA。地球以外透明。
    night = ([geocolor] 設定, NightLights 或 None)：夜光底圖 (只用於紅外)，改用 GeoColor 夜間圖層的做法：
      底色 = [geocolor] land_color / water_color 依陸地比例混合，城市燈光處改為燈光的顏色；
      輸出 = N_IR x 白 + (1 - N_IR) x 底色，N_IR 為 GeoColor 高雲圖層的不透明度 ([geocolor] ir_range)。"""
    if night is not None:
        from .truecolor import geocolor
        geo, lights = night
        land = np.asarray(geo['land_color'], dtype=np.float32)
        ocean = np.asarray(geo['water_color'], dtype=np.float32)
    else:
        land = np.asarray(c.land_color, dtype=np.float32) / 255.
        ocean = np.asarray(c.ocean_color, dtype=np.float32) / 255.
    power = float(c.alpha_power)
    data = np.asarray(data, dtype=np.float32)
    rows = max(1, (1 << 21) // max(1, data.shape[-1]))
    for i0 in range(0, data.shape[0], rows):
        i1 = min(data.shape[0], i0 + rows)
        v, lo, la = data[i0:i1], lon[i0:i1], lat[i0:i1]
        f = mask.fraction(lo, la)[..., None]
        bg = ocean + f * (land - ocean)
        if night is not None:
            if lights is not None:
                geocolor.add_lights(bg, lights.sample(lo, la))
            a = geocolor.high_cloud(v.astype(np.float64) + 273.15, la, geo)[..., None]
            out[i0:i1, :, :3] = (np.clip(a + (1. - a) * bg, 0., 1.) * 255. + 0.5).astype(np.uint8)
            out[i0:i1, :, 3] = np.where(np.isfinite(lo) & np.isfinite(la), 255, 0)
            continue
        if reflective:
            r0, r1 = (float(x) for x in c.vis_alpha)
            gray = np.clip(v, 0., 1.)
            a = np.clip((v - r0) / (r1 - r0), 0., 1.)
        else:
            nodes = np.asarray(c.ir_gray, dtype=np.float64)
            order = np.argsort(nodes[:, 0])
            gray = (np.interp(v, nodes[order, 0], nodes[order, 1]) / 255.).astype(np.float32)
            t0, t1 = (float(x) for x in c.ir_alpha)
            a = np.clip((t0 - v) / (t0 - t1), 0., 1.)
        a = np.where(np.isfinite(v), a ** power, 0.).astype(np.float32)[..., None]
        gray = np.nan_to_num(gray)[..., None]
        rgb = a * gray + (1. - a) * bg
        out[i0:i1, :, :3] = (np.clip(rgb, 0., 1.) * 255. + 0.5).astype(np.uint8)
        out[i0:i1, :, 3] = np.where(np.isfinite(lo) & np.isfinite(la), 255, 0)
    return out
