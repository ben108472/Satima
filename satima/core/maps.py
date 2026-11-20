"""地圖投影與底圖：basemap 或 cartopy，由設定檔 [map] backend 決定

投影 ([map] projection；執行時詢問，直接 Enter = 設定檔的值)：
  merc  麥卡托 (預設)
  cyl   等距圓柱 (經緯度等間距)
  lcc   蘭伯特正形圓錐：標準緯線、中央經線見 [map] lcc_parallels、lcc_central_longitude
        (預設依地圖範圍自動決定)；圖框為涵蓋整個經緯度範圍的最小矩形
  fd    全圓盤 (靜止衛星投影，衛星視角)：星下點與衛星高度取自資料本身，不需輸入經緯度範圍
merc、cyl、lcc 由地圖範圍 (南 北 西 東) 決定。不論哪一種投影，衛星影像都是每個輸出像素以經緯度
經精確導航直接取樣 (resample.py)，與投影無關。

兩者介面相同，Canvas 與輸出格點不需要知道用的是哪一個：
  extent            投影座標範圍 (x0, x1, y0, y1)
  separable         經度只與 x、緯度只與 y 有關 (merc、cyl：輸出格點只需反算一列與一行)
  geos, origin      全圓盤的 (星下點經度, 高度, a, b, sweep) 與星下點的投影座標 (其他投影 geos 為 None)
  inverse(x, y)     投影座標 -> 經緯度 (度)；地球以外 (全圓盤的太空) 為 NaN
  forward(lon, lat) 經緯度 -> 投影座標
  axes(fig)         建立填滿整張圖的座標軸
  imshow(ax, image) 影像 (輸出格點，第一列為北邊) 貼滿地圖範圍，不再重新投影
  lines(ax, cfg)    海岸線、國界、經緯線
海岸線資料：basemap 為 GSHHS (c / l / i / h / f)，cartopy 為 Natural Earth (110m / 50m / 10m)；
[map] resolution 兩種寫法皆可，自動對應 (同 Satima-PO)。
"""
from dataclasses import dataclass

import numpy as np

BACKENDS = ('basemap', 'cartopy')
RESOLUTIONS = ('c', 'l', 'i', 'h', 'f', '110m', '50m', '10m')
# 投影代號 -> 選單顯示名稱
PROJECTIONS = {'merc': 'MERC', 'cyl': 'CYL', 'lcc': 'LCC', 'fd': 'FD'}
WGS84 = (6378137.0, 6356752.31414)
# merc / lcc 的地球模型：WGS84 橢球 (同 PROJ 的 ellps=WGS84，即 cartopy 的預設)。basemap 預設是半徑 6370997 m 的圓球，
# 這裡改成相同的橢球，兩種底圖的格點與影像才會逐像素相同 (衛星導航的經緯度本來就是 WGS84 大地緯度)
WGS84_ELLPS = (6378137.0, 6378137.0 * (1. - 1. / 298.257223563))
FD_MARGIN = 0.01            # 全圓盤：地球邊緣外多留的比例 (圓盤不貼齊圖框)
_TO_CARTOPY = {'c': '110m', 'l': '110m', 'i': '50m', 'h': '10m', 'f': '10m'}
_TO_BASEMAP = {'110m': 'l', '50m': 'i', '10m': 'h'}
_CACHE = {}
_GEOS = {}


@dataclass(frozen=True)
class View:
    """一張地圖的投影與範圍 (可當 dict key：相同 View 共用底圖與輸出格點)。"""
    projection: str                 # merc / cyl / lcc / fd
    georange: tuple = None          # merc / cyl / lcc：(南, 北, 西, 東) 度
    lcc: tuple = None               # lcc：(標準緯線1, 標準緯線2, 中央緯度, 中央經線)
    geos: tuple = None              # fd：(星下點經度, 衛星離地高度 m, 地球長半軸 a, 短半軸 b, sweep)
    extent: tuple = None            # fd：靜止衛星投影座標範圍 (x0, x1, y0, y1)，公尺 (= 掃描角 x 高度)
    width_km: float = None          # fd：圖寬相當的星下點距離 (km；決定原始解析度的圖寬)
    aligned: bool = False           # fd：圖框對齊原始像素 (原始解析度時輸出像素與原始像素一一對應)

    @property
    def tag(self):
        """檔名欄位 {projection}：麥卡托為空 (與舊版檔名相同)，其餘為 CYL / LCC / FD。"""
        return '' if self.projection == 'merc' else PROJECTIONS[self.projection]


def resolution(value, backend):
    res = str(value).strip().lower()
    return _TO_CARTOPY.get(res, res) if backend == 'cartopy' else _TO_BASEMAP.get(res, res)


def make(view, cfg):
    """依設定建立底圖；同範圍重複使用 (Basemap 建立時讀海岸線很慢)。"""
    backend = str(cfg.map.backend).lower()
    need_data = cfg.coastline.draw or cfg.borderline.draw          # 經緯線不需要海岸線資料
    res = resolution(cfg.map.resolution, backend) if need_data else None
    key = (backend, view, res)
    if key not in _CACHE:
        _CACHE[key] = (CartopyMap if backend == 'cartopy' else BasemapMap)(view, res)
    return _CACHE[key]


# ---- 投影參數 --------------------------------------------------------------
def lcc_params(georange, mcfg):
    """LCC 的 (標準緯線1, 標準緯線2, 中央緯度, 中央經線)。
    標準緯線預設取南北範圍的 1/6 與 5/6 處 (使整個範圍的長度誤差最小的常用取法)；中央經線預設為範圍中央。"""
    latmin, latmax, lonmin, lonmax = (float(v) for v in georange)
    par = [float(v) for v in (mcfg.get('lcc_parallels') or [])]
    if par:
        lat1, lat2 = par[0], par[-1]
    else:
        lat1, lat2 = latmin + (latmax - latmin) / 6., latmax - (latmax - latmin) / 6.
    if abs(lat1 + lat2) < 1.:
        # 兩條標準緯線對稱於赤道時圓錐退化 (proj 無法計算)：整體平移，使圓錐接近麥卡托
        shift = ((1. if lat1 + lat2 >= 0 else -1.) - (lat1 + lat2)) / 2.
        lat1, lat2 = lat1 + shift, lat2 + shift
    lon0 = mcfg.get('lcc_central_longitude', 'auto')
    lon0 = 0.5 * (lonmin + lonmax) if isinstance(lon0, str) else float(lon0)
    return lat1, lat2, 0.5 * (latmin + latmax), lon0


def fd_view(reader, scene, ch, backend='basemap'):
    """全圓盤的 View (依產品的波段 ch (Channel)；真實色用紅光)。
    該波段是 CGMS 標準靜止衛星投影 (GeosNav) 時：
      * 全圓盤觀測：圖框 = 原始影像範圍 (例如 MTSAT VIS 11000x11000)，原始解析度輸出時
        輸出像素與原始像素一一對應，等於直接畫原始影像
      * 區域觀測 (半圓盤、RSS、mesoscale 等)：圖框為地球邊緣外加 FD_MARGIN，對齊原始像素邊界
    其他導航 (GMS、GOES 舊衛星等) 以星下點與高度建立理想的靜止衛星投影 (圖框為地球邊緣外加 FD_MARGIN)，
    每個輸出像素仍以精確導航取樣。"""
    from ..geo.nav import GeosNav
    nav = reader.nav(scene, ch.name)
    if isinstance(nav, GeosNav):
        lon0, h, a, b, sweep = nav.sub_lon, nav.height, nav.a, nav.b, nav.sweep
    else:
        lon0, alt = reader.satellite_position(scene)
        # 每張影像的位置略有不同，取整數 (0.1 度、10 km) 讓多個時次共用底圖；取樣仍用各自的精確導航
        lon0, h, (a, b), sweep = round(float(lon0), 1), round(float(alt), -1) * 1000., WGS84, 'y'
    lon0 = (float(lon0) + 180.) % 360. - 180.
    aligned = isinstance(nav, GeosNav)
    if backend == 'basemap' and sweep != 'y':
        sweep, aligned = 'y', False          # basemap 不支援 sweep = x (GOES-R)：改用 y，不對齊原始像素
    limb = float(np.arcsin(a / (a + h)))     # 地球邊緣的掃描角 (弧度)
    half = limb * (1. + FD_MARGIN)
    if aligned:
        # 原始影像範圍 (像素外緣)
        fx = sorted((nav.x0 - nav.dx / 2., nav.x0 + (nav.shape[1] - 0.5) * nav.dx))
        fy = sorted((nav.y0 - nav.dy / 2., nav.y0 + (nav.shape[0] - 0.5) * nav.dy))
        if min(-fx[0], fx[1], -fy[0], fy[1]) >= 0.99 * limb:        # 全圓盤觀測 (GOES-R 的圖框略小於地球)
            (x0, x1), (y0, y1), nx = fx, fy, nav.shape[1]
        else:
            x0, x1, nx = _edges(nav.x0, nav.dx, half)
            y0, y1, _ = _edges(nav.y0, nav.dy, half)
        width_km = nx * ch.res
    else:
        x0, x1, y0, y1 = -half, half, -half, half
        width_km = 2. * half * h / 1000.
    extent = tuple(float(v) * h for v in (x0, x1, y0, y1))
    return View('fd', geos=(float(lon0), float(h), float(a), float(b), sweep), extent=extent, width_km=width_km,
                aligned=aligned)


def source_pixel_km(reader, scene, ch):
    """原始像素在星下點相當的大小 (東西, 南北) km；CGMS 標準投影 (GeosNav，像素是方的) 或無法判斷時回傳 None。

    GMS、GOES-1~15 等掃描式儀器，東西向的取樣間距通常比南北向的掃描線間距密 (GOES Imager 紅外 2.3 x 4 km、
    可見光 0.57 x 1 km；GMS-5 紅外 3.4 x 5 km)，Channel.res 只記南北向的標稱值。
    做法：影像中央附近幾個點，相鄰半個像素的經緯度投影到理想的靜止衛星投影 (座標 = 掃描角 x 高度)，
    量出每個像素的掃描角間距；掃描角間距在整張影像上固定，所以區域觀測 (星下點不在影像內) 也適用。"""
    from pyproj import Proj

    from ..geo.nav import GeosNav
    try:
        nav = reader.nav(scene, ch.name)
        if isinstance(nav, GeosNav):
            return None
        lon0, alt = reader.satellite_position(scene)
        proj = Proj(proj='geos', h=float(alt) * 1000., lon_0=float(lon0), a=WGS84[0], b=WGS84[1], sweep='y')
        ny, nx = nav.shape
        rows = ny * np.array([0.5, 0.35, 0.65, 0.5, 0.5])
        cols = nx * np.array([0.5, 0.5, 0.5, 0.35, 0.65])

        def step(dr, dc):
            xy = []
            for sign in (-1., 1.):
                lon, lat = nav.to_lonlat(rows + sign * dr, cols + sign * dc)
                with np.errstate(invalid='ignore'):
                    x, y = proj(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
                xy.append(_finite(x, y))
            d = np.hypot(xy[1][0] - xy[0][0], xy[1][1] - xy[0][1]) / 1000.
            return float(np.nanmedian(d)) if np.isfinite(d).any() else np.nan

        ew, ns = step(0., 0.5), step(0.5, 0.)
    except Exception:
        return None
    if not (np.isfinite(ew) and np.isfinite(ns) and ew > 0. and ns > 0.):
        return None
    return ew, ns


def _edges(v0, dv, half):
    """涵蓋 -half ~ half 且落在像素邊界 (v0 - dv/2 + k|dv|) 上的範圍；回傳 (下限, 上限, 像素數)。"""
    step = abs(float(dv))
    off = (float(v0) - float(dv) / 2.) % step
    lo = off + np.floor((-half - off) / step) * step
    hi = off + np.ceil((half - off) / step) * step
    return lo, hi, int(round((hi - lo) / step))


def _box(georange, n=361):
    """經緯度範圍的邊框 (密集取點，求 LCC 的投影範圍用)。"""
    latmin, latmax, lonmin, lonmax = georange
    lons, lats = np.linspace(lonmin, lonmax, n), np.linspace(latmin, latmax, n)
    lon = np.concatenate([lons, np.full(n, lonmax), lons[::-1], np.full(n, lonmin)])
    lat = np.concatenate([np.full(n, latmin), lats, np.full(n, latmax), lats[::-1]])
    return lon, lat


def _finite(lon, lat):
    """看不到的位置 (basemap 為 1e30、cartopy 為 inf) 改為 NaN。"""
    lon, lat = np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64)
    bad = ~(np.isfinite(lon) & np.isfinite(lat) & (np.abs(lon) < 1e20) & (np.abs(lat) < 1e20))
    if bad.any():
        lon, lat = np.where(bad, np.nan, lon), np.where(bad, np.nan, lat)
    return lon, lat


# ---- basemap ---------------------------------------------------------------
class BasemapMap:
    # 存檔的快速路徑 (canvas.save)：影像不經 matplotlib，直接與線條層疊合 (已驗證與 matplotlib 繪製相同)
    overlay_ok = True

    def __init__(self, view, res):
        from mpl_toolkits.basemap import Basemap
        p = view.projection
        self.view = view
        self.separable = p in ('merc', 'cyl')
        self.geos = view.geos                  # 全圓盤：(星下點經度, 高度, a, b, sweep)，resample 的快速路徑用
        if p == 'fd':
            lon0, h, a, b, _ = view.geos
            key = (view.geos, res)             # 同一衛星不同波段 (圖框不同) 共用 Basemap (讀海岸線很慢)
            if key not in _GEOS:
                _GEOS[key] = Basemap(projection='geos', lon_0=lon0, satellite_height=h, rsphere=(a, b),
                                     resolution=res)
            self.m = _GEOS[key]
            ox, oy = self.m(lon0, 0.)          # basemap 的座標原點在左下角：星下點 = 投影座標 (0, 0)
            self.origin = (float(ox), float(oy))
            x0, x1, y0, y1 = view.extent
            self.extent = (x0 + ox, x1 + ox, y0 + oy, y1 + oy)
            return
        latmin, latmax, lonmin, lonmax = view.georange
        if p == 'lcc':
            lat1, lat2, lat0, lon0 = view.lcc
            kw = dict(projection='lcc', lat_1=lat1, lat_2=lat2, lat_0=lat0, lon_0=lon0, rsphere=WGS84_ELLPS)
            tmp = Basemap(llcrnrlon=lonmin, llcrnrlat=latmin, urcrnrlon=lonmax, urcrnrlat=latmax,
                          resolution=None, **kw)
            bx, by = tmp(*_box(view.georange))
            lonmin, latmin = tmp(float(np.min(bx)), float(np.min(by)), inverse=True)
            lonmax, latmax = tmp(float(np.max(bx)), float(np.max(by)), inverse=True)
        else:
            kw = dict(projection=p)
            if p == 'merc':
                kw['rsphere'] = WGS84_ELLPS
        self.m = Basemap(llcrnrlat=latmin, urcrnrlat=latmax, llcrnrlon=lonmin, urcrnrlon=lonmax,
                         resolution=res, **kw)
        self.extent = (self.m.llcrnrx, self.m.urcrnrx, self.m.llcrnry, self.m.urcrnry)

    def inverse(self, x, y):
        lon, lat = self.m(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64), inverse=True)
        return _finite(lon, lat)

    def forward(self, lon, lat):
        x, y = self.m(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        return np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64)

    @staticmethod
    def axes(fig):
        return fig.add_axes([0, 0, 1, 1])

    def imshow(self, ax, image, **kwargs):
        if self.view.projection == 'fd':
            from matplotlib.axes import Axes
            return Axes.imshow(ax, image, origin='upper', interpolation='nearest', extent=self.extent, **kwargs)
        return self.m.imshow(image, origin='upper', interpolation='nearest', ax=ax, **kwargs)

    def lines(self, ax, cfg):
        c, b, l = cfg.coastline, cfg.borderline, cfg.latlon
        if self.view.projection == 'fd':
            # 全圓盤時 basemap 預設會沿地球邊緣畫一圈黑線：先放一個不畫線的邊界 (仍作為海岸線的裁切範圍)
            self.m.drawmapboundary(fill_color='none', linewidth=0, ax=ax)
        if c.draw:
            self.m.drawcoastlines(color=c.color, linewidth=c.linewidth, ax=ax)
        if b.draw:
            self.m.drawcountries(color=b.color, linewidth=b.linewidth, ax=ax)
        if l.draw:
            self.m.drawparallels(np.arange(-90, 90, l.step), linewidth=l.linewidth, dashes=(None, None),
                                 color=l.color, ax=ax)
            self.m.drawmeridians(np.arange(0, 360, l.step), linewidth=l.linewidth, dashes=(None, None),
                                 color=l.color, ax=ax)
        if self.view.projection == 'fd':
            # 全圓盤的圖框 (對齊原始像素) 與 basemap 的圓盤範圍略有不同：畫線後還原
            x0, x1, y0, y1 = self.extent
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)


# ---- cartopy ---------------------------------------------------------------
class CartopyMap:
    """中央經線設在地圖範圍中央 (跨換日線時經度可用 0~360，例如 165 195)。"""
    # 存檔同樣走疊合路徑 (canvas.save)：GeoAxes 的座標轉換是非仿射的，matplotlib 繪製影像時會再重新取樣一次，
    # 即使畫布與影像一樣大，也有約一半的像素偏移 1 格 (2026/11 實測 GOES-19 麥卡托 54%)。
    # 疊合時影像逐像素不變；線條與 matplotlib 繪製相同 (最多差 1 個色階)。全圓盤地球外的影像本來就是透明的
    # (背景色)，地球邊緣則依 Satima-GEO 逐像素的可見性判斷，不受 cartopy 圓盤多邊形裁切。
    overlay_ok = True

    def __init__(self, view, res):
        import cartopy.crs as ccrs
        p = view.projection
        self.view = view
        self.res = res
        self.separable = p in ('merc', 'cyl')
        self.geos = view.geos
        if p == 'fd':
            lon0, h, a, b, sweep = view.geos
            globe = ccrs.Globe(semimajor_axis=a, semiminor_axis=b, ellipse=None)
            self.crs = ccrs.Geostationary(central_longitude=lon0, satellite_height=h, sweep_axis=sweep, globe=globe)
            self.clon = 0.
            self._src = ccrs.PlateCarree(globe=globe)
            self.origin = (0., 0.)
            self.extent = tuple(view.extent)
            return
        latmin, latmax, lonmin, lonmax = (float(v) for v in view.georange)
        if p == 'lcc':
            lat1, lat2, lat0, lon0 = view.lcc
            north = lat1 + lat2 >= 0
            cutoff = min(-30., latmin - 10.) if north else max(30., latmax + 10.)
            self.clon = lon0
            self.crs = ccrs.LambertConformal(central_longitude=lon0, central_latitude=lat0,
                                             standard_parallels=(lat1, lat2), cutoff=cutoff)
        else:
            self.clon = 0.5 * (lonmin + lonmax)
            self.crs = (ccrs.Mercator if p == 'merc' else ccrs.PlateCarree)(central_longitude=self.clon)
        self._src = ccrs.PlateCarree(central_longitude=self.clon)
        if p == 'lcc':
            lon, lat = _box(view.georange)
            pts = self.crs.transform_points(self._src, lon - self.clon, lat)
            self.extent = (float(pts[:, 0].min()), float(pts[:, 0].max()),
                           float(pts[:, 1].min()), float(pts[:, 1].max()))
        else:
            pts = self.crs.transform_points(self._src, np.array([lonmin, lonmax]) - self.clon,
                                            np.array([latmin, latmax]))
            self.extent = (float(pts[0, 0]), float(pts[1, 0]), float(pts[0, 1]), float(pts[1, 1]))

    def inverse(self, x, y):
        x, y = np.broadcast_arrays(np.asarray(x, dtype=np.float64), np.asarray(y, dtype=np.float64))
        pts = self._src.transform_points(self.crs, x, y)
        return _finite(pts[..., 0] + self.clon, pts[..., 1])

    def forward(self, lon, lat):
        lon, lat = np.broadcast_arrays(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        pts = self.crs.transform_points(self._src, lon - self.clon, lat)
        return pts[..., 0], pts[..., 1]

    def axes(self, fig):
        return fig.add_axes([0, 0, 1, 1], projection=self.crs)

    def imshow(self, ax, image, **kwargs):
        from matplotlib.axes import Axes
        # 影像已經是本投影的規則格點：直接用 matplotlib 的 imshow，不經 cartopy 重新投影
        return Axes.imshow(ax, image, origin='upper', interpolation='nearest', extent=self.extent, **kwargs)

    def lines(self, ax, cfg):
        import cartopy.crs as ccrs
        import cartopy.feature as cfeature
        c, b, l = cfg.coastline, cfg.borderline, cfg.latlon
        if c.draw:
            ax.coastlines(resolution=self.res, color=c.color, linewidth=c.linewidth)
        if b.draw:
            ax.add_feature(cfeature.BORDERS.with_scale(self.res), edgecolor=b.color, linewidth=b.linewidth,
                           facecolor='none')
        if l.draw:
            ax.gridlines(crs=ccrs.PlateCarree(), draw_labels=False, color=l.color, linewidth=l.linewidth,
                         linestyle='-', xlocs=np.arange(-180., 180. + l.step / 2, l.step),
                         ylocs=np.arange(-90., 90. + l.step / 2, l.step))
        if self.view.projection == 'fd':
            # 全圓盤的圖框四角在地球之外：set_extent 會與投影範圍取交集而失敗，直接設定座標軸範圍
            x0, x1, y0, y1 = self.extent
            ax.set_xlim(x0, x1)
            ax.set_ylim(y0, y1)
        else:
            ax.set_extent(self.extent, crs=self.crs)
