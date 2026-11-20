"""地圖繪圖 (GEO / MW / PO 共用；投影 merc / cyl / lcc / fd 與底圖 basemap / cartopy 見 maps.py)

影像一律先取樣到輸出地圖的格點 (各類衛星自己的取樣程式)，這裡只負責決定圖的大小、畫海岸線與標註、存檔；
影像與畫布像素一對一，不經 matplotlib 重新取樣。圖的長寬比依地圖本身，地圖剛好填滿整張圖。
兩種決定大小的方式 ([image] mode)：
  native  原始解析度：地圖中央的像素大小 = 產品最細波段的解析度 / 取樣間距 (width_for)
          (merc、cyl 為範圍中央緯度的東西向距離；lcc 為範圍中央的東西向距離；fd 為星下點)。
          圖寬固定當成 WIDTH_INCH 英吋，所以以點數指定的字級、線寬占圖片的比例固定，不論輸出多少像素看起來都一樣。
  dpi     依 dpi 輸出 (dpi_size)：圖的英吋數照 matplotlib 預設 (rcParams 的 figure.figsize 與 subplot 邊界，
          預設 6.4 x 4.8 英吋的圖、座標軸約 4.96 x 3.70 英吋)，地圖等比例放進座標軸的範圍後裁掉白邊，
          像素數 = 英吋 x dpi。字級、線寬以實際的英吋數為準，和以前 pcolormesh 依 dpi 存檔的圖相同。
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

from . import maps
from .i18n import tr

WIDTH_INCH = 5.0
DEFAULT_WIDTH = 2000
EARTH_KM = 6371.0


def _is_ascii(text):
    return all(ord(c) < 128 for c in text)


def _distance_km(lon1, lat1, lon2, lat2):
    """大圓距離 (km，球面)。"""
    p1, p2 = np.deg2rad(lat1), np.deg2rad(lat2)
    dl = np.deg2rad(lon2 - lon1)
    h = np.sin((p2 - p1) / 2.) ** 2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2.) ** 2
    return 2. * EARTH_KM * np.arcsin(np.sqrt(h))


class Canvas:
    GRID = None             # 輸出格點的類別 (各類衛星的取樣程式提供；見 satima/<類別>/canvas.py)

    def __init__(self, view, cfg):
        self.cfg = cfg
        self.view = view
        self.map = maps.make(view, cfg)
        x0, x1, y0, y1 = self.map.extent
        self.aspect = (y1 - y0) / (x1 - x0)
        self.width_km = self._width_km()
        self._grids = {}

    def _width_km(self):
        """圖寬相當的距離 (km)，決定原始解析度的圖寬。"""
        v = self.view
        x0, x1, _, _ = self.map.extent
        if v.width_km:
            return float(v.width_km)
        latmin, latmax, lonmin, lonmax = v.georange
        lat_c = 0.5 * (latmin + latmax)
        if v.projection == 'merc':          # 麥卡托座標為公尺 (赤道)：乘上中央緯度的縮放 (同舊版)
            return (x1 - x0) * np.cos(np.deg2rad(lat_c)) / 1000.
        if v.projection == 'cyl':           # 等距圓柱座標為經度 (度)
            return (x1 - x0) * np.deg2rad(1.) * EARTH_KM * np.cos(np.deg2rad(lat_c))
        # 其他投影：在範圍中央往東移一小段，量實際的大圓距離
        lon_c = 0.5 * (lonmin + lonmax)
        xc, yc = self.map.forward(lon_c, lat_c)
        d = (x1 - x0) * 1e-4
        lon2, lat2 = self.map.inverse(float(xc) + d, float(yc))
        return (x1 - x0) / d * _distance_km(lon_c, lat_c, float(lon2), float(lat2))

    # ---- 尺寸與格點 --------------------------------------------------------
    def width_for(self, res_km, fixed_width=None, scale=1):
        """輸出寬度 (像素)：指定寬度，或原始解析度 res_km 的整數倍 (每個原始像素的大小才一致)。
        * scale：[image] size = "xN" 或執行時輸入的倍數。全圓盤與原始像素對齊時 (view.aligned) 不套用：
          原始解析度已經一格對一個原始像素，放大只是把每個原始像素複製成 N x N，
          高解析波段還會碰到長邊上限而變成非整數倍。
        * 寬度不足 [image] min_size 時 (小範圍、解析度粗)，取達到下限的最小整數倍，且至少 [image] min_scale 倍
          (預設 4：1500 -> x4 = 6000 而不是 x2 = 3000)。
        * 長邊超過 max_size 時改用放得下的最大整數倍；原始解析度也放不下時才縮小到上限。"""
        if fixed_width:
            return int(fixed_width)
        if not res_km:
            return DEFAULT_WIDTH
        img = self.cfg.image
        native = max(1, int(round(self.width_km / res_km)))
        scale = max(1, int(scale))
        if scale > 1 and getattr(self.view, 'aligned', False):
            print(tr('  全圓盤與原始像素一一對應，維持原始解析度 (不套用 x%d)。',
                     '  The full disk maps one-to-one onto the original pixels; keeping native resolution '
                     '(x%d not applied).') % scale)
            scale = 1
        k = scale
        if native < int(img.min_size):                         # 達到寬度下限的最小整數倍，且至少 min_scale 倍
            k = max(scale, -(-int(img.min_size) // native), int(img.min_scale))
        limit = int(img.max_size / max(1., self.aspect))       # 長邊 = max_size 時的寬度
        if native * k > limit:
            fit = limit // native
            if fit < 1:
                print(tr('  原始解析度需 %d 像素，超過長邊上限 %d，改為縮小輸出。',
                         '  Native resolution needs %d px, over the limit [image] max_size = %d; shrunk to fit.') % (
                    max(native, native * self.aspect), img.max_size))
                return max(1, limit)
            print(tr('  原始解析度 x%d 需 %d 像素寬，超過長邊上限 %d，改為 x%d。',
                     '  Native resolution x%d needs %d px width, over the limit [image] max_size = %d; using x%d.')
                  % (k, native * k, img.max_size, fit))
            return native * fit
        if k > scale:
            print(tr('  原始解析度只有 %d 像素寬，放大 %d 倍為 %d (寬度下限 %d，不足時至少 x%d)。',
                     '  Native width is only %d px; enlarged x%d to %d ([image] min_size = %d, at least x%d).') % (
                native, k, native * k, img.min_size, img.min_scale))
        return native * k

    def dpi_size(self, dpi):
        """dpi 模式的 (寬度像素, 圖寬英吋)：matplotlib 預設的圖與座標軸範圍，地圖等比例放進去後裁掉白邊。"""
        rc = plt.rcParams
        fw, fh = rc['figure.figsize']
        aw = fw * (rc['figure.subplot.right'] - rc['figure.subplot.left'])
        ah = fh * (rc['figure.subplot.top'] - rc['figure.subplot.bottom'])
        inch = min(aw, ah / self.aspect)
        return max(1, int(round(inch * float(dpi)))), inch

    def grid(self, width):
        if width not in self._grids:
            shape = (max(1, int(round(width * self.aspect))), width)
            self._grids[width] = self.GRID(self.map, shape, self.width_km / width)
        return self._grids[width]

    # ---- 繪圖 --------------------------------------------------------------
    def draw(self, image, width, title='', info='', lines=True, inch=None, **imshow_kwargs):
        """畫一張圖並回傳 Figure；image 為已上色的 8 位元 RGBA / RGB (products.render)，第一列為北邊
        (也可傳 2-D 陣列配合 cmap 等 imshow 參數)。lines = False 時不畫海岸線、國界與經緯線。
        inch：圖寬的英吋數 (dpi 模式，見 dpi_size)；None = WIDTH_INCH (原始解析度模式)。"""
        # 圖高取格點的整數列數 (與 grid() 相同)：存檔時畫布剛好 width x 列數 像素，影像與線條都一對一對齊
        # (以前用浮點的長寬比，畫布高度不是整數像素，matplotlib 會把影像縮放約 0.4 像素，最近鄰取樣時漏掉一列)
        # 圖的英吋 x dpi 要剛好是整數像素：浮點誤差 (例如 1691.9999999999998) 會被 Agg 截成少一列，
        # 所以多給 0.001 像素 (Agg 取整數部分，畫布正好 width x ny)
        ny = max(1, int(round(width * self.aspect)))
        dpi = width / float(inch or WIDTH_INCH)
        fig = plt.figure(figsize=((width + 1e-3) / dpi, (ny + 1e-3) / dpi))
        fig._satima_dpi = dpi
        ax = self.map.axes(fig)
        fig._satima_image = None
        if (not imshow_kwargs and isinstance(image, np.ndarray) and image.dtype == np.uint8 and image.ndim == 3
                and image.shape[2] in (3, 4)):
            # 8 位元影像在存檔時才疊上 (save 的快速路徑)：先放一張不顯示的 1x1 佔位影像，座標範圍與 imshow 完全相同
            fig._satima_im = self.map.imshow(ax, np.zeros((1, 1, image.shape[2]), np.uint8))
            fig._satima_im.set_visible(False)
            fig._satima_image = image
        else:
            self.map.imshow(ax, image, **imshow_kwargs)
        fig._satima_ax = ax
        if lines:
            self.map.lines(ax, self.cfg)
        # 座標軸填滿整張圖、範圍就是地圖範圍 (格點的長寬比與地圖只差不到一個像素，不再為「等比例」縮小座標軸)
        x0, x1, y0, y1 = self.map.extent
        ax.set_aspect('auto')
        ax.set_xlim(x0, x1)
        ax.set_ylim(y0, y1)
        if title:
            t = self.cfg.title
            ascii_ = _is_ascii(title)
            self._annotate(fig, ax, title, (0.04, 0.04), 'left', t.fontsize,
                           t.family if ascii_ else t.family_nonascii, t.weight if ascii_ else t.weight_nonascii,
                           'black', 'w', t.box_alpha, t.box_pad)
        if info:                # 是否加標註由呼叫端決定 ([annotation] draw、[[output.image]] 的 annotation)
            a = self.cfg.annotation
            if a.upper:     # 轉大寫時 µm 保持原樣 (µ.upper() 會變成希臘字母 Μ)
                info = info.upper().replace('ΜM', 'µm')
            self._annotate(fig, ax, info, (0.5, 0), 'center', a.fontsize,
                           a.family, a.weight, a.color, a.box_color, a.box_alpha, a.box_pad)
        ax.axis('off')
        return fig

    @staticmethod
    def _annotate(fig, ax, text, xy, ha, fontsize, family, weight, color, box_color, box_alpha, box_pad=4.):
        """標註框 (xy 為 axes fraction 錨點)；box_pad 為文字到框邊的留白 (點)。
        只有超出圖寬時 (例如很長的標題) 才縮小字級，避免存檔時把圖撐大。"""
        # 下方標註：文字底部貼齊圖片下緣，框的下方留白落在圖外 (同舊版外觀)
        ann = ax.annotate(text, xy=xy, xycoords='axes fraction', ha=ha, va='bottom', color=color,
                          bbox=dict(facecolor=box_color, edgecolor='none', alpha=box_alpha, pad=float(box_pad)),
                          fontsize=fontsize, family=family, weight=weight)
        # 外框要到繪製時才定位，所以用文字範圍加上外框留白來估算
        box = ann.get_window_extent(fig.canvas.get_renderer())
        pad = float(box_pad) * fig.dpi / 72.
        lo, hi = 0.01 * fig.bbox.width, 0.99 * fig.bbox.width
        if box.x0 - pad < lo or box.x1 + pad > hi:
            scale = (hi - lo) / (box.width + 2 * pad) if ha == 'center' else (hi - box.x0) / (box.width + pad)
            ann.set_fontsize(ann.get_fontsize() * scale)
        return ann

    def save(self, fig, path, fast=True):
        """存檔 (PNG 時依序嘗試，條件不符才退到下一種)：
        1. 座標軸剛好填滿整張圖、影像與畫布像素一對一 (全圓盤等)：matplotlib 只在透明畫布上畫一次線條與標註，
           影像與這層以 Agg 相同的 8 位元混色公式疊合 (不經 matplotlib 對整張影像重新取樣)。線條重疊處偶有 1~3 個色階差。
        2. 裁切範圍 (bbox_inches='tight') 就是整張圖：連同影像繪製一次，直接取繪好的影像 (與 savefig 逐像素相同)。
        3. 其他情況：matplotlib 的 savefig (繪製兩次)。
        1、2 以多核心寫出 PNG (satima/pngio.py)，不含 matplotlib 的中繼資料，背景不透明時存成 RGB。"""
        from matplotlib.colors import to_rgba
        from . import parallel, pngio
        dpi, bg = fig._satima_dpi, self.cfg.image.background
        image = fig._satima_image
        out = None
        if fast and str(path).lower().endswith('.png'):
            if (image is not None and to_rgba(bg)[3] == 1. and getattr(self.map, 'overlay_ok', False)
                    and self._fills(fig)):
                out = self._overlay(fig, dpi, image, bg)
            if out is None:
                self._restore(fig)
                out = self._render(fig, dpi, bg)
        if out is None:
            self._restore(fig)
            fig.savefig(path, dpi=dpi, facecolor=bg, bbox_inches='tight', pad_inches=0)
        else:
            pngio.write(path, out, parallel.workers(self.cfg))
        plt.close(fig)

    @staticmethod
    def _restore(fig):
        """影像放回佔位的 imshow (由 matplotlib 繪製)。"""
        if fig._satima_image is not None:
            fig._satima_im.set_data(fig._satima_image)
            fig._satima_im.set_visible(True)
            fig._satima_image = None

    @staticmethod
    def _fills(fig):
        """座標軸 (調整長寬比後) 是否剛好是整張圖。"""
        ax = fig._satima_ax
        ax.apply_aspect()
        return all(abs(v - e) < 1e-9 for v, e in zip(ax.get_position().bounds, (0., 0., 1., 1.)))

    @staticmethod
    def _draw(fig, dpi, facecolor):
        """以存檔的 dpi 繪製一次；裁切範圍與整張圖相同時回傳畫布 (RGBA 陣列)，否則 None。"""
        from matplotlib.backends.backend_agg import FigureCanvasAgg
        fig.set_dpi(dpi)
        fig.patch.set_facecolor(facecolor)
        canvas = FigureCanvasAgg(fig)
        canvas.draw()
        tight = fig.get_tightbbox(canvas.get_renderer())
        tol = 0.5 / dpi                                  # 半個像素以內
        if not all(abs(u - v) <= tol for u, v in zip(tight.extents, fig.bbox_inches.extents)):
            return None
        return np.asarray(canvas.buffer_rgba())

    def _render(self, fig, dpi, bg):
        buf = self._draw(fig, dpi, bg)
        if buf is None:
            return None
        return buf[..., :3] if (buf[..., 3] == 255).all() else buf

    def _overlay(self, fig, dpi, image, bg):
        """線條層 (透明畫布) 與影像疊合，回傳 (列, 行, 3) uint8；尺寸不合時 None。"""
        from matplotlib.colors import to_rgba
        from . import parallel
        over = self._draw(fig, dpi, (0., 0., 0., 0.))
        if over is None or over.shape[:2] != image.shape[:2]:
            return None
        back = np.array([round(v * 255) for v in to_rgba(bg)[:3]], np.uint8)
        out = np.empty(image.shape[:2] + (3,), np.uint8)
        rows = [(i, min(image.shape[0], i + 512)) for i in range(0, image.shape[0], 512)]
        parallel.run([lambda a=a, b=b: _composite(image[a:b], over[a:b], back, out[a:b]) for a, b in rows],
                     parallel.workers(self.cfg))
        return out


def _composite(image, over, back, out):
    """背景 -> 影像 (RGBA 時依 alpha) -> 線條層 (透明畫布)，混色同 matplotlib Agg (不透明目標的 plain alpha)。"""
    try:
        from ._composite_nb import composite
    except ImportError:
        composite = _composite_numpy
    composite(image, over, back, out)


def _composite_numpy(image, over, back, out):
    def blend(dst, src, alpha):
        dst, src, alpha = dst.astype(np.int64), src.astype(np.int64), alpha.astype(np.int64)
        mix = (256 * src * alpha - 255 * dst * alpha + 65280 * dst) // (65280 + alpha)
        return np.where(alpha == 0, dst, np.where(alpha == 255, src, mix))
    a = np.full(image.shape[:2] + (1,), 255) if image.shape[2] == 3 else image[..., 3:4]
    base = blend(np.broadcast_to(back, image.shape[:2] + (3,)), image[..., :3], a)
    out[...] = blend(base, over[..., :3], over[..., 3:4])
