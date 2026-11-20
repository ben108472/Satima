"""Blue Marble 底圖 (GRAY 的第 2 層)：雲依亮溫 / 反射率半透明，疊在 NASA Blue Marble 上
(仿 NICT「Himawari-8 Real-time Web」https://himawari8.nict.go.jp/ 各波段疊 Blue Marble 的畫面)
兩種底圖，雲的合成相同 (設定都在 [bluemarble])：
  NASA ([bluemarble])  Blue Marble Next Generation 2004 (NASA GIBS)，較暗、高原有積雪，約 1 km
  NICT ([nictmap])     Blue Marble 2002 (land_shallow_topo)，即 NICT 網站所用的版本，約 1.85 km

合成：輸出 = a x 雲色 (白) + (1 - a) x Blue Marble
  紅外    a 與 10.4 µm 的 Planck 輻射量成線性：ir_alpha 的兩點 [亮溫 °C, 不透明度 0~255] 決定直線，
          冷端外插，限制在 0~255 (所有紅外波段都用同一條 亮溫 -> 不透明度 的曲線)
  可見光  a 與反射率成線性：vis_alpha 的兩點 [反射率, 不透明度 0~255] ([visible] 的處理之後)
  無資料的像素 a = 0 (只有底圖)；地球以外透明。

預設值是以 NICT 網站 (FULL_24h 模式) 的 Himawari-9 雲圖 (白色 + 透明度的 PNG，10d 層 = 5500 像素全圓盤)
對照同時次的 HSD 逐像素反推：
  B13  2026-10-20 21:00 UTC 全圓盤、2026-11-09 05:00 UTC 北半部：10.5°C 全透明、-80°C 為 230.5，
       兩點之間 (與更冷處) 與 10.4 µm 輻射量成線性 (誤差 0.4 個色階，即 8 位元量化；與亮溫線性則差到 28 個色階)
  B03  2026-11-09 05:00 UTC：不透明度 = 反射率 x 255 (2 km 的 4x4 平均，沒有太陽天頂角校正；誤差 0.3 個色階)
       = [visible] 預設值 (不校正、倍數 1、無 gamma) 的單波段灰階

底圖 (NASA 公有領域；引用時註明 NASA Earth Observatory / Goddard Space Flight Center)，第一次使用時自動下載，
存成分塊 JPEG GeoTIFF (含縮圖層，內容相同的分塊 (大洋) 共用同一份資料)；取樣時選不比輸出像素粗的最粗一層，雙線性內插。
  NASA  Blue Marble Next Generation (MODIS 2004 年合成)：NASA GIBS 的 WMTS 圖層 BlueMarble_NextGeneration
        (EPSG:4326，第 6 層 0.0088°)。GIBS 的 512x512 JPEG 圖磚原封不動當成 TIFF 的分塊 (不重新壓縮)，
        縮圖層也直接取 GIBS 較粗的層級。
  NICT  Blue Marble 2002 land_shallow_topo (陸地、淺海、地形陰影；海洋為均勻的 RGB 10, 10, 51)：
        Wikimedia Commons 的 21600x10800 JPEG (1/60°)，切成 512x512 分塊重新以 JPEG 壓縮。
        與 NICT 網站的底圖 (FULL_24h/BlueMarble) 逐像素比對，陸地與海洋的平均色都在 1~2 個色階內；
        GIBS 的 2004 版陸地暗約 25%、海洋近黑，NICT 網站用的是 2002 版。
"""
import hashlib
import io
import os
import struct
import threading
import time
import urllib.request
from concurrent.futures import ThreadPoolExecutor

import numpy as np

from ..core.i18n import tr
from .landsea import LandSea

GIBS = ('https://gibs.earthdata.nasa.gov/wmts/epsg4326/best/BlueMarble_NextGeneration/default/500m/'
        '{z}/{row}/{col}.jpeg')
LEVEL = 6                   # 下載的最細層級：GIBS EPSG:4326 第 6 層 = 0.0087890625° (約 0.98 km)，67.5 MB
DEG0 = 0.5625               # GIBS EPSG:4326 第 0 層每像素的度數 (每一層減半)
TILE = 512
# Blue Marble 2002 land_shallow_topo 21600x10800 (Commons 的檔名雖為 2048，實際為全解析度，57.6 MB)
BM2002 = 'https://upload.wikimedia.org/wikipedia/commons/9/91/Land_shallow_topo_2048.jpg'
JPEG_QUALITY = 92           # BM2002 重新壓縮的品質
WAVELENGTH = 10.4e-6        # 紅外不透明度的 Planck 輻射量所用波長 (m)
_LOCK = threading.Lock()
_OPEN = {}
_FAILED = {}


# ---- 下載與存檔 ----------------------------------------------------------------
def _fetch(url, tries=3):
    for i in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=30) as r:
                data = r.read()
            if data[:2] != b'\xff\xd8':
                raise ValueError(tr('回應不是 JPEG', 'the response is not a JPEG'))
            return data
        except Exception:
            if i == tries - 1:
                raise
            time.sleep(2. * (i + 1))


def _subsampling(jpeg):
    """JPEG 的 SOF 中第 1 個成分 (Y) 的取樣倍數 (水平, 垂直)，例如 4:2:0 為 (2, 2)。"""
    i = 2
    while i + 4 <= len(jpeg):
        if jpeg[i] != 0xFF:
            break
        marker = jpeg[i + 1]
        size = struct.unpack('>H', jpeg[i + 2:i + 4])[0]
        if marker in (0xC0, 0xC1, 0xC2):
            hv = jpeg[i + 11]
            return hv >> 4, hv & 15
        i += 2 + size
    return 2, 2


def build(path, level=LEVEL, workers=8):
    """由 NASA GIBS 下載 Blue Marble (第 level 層到第 0 層) 並存成分塊 GeoTIFF (需要網路，第 6 層實測約 6 分鐘)。"""
    t0 = time.time()
    print(tr('  第一次使用 NASA 底圖：由 NASA GIBS 下載 Blue Marble Next Generation (只需一次) ...',
             '  First use of the NASA basemap: downloading Blue Marble Next Generation from NASA GIBS (one time only) ...'))
    levels = []
    for z in range(level, -1, -1):
        res = DEG0 / 2 ** z
        w, h = int(round(360. / res)), int(round(180. / res))
        levels.append((z, w, h, -(-w // TILE), -(-h // TILE)))
    jobs = [(z, r, c) for z, _, _, nc, nr in levels for r in range(nr) for c in range(nc)]
    tiles, done, errors = {}, [0], []

    def get(job):
        if errors:                  # 已有圖磚失敗 (例如沒有網路)：其餘不再嘗試
            return
        z, r, c = job
        try:
            tiles[job] = _fetch(GIBS.format(z=z, row=r, col=c))
        except Exception as err:
            errors.append(err)
            return
        done[0] += 1
        if done[0] % 50 == 0 or done[0] == len(jobs):
            print(tr('\r    %d / %d 個圖磚', '\r    %d / %d tiles') % (done[0], len(jobs)), end='', flush=True)

    with ThreadPoolExecutor(workers) as ex:
        list(ex.map(get, jobs))
    print()
    if errors:
        raise errors[0]
    out = [(w, h, [tiles[(z, r, c)] for r in range(nr) for c in range(nc)]) for z, w, h, nc, nr in levels]
    _write(path, out, DEG0 / 2 ** level, _subsampling(tiles[jobs[0]]))
    print(tr('  Blue Marble 已存檔 %s (%dx%d，%.1f MB，%.0f 秒)', '  Blue Marble saved %s (%dx%d, %.1f MB, %.0f s)') % (
        path, levels[0][1], levels[0][2], os.path.getsize(path) / 1e6, time.time() - t0))


def _download(url, label):
    """下載整個檔案 (印出進度)；Wikimedia 不接受沒有 User-Agent 的程式，給一個一般的名稱。"""
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (compatible; satellite-imagery-basemap)'})
    chunks, got = [], 0
    with urllib.request.urlopen(req, timeout=60) as r:
        total = int(r.headers.get('Content-Length') or 0)
        while True:
            b = r.read(1 << 20)
            if not b:
                break
            chunks.append(b)
            got += len(b)
            print('\r    %s %.1f / %.1f MB' % (label, got / 1e6, total / 1e6), end='', flush=True)
    print()
    if total and got != total:
        raise IOError(tr('下載不完整 (%d / %d bytes)', 'incomplete download (%d / %d bytes)') % (got, total))
    return b''.join(chunks)


def _half(a):
    """2x2 區塊平均 (分段計算，不需整張的浮點陣列)。"""
    h, w = a.shape[0] // 2, a.shape[1] // 2
    out = np.empty((h, w) + a.shape[2:], np.uint8)
    for r0 in range(0, h, 512):
        r1 = min(h, r0 + 512)
        s = a[2 * r0:2 * r1, :2 * w].astype(np.uint16)
        out[r0:r1] = (s[0::2, 0::2] + s[1::2, 0::2] + s[0::2, 1::2] + s[1::2, 1::2] + 2) // 4
    return out


def _jpeg_tiles(a, workers=8):
    """影像切成 512x512 分塊 (邊緣以最後一列 / 行補滿) 並以 JPEG 壓縮 (4:2:0)；列優先。"""
    from PIL import Image
    h, w = a.shape[:2]
    jobs = [(r, c) for r in range(0, h, TILE) for c in range(0, w, TILE)]

    def enc(rc):
        r, c = rc
        t = a[r:r + TILE, c:c + TILE]
        if t.shape[:2] != (TILE, TILE):
            t = np.pad(t, ((0, TILE - t.shape[0]), (0, TILE - t.shape[1]), (0, 0)), mode='edge')
        buf = io.BytesIO()
        Image.fromarray(np.ascontiguousarray(t)).save(buf, 'JPEG', quality=JPEG_QUALITY, subsampling=2)
        return buf.getvalue()

    with ThreadPoolExecutor(workers) as ex:
        return list(ex.map(enc, jobs))


def build_bm2002(path, url=BM2002):
    """由 Wikimedia Commons 下載 Blue Marble 2002 (21600x10800 JPEG，57.6 MB) 並存成分塊 GeoTIFF。
    url 也可以是本機檔案 (例如被 Wikimedia 限流 (HTTP 429) 時，先用瀏覽器下載)。
    解碼時暫時需要約 1.5 GB 記憶體 (只有產生時)。"""
    from PIL import Image
    t0 = time.time()
    if os.path.isfile(url):
        print(tr('  NICT 底圖：由 %s 轉檔 ...', '  NICT basemap: converting from %s ...') % url)
        with open(url, 'rb') as f:
            data = f.read()
    else:
        print(tr('  第一次使用 NICT 底圖：下載 NASA Blue Marble 2002 (只需一次) ...',
                 '  First use of the NICT basemap: downloading NASA Blue Marble 2002 (one time only) ...'))
        data = _download(url, 'Blue Marble 2002')
    limit = Image.MAX_IMAGE_PIXELS
    Image.MAX_IMAGE_PIXELS = None            # 2.3 億像素，超過 Pillow 預設的上限
    try:
        with Image.open(io.BytesIO(data)) as im:
            a = np.asarray(im if im.mode == 'RGB' else im.convert('RGB'))
    finally:
        Image.MAX_IMAGE_PIXELS = limit
    del data
    if a.shape[1] != 2 * a.shape[0]:
        raise ValueError(tr('Blue Marble 2002 須為 2:1 的全球圖 (目前為 %dx%d)',
                            'Blue Marble 2002 must be a 2:1 global image (got %dx%d)') % (a.shape[1], a.shape[0]))
    print(tr('    切成分塊並壓縮 ...', '    Tiling and compressing ...'))
    levels, level = [], a
    while True:
        levels.append((level.shape[1], level.shape[0], _jpeg_tiles(level)))
        if level.shape[1] <= 2 * TILE or level.shape[0] % 2 or level.shape[1] % 2:
            break
        level = _half(level)
    _write(path, levels, 360. / a.shape[1], (2, 2))
    print(tr('  Blue Marble 2002 已存檔 %s (%dx%d，%.1f MB，%.0f 秒)',
             '  Blue Marble 2002 saved %s (%dx%d, %.1f MB, %.0f s)') % (
        path, a.shape[1], a.shape[0], os.path.getsize(path) / 1e6, time.time() - t0))


SOURCES = {'gibs': build, 'bm2002': build_bm2002}       # 底圖檔不存在時的產生方式

_FMT = {3: 'H', 4: 'I', 5: 'I', 12: 'd'}


def _ifd(f, entries, link):
    """寫一個 IFD (陣列值先寫在前面)，並把上一個 IFD 的「下一個 IFD」指到這裡；回傳本 IFD 的該欄位位置。"""
    packed = []
    for tag, typ, vals in entries:
        data = struct.pack('<%d%s' % (len(vals), _FMT[typ]), *vals)
        count = len(vals) // 2 if typ == 5 else len(vals)
        if len(data) > 4:
            f.write(b'\0' * (f.tell() % 2))
            off = f.tell()
            f.write(data)
            data = struct.pack('<I', off)
        packed.append((tag, typ, count, data.ljust(4, b'\0')))
    f.write(b'\0' * (f.tell() % 2))
    pos = f.tell()
    f.write(struct.pack('<H', len(packed)))
    for tag, typ, count, data in sorted(packed):
        f.write(struct.pack('<HHI', tag, typ, count) + data)
    nxt = f.tell()
    f.write(b'\0\0\0\0')
    f.seek(link)
    f.write(struct.pack('<I', pos))
    f.seek(0, 2)
    return nxt


def _write(path, levels, deg, sub):
    """levels = [(寬, 高, [JPEG 分塊，列優先])]，第 1 個為原圖，其後為縮圖層 (NewSubfileType = 1)。
    JPEG 壓縮 (7)、YCbCr、像素交錯的 RGB 分塊 TIFF；經緯度 (WGS 84)、PixelIsArea、左上角 (180°W, 90°N)。"""
    os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
    tmp = path + '.tmp'
    geokeys = (1, 1, 0, 3, 1024, 0, 1, 2, 1025, 0, 1, 1, 2048, 0, 1, 4326)
    with open(tmp, 'wb') as f:
        f.write(b'II*\x00\x00\x00\x00\x00')
        link = 4
        for k, (w, h, tiles) in enumerate(levels):
            seen, offsets, counts = {}, [], []
            for t in tiles:
                key = hashlib.sha1(t).digest()
                if key not in seen:                 # 內容相同的分塊 (大洋) 共用同一份資料
                    f.write(b'\0' * (f.tell() % 2))
                    seen[key] = f.tell()
                    f.write(t)
                offsets.append(seen[key])
                counts.append(len(t))
            entries = [(254, 4, [1 if k else 0]), (256, 4, [w]), (257, 4, [h]), (258, 3, [8, 8, 8]),
                       (259, 3, [7]), (262, 3, [6]), (277, 3, [3]), (284, 3, [1]), (322, 3, [TILE]),
                       (323, 3, [TILE]), (324, 4, offsets), (325, 4, counts), (339, 3, [1, 1, 1]),
                       (530, 3, list(sub)), (532, 5, [0, 1, 255, 1, 128, 1, 255, 1, 128, 1, 255, 1])]
            if k == 0:
                entries += [(33550, 12, [deg, deg, 0.]), (33922, 12, [0., 0., 0., -180., 90., 0.]),
                            (34735, 3, list(geokeys))]
            link = _ifd(f, entries, link)
    os.replace(tmp, path)


# ---- 取樣 ----------------------------------------------------------------------
class BlueMarble(LandSea):
    """依經緯度取 Blue Marble 的 RGB (0~1)，雙線性內插。"""

    LABEL = LABEL_EN = 'Blue Marble'

    def __init__(self, path, pixel_km):
        super().__init__(path, pixel_km, cache_tiles=96)     # RGB 分塊每塊 768 KB
        if self.bands != 3:
            raise ValueError(tr('Blue Marble 底圖須為 RGB 的分塊 GeoTIFF (%s)',
                                'the Blue Marble basemap must be an RGB tiled GeoTIFF (%s)') % path)

    def rgb(self, lon, lat):
        return self.bilinear(lon, lat) / 255.


def open_image(path, pixel_km, source='gibs'):
    """取得底圖 (檔案不存在時先依 source 下載產生：'gibs' = NASA、'bm2002' = NICT)；同一檔案與解析度共用。
    無法下載或讀取時引發 ValueError。"""
    with _LOCK:
        if path in _FAILED:
            raise ValueError(_FAILED[path])
        try:
            if not os.path.isfile(path):
                SOURCES[source](path)
            key = (path, round(float(pixel_km), 6))
            if key not in _OPEN:
                _OPEN[key] = BlueMarble(path, pixel_km)
            return _OPEN[key]
        except Exception as err:
            if isinstance(err, ValueError) and os.path.isfile(path):
                raise
            _FAILED[path] = tr('Blue Marble 底圖無法取得 (%s: %s)；可稍後再試，或把檔案放在 %s',
                               'Blue Marble basemap unavailable (%s: %s); try again later, or put the file at %s') % (
                type(err).__name__, err, path)
            raise ValueError(_FAILED[path]) from None


# ---- 合成 ----------------------------------------------------------------------
def _radiance(t_c):
    """10.4 µm 的 Planck 輻射量 (相對值；常數在兩點的比例中消去)。"""
    with np.errstate(divide='ignore', over='ignore', invalid='ignore'):
        return 1. / np.expm1(1.4387769e-2 / WAVELENGTH / (np.asarray(t_c, dtype=np.float64) + 273.15))


def opacity(v, reflective, c):
    """不透明度 0~1 (float32)：紅外 v 為亮溫 °C (與 10.4 µm 輻射量線性)，可見光 v 為反射率 (線性)；無資料為 0。"""
    v = np.asarray(v, dtype=np.float64)
    if reflective:
        (x0, a0), (x1, a1) = c.vis_alpha
        x, x0, x1 = v, float(x0), float(x1)
    else:
        (t0, a0), (t1, a1) = c.ir_alpha
        x, x0, x1 = _radiance(v), _radiance(float(t0)), _radiance(float(t1))
    with np.errstate(invalid='ignore'):
        a = np.clip(float(a0) + (float(a1) - float(a0)) * (x - x0) / (x1 - x0), 0., 255.) / 255.
    return np.where(np.isfinite(a), a, 0.).astype(np.float32)


def compose(data, lon, lat, image, reflective, c, out):
    """data：亮溫 (°C) 或反射率 (已依 [visible] 處理)；寫入 out (列, 行, 4) uint8 RGBA。地球以外透明。"""
    cloud = np.asarray(c.cloud_color, dtype=np.float32) / 255.
    data = np.asarray(data, dtype=np.float32)
    rows = max(1, (1 << 20) // max(1, data.shape[-1]))
    for i0 in range(0, data.shape[0], rows):
        i1 = min(data.shape[0], i0 + rows)
        lo, la = lon[i0:i1], lat[i0:i1]
        a = opacity(data[i0:i1], reflective, c)[..., None]
        rgb = a * cloud + (1. - a) * image.rgb(lo, la)
        out[i0:i1, :, :3] = (np.clip(rgb, 0., 1.) * 255. + 0.5).astype(np.uint8)
        out[i0:i1, :, 3] = np.where(np.isfinite(lo) & np.isfinite(la), 255, 0)
    return out
