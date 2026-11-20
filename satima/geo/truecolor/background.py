"""GeoColor 夜間地表背景：城市燈光與海陸遮罩

城市燈光：VIIRS DNB 年平均夜光 (EOG VNL v2，單位 nW/cm²/sr) 的 Cloud-Optimized GeoTIFF。
  以內建的分塊 TIFF 讀取器 (TiledTiff：numpy + zlib；JPEG 分塊用 Pillow) 讀取，不需要 GDAL
  (GDAL 與 basemap 自帶的 GEOS DLL 版本衝突，basemap 先載入時 GDAL 無法匯入)。
  依輸出像素大小選用檔內的縮圖層 (500 m、1 km、2 km…)，每一條只解壓涵蓋到的 512x512 分塊。
  找不到檔案或格式不支援時不畫燈光 (只印一次提示)。
海陸遮罩：basemap 內建的 1.25 角分 (約 2.3 km) 海陸資料 (湖泊算水)，不需另外下載；沒有 basemap 時視為陸地。
"""
import glob
import os
import struct
import threading
import zlib
from collections import OrderedDict

import numpy as np

from ...core.i18n import tr

_NOTES = set()
_LANDMASK = None


def _note(key, text):
    if key not in _NOTES:
        _NOTES.add(key)
        print(tr('  注意：', '  NOTE: ') + text)


def land_mask(lon, lat):
    """True 為陸地 (與 lon、lat 同形狀)；無法取得遮罩時回傳 None。"""
    global _LANDMASK
    if _LANDMASK is None:       # GeoColor 在平行計算前先呼叫一次載入
        try:
            from mpl_toolkits.basemap import _readlsmask
            lons, lats, mask = _readlsmask(lakes=True, resolution='i', grid=1.25)
            _LANDMASK = (mask, float(lons[0]), float(lats[0]), float(lons[1] - lons[0]))
        except Exception as err:     # 只用 cartopy 的環境
            _LANDMASK = False
            _note('landmask', tr('無法讀取 basemap 的海陸遮罩 (%s: %s)，夜間低雲一律使用陸地門檻',
                                 'cannot read the basemap land/sea mask (%s: %s); the land threshold is used for all '
                                 'night-time low clouds') % (type(err).__name__, err))
    if _LANDMASK is False:
        return None
    mask, lon0, lat0, step = _LANDMASK
    ny, nx = mask.shape
    col = np.rint((np.mod(lon - lon0, 360.)) / step).astype(np.int64) % nx
    row = np.clip(np.rint((lat - lat0) / step).astype(np.int64), 0, ny - 1)
    return mask[row, col] == 1


def find_lights(setting, home):
    """設定值 (檔案或資料夾，相對路徑以 Satima-GEO.py 所在目錄為準) -> GeoTIFF 路徑；空字串或找不到時回傳 None。"""
    setting = str(setting or '').strip()
    if not setting:
        return None
    path = os.path.normpath(os.path.join(home, os.path.expanduser(setting)))
    if os.path.isdir(path):
        files = sorted(glob.glob(os.path.join(path, '*.tif')) + glob.glob(os.path.join(path, '*.tiff')))
        path = files[-1] if files else None
    if not path or not os.path.isfile(path):
        _note('lights-missing', tr('找不到城市燈光檔 (%s)，不畫城市燈光',
                                   'city lights file not found (%s); city lights are not drawn') % setting)
        return None
    return path


_TYPES = {1: 'B', 2: 'c', 3: 'H', 4: 'I', 5: 'II', 6: 'b', 7: 'B', 8: 'h', 9: 'i', 10: 'ii', 11: 'f', 12: 'd',
          16: 'Q', 17: 'q', 18: 'Q'}
_DTYPES = {(1, 8): 'u1', (1, 16): 'u2', (1, 32): 'u4', (2, 8): 'i1', (2, 16): 'i2', (2, 32): 'i4',
           (3, 32): 'f4', (3, 64): 'f8'}


class TiledTiff:
    """分塊 (Big)TIFF 讀取器；levels[0] 為原圖，其後為縮圖層。
    單波段：未壓縮或 DEFLATE；多波段 (像素交錯，例如 RGB)：另支援 JPEG 壓縮 (以 Pillow 解碼)。
    分塊為 (th, tw) 或 (th, tw, 波段數)。"""

    def __init__(self, path, cache_tiles=256):
        self.f = open(path, 'rb')
        head = self.f.read(16)
        self.bo = {b'II': '<', b'MM': '>'}[head[:2]]
        magic = struct.unpack(self.bo + 'H', head[2:4])[0]
        self.big = magic == 43
        offset = struct.unpack(self.bo + ('Q' if self.big else 'I'), head[8:16] if self.big else head[4:8])[0]
        self.levels = []
        while offset:
            tags, offset = self._ifd(offset)
            if 322 not in tags:
                raise ValueError(tr('不是分塊 (tiled) TIFF', 'not a tiled TIFF'))
            comp = tags.get(259, (1,))[0]
            if comp not in (1, 7, 8, 32946):
                raise ValueError(tr('不支援的 TIFF 壓縮方式 %d', 'unsupported TIFF compression %d') % comp)
            spp = tags.get(277, (1,))[0]
            if spp > 1 and tags.get(284, (1,))[0] != 1:
                raise ValueError(tr('多波段 TIFF 只支援像素交錯 (PlanarConfiguration = 1)',
                                    'multi-band TIFF must be pixel-interleaved (PlanarConfiguration = 1)'))
            dt = np.dtype(self.bo + _DTYPES[(tags.get(339, (1,))[0], tags[258][0])])
            if comp == 7 and (dt != np.uint8 or spp not in (1, 3)):
                raise ValueError(tr('JPEG 壓縮只支援 8 位元的灰階或 RGB',
                                    'JPEG compression is supported only for 8-bit grayscale or RGB'))
            tables = tags.get(347)
            self.levels.append(dict(width=tags[256][0], height=tags[257][0], tw=tags[322][0], th=tags[323][0],
                                    offsets=tags[324], counts=tags[325], comp=comp, spp=spp,
                                    jpeg_tables=bytes(tables) if tables else b'',
                                    predictor=tags.get(317, (1,))[0], dtype=dt, tags=tags))
        t0 = self.levels[0]['tags']
        scale, tie = t0[33550], t0[33922]
        self.lon0, self.lat0 = tie[3] - tie[0] * scale[0], tie[4] + tie[1] * scale[1]
        self.dlon, self.dlat = scale[0], scale[1]
        nod = t0.get(42113)
        self.nodata = float(nod.strip('\x00 ')) if nod else None
        self._cache = OrderedDict()
        self._cache_n = cache_tiles
        self._lock = threading.Lock()

    def _ifd(self, offset):
        bo, f = self.bo, self.f
        f.seek(offset)
        if self.big:
            n = struct.unpack(bo + 'Q', f.read(8))[0]
            size, entry = 20, bo + 'HHQ8s'
        else:
            n = struct.unpack(bo + 'H', f.read(2))[0]
            size, entry = 12, bo + 'HHI4s'
        raw = f.read(n * size)
        nxt = struct.unpack(bo + ('Q' if self.big else 'I'), f.read(8 if self.big else 4))[0]
        tags = {}
        for i in range(n):
            tag, typ, count, val = struct.unpack(entry, raw[i * size:(i + 1) * size])
            fmt = _TYPES.get(typ)
            if fmt is None:
                continue
            nbytes = struct.calcsize(bo + fmt) * count
            if nbytes > len(val):
                f.seek(struct.unpack(bo + ('Q' if self.big else 'I'), val)[0])
                data = f.read(nbytes)
            else:
                data = val[:nbytes]
            if typ == 2:
                tags[tag] = data.decode('latin-1')
            else:
                tags[tag] = struct.unpack(bo + fmt * count, data)
        return tags, nxt

    def _tile(self, level, ty, tx):
        key = (level, ty, tx)
        with self._lock:
            if key in self._cache:
                self._cache.move_to_end(key)
                return self._cache[key]
        lv = self.levels[level]
        across = -(-lv['width'] // lv['tw'])
        k = ty * across + tx
        off, cnt = lv['offsets'][k], lv['counts'][k]
        shape = (lv['th'], lv['tw']) + ((lv['spp'],) if lv['spp'] > 1 else ())
        if cnt == 0:          # 稀疏分塊 (全為無資料)
            tile = np.zeros(shape, dtype=lv['dtype'])
        else:
            with self._lock:                 # 共用檔案指標 (平行計算時)
                self.f.seek(off)
                raw = self.f.read(cnt)
            if lv['comp'] == 7:
                tile = self._jpeg(raw, lv, shape)
            else:
                if lv['comp'] != 1:
                    raw = zlib.decompress(raw)
                tile = np.frombuffer(raw, dtype=lv['dtype']).reshape(shape).astype(lv['dtype'].newbyteorder('='))
                if lv['predictor'] == 2 and tile.dtype.kind in 'iu':     # 水平差分：以無號整數循環相加還原
                    u = tile.view(tile.dtype.str.replace('i', 'u'))
                    tile = np.cumsum(u, axis=1, dtype=u.dtype).view(tile.dtype)
        with self._lock:
            self._cache[key] = tile
            if len(self._cache) > self._cache_n:
                self._cache.popitem(last=False)
        return tile

    @staticmethod
    def _jpeg(raw, lv, shape):
        """JPEG 分塊 (完整的 JPEG；共用量化/霍夫曼表放在 JPEGTables 時先接回去) -> (th, tw[, 3]) uint8。
        YCbCr 由 Pillow 轉成 RGB。"""
        import io
        from PIL import Image
        tables = lv['jpeg_tables']
        if tables and raw[:2] == b'\xff\xd8':
            raw = tables[:-2] + raw[2:]           # 表格的 SOI...(EOI 去掉) + 分塊的 (SOI 去掉)...EOI
        with Image.open(io.BytesIO(raw)) as im:
            im = im.convert('RGB' if lv['spp'] == 3 else 'L')
            tile = np.asarray(im)
        if tile.shape != shape:                   # 不足一個分塊時補滿 (只會出現在格式不標準的檔案)
            full = np.zeros(shape, np.uint8)
            full[:tile.shape[0], :tile.shape[1]] = tile[:shape[0], :shape[1]]
            tile = full
        return tile

    def read(self, level, r0, r1, c0, c1):
        """讀取 [r0:r1, c0:c1] (縮圖層的列行)；多波段時多一維 (波段)。"""
        lv = self.levels[level]
        extra = (lv['spp'],) if lv['spp'] > 1 else ()
        out = np.empty((r1 - r0, c1 - c0) + extra, dtype=lv['dtype'].newbyteorder('='))
        th, tw = lv['th'], lv['tw']
        for ty in range(r0 // th, (r1 - 1) // th + 1):
            for tx in range(c0 // tw, (c1 - 1) // tw + 1):
                tile = self._tile(level, ty, tx)
                a0, a1 = max(r0, ty * th), min(r1, (ty + 1) * th)
                b0, b1 = max(c0, tx * tw), min(c1, (tx + 1) * tw)
                out[a0 - r0:a1 - r0, b0 - c0:b1 - c0] = tile[a0 - ty * th:a1 - ty * th, b0 - tx * tw:b1 - tx * tw]
        return out


class NightLights:
    """依經緯度取城市燈光輻射量 (nW/cm²/sr)。"""

    def __init__(self, path, pixel_km):
        self.path = path
        self.tif = TiledTiff(path)
        base = self.tif.levels[0]
        base_km = self.tif.dlon * 111.32
        # 不比輸出像素粗的最粗一層 (例如輸出 4 km 用 2 km 層；輸出 0.5 km 用原始 500 m)
        pick = 0
        for i, lv in enumerate(self.tif.levels):
            if base_km * base['width'] / lv['width'] <= pixel_km * 1.001:
                pick = i
        self.level = pick
        lv = self.tif.levels[pick]
        self.nx, self.ny = lv['width'], lv['height']
        self.lon0, self.lat0 = self.tif.lon0, self.tif.lat0
        self.dlon = self.tif.dlon * base['width'] / self.nx
        self.dlat = self.tif.dlat * base['height'] / self.ny
        print(tr('  城市燈光：%s (%dx%d，約 %.2f km)', '  City lights: %s (%dx%d, about %.2f km)')
              % (os.path.basename(path), self.nx, self.ny, self.dlon * 111.32))

    def sample(self, lon, lat):
        """最近格點的輻射量；資料範圍外與無資料為 0。"""
        out = np.zeros(lon.shape, dtype=np.float32)
        ok = np.isfinite(lon) & np.isfinite(lat)
        col = np.zeros(lon.shape, dtype=np.int64)
        row = np.zeros(lon.shape, dtype=np.int64)
        col[ok] = (np.floor(np.mod(lon[ok] - self.lon0, 360.) / self.dlon)).astype(np.int64) % self.nx
        row[ok] = np.floor((self.lat0 - lat[ok]) / self.dlat).astype(np.int64)
        ok &= (row >= 0) & (row < self.ny)
        if not ok.any():
            return out
        r0, r1 = int(row[ok].min()), int(row[ok].max()) + 1
        c0, c1 = int(col[ok].min()), int(col[ok].max()) + 1      # 跨換日線時 c0≈0、c1≈nx (整列讀取)
        data = self.tif.read(self.level, r0, r1, c0, c1)
        vals = data[row[ok] - r0, col[ok] - c0].astype(np.float32)
        out[ok] = np.where(vals > 0, vals, 0.)                     # 無資料 (-32768) 與負值為 0
        return out


def night_lights(setting, home, pixel_km):
    """建立 NightLights；找不到檔案或無法讀取時回傳 None。"""
    path = find_lights(setting, home)
    if path is None:
        return None
    try:
        return NightLights(path, pixel_km)
    except Exception as err:
        _note('lights-read', tr('城市燈光檔讀取失敗 (%s: %s)，不畫城市燈光',
                                'failed to read the city lights file (%s: %s); city lights are not drawn')
              % (type(err).__name__, err))
    return None
