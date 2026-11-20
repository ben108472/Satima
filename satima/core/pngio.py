"""多核心 PNG 寫入

PIL / matplotlib 以單核心壓縮 PNG，22000x22000 的全圓盤要 1 分鐘以上。這裡的做法同 pigz：
  1. 每列依 libpng 的啟發法自動選擇濾波 (None/Sub/Up/Average/Paeth 中，絕對值和最小者；numba 計算)
  2. 影像分成多段，各段以 zlib 平行壓縮成 raw deflate (段尾 Z_SYNC_FLUSH 對齊位元組，最後一段 Z_FINISH)，
     依序接成一個 deflate 串流；Adler-32 由各段的值合併 (zlib 的 adler32_combine)
  3. 寫出標準 PNG (IHDR / IDAT / IEND)
解出來的像素與原圖完全相同；檔案大小與 PIL 相近。沒有 numba 時濾波改用 numpy (只用 Paeth 以外的固定濾波 Up)。
"""
import struct
import zlib

import numpy as np

from . import parallel

ROWS_PER_JOB = 256
IDAT_SIZE = 1 << 23

try:
    from numba import njit
except ImportError:
    njit = None


def _make_filter():
    if njit is None:
        return None

    @njit(nogil=True, cache=True)
    def filt(img, prev, bpp, out):
        """img (n, w*c) uint8、prev 為上一列 (第一列為全 0) -> out (n, 1 + w*c)：每列濾波類型 + 濾波後位元組"""
        n, m = img.shape
        cand = np.empty((5, m), np.uint8)
        for i in range(n):
            row = img[i]
            up = prev if i == 0 else img[i - 1]
            best, best_sum = 0, -1
            for f in range(5):
                s = 0
                for x in range(m):
                    a = np.int32(row[x - bpp]) if x >= bpp else 0
                    b = np.int32(up[x])
                    c = np.int32(up[x - bpp]) if x >= bpp else 0
                    r = np.int32(row[x])
                    if f == 0:
                        v = r
                    elif f == 1:
                        v = r - a
                    elif f == 2:
                        v = r - b
                    elif f == 3:
                        v = r - ((a + b) >> 1)
                    else:
                        p = a + b - c
                        pa, pb, pc = abs(p - a), abs(p - b), abs(p - c)
                        if pa <= pb and pa <= pc:
                            pr = a
                        elif pb <= pc:
                            pr = b
                        else:
                            pr = c
                        v = r - pr
                    v &= 0xFF
                    cand[f, x] = v
                    s += v if v < 128 else 256 - v
                if best_sum < 0 or s < best_sum:
                    best, best_sum = f, s
                    out[i, 1:] = cand[f]
            out[i, 0] = best            # 找到更好的濾波時才複製候選列，最後選定的就是 out[i, 1:]
        return out
    return filt


_FILTER = _make_filter()


def _filtered(img2d, prev, bpp):
    n, m = img2d.shape
    out = np.empty((n, m + 1), np.uint8)
    if _FILTER is not None:
        return _FILTER(img2d, prev, bpp, out)
    out[:, 0] = 2                                           # Up
    up = np.vstack([prev[None, :], img2d[:-1]])
    out[:, 1:] = img2d - up
    return out


def _adler32_combine(a1, a2, len2):
    """zlib adler32_combine：a1 = adler(A)、a2 = adler(B)、len2 = len(B) -> adler(A + B)"""
    base = 65521
    rem = len2 % base
    sum1 = a1 & 0xFFFF
    sum2 = (rem * sum1) % base
    sum1 += (a2 & 0xFFFF) + base - 1
    sum2 += ((a1 >> 16) & 0xFFFF) + ((a2 >> 16) & 0xFFFF) + base - rem
    if sum1 >= base:
        sum1 -= base
    if sum1 >= base:
        sum1 -= base
    if sum2 >= (base << 1):
        sum2 -= (base << 1)
    if sum2 >= base:
        sum2 -= base
    return sum1 | (sum2 << 16)


def _chunk(tag, data):
    return struct.pack('>I', len(data)) + tag + data + struct.pack('>I', zlib.crc32(data, zlib.crc32(tag)) & 0xFFFFFFFF)


def write(path, image, workers=1, level=6):
    """image：(列, 行, 3 或 4) uint8，寫成 RGB / RGBA PNG。"""
    image = np.ascontiguousarray(image, dtype=np.uint8)
    h, w, c = image.shape
    if c not in (3, 4):
        raise ValueError('只支援 RGB / RGBA')
    flat = image.reshape(h, w * c)
    zero = np.zeros(w * c, np.uint8)
    ranges = [(i, min(h, i + ROWS_PER_JOB)) for i in range(0, h, ROWS_PER_JOB)]
    last = len(ranges) - 1

    def job(k):
        a, b = ranges[k]
        raw = _filtered(flat[a:b], flat[a - 1] if a > 0 else zero, c).tobytes()
        comp = zlib.compressobj(level, zlib.DEFLATED, -15)
        data = comp.compress(raw) + comp.flush(zlib.Z_FINISH if k == last else zlib.Z_SYNC_FLUSH)
        return data, zlib.adler32(raw), len(raw)

    parts = parallel.run([lambda k=k: job(k) for k in range(len(ranges))], workers)
    adler = 1
    for _, a, n in parts:
        adler = _adler32_combine(adler, a, n)
    with open(path, 'wb') as f:
        f.write(b'\x89PNG\r\n\x1a\n')
        f.write(_chunk(b'IHDR', struct.pack('>IIBBBBB', w, h, 8, 2 if c == 3 else 6, 0, 0, 0)))
        stream = [b'\x78\x9c'] + [p[0] for p in parts] + [struct.pack('>I', adler)]
        buf = bytearray()
        for piece in stream:
            buf += piece
            while len(buf) >= IDAT_SIZE:
                f.write(_chunk(b'IDAT', bytes(buf[:IDAT_SIZE])))
                del buf[:IDAT_SIZE]
        if buf:
            f.write(_chunk(b'IDAT', bytes(buf)))
        f.write(_chunk(b'IEND', b''))
