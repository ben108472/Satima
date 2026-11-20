"""resample 的 numba 核心：區塊平均與取值 (nogil，各條可由多個執行緒同時計算)

與 numpy 版的結果相同 (區塊平均只差 float32 加總順序)；沒有 numba 時 import 失敗，resample 自動改用 numpy。
"""
import math

import numpy as np
from numba import njit

F32 = np.float32


@njit(nogil=True, cache=True)
def block_mean(arr, k):
    """k x k 區塊平均，忽略 NaN (全為 NaN 的區塊為 NaN)；不足一個區塊的邊緣捨去 (同 numpy 版)。"""
    ny, nx = arr.shape[0] // k, arr.shape[1] // k
    out = np.empty((ny, nx), F32)
    for by in range(ny):
        for bx in range(nx):
            s = F32(0.)
            c = 0
            for y in range(by * k, by * k + k):
                for x in range(bx * k, bx * k + k):
                    v = arr[y, x]
                    if not math.isnan(v) and not math.isinf(v):
                        s += v
                        c += 1
            out[by, bx] = s / F32(c) if c > 0 else F32(np.nan)
    return out


@njit(nogil=True, cache=True)
def geos_rowcol(tx, ty, r1, c1, rg, ab, sweep_x, missing, rows, cols):
    """全圓盤快速路徑：每列的 tan(掃描角 y)、每行的 tan(掃描角 x) 與最近索引 -> 逐像素可見性與行列號
    (可見性判斷同 resample._geos_visible，看不到的位置行列號都設為 missing)。"""
    rg2 = rg * rg
    for i in range(ty.size):
        for j in range(tx.size):
            if sweep_x:
                vz = ty[i]
                vy = tx[j] * math.sqrt(1. + vz * vz)
            else:
                vy = tx[j]
                vz = ty[i] * math.sqrt(1. + vy * vy)
            q = vy * vy + (vz * ab) ** 2 + 1.
            if rg2 >= q * (rg2 - 1.) and r1[i] != missing and c1[j] != missing:
                rows[i, j] = r1[i]
                cols[i, j] = c1[j]
            else:
                rows[i, j] = missing
                cols[i, j] = missing


@njit(nogil=True, cache=True)
def geos_inverse(sx, sy, lon0, rg, ab, sweep_x, lon, lat):
    """靜止衛星投影的反算 (同 PROJ geos 的橢球公式)：掃描角 (弧度) -> 經緯度 (度)，看不到地球為 NaN。
    rg = 1 + h / a、ab = a / b。"""
    radius_p_inv2 = ab * ab
    c = rg * rg - 1.
    for i in range(sx.size):
        if sweep_x:
            vz = math.tan(sy[i])
            vy = math.tan(sx[i]) * math.sqrt(1. + vz * vz)
        else:
            vy = math.tan(sx[i])
            vz = math.tan(sy[i]) * math.sqrt(1. + vy * vy)
        a = vy * vy + (vz * ab) ** 2 + 1.
        det = 4. * rg * rg - 4. * a * c
        if det < 0.:
            lon[i] = np.nan
            lat[i] = np.nan
            continue
        k = (2. * rg - math.sqrt(det)) / (2. * a)
        vx = rg - k
        vy *= k
        vz *= k
        lam = math.atan2(vy, vx)
        phi = math.atan(vz * math.cos(lam) / vx)
        phi = math.atan(radius_p_inv2 * math.tan(phi))
        lon[i] = math.degrees(lam) + lon0
        lat[i] = math.degrees(phi)


@njit(nogil=True, cache=True)
def bounds(rows, cols, missing):
    """有效索引的 (列最小, 列最大, 行最小, 行最大, 個數)。"""
    rmin, rmax, cmin, cmax, n = 0, -1, 0, -1, 0
    for i in range(rows.size):
        r = rows[i]
        if r == missing:
            continue
        c = cols[i]
        if n == 0:
            rmin = rmax = r
            cmin = cmax = c
        else:
            rmin, rmax = min(rmin, r), max(rmax, r)
            cmin, cmax = min(cmin, c), max(cmax, c)
        n += 1
    return rmin, rmax, cmin, cmax, n


@njit(nogil=True, cache=True)
def _blend(dst, src, alpha):
    """matplotlib Agg (fixed_blender_rgba_plain) 對不透明目標的混色：alpha 255 直接複製、0 不變"""
    if alpha == 0:
        return dst
    if alpha == 255:
        return src
    return (256 * src * alpha - 255 * dst * alpha + 65280 * dst) // (65280 + alpha)


@njit(nogil=True, cache=True)
def composite(image, over, back, out):
    """背景 back (3,) -> 影像 image (RGB 或 RGBA) -> 線條層 over (RGBA，透明畫布上畫的線條與標註)，結果寫入 out (RGB)"""
    h, w, c = image.shape
    for i in range(h):
        for j in range(w):
            a = 255 if c == 3 else np.int64(image[i, j, 3])
            oa = np.int64(over[i, j, 3])
            for k in range(3):
                v = _blend(np.int64(back[k]), np.int64(image[i, j, k]), a)
                out[i, j, k] = _blend(v, np.int64(over[i, j, k]), oa)


@njit(nogil=True, cache=True)
def gather(data, rows, cols, r0, c0, factor, out):
    """out[i] = data[(rows[i] - r0) // factor, (cols[i] - c0) // factor]；範圍外 (含 MISSING) 不變 (NaN)。"""
    ny, nx = data.shape
    for i in range(rows.size):
        br = (rows[i] - r0) // factor
        bc = (cols[i] - c0) // factor
        if br >= 0 and bc >= 0 and br < ny and bc < nx:
            out[i] = data[br, bc]


@njit(nogil=True, cache=True)
def _pair(b, q, n):
    """最近索引 b 與量化偏移 q (偏移 = q / 254，-0.5~0.5 像素) -> (上/左索引, 下/右的權重, 下/右索引)；
    偏移為 0 或鄰點超出範圍時 (b, 0, b)。"""
    if q > 0:
        a, t = b, q / 254.
    elif q < 0:
        a, t = b - 1, 1. + q / 254.
    else:
        return b, 0., b
    if a < 0 or a + 1 >= n:
        return b, 0., b
    return a, t, a + 1


@njit(nogil=True, cache=True)
def bilinear(data, rows, cols, fr, fc, r0, c0, out):
    """雙線性內插 (factor = 1)：最近像素 (rows, cols) 加上量化偏移 (fr, fc)；
    四個鄰點有 NaN 或在範圍外時用最近像素的值 (與最近鄰相同)，偏移為 0 時也直接取最近像素。"""
    ny, nx = data.shape
    for i in range(rows.size):
        br = rows[i] - r0
        bc = cols[i] - c0
        if br < 0 or bc < 0 or br >= ny or bc >= nx:
            continue
        v = data[br, bc]
        ra, t, rb = _pair(br, fr[i], ny)
        ca, u, cb = _pair(bc, fc[i], nx)
        if t == 0. and u == 0.:
            out[i] = v
            continue
        a, b, c, d = data[ra, ca], data[ra, cb], data[rb, ca], data[rb, cb]
        if not (math.isfinite(a) and math.isfinite(b) and math.isfinite(c) and math.isfinite(d)):
            out[i] = v
            continue
        out[i] = (1. - t) * ((1. - u) * a + u * b) + t * ((1. - u) * c + u * d)
