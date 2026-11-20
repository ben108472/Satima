"""存檔疊合的 numba 核心 (nogil，可多執行緒)；沒有 numba 時 canvas.py 改用 numpy。"""
import numpy as np
from numba import njit


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
