"""掃描帶內插的 numba 核心 (nogil，可多執行緒)；沒有 numba 時 swath.py 改用 numpy 的版本。

與 swath.py 的 _lines / _combine / _hermite 是同一套算法 (逐像素的寫法)，結果相同 (只差浮點捨入)。
"""
import math

import numpy as np
from numba import njit

EPS = 1e-6


@njit(nogil=True, cache=True)
def _slope(sa, sb, ha, hb):
    if sa * sb <= 0.:
        return 0.
    central = (hb * sa + ha * sb) / (ha + hb)
    m = min(abs(central), 2. * min(abs(sa), abs(sb)))
    return m if central > 0. else -m


@njit(nogil=True, cache=True)
def _hermite(y0, y1, y2, y3, x0, x1, x2, x3, x):
    h = max(x2 - x1, EPS)
    t = min(max((x - x1) / h, 0.), 1.)
    s1 = (y2 - y1) / h
    if math.isfinite(y0):
        hl = max(x1 - x0, EPS)
        m1 = _slope((y1 - y0) / hl, s1, hl, h)
    else:
        m1 = s1
    if math.isfinite(y3):
        hr = max(x3 - x2, EPS)
        m2 = _slope(s1, (y3 - y2) / hr, h, hr)
    else:
        m2 = s1
    t2 = t * t
    t3 = t2 * t
    return ((2. * t3 - 3. * t2 + 1.) * y1 + (t3 - 2. * t2 + t) * h * m1 +
            (3. * t2 - 2. * t3) * y2 + (t3 - t2) * h * m2)


@njit(nogil=True, cache=True)
def _lines(tb, r, c, m, t, v, cnt):
    """一組掃描在像素附近的掃描線 (沿掃描方向內插後的值) 加入 t / v，回傳新的個數。"""
    ny, nx = tb.shape
    rc = min(max(r, 0.), ny - 1.)
    cc = min(max(c, 0.), nx - 1.)
    if m == 0:
        i = int(np.rint(rc))
        j = int(np.rint(cc))
        t[cnt] = i - rc
        v[cnt] = tb[i, j]
        return cnt + 1
    i0 = min(max(int(math.floor(rc)), 0), ny - 2)
    j0 = min(max(int(math.floor(cc)), 0), nx - 2)
    g = cc - j0
    lo, hi = (0, 1) if m == 1 else (-1, 2)
    for off in range(lo, hi + 1):
        i = i0 + off
        if i < 0 or i >= ny:
            continue
        a = tb[i, j0]
        b = tb[i, j0 + 1]
        if m == 1:
            val = a + (b - a) * g
        else:
            yl = tb[i, j0 - 1] if j0 >= 1 else np.nan
            yr = tb[i, j0 + 2] if j0 + 2 < nx else np.nan
            val = _hermite(yl, a, b, yr, -1., 0., 1., 2., g)
        t[cnt] = i - rc
        v[cnt] = val
        cnt += 1
    return cnt


@njit(nogil=True, cache=True)
def _combine(t, v, cnt, m):
    """沿軌方向內插：掃描線聯集中，像素之前 / 之後最近的各 1 (bilinear) 或 2 (cubic) 條。"""
    if m == 0:
        best, val = np.inf, np.nan
        for k in range(cnt):
            if math.isfinite(v[k]) and abs(t[k]) < best:
                best, val = abs(t[k]), v[k]
        return val
    x0 = x1 = -np.inf
    x2 = x3 = np.inf
    y0 = y1 = y2 = y3 = np.nan
    for k in range(cnt):
        if not math.isfinite(v[k]):
            continue
        if t[k] <= 0.:
            if t[k] > x1:
                x0, y0 = x1, y1
                x1, y1 = t[k], v[k]
            elif t[k] > x0:
                x0, y0 = t[k], v[k]
        else:
            if t[k] < x2:
                x3, y3 = x2, y2
                x2, y2 = t[k], v[k]
            elif t[k] < x3:
                x3, y3 = t[k], v[k]
    out = np.nan
    if math.isfinite(y1) and math.isfinite(y2):
        if m == 1:
            w = min(max(-x1 / max(x2 - x1, EPS), 0.), 1.)
            out = y1 + (y2 - y1) * w
        else:
            if not math.isfinite(y0):
                x0 = x1 - 1.
            if not math.isfinite(y3):
                x3 = x2 + 1.
            out = _hermite(y0, y1, y2, y3, x0, x1, x2, x3, 0.)
    if not math.isfinite(out) and x1 == 0.:
        out = y1                    # 剛好在掃描線上 (含掃描帶最後一條)：不需要另一側
    return out


@njit(nogil=True, cache=True)
def interp(tb0, r0, c0, ok0, tb1, r1, c1, ok1, nsw, order, out):
    """一或兩組掃描 (nsw) 的內插；order：0 nearest、1 bilinear、2 cubic。
    算不出來 (周圍有無效值) 的像素改用較低階的方法。"""
    t = np.empty(8)
    v = np.empty(8)
    for p in range(out.size):
        val = np.nan
        a = ok0[p]
        b = nsw > 1 and ok1[p]
        if a or b:
            m = order
            while m >= 0:
                cnt = 0
                if a:
                    cnt = _lines(tb0, r0[p], c0[p], m, t, v, cnt)
                if b:
                    cnt = _lines(tb1, r1[p], c1[p], m, t, v, cnt)
                val = _combine(t, v, cnt, m)
                if math.isfinite(val):
                    break
                m -= 1
        out[p] = val


@njit(nogil=True, cache=True)
def gauss(P, sx, sy, tb, r, c, ok, sigma, factor, nsw, cut, fill, maxw, num, den, near):
    """高斯距離加權 (一組掃描)：權重 exp(-d² / σ²)，d 為取樣點到像素的實際距離 (單位球上的弦長)。
    σ = sigma (> 0，固定) 或 factor x 該處較粗方向的取樣間距 (sigma <= 0)。結果累加到 num / den；
    near：像素附近 (fill x 較粗的取樣間距以內) 有有效的取樣點。"""
    ny, nx = tb.shape
    for p in range(r.size):
        if not ok[p]:
            continue
        rc = min(max(r[p], 0.), ny - 1.)
        cc = min(max(c[p], 0.), nx - 1.)
        i_n = int(np.rint(rc))
        j_n = int(np.rint(cc))
        i0 = min(max(int(math.floor(rc)), 0), ny - 2)
        j0 = min(max(int(math.floor(cc)), 0), nx - 2)
        fr = r[p] - i0                              # 掃描帶邊緣外半格以內：由邊緣外插
        fc = c[p] - j0
        w00 = (1. - fr) * (1. - fc)
        w01 = (1. - fr) * fc
        w10 = fr * (1. - fc)
        w11 = fr * fc
        qx = w00 * P[i0, j0, 0] + w01 * P[i0, j0 + 1, 0] + w10 * P[i0 + 1, j0, 0] + w11 * P[i0 + 1, j0 + 1, 0]
        qy = w00 * P[i0, j0, 1] + w01 * P[i0, j0 + 1, 1] + w10 * P[i0 + 1, j0, 1] + w11 * P[i0 + 1, j0 + 1, 1]
        qz = w00 * P[i0, j0, 2] + w01 * P[i0, j0 + 1, 2] + w10 * P[i0 + 1, j0, 2] + w11 * P[i0 + 1, j0 + 1, 2]
        if not math.isfinite(qx + qy + qz):         # 所在格子有無效的經緯度：以最近的取樣點為準
            qx, qy, qz = P[i_n, j_n, 0], P[i_n, j_n, 1], P[i_n, j_n, 2]
            if not math.isfinite(qx + qy + qz):
                continue
        sc = sx[i_n, j_n]
        sr = sy[i_n, j_n]
        coarse = max(sc, sr)
        sig = sigma if sigma > 0. else factor * max(sc, sr / nsw)
        reach = max(cut * sig, fill * coarse)
        nr = min(maxw, int(reach / sr) + 1)
        nc = min(maxw, int(reach / sc) + 1)
        r2 = reach * reach
        inv = 1. / (sig * sig)
        lim = (fill * coarse) ** 2
        sw = 0.
        sv = 0.
        hit = False
        for i in range(max(0, i_n - nr), min(ny, i_n + nr + 1)):
            for j in range(max(0, j_n - nc), min(nx, j_n + nc + 1)):
                val = tb[i, j]
                if not math.isfinite(val):
                    continue
                dx = P[i, j, 0] - qx
                dy = P[i, j, 1] - qy
                dz = P[i, j, 2] - qz
                d2 = dx * dx + dy * dy + dz * dz
                if not d2 <= r2:
                    continue
                if d2 <= lim:
                    hit = True
                w = math.exp(-min(d2 * inv, 700.))
                sw += w
                sv += w * val
        num[p] += sv
        den[p] += sw
        if hit:
            near[p] = True
