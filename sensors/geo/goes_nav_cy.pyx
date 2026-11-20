# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""GOES Imager (GVAR) 導航的 Cython 版 ([navigation] engine = "cython" 時使用；計算與 goes_imager.py 的 numba 版完全相同)

狀態表 S (K x 28)：[elvmax, scnmax, doff, ff, rma, pma, aebe2c, xs(3), bt(9), binv(9)]，依時間線性內插 (K = 1 時固定)。
編譯：python sensors/geo/build_vissr_nav_cy.py (同時編譯 VISSR 與 GOES 兩個模組)。
"""
from cython.parallel import parallel, prange
from libc.math cimport sin, cos, tan, sqrt, atan, atan2, asin, floor, fabs, isfinite, M_PI, NAN
from libc.stdlib cimport malloc, free

cdef enum:
    NSTATE = 28
cdef double ELVLN = 28.e-6
cdef double SCNPX = 16.e-6


cdef inline void state_at(const double[:, ::1] S, double g0, double dg, double ts, double* out) noexcept nogil:
    cdef Py_ssize_t K = S.shape[0], i, j
    cdef double f, w
    if K == 1:
        for j in range(NSTATE):
            out[j] = S[0, j]
        return
    f = (ts - g0) / dg
    i = <Py_ssize_t> floor(f)
    if i < 0:
        i = 0
    elif i > K - 2:
        i = K - 2
    w = f - i
    for j in range(NSTATE):
        out[j] = S[i, j] + w * (S[i + 1, j] - S[i, j])


cdef inline void fwd1(double L, double P, const double* s, double* lon, double* lat) noexcept nogil:
    cdef double elvmax = s[0], scnmax = s[1], doff = s[2], ff = s[3], rma = s[4], pma = s[5], k = s[6]
    cdef double x0 = s[7], x1 = s[8], x2 = s[9]
    cdef double a0, z0, al, ze, ca, sa, cz, da, dz, g0, g1, g2, h0, h1, h2, q1, q2, q3, d, hh, u0, u1, u2, d1
    a0 = elvmax - (L - 4.5) * ELVLN
    z0 = (P - 1.) * SCNPX - scnmax
    al = a0 - a0 * z0 * doff
    ze = z0 + 0.5 * a0 * a0 * doff
    ca = cos(al)
    sa = sin(al)
    cz = cos(ze)
    da = al - pma * sa * (ff / cz + tan(ze)) - rma * (1. - ca / cz)
    dz = ze + ff * rma * sa
    g0 = sin(dz)
    g1 = -cos(dz) * sin(da)
    g2 = cos(dz) * cos(da)
    h0 = s[10] * g0 + s[11] * g1 + s[12] * g2
    h1 = s[13] * g0 + s[14] * g1 + s[15] * g2
    h2 = s[16] * g0 + s[17] * g1 + s[18] * g2
    q1 = h0 * h0 + h1 * h1 + k * h2 * h2
    q2 = x0 * h0 + x1 * h1 + k * x2 * h2
    q3 = x0 * x0 + x1 * x1 + k * x2 * x2 - 1.
    d = q2 * q2 - q1 * q3
    if fabs(d) < 1e-9:
        d = 0.
    if d < 0.:
        lon[0] = NAN
        lat[0] = NAN
        return
    hh = -(q2 + sqrt(d)) / q1
    u0 = x0 + hh * h0
    u1 = x1 + hh * h1
    u2 = x2 + hh * h2
    d1 = u2 / sqrt(u0 * u0 + u1 * u1 + u2 * u2)
    lon[0] = atan2(u1, u0) * 180. / M_PI
    lat[0] = atan(k * d1 / sqrt(1. - d1 * d1)) * 180. / M_PI


cdef inline void inv1(double u0, double u1, double u2, const double* s, int iters, double* L, double* P) noexcept nogil:
    cdef double elvmax = s[0], scnmax = s[1], doff = s[2], ff = s[3], rma = s[4], pma = s[5], k = s[6]
    cdef double x0 = s[7], x1 = s[8], x2 = s[9]
    cdef double v0, v1, v2, g0, g1, g2, n, dz, da, al, ze, a0, z0
    cdef int it
    L[0] = NAN
    P[0] = NAN
    if (x0 - u0) * u0 + (x1 - u1) * u1 + k * (x2 - u2) * u2 <= 0.:     # 地球背面
        return
    v0 = u0 - x0
    v1 = u1 - x1
    v2 = u2 - x2
    g0 = s[19] * v0 + s[20] * v1 + s[21] * v2
    g1 = s[22] * v0 + s[23] * v1 + s[24] * v2
    g2 = s[25] * v0 + s[26] * v1 + s[27] * v2
    n = sqrt(g0 * g0 + g1 * g1 + g2 * g2)
    dz = asin(g0 / n)
    da = atan2(-g1, g2)
    al = da
    ze = dz
    for it in range(iters):                                            # 失準修正的反函數
        ze = dz - ff * rma * sin(al)
        al = da + pma * sin(al) * (ff / cos(ze) + tan(ze)) + rma * (1. - cos(al) / cos(ze))
    a0 = al
    z0 = ze
    for it in range(iters):                                            # doff 修正的反函數
        a0 = al / (1. - z0 * doff)
        z0 = ze - 0.5 * a0 * a0 * doff
    L[0] = 4.5 + (elvmax - a0) / ELVLN
    P[0] = 1. + (z0 + scnmax) / SCNPX


def forward(const double[::1] line, const double[::1] pixel, const double[:, ::1] S, double g0, double dg, double ta,
            double tb, double L0, double[::1] lon, double[::1] lat):
    """可見光行號、像素號 -> 經緯度。"""
    cdef Py_ssize_t i, n = line.shape[0]
    cdef double* s
    with nogil, parallel():
        s = <double*> malloc(NSTATE * sizeof(double))
        for i in prange(n, schedule='static'):
            if not (isfinite(line[i]) and isfinite(pixel[i])):
                lon[i] = NAN
                lat[i] = NAN
                continue
            state_at(S, g0, dg, ta + tb * (line[i] - L0), s)
            fwd1(line[i], pixel[i], s, &lon[i], &lat[i])
        free(s)


def inverse(const double[::1] lon, const double[::1] lat, const double[:, ::1] S, double g0, double dg, double ta,
            double tb, double L0, int iters, double[::1] line, double[::1] pixel):
    """經緯度 -> 可見光行號、像素號 (精確反算；IMC 關閉時依行號的掃描時間反覆修正；看不到的點為 NaN)。"""
    cdef Py_ssize_t i, n = lon.shape[0]
    cdef double k = 1. / (1. - (1. - 6356.7533 / 6378.137)) ** 2
    cdef double psi, lam, r, u0, u1, u2, ts, L, P
    cdef int rounds = 4 if S.shape[0] > 1 else 1, it
    cdef double* s
    with nogil, parallel():
        s = <double*> malloc(NSTATE * sizeof(double))
        for i in prange(n, schedule='static'):
            line[i] = NAN
            pixel[i] = NAN
            if not (isfinite(lon[i]) and isfinite(lat[i])):
                continue
            psi = atan(tan(lat[i] * M_PI / 180.) / k)                 # 大地緯度 -> 地心緯度
            lam = lon[i] * M_PI / 180.
            r = 1. / sqrt(cos(psi) * cos(psi) + k * sin(psi) * sin(psi))
            u0 = r * cos(psi) * cos(lam)
            u1 = r * cos(psi) * sin(lam)
            u2 = r * sin(psi)
            ts = ta
            L = NAN
            P = NAN
            for it in range(rounds):
                state_at(S, g0, dg, ts, s)
                inv1(u0, u1, u2, s, iters, &L, &P)
                if not isfinite(L):
                    break
                ts = ta + tb * (L - L0)
            line[i] = L
            pixel[i] = P
        free(s)
