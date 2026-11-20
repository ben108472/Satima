# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""McIDAS 'GOES' 導航 (SMS、GOES-1~7) 的 Cython 版 ([navigation] engine = "cython" 時使用；
計算與 mcidas_goes_nav.py 的 numba 版完全相同)

常數 C (26)：numsen piclin radlin picele radele gamma gamdot rotm11 13 21 23 31 33 b11..b33 sqrt(rfact) roasin pictim tmpscl
每圈狀態 S (K x 8)：samtim、衛星 x y z、rot、cos / sin(地球自轉角)、gamma-dot x samtim (由掃描序號 k0 起)。
編譯：python sensors/geo/build_vissr_nav_cy.py
"""
from cython.parallel import prange
from libc.math cimport sin, cos, sqrt, atan, atan2, asin, floor, isfinite, M_PI, NAN

cdef double ASQ = 40683833.48
cdef double BSQ = 40410330.18
cdef double AB = 40546851.22
cdef double RAD2DEG = 180. / M_PI


cdef inline long scan_index(double L, double numsen) noexcept nogil:
    cdef double ilin = floor(L + 0.5)
    return <long> ((ilin - 1.) / numsen) + 1          # C 的整數轉換向 0 截斷，同 Java / Fortran


cdef inline Py_ssize_t clamp(long i, Py_ssize_t K) noexcept nogil:
    if i < 0:
        return 0
    if i > K - 1:
        return K - 1
    return i


cdef inline void fwd1(double L, double E, const double[::1] C, const double[:, ::1] S, long k0,
                      double* lon, double* lat) noexcept nogil:
    cdef Py_ssize_t i = clamp(scan_index(L, C[0]) - k0, S.shape[0])
    cdef double xs = S[i, 1], ys = S[i, 2], zs = S[i, 3], rot = S[i, 4], ct = S[i, 5], st = S[i, 6], gd = S[i, 7]
    cdef double ylin, yele, cl, sl, ce, se, e1, e2, eni, eli, emi, elo, emo, eno, basq, aq, bq, cq, rad, s
    cdef double x, y, z, x1, y1, a
    ylin = (L - C[1]) * C[2]
    yele = (E - C[3] + C[5] + gd) * C[4] - rot
    cl = cos(ylin)
    sl = sin(ylin)
    ce = cos(yele)
    se = sin(yele)
    e1 = C[7] * cl - C[8] * sl
    e2 = C[9] * cl - C[10] * sl
    eni = C[11] * cl - C[12] * sl
    eli = ce * e1 + se * e2
    emi = -se * e1 + ce * e2
    elo = C[13] * eli + C[16] * emi + C[19] * eni
    emo = C[14] * eli + C[17] * emi + C[20] * eni
    eno = C[15] * eli + C[18] * emi + C[21] * eni
    basq = BSQ / ASQ
    aq = basq + (1. - basq) * eno * eno
    bq = 2. * ((elo * xs + emo * ys) * basq + eno * zs)
    cq = (xs * xs + ys * ys) * basq + zs * zs - BSQ
    rad = bq * bq - 4. * aq * cq
    if rad < 1.:
        lon[0] = NAN
        lat[0] = NAN
        return
    s = -(bq + sqrt(rad)) / (2. * aq)
    x = xs + elo * s
    y = ys + emo * s
    z = zs + eno * s
    x1 = ct * x + st * y
    y1 = -st * x + ct * y
    a = atan(z / sqrt(x1 * x1 + y1 * y1))
    lat[0] = atan2(ASQ * sin(a), BSQ * cos(a)) * RAD2DEG
    lon[0] = atan2(y1, x1) * RAD2DEG


cdef inline void inv1(double lon, double lat, const double[::1] C, const double[:, ::1] S, long k0, long kstart,
                      double* Lout, double* Eout) noexcept nogil:
    cdef Py_ssize_t K = S.shape[0], i, j
    cdef int it
    cdef double ylat, ylon, sn, cs, tn, r, xe, ye, ze, xs, ys, zs, rot, ct, st, gd, x, y, vx, vy, vz, n
    cdef double elo, emo, eno, eli, emi, eni, q, ylin, cl, sl, e1, e2, ang, L = NAN, E = NAN
    ylat = atan2(BSQ * sin(lat / RAD2DEG), ASQ * cos(lat / RAD2DEG))
    ylon = lon / RAD2DEG
    sn = sin(ylat)
    cs = cos(ylat)
    tn = (sn / cs) * (sn / cs)
    r = AB * sqrt((1. + tn) / (BSQ + ASQ * tn))
    xe = r * cs * cos(ylon)
    ye = r * cs * sin(ylon)
    ze = r * sn
    i = clamp(kstart - k0, K)
    for it in range(6):
        xs = S[i, 1]
        ys = S[i, 2]
        zs = S[i, 3]
        rot = S[i, 4]
        ct = S[i, 5]
        st = S[i, 6]
        gd = S[i, 7]
        x = ct * xe - st * ye
        y = st * xe + ct * ye
        vx = x - xs
        vy = y - ys
        vz = ze - zs
        if vx * x / ASQ + vy * y / ASQ + vz * ze / BSQ >= 0.:
            Lout[0] = NAN
            Eout[0] = NAN
            return
        n = sqrt(vx * vx + vy * vy + vz * vz)
        elo = vx / n
        emo = vy / n
        eno = vz / n
        eli = C[13] * elo + C[14] * emo + C[15] * eno
        emi = C[16] * elo + C[17] * emo + C[18] * eno
        eni = C[19] * elo + C[20] * emo + C[21] * eno
        q = eni / C[22]
        if q > 1. or q < -1.:
            Lout[0] = NAN
            Eout[0] = NAN
            return
        ylin = C[23] - asin(q)
        cl = cos(ylin)
        sl = sin(ylin)
        e1 = C[7] * cl - C[8] * sl
        e2 = C[9] * cl - C[10] * sl
        ang = atan2(e2, e1) - atan2(emi, eli) + rot - (C[5] + gd) * C[4]
        ang = ang - 2. * M_PI * floor((ang + M_PI) / (2. * M_PI))
        L = ylin / C[2] + C[1]
        E = ang / C[4] + C[3]
        j = clamp(scan_index(L, C[0]) - k0, K)
        if j == i:
            break
        i = j
    Lout[0] = L
    Eout[0] = E


def forward(const double[::1] line, const double[::1] ele, const double[::1] C, const double[:, ::1] S, long k0,
            double[::1] lon, double[::1] lat):
    cdef Py_ssize_t i, n = line.shape[0]
    for i in prange(n, nogil=True, schedule='static'):
        if not (isfinite(line[i]) and isfinite(ele[i])):
            lon[i] = NAN
            lat[i] = NAN
        else:
            fwd1(line[i], ele[i], C, S, k0, &lon[i], &lat[i])


def inverse(const double[::1] lon, const double[::1] lat, const double[::1] C, const double[:, ::1] S, long k0,
            long kstart, double[::1] line, double[::1] ele):
    cdef Py_ssize_t i, n = lon.shape[0]
    for i in prange(n, nogil=True, schedule='static'):
        if not (isfinite(lon[i]) and isfinite(lat[i])):
            line[i] = NAN
            ele[i] = NAN
        else:
            inv1(lon[i], lat[i], C, S, k0, kstart, &line[i], &ele[i])
