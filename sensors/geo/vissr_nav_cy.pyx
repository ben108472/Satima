# cython: language_level=3, boundscheck=False, wraparound=False, cdivision=True, initializedcheck=False
"""JMA VISSR 導航的 Cython 版 ([navigation] engine = "cython" 時使用；計算與 vissr_nav.py 的 numba 版完全相同)

結構沿用舊版 gms_nav.pyx (MSC Tech. Note No.23 的 MGIVSR)：MGI100 (姿態內插)、MGI110 (軌道內插)、
MGI200 (單位向量)、MGI220 (外積單位向量)；另加精確反算 (直接解掃描角，觀測時刻反覆修正到收斂)。
以 OpenMP 平行計算。編譯：python sensors/geo/build_vissr_nav_cy.py (需要 Cython 與 C 編譯器)。
"""
from cython.parallel import parallel, prange
from libc.math cimport sin, cos, sqrt, atan, atan2, acos, floor, fabs, fmod, isfinite, M_PI, NAN
from libc.stdlib cimport malloc, free

import numpy as np


cdef inline Py_ssize_t interval(double x, const double[::1] xp) noexcept nogil:
    """區間 xp[i] <= x < xp[i+1] (超出範圍時用第一或最後一個區間外插)。"""
    cdef Py_ssize_t i, n = xp.shape[0]
    for i in range(n - 2):
        if x < xp[i + 1]:
            return i
    return n - 2


cdef inline void mgi200(double* v) noexcept nogil:
    """單位向量 (MGI200)。"""
    cdef double r = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
    v[0] /= r
    v[1] /= r
    v[2] /= r


cdef inline void mgi220(const double* a, const double* b, double* c) noexcept nogil:
    """外積的單位向量 (MGI220)。"""
    c[0] = a[1] * b[2] - a[2] * b[1]
    c[1] = a[2] * b[0] - a[0] * b[2]
    c[2] = a[0] * b[1] - a[1] * b[0]
    mgi200(c)


cdef void frame(double t, const double[::1] att_t, const double[:, ::1] att, const double[::1] orb_t,
                const double[:, ::1] orb, const double[:, :, ::1] npa, double* f) noexcept nogil:
    """MGI100 + MGI110 + 自旋座標系：時刻 t -> f[0:3] 衛星位置、f[3:6] SP、f[6:9] SX、f[9:12] SY。"""
    cdef Py_ssize_t i = interval(t, att_t), j = interval(t, orb_t), m
    cdef double w = (t - att_t[i]) / (att_t[i + 1] - att_t[i])
    cdef double v = (t - orb_t[j]) / (orb_t[j + 1] - orb_t[j])
    cdef double alp = att[i, 0] + w * (att[i + 1, 0] - att[i, 0])
    cdef double dlt = att[i, 1] + w * (att[i + 1, 1] - att[i, 1])
    cdef double bet = att[i, 2] + w * (att[i + 1, 2] - att[i, 2])
    cdef double gst = orb[j, 3] + v * (orb[j + 1, 3] - orb[j, 3])
    cdef double az = orb[j, 4] + v * (orb[j + 1, 4] - orb[j, 4])
    cdef double el = orb[j, 5] + v * (orb[j + 1, 5] - orb[j, 5])
    cdef double a1[3]
    cdef double a2[3]
    cdef double ss[3]
    cdef double w1[3]
    cdef double w2[3]
    for m in range(3):
        f[m] = orb[j, m] + v * (orb[j + 1, m] - orb[j, m])
    a1[0], a1[1], a1[2] = sin(dlt), -cos(dlt) * sin(alp), cos(dlt) * cos(alp)
    for m in range(3):
        a2[m] = npa[j, m, 0] * a1[0] + npa[j, m, 1] * a1[1] + npa[j, m, 2] * a1[2]
    f[3] = cos(gst) * a2[0] + sin(gst) * a2[1]
    f[4] = -sin(gst) * a2[0] + cos(gst) * a2[1]
    f[5] = a2[2]
    mgi200(&f[3])
    ss[0], ss[1], ss[2] = cos(el) * cos(az), cos(el) * sin(az), sin(el)
    mgi220(&f[3], ss, w1)
    mgi220(w1, &f[3], w2)
    for m in range(3):
        f[6 + m] = w1[m] * sin(bet) + w2[m] * cos(bet)
    mgi200(&f[6])
    mgi220(&f[3], &f[6], &f[9])


cdef inline double obs_time(double line, double pixel, const double[::1] p) noexcept nogil:
    return p[5] + (floor((line - 1.) / p[4]) + pixel * p[1] / (2. * M_PI)) / (p[6] * 1440.)


def forward(const double[::1] line, const double[::1] pixel, const double[::1] p, const double[:, ::1] e,
            const double[::1] att_t, const double[:, ::1] att, const double[::1] orb_t, const double[:, ::1] orb,
            const double[:, :, ::1] npa, double[::1] lon, double[::1] lat):
    """行號、像素號 -> 經緯度 (MGIVSR 正算)。p = [step, samp, line0, pixel0, sensors, t0, spin, 扁率, 赤道半徑]。"""
    cdef Py_ssize_t k, n = line.shape[0]
    cdef double f2 = (1. - p[7]) * (1. - p[7]), ea = p[8]
    cdef double L, P, ap, aq, pc, ps, qc, qs, w1x, w1y, w1z, w2x, w2y, a, b, c, disc, d, px, py, pz
    cdef double* fr                      # 每個執行緒自己的工作區：[0:12] 衛星座標系、[12:15] 視線
    with nogil, parallel():
        fr = <double*> malloc(15 * sizeof(double))
        for k in prange(n, schedule='static'):
            L = line[k]
            P = pixel[k]
            if not (isfinite(L) and isfinite(P)):
                lon[k] = NAN
                lat[k] = NAN
                continue
            frame(obs_time(L, P, p), att_t, att, orb_t, orb, npa, fr)
            ap = p[0] * (L - p[2])
            aq = p[1] * (P - p[3])
            pc = cos(ap)
            ps = sin(ap)
            qc = cos(aq)
            qs = sin(aq)
            w1x = e[0, 0] * pc + e[0, 2] * ps
            w1y = e[1, 0] * pc + e[1, 2] * ps
            w1z = e[2, 0] * pc + e[2, 2] * ps
            w2x = qc * w1x - qs * w1y
            w2y = qs * w1x + qc * w1y
            fr[12] = fr[6] * w2x + fr[9] * w2y + fr[3] * w1z
            fr[13] = fr[7] * w2x + fr[10] * w2y + fr[4] * w1z
            fr[14] = fr[8] * w2x + fr[11] * w2y + fr[5] * w1z
            mgi200(&fr[12])
            a = f2 * (fr[12] * fr[12] + fr[13] * fr[13]) + fr[14] * fr[14]
            b = f2 * (fr[0] * fr[12] + fr[1] * fr[13]) + fr[2] * fr[14]
            c = f2 * (fr[0] * fr[0] + fr[1] * fr[1] - ea * ea) + fr[2] * fr[2]
            disc = b * b - a * c
            if disc < 0.:
                lon[k] = NAN
                lat[k] = NAN
                continue
            d = (-b - sqrt(disc)) / a
            px = fr[0] + d * fr[12]
            py = fr[1] + d * fr[13]
            pz = fr[2] + d * fr[14]
            lat[k] = atan2(pz, f2 * sqrt(px * px + py * py)) * 180. / M_PI
            lon[k] = atan2(py, px) * 180. / M_PI
        free(fr)


def inverse(const double[::1] lon, const double[::1] lat, const double[::1] p, const double[:, ::1] e,
            const double[::1] att_t, const double[:, ::1] att, const double[::1] orb_t, const double[:, ::1] orb,
            const double[:, :, ::1] npa, int iters, double[::1] out_line, double[::1] out_pixel):
    """經緯度 -> 行號、像素號 (精確反算；看不到的點為 NaN)。"""
    cdef Py_ssize_t k, n = lon.shape[0]
    cdef int it
    cdef double ea = p[8], ee = 2. * p[7] - p[7] * p[7], eb = p[8] * (1. - p[7])
    cdef double rr = sqrt(e[2, 0] * e[2, 0] + e[2, 2] * e[2, 2]), base = atan2(e[2, 2], e[2, 0])
    cdef double phi, lam, en, stx, sty, stz, L, P, Ln, Pn, wx, wy, wz, ap, aq, ux, uy, ratio
    cdef double* fr                      # 每個執行緒自己的工作區：[0:12] 衛星座標系、[12:15] 視線
    with nogil, parallel():
        fr = <double*> malloc(15 * sizeof(double))
        for k in prange(n, schedule='static'):
            out_line[k] = NAN
            out_pixel[k] = NAN
            if not (isfinite(lon[k]) and isfinite(lat[k])):
                continue
            phi = lat[k] * M_PI / 180.
            lam = lon[k] * M_PI / 180.
            en = ea / sqrt(1. - ee * sin(phi) * sin(phi))
            stx = en * cos(phi) * cos(lam)
            sty = en * cos(phi) * sin(lam)
            stz = en * (1. - ee) * sin(phi)
            # 初值同 MGIVSR：行號由緯度粗估、像素號取中心
            L = p[2] - atan(sin(phi) / (6.610689 - cos(phi))) / p[0]
            P = p[3]
            for it in range(iters):
                frame(obs_time(L, P, p), att_t, att, orb_t, orb, npa, fr)
                fr[12] = stx - fr[0]
                fr[13] = sty - fr[1]
                fr[14] = stz - fr[2]
                mgi200(&fr[12])
                wx = fr[6] * fr[12] + fr[7] * fr[13] + fr[8] * fr[14]          # 自旋座標系的視線
                wy = fr[9] * fr[12] + fr[10] * fr[13] + fr[11] * fr[14]
                wz = fr[3] * fr[12] + fr[4] * fr[13] + fr[5] * fr[14]
                # 解 Rz(q) · ELMIS · (cos p, 0, sin p) = w：z 分量與 q 無關，先解 p，再由 x、y 分量解 q
                ratio = wz / rr
                if ratio > 1.:
                    ratio = 1.
                elif ratio < -1.:
                    ratio = -1.
                ap = base - acos(ratio)
                ux = e[0, 0] * cos(ap) + e[0, 2] * sin(ap)
                uy = e[1, 0] * cos(ap) + e[1, 2] * sin(ap)
                aq = fmod(atan2(wy, wx) - atan2(uy, ux) + M_PI, 2. * M_PI)
                if aq < 0.:
                    aq = aq + 2. * M_PI
                aq = aq - M_PI
                Ln = p[2] + ap / p[0]
                Pn = p[3] + aq / p[1]
                if fabs(Ln - L) < 1e-6 and fabs(Pn - P) < 1e-6:
                    L = Ln
                    P = Pn
                    break
                L = Ln
                P = Pn
            # 地球背面看不到：衛星在該點的地平面之下
            if (fr[0] - stx) * stx / (ea * ea) + (fr[1] - sty) * sty / (ea * ea) + (fr[2] - stz) * stz / (eb * eb) > 0.:
                out_line[k] = L
                out_pixel[k] = P
        free(fr)
