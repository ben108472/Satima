"""真實色逐像素合成的 numba 版本 ([truecolor] engine = "numba"，預設)

與 composite._compose (numpy 版) 是同一套計算，只是融合成單一迴圈：
  太陽與衛星觀測幾何 (pyorbital 的公式)、太陽天頂角校正、Rayleigh 路徑反射 (pyspectral 查表的多線性內插)、
  雲頂削減、穿透率與球面反照率、綠光、紅光銳化、晨昏線淡出、色調曲線 (CIRA 對數拉伸或 enhanced)、飽和度。
沒有大量的暫存陣列，資料留在 CPU 快取，單核心就快數倍；nogil = True，所以各條可由多個執行緒同時計算
(satima/parallel.py)，numpy 版受記憶體頻寬與 GIL 限制，多核心時只快約 2 倍。
與 numpy 版的差異只來自浮點運算順序 (色階差 0~1，見 README)。
"""
import math

import numpy as np

try:
    from numba import njit
    HAVE_NUMBA = True
except ImportError:          # 沒有 numba 時改用 numpy 版
    HAVE_NUMBA = False

    def njit(*args, **kwargs):
        return lambda f: f

from .correction import OZONE_HEIGHT

F32 = np.float32
GREEN_REAL, GREEN_AHI, GREEN_LINEAR = 0, 1, 2
# 高雲的球面反照率 S(tau x 削減係數)：削減係數 0.3~1 之間每 0.005 一點的查表，線性內插 (S 幾乎是 tau 的線性函數)
SPH_F0, SPH_N = 0.3, 141
SPH_INV = (SPH_N - 1) / (1. - SPH_F0)


def prepare(tc):
    """由 TrueColor 物件整理核心需要的常數 (每個產品一次)。"""
    from pyorbital import astronomy as pa
    from pyspectral.rayleigh import _get_wavelength_adjusted_lut_rayleigh_reflectance, get_reflectance_lut_from_file

    from .composite import (AHI_GREEN, AHI_GREEN_FLOOR, AHI_GREEN_MIN_BR, LOG_RANGE, LUMA, REFL_RANGE, VIBRANCE_HIGH,
                            VIBRANCE_SPAN)
    from .correction import LUT_MAX_WAVELENGTH, rayleigh_optical_depth, rayleigh_spherical_albedo

    t = tc.scene.time
    gmst = float(pa.gmst(t))
    ra, dec = (float(v) for v in pa.sun_ra_dec(t))
    (sx, sy, sz), _ = pa.observer_position(t, tc.sat_lon, 0., tc.sat_alt)
    ray = tc.rayleigh
    fname = ray._rayleigh.reflectance_lut_filename
    azid, satsec, sunsec = get_reflectance_lut_from_file(fname)
    smin = np.array([sunsec[0], azid[0], satsec[0]], F32)
    smax = np.array([sunsec[-1], azid[-1], satsec[-1]], F32)
    orders = np.array([len(sunsec), len(azid), len(satsec)], np.int64)
    clip_sun = F32(np.nan_to_num(np.rad2deg(np.arccos(F32(1.) / sunsec.astype(F32).max()))))
    clip_sat = F32(np.nan_to_num(np.rad2deg(np.arccos(F32(1.) / satsec.astype(F32).max()))))

    roles = ('blue', 'green', 'red', 'nir')
    luts = np.zeros((4, int(orders.prod())), F32)
    tau = np.zeros(4, F32)
    sph = np.zeros(4, F32)
    scale = np.ones(4, F32)          # 近紅外：查表範圍外，依光學厚度比例由 0.78 µm 換算
    sphlut = np.zeros((4, SPH_N), F32)
    o3 = np.zeros(4, np.float64)     # 臭氧量 x 吸收係數 (0 = 不校正)
    for i, role in enumerate(roles):
        if role not in tc.wave:
            continue
        wvl = ray.effective_wavelength(*tc.wave[role])
        lut_wvl = wvl
        if wvl > LUT_MAX_WAVELENGTH:
            lut_wvl = 0.78
            scale[i] = rayleigh_optical_depth(wvl) / rayleigh_optical_depth(0.78)
        luts[i] = _get_wavelength_adjusted_lut_rayleigh_reflectance(fname, lut_wvl * 1000.).astype(F32).ravel()
        tau[i] = rayleigh_optical_depth(wvl)
        sph[i] = rayleigh_spherical_albedo(float(tau[i]))
        sphlut[i] = rayleigh_spherical_albedo(float(tau[i]) * np.linspace(SPH_F0, 1., SPH_N))
        o3[i] = tc.ozone.get(role, 0.)

    o = tc.opts
    green = o['abi_green']
    if tc.recipe.get('green'):
        mode, lin = GREEN_REAL, (0., 0., 0.)
    elif isinstance(green, str):
        mode, lin = GREEN_AHI, (0., 0., 0.)
    else:
        mode, lin = GREEN_LINEAR, tuple(float(v) for v in green)
    fade0, fade1 = float(o['terminator_fade_start']), float(o['terminator_fade_end'])
    return dict(
        geo=np.array([gmst, ra, dec, sx, sy, sz], np.float64),
        luts=luts, smin=smin, smax=smax, orders=orders, clip=np.array([clip_sun, clip_sat], F32),
        tau=tau, sph=sph, scale=scale, sphlut=sphlut, o3=o3,
        par=np.array([o['sza_limit'], o['sza_max'], fade0, fade1, o['hybrid_green_fraction'],
                      lin[0], lin[1], lin[2], AHI_GREEN[0], AHI_GREEN[1], AHI_GREEN[2], AHI_GREEN_FLOOR,
                      AHI_GREEN_MIN_BR, REFL_RANGE[0], REFL_RANGE[1], LOG_RANGE[0], LOG_RANGE[1],
                      o.get('white_point', 1.25), o.get('shoulder', 0.95), o.get('shadow_toe', 0.006),
                      o.get('saturation', 1.), VIBRANCE_HIGH, VIBRANCE_SPAN, LUMA[0], LUMA[1], LUMA[2]], np.float64),
        flags=np.array([mode, 1 if ray.mode == 'full' else 0, 1 if fade1 > fade0 else 0,
                        1 if str(o.get('tone', 'cira')).lower() == 'enhanced' else 0,
                        1 if tc.cloud_tau else 0], np.int64))


def compose(values, lon, lat, k):
    """與 composite._compose 相同的介面與結果：回傳 (n, 3) float32 RGB。"""
    n = lon.size
    empty = np.zeros(0, F32)
    out = np.empty((n, 3), F32)
    _compose_kernel(values['blue'], values.get('green', empty), values['red_lo'], values['nir'],
                    values.get('red_hi', empty), values.get('bt', empty),
                    np.ascontiguousarray(lon, np.float64), np.ascontiguousarray(lat, np.float64),
                    k['geo'], k['luts'], k['smin'], k['smax'], k['orders'], k['clip'], k['tau'], k['sph'],
                    k['scale'], k['sphlut'], k['o3'], k['par'], k['flags'], out)
    return out


@njit(nogil=True, cache=True)
def _interp(lut, smin, smax, orders, s0, s1, s2):
    """pyspectral (geotiepoints MultilinearInterpolator) 的 3 維規則格點多線性內插 (float32)"""
    o0, o1, o2 = orders[0], orders[1], orders[2]
    n0 = (s0 - smin[0]) / (smax[0] - smin[0])
    n1 = (s1 - smin[1]) / (smax[1] - smin[1])
    n2 = (s2 - smin[2]) / (smax[2] - smin[2])
    q0 = max(min(int(n0 * F32(o0 - 1)), o0 - 2), 0)
    q1 = max(min(int(n1 * F32(o1 - 1)), o1 - 2), 0)
    q2 = max(min(int(n2 * F32(o2 - 1)), o2 - 2), 0)
    l0 = n0 * F32(o0 - 1) - F32(q0)
    l1 = n1 * F32(o1 - 1) - F32(q1)
    l2 = n2 * F32(o2 - 1) - F32(q2)
    m0, m1 = o1 * o2, o2
    b = m0 * q0 + m1 * q1 + q2
    v000, v001 = lut[b], lut[b + 1]
    v010, v011 = lut[b + m1], lut[b + m1 + 1]
    v100, v101 = lut[b + m0], lut[b + m0 + 1]
    v110, v111 = lut[b + m0 + m1], lut[b + m0 + m1 + 1]
    one = F32(1.)
    return ((one - l0) * ((one - l1) * ((one - l2) * v000 + l2 * v001) + l1 * ((one - l2) * v010 + l2 * v011)) +
            l0 * ((one - l1) * ((one - l2) * v100 + l2 * v101) + l1 * ((one - l2) * v110 + l2 * v111)))


@njit(nogil=True, cache=True)
def _path(lut, smin, smax, orders, sunsec, azi, satsec, red_pct, use_red):
    """Rayleigh 路徑反射 (0~1)：查表 x 100、(紅光削減)、截斷 0~100、/ 100 (同 pyspectral 與 correction.path_reflectance)"""
    res = _interp(lut, smin, smax, orders, sunsec, azi, satsec) * F32(100.)
    if use_red and not (red_pct < F32(20.)):
        res = (F32(1.) - (red_pct - F32(20.)) / F32(80.)) * res
    res = min(max(res, F32(0.)), F32(100.))
    return res / F32(100.)


@njit(nogil=True, cache=True)
def _trans(mu, tau):
    return ((F32(2. / 3.) + mu) + (F32(2. / 3.) - mu) * F32(math.exp(-tau / mu))) / (F32(4. / 3.) + tau)


@njit(nogil=True, cache=True)
def _sunz_corr(x, cosz, limit, smax):
    """correction.sunzen_corr_cos (float64 校正係數)"""
    if math.isnan(cosz):
        return F32(0.) * x
    limit_rad = math.radians(limit)
    limit_cos = math.cos(limit_rad)
    if cosz > limit_cos:
        corr = 1. / cosz
    else:
        grad = (math.acos(cosz) - limit_rad) / (math.radians(smax) - limit_rad)
        grad = 1. - math.log(grad + 1.) / math.log(2.)
        corr = max(grad, 0.) / limit_cos
    return F32(np.float64(x) * corr)


@njit(nogil=True, cache=True)
def _fade(sunz, start, end):
    """correction.sunzen_reduction (strength 1.3)"""
    if math.isnan(sunz):
        return F32(0.)
    if sunz < start:
        return F32(1.)
    f = min(max((sunz - F32(start)) / F32(end - start), F32(0.)), F32(1.))
    f = F32(1.) - F32(math.log2(f + F32(1.)))
    a = f ** F32(1.3)
    return a / (a + (F32(1.) - f) ** F32(1.3))


@njit(nogil=True, cache=True)
def _stretch(x, lo, hi, l0, l1):
    if math.isnan(x):
        return F32(0.)
    x = min(max(x, F32(lo)), F32(hi))
    v = (F32(math.log10(x)) - F32(l0)) / F32(l1 - l0)
    return min(max(v, F32(0.)), F32(1.))


@njit(nogil=True, cache=True)
def _stretch_enhanced(x, black, white, knee, toe):
    """composite.enhanced_stretch：平滑暗部下限、白點 white、knee 以上平滑壓縮"""
    if math.isnan(x):
        return F32(0.)
    d = np.float64(x) - black
    if toe > 0.:
        d = 0.5 * (d + math.sqrt(d * d + 4. * toe * toe))
    elif d < 0.:
        d = 0.
    v = math.log10((black + d) / black) / math.log10(white / black)
    if knee < 1. and v > knee:
        v = knee + (1. - knee) * (1. - math.exp(-(v - knee) / (1. - knee)))
    return F32(min(max(v, 0.), 1.))


@njit(nogil=True, cache=True)
def _compose_kernel(blue, green, red_lo, nir, red_hi, bt, lon, lat, geo, luts, smin, smax, orders, clipa,
                    tau, sph, scale, sphlut, o3, par, flags, out):
    gmst, ra, dec, sx, sy, sz = geo[0], geo[1], geo[2], geo[3], geo[4], geo[5]
    sin_dec, cos_dec, tan_dec = math.sin(dec), math.cos(dec), math.tan(dec)
    A, F = 6378.137, 1. / 298.257223563
    two_pi = 2. * math.pi
    sza_limit, sza_max, fade0, fade1, frac = par[0], par[1], par[2], par[3], F32(par[4])
    lin0, lin1, lin2 = F32(par[5]), F32(par[6]), F32(par[7])
    gk, ga, gb, gfloor, gmin = F32(par[8]), F32(par[9]), F32(par[10]), F32(par[11]), F32(par[12])
    rlo, rhi, llo, lhi = par[13], par[14], par[15], par[16]
    white, knee, toe, sat = par[17], par[18], par[19], F32(par[20])
    vib_hi, vib_span, lr, lg, lb = F32(par[21]), F32(par[22]), F32(par[23]), F32(par[24]), F32(par[25])
    mode, full, fade, enhanced = flags[0], flags[1] == 1, flags[2] == 1, flags[3] == 1
    has_green, has_hi, has_bt = green.size > 0, red_hi.size > 0, bt.size > 0
    cloud_tau = flags[4] == 1 and full and has_bt
    has_o3 = o3[0] + o3[1] + o3[2] + o3[3] > 0.
    oz_num, oz_2h = 1. + OZONE_HEIGHT, 2. * OZONE_HEIGHT
    for i in range(lon.size):
        lonr, latr = math.radians(lon[i]), math.radians(lat[i])
        sl, cl = math.sin(latr), math.cos(latr)
        # ---- 太陽 (pyorbital cos_zen / get_alt_az) ----
        h = gmst + lonr - ra
        cosz = sl * sin_dec + cl * cos_dec * math.cos(h)
        sunz64 = math.degrees(math.acos(min(max(cosz, -1.), 1.)))
        suna = math.degrees(math.atan2(-math.sin(h), cl * tan_dec - sl * math.cos(h))) % 360.
        # ---- 衛星 (pyorbital get_observer_look，觀測點高度 0) ----
        theta = (gmst + lonr) % two_pi
        c = 1. / math.sqrt(1. + F * (F - 2.) * sl * sl)
        sq = c * (1. - F) ** 2
        achcp = A * c * cl
        st, ct = math.sin(theta), math.cos(theta)
        rx, ry, rz = sx - achcp * ct, sy - achcp * st, sz - A * sq * sl
        top_s = sl * ct * rx + sl * st * ry - cl * rz
        top_e = -st * rx + ct * ry
        top_z = cl * ct * rx + cl * st * ry + sl * rz
        sata = math.degrees((math.atan2(-top_e, top_s) + math.pi) % two_pi)
        rg = math.sqrt(rx * rx + ry * ry + rz * rz)
        satz64 = 90. - math.degrees(math.asin(min(top_z / rg, 1.)))
        d = abs(suna - sata)
        ssa = F32(min(d, 360. - d))
        sunz, satz = F32(sunz64), F32(satz64)
        # ---- Rayleigh 查表座標 (同 pyspectral：NaN 為 0、截斷在查表範圍內、取正割) ----
        zs = min(max(sunz if not math.isnan(sunz) else F32(0.), F32(0.)), clipa[0])
        zv = min(max(satz if not math.isnan(satz) else F32(0.), F32(0.)), clipa[1])
        sunsec = F32(1.) / F32(math.cos(F32(math.radians(zs))))
        satsec = F32(1.) / F32(math.cos(F32(math.radians(zv))))
        azi = F32(180.) - ssa
        # ---- 太陽天頂角校正 ----
        b = _sunz_corr(blue[i], cosz, sza_limit, sza_max)
        r = _sunz_corr(red_lo[i], cosz, sza_limit, sza_max)
        n = _sunz_corr(nir[i], cosz, sza_limit, sza_max)
        g = _sunz_corr(green[i], cosz, sza_limit, sza_max) if has_green else F32(0.)
        rh = _sunz_corr(red_hi[i], cosz, sza_limit, sza_max) if has_hi else F32(0.)
        # ---- 臭氧吸收 (correction.ozone_airmass；臭氧在雲與 Rayleigh 散射之上，先除以臭氧穿透率) ----
        if has_o3:
            cs = math.cos(math.radians(min(max(sunz64, 0.), 90.)))
            cv = math.cos(math.radians(min(max(satz64, 0.), 90.)))
            am = oz_num / math.sqrt(cs * cs + oz_2h) + oz_num / math.sqrt(cv * cv + oz_2h)
            b = b * F32(math.exp(o3[0] * am))
            g = g * F32(math.exp(o3[1] * am))
            r = r * F32(math.exp(o3[2] * am))
            rh = rh * F32(math.exp(o3[2] * am))
            n = n * F32(math.exp(o3[3] * am))
        # ---- 雲頂削減 (CIRA) ----
        fac = F32(1.)
        if has_bt:
            red_ = min(max((F32(283.) - (bt[i] + F32(273.15))) / F32(50.), F32(0.)), F32(1.))
            if not math.isnan(red_):
                fac = F32(1.) - F32(0.7) * red_
        mus = min(max(F32(math.cos(F32(math.radians(sunz)))), F32(0.02)), F32(1.))
        muv = min(max(F32(math.cos(F32(math.radians(satz)))), F32(0.02)), F32(1.))
        # 0 藍、1 綠、2 紅、3 近紅外；紅光削減依據：沒有亮溫時用 pyspectral 的紅光反射率規則
        rp = r * F32(100.)
        corr = np.empty(4, F32)
        rays = np.empty(4, F32)
        trans = np.empty(4, F32)
        sphs = np.empty(4, F32)
        for j in range(4):
            if j == 1 and not has_green:
                continue
            ray = _path(luts[j], smin, smax, orders, sunsec, azi, satsec, rp, not has_bt) * scale[j]
            rays[j] = ray * fac
            if cloud_tau and fac < F32(1.):          # 高雲：穿透率與球面反照率的光學厚度同樣乘上削減係數
                te = tau[j] * fac
                trans[j] = _trans(mus, te) * _trans(muv, te)
                p = (fac - F32(SPH_F0)) * F32(SPH_INV)
                q = min(max(int(p), 0), SPH_N - 2)
                sphs[j] = sphlut[j, q] + (p - F32(q)) * (sphlut[j, q + 1] - sphlut[j, q])
            else:
                trans[j] = _trans(mus, tau[j]) * _trans(muv, tau[j]) if full else F32(1.)
                sphs[j] = sph[j]
        vals = (b, g, r, n)
        for j in range(4):
            if j == 1 and not has_green:
                corr[j] = F32(0.)
                continue
            if j == 3 and not full:              # 舊版：近紅外不校正
                corr[j] = n
                continue
            x = vals[j] - rays[j]
            if full:
                x = x / trans[j]
                x = x / (F32(1.) + sphs[j] * x)
            corr[j] = x
        cb, cg, cr, cn = corr[0], corr[1], corr[2], corr[3]
        # ---- 綠光 ----
        if mode == 0:
            gg = (F32(1.) - frac) * cg + frac * cn
        elif mode == 1:
            rr = max(cr, gfloor)
            bb = max(max(cb, gfloor), gmin * rr)
            gg = (F32(1.) - frac) * (gk * bb ** ga * rr ** gb) + frac * cn
        else:
            gg = lin0 * cb + lin1 * cr + lin2 * cn
        # ---- 紅光銳化 (高解析紅光與低解析紅光共用路徑反射等各項) ----
        if has_hi:
            ray_hi = rays[2]
            if not has_bt:                       # 沒有亮溫時，削減依據是高解析紅光自己的反射率 (同 numpy 版)
                ray_hi = _path(luts[2], smin, smax, orders, sunsec, azi, satsec, rh * F32(100.), True) * scale[2]
            x = rh - ray_hi
            if full:
                x = x / trans[2]
                x = x / (F32(1.) + sphs[2] * x)
            # 低解析紅光校正後剛好為 0 時 (晨昏線附近，反射率與路徑反射相等) 比例為 1：numba 的浮點除以 0 會丟出
            # ZeroDivisionError，不像 numpy 得到 inf / nan (composite.sharpening_ratio 把非有限值設為 1)
            ratio = x / cr if cr != F32(0.) else F32(1.)
            if not (ratio >= F32(0.)) or math.isinf(ratio):
                ratio = F32(1.)
            ratio = min(ratio, F32(1.5))
            cr, gg, cb = x, gg * ratio, cb * ratio
        # ---- 晨昏線淡出 ----
        if fade:
            f = _fade(sunz, fade0, fade1)
            cr, gg, cb = cr * f, gg * f, cb * f
        # ---- 色調曲線 ----
        if enhanced:
            o0 = _stretch_enhanced(cr, rlo, white, knee, toe)
            o1 = _stretch_enhanced(gg, rlo, white, knee, toe)
            o2 = _stretch_enhanced(cb, rlo, white, knee, toe)
        else:
            o0 = _stretch(cr, rlo, rhi, llo, lhi)
            o1 = _stretch(gg, rlo, rhi, llo, lhi)
            o2 = _stretch(cb, rlo, rhi, llo, lhi)
        # ---- 飽和度 (composite.vibrance) ----
        if sat != F32(1.):
            lum = lr * o0 + lg * o1 + lb * o2
            k = F32(1.) + (sat - F32(1.)) * min(max((vib_hi - lum) / vib_span, F32(0.)), F32(1.))
            o0 = min(max(lum + k * (o0 - lum), F32(0.)), F32(1.))
            o1 = min(max(lum + k * (o1 - lum), F32(0.)), F32(1.))
            o2 = min(max(lum + k * (o2 - lum), F32(0.)), F32(1.))
        out[i, 0] = o0
        out[i, 1] = o1
        out[i, 2] = o2
