"""McIDAS 'GOES' 導航 (自旋掃描衛星：SMS-1/2、GOES-1~7 的 VISSR / VAS，McIDAS AREA 導航區塊類型 'GOES')
— GOES-PREGVAR 讀取器用 (檔名沒有 '-'，不會被當成讀取器載入)

正算 (行、像素 -> 經緯度) 與 McIDAS nvxgoes / VisAD edu.wisc.ssec.mcidas.GOESnav 逐式相同：
  導航區塊 (128 個 int32)：Kepler 軌道 (epoch、半長軸、離心率、傾角、近地點、升交點、平近點角)、
  自旋軸 (赤緯、赤經)、自旋週期、掃描幾何 (掃描線數 LINTOT = 感測器數 x 100000 + 掃描數、
  南北 / 東西掃描角)、姿態 (pitch、yaw、roll、skew)、元素偏移 gamma / gamma-dot。
  掃描時間：同一圈 (自旋) 的 NUMSEN 條線同時取樣，
    samtim = pictim + 自旋週期 x ((NINT(行號) - 1) // NUMSEN + 1) (小時，// 為整數除法)
  -> 該時刻衛星位置 (Kepler 方程) -> 掃描方向 -> 與地球橢球 (a = 6378.2064、b = 6356.5838 km) 的交點
  -> 依地球自轉轉為地固座標 -> 大地經緯度。
  (舊版 Python 移植 nav.py 的 ietimh、NUMSEN、掃描序號用了浮點除法，與 McIDAS 不同。)
反算 (經緯度 -> 行、像素)：直接解上述旋轉 (不是 McIDAS 以 beta 紀錄內插時間的近似反算)，
  每次以行號求掃描時間、重算到掃描序號不變為止，與正算互為精確反函數。
每一圈的衛星位置、地球自轉角等只與掃描序號有關，先逐圈精確計算成表 (不內插)，逐點計算時查表。
計算方式依 [navigation] engine：numba (預設) / cython (mcgoes_nav_cy) / numpy，結果相同。
"""
import math

import numpy as np

from satima.core.i18n import tr
from satima.geo.nav import Nav

try:
    import numba
except ImportError:
    numba = None

MCMISSING = np.int32(np.uint32(0x80808080).view(np.int32))
RDPDG = 1.745329252e-02
ASQ, BSQ, AB = 40683833.48, 40410330.18, 40546851.22       # 橢球 a^2、b^2、ab (km^2)
EMEGA = 0.26251617                                        # 地球自轉 (rad / 小時)
SOLSID, SHA = 1.00273791, 100.26467
NSTATE = 8                                                # 每圈：samtim、衛星 x y z、rot、cos、sin (地球自轉角)、gamma 項
NCONST = 26


# ---- McIDAS 時間 / 角度工具 (同 nvxgoes 與 McIDASUtil) -------------------------------------------
def flalo(m):
    """DDDMMSS / HHMMSS 整數 -> 度 / 小時。"""
    n = abs(int(m))
    v = n // 10000 + (n // 100 % 100) / 60. + (n % 100) / 3600.
    return -v if m < 0 else v


def _itime(t):
    """小時 -> HHMMSS 整數 (同 McIDASUtil.mcDoubleToPackedInteger：四捨五入到秒)。"""
    j = int(3600. * abs(t) + 0.5)
    v = 10000 * (j // 3600) + 100 * (j // 60 % 60) + j % 60
    return -v if t < 0 else v


def _icon1(yymmdd):
    num = (0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334)
    year, month, day = yymmdd // 10000 % 100, yymmdd // 100 % 100, yymmdd % 100
    if month < 1 or month > 12:
        month = 1
    jday = day + num[month - 1]
    if year % 4 == 0 and month > 2:
        jday += 1
    return 1000 * year + jday


def _leapyr(iy):
    return 366 - (iy % 4 + 3) // 4


def _days(yyddd):
    iy, idd = yyddd // 1000 % 100, yyddd % 1000
    return 365 * (iy - 1) + (iy - 1) // 4 + 1 + idd - 1


def _timdif(yd1, hms1, yd2, hms2):
    """兩時刻相差的分鐘數 (McIDAS TIMDIF)。"""
    return 1440. * (_days(yd2) - _days(yd1)) + 60. * (flalo(hms2) - flalo(hms1))


def _epoch(ietimy, ietimh, semima, oeccen, xmeana):
    """Kepler epoch -> 近地點時刻 (McIDAS EPOCH)。"""
    xmmc = 0.07436574 * math.sqrt(6378.388 / semima) ** 3
    xmanom = RDPDG * xmeana
    t = flalo(ietimh) - (xmanom - oeccen * math.sin(xmanom)) / (60. * xmmc)
    iday = 0
    if t > 48.:
        t, iday = t - 48., 2
    elif t > 24.:
        t, iday = t - 24., 1
    elif t < -24.:
        t, iday = t + 48., -2
    elif t < 0.:
        t, iday = t + 24., -1
    ietimh = _itime(t)
    if iday == 0:
        return ietimy, ietimh
    jyear, jday = ietimy // 1000 % 100, ietimy % 1000 + iday
    if jday < 1:
        jyear -= 1
        jday += _leapyr(jyear)
    elif jday > _leapyr(jyear):
        jday -= _leapyr(jyear)
        jyear += 1
    return 1000 * jyear + jday, ietimh


# ---- 導航區塊 -> 參數 -------------------------------------------------------------------------
def params(nav):
    """McIDAS AREA 'GOES' 導航區塊 (128 個 int32，已處理位元組順序) -> 參數 dict；無法導航時 ValueError。"""
    a = [int(x) for x in np.asarray(nav, dtype=np.int64)]
    if a[6] <= 0 or a[7] < 0 or (a[6] == 0 and a[7] == 0 and a[8] == 0):
        raise ValueError(tr('導航區塊沒有軌道參數 (半長軸 = %d)',
                            'the navigation block has no orbit parameters (semi-major axis = %d)') % a[6])
    if a[12] == 0 and a[13] == 0 and a[14] == 0:
        raise ValueError(tr('導航區塊沒有自旋軸 (赤緯 / 赤經)',
                            'the navigation block has no spin axis (declination / right ascension)'))
    if a[15] == 0:
        raise ValueError(tr('導航區塊沒有自旋週期', 'the navigation block has no spin period'))
    jday, jtime = a[1], a[2]
    p = dict(navday=jday % 100000, raw=a)
    ietimy = _icon1(a[4])
    ietimh = 100 * (a[5] // 100) + int(math.floor(.6 * (a[5] % 100) + .5))
    p['semima'], p['oeccen'] = a[6] / 100., a[7] / 1e6
    p['orbinc'], p['xmeana'], p['perhel'], p['asnode'] = a[8] / 1000., a[9] / 1000., a[10] / 1000., a[11] / 1000.
    p['ietimy'], p['ietimh'] = _epoch(ietimy, ietimh, p['semima'], p['oeccen'], p['xmeana'])
    p['declin'], p['rascen'] = flalo(a[12]), flalo(a[13])
    p['piclin'] = a[14] / 10000. if a[14] > 1000000 else float(a[14])
    spinra = a[15] / 1000.
    if spinra < 300.:
        spinra = 60000. / spinra
    p['spinra'] = spinra                                  # 自旋週期 (毫秒)
    p['deglin'], p['lintot'], p['degele'], p['eletot'] = flalo(a[16]), a[17], flalo(a[18]), a[19]
    p['pitch'], p['yaw'], p['roll'] = flalo(a[20]), flalo(a[21]), flalo(a[22])
    p['skew'] = 0. if a[28] == MCMISSING else a[28] / 100000.
    p['numsen'] = max(1, a[17] // 100000 % 100)
    p['scans'] = a[17] % 100000
    p['pictim'] = flalo(jtime)
    p['gamma'], p['gamdot'] = a[38] / 100., a[39] / 100.
    p['tmpscl'] = spinra / 3600000.                       # 每圈的小時數
    p['time'] = (jday // 100000, jday % 100000, jtime)
    # 衛星編號 > 25 或 12 (VAS、GMS)：beta 紀錄的掃描序號與時間 (McIDAS 反算與子衛星點用；正算不用)
    if (jday // 100000 > 25 or jday // 100000 == 12) and a[30] > 0:
        p['beta'] = (float(a[30]), flalo(a[31]), float(a[34]), flalo(a[35]))
    else:
        p['beta'] = (1., p['pictim'], float(p['scans']), p['pictim'] + p['scans'] * p['tmpscl'])
    return p


def constants(p):
    """逐點計算用的常數 (NCONST 個)。"""
    totlin = p['numsen'] * p['scans']
    radlin = RDPDG * p['deglin'] / (totlin - 1.)
    radele = RDPDG * p['degele'] / (p['eletot'] - 1.)
    picele = (1. + p['eletot']) / 2.
    cpitch, cyaw, croll = RDPDG * p['pitch'], RDPDG * p['yaw'], RDPDG * p['roll']
    pskew = math.atan2(p['skew'], radlin / radele)
    stp, ctp = math.sin(cpitch), math.cos(cpitch)
    sty, cty = math.sin(cyaw - pskew), math.cos(cyaw - pskew)
    s_r, c_r = math.sin(croll), math.cos(croll)
    rotm = (c_r * ctp, sty * s_r * ctp + cty * stp, -s_r, sty * c_r, -c_r * stp, cty * ctp - sty * s_r * stp)
    dec, ras = p['declin'] * RDPDG, p['rascen'] * RDPDG
    sd, cd, sr, cr = math.sin(dec), math.cos(dec), math.sin(ras), math.cos(ras)
    b = (-sr, cr, 0., -sd * cr, -sd * sr, cd, cd * cr, cd * sr, sd)
    c = [float(p['numsen']), p['piclin'], radlin, picele, radele, p['gamma'], p['gamdot'], *rotm, *b,
         math.sqrt(rotm[4] ** 2 + rotm[5] ** 2), math.atan2(rotm[4], rotm[5]), p['pictim'], p['tmpscl']]
    assert len(c) == NCONST
    return np.array(c, dtype=np.float64)


def _xref(p):
    """導航日 0 時的格林威治赤經 (rad) (McIDAS RAERAC(navday, 0, 0))。"""
    raha = 0. + _timdif(74001, 0, p['navday'], 0) * SOLSID / 4. + SHA
    return (raha % 360.) * RDPDG


def satvec(p, samtim):
    """時刻 samtim (導航日的小時) 的衛星位置 (慣性座標，km) (McIDAS SATVEC)。"""
    o, pp, a = RDPDG * p['orbinc'], RDPDG * p['perhel'], RDPDG * p['asnode']
    so, co = math.sin(o), math.cos(o)
    sp, cp = math.sin(pp) * p['semima'], math.cos(pp) * p['semima']
    sa, ca = math.sin(a), math.cos(a)
    px, py, pz = cp * ca - sp * sa * co, cp * sa + sp * ca * co, sp * so
    qx, qy, qz = -sp * ca - cp * sa * co, -sp * sa + cp * ca * co, cp * so
    e = p['oeccen']
    srome2 = math.sqrt(1. - e) * math.sqrt(1. + e)
    xmmc = 0.07436574 * 6378.388 * math.sqrt(6378.388 / p['semima']) / p['semima']
    te = 1440. * _days(p['ietimy']) + 60. * flalo(p['ietimh'])
    tdife = _days(p['navday']) * 1440. - te
    xmanom = xmmc * (tdife + samtim * 60.)
    ecanm1 = xmanom
    ecanom = xmanom
    for _ in range(20):
        ecanom = xmanom + e * math.sin(ecanm1)
        if abs(ecanom - ecanm1) < 1e-8:
            break
        ecanm1 = ecanom
    xo, yo = math.cos(ecanom) - e, srome2 * math.sin(ecanom)
    return xo * px + yo * qx, xo * py + yo * qy, xo * pz + yo * qz


def scan_table(p, c, k0, k1):
    """掃描序號 k0..k1 每圈的狀態 (K x NSTATE)：samtim、衛星 x y z、rot、cos / sin(地球自轉角)、gamma-dot 項。"""
    xref = _xref(p)
    b11, b12, b13, b21, b22, b23 = c[13:19]
    rows = []
    for k in range(k0, k1 + 1):
        samtim = p['tmpscl'] * k + p['pictim']
        xs, ys, zs = satvec(p, samtim)
        rot = math.atan2(b21 * xs + b22 * ys + b23 * zs, b11 * xs + b12 * ys + b13 * zs) + math.pi
        w = EMEGA * samtim + xref
        rows.append((samtim, xs, ys, zs, rot, math.cos(w), math.sin(w), p['gamdot'] * samtim))
    return np.ascontiguousarray(rows, dtype=np.float64)


# ---- 逐點計算 ---------------------------------------------------------------------------------
# 常數 C：0 numsen 1 piclin 2 radlin 3 picele 4 radele 5 gamma 6 gamdot 7-12 rotm11 13 21 23 31 33
#         13-21 b11..b33 22 sqrt(rfact) 23 roasin 24 pictim 25 tmpscl
def _scan_index(L, numsen):
    """行號 -> 掃描序號 (GOESnav：Math.round 後整數除法，整數除法向 0 截斷)。"""
    ilin = math.floor(L + 0.5)
    q = (ilin - 1.) / numsen
    return int(q) + 1                                   # int() 向 0 截斷


if numba is not None:
    _scan_index_nb = numba.njit(cache=True)(_scan_index)

    @numba.njit(cache=True)
    def _fwd1(L, E, C, S, k0):
        K = S.shape[0]
        i = _scan_index_nb(L, C[0]) - k0
        i = 0 if i < 0 else (K - 1 if i > K - 1 else i)
        xs, ys, zs, rot, ct, st, gd = S[i, 1], S[i, 2], S[i, 3], S[i, 4], S[i, 5], S[i, 6], S[i, 7]
        ylin = (L - C[1]) * C[2]
        yele = (E - C[3] + C[5] + gd) * C[4] - rot
        cl, sl, ce, se = math.cos(ylin), math.sin(ylin), math.cos(yele), math.sin(yele)
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
            return np.nan, np.nan
        s = -(bq + math.sqrt(rad)) / (2. * aq)
        x, y, z = xs + elo * s, ys + emo * s, zs + eno * s
        x1, y1 = ct * x + st * y, -st * x + ct * y
        a = math.atan(z / math.sqrt(x1 * x1 + y1 * y1))
        lat = math.degrees(math.atan2(ASQ * math.sin(a), BSQ * math.cos(a)))
        lon = math.degrees(math.atan2(y1, x1))
        return lon, lat

    @numba.njit(cache=True)
    def _inv1(lon, lat, C, S, k0, kstart):
        K = S.shape[0]
        ylat = math.atan2(BSQ * math.sin(math.radians(lat)), ASQ * math.cos(math.radians(lat)))
        ylon = math.radians(lon)
        sn, cs = math.sin(ylat), math.cos(ylat)
        tn = (sn / cs) ** 2
        r = AB * math.sqrt((1. + tn) / (BSQ + ASQ * tn))
        xe, ye, ze = r * cs * math.cos(ylon), r * cs * math.sin(ylon), r * sn
        i = kstart - k0
        L, E = np.nan, np.nan
        for _ in range(6):
            i = 0 if i < 0 else (K - 1 if i > K - 1 else i)
            xs, ys, zs, rot, ct, st, gd = S[i, 1], S[i, 2], S[i, 3], S[i, 4], S[i, 5], S[i, 6], S[i, 7]
            x, y = ct * xe - st * ye, st * xe + ct * ye                     # 地固 -> 慣性
            vx, vy, vz = x - xs, y - ys, ze - zs
            # 看得到 (視線與該點的橢球法向量夾角 > 90 度)
            if vx * x / ASQ + vy * y / ASQ + vz * ze / BSQ >= 0.:
                return np.nan, np.nan
            n = math.sqrt(vx * vx + vy * vy + vz * vz)
            elo, emo, eno = vx / n, vy / n, vz / n
            eli = C[13] * elo + C[14] * emo + C[15] * eno
            emi = C[16] * elo + C[17] * emo + C[18] * eno
            eni = C[19] * elo + C[20] * emo + C[21] * eno
            q = eni / C[22]
            if q > 1. or q < -1.:
                return np.nan, np.nan
            ylin = C[23] - math.asin(q)
            cl, sl = math.cos(ylin), math.sin(ylin)
            e1 = C[7] * cl - C[8] * sl
            e2 = C[9] * cl - C[10] * sl
            # 元素角 (自旋角) 取最接近影像中心的一圈
            ang = math.atan2(e2, e1) - math.atan2(emi, eli) + rot - (C[5] + gd) * C[4]
            ang = (ang + math.pi) % (2. * math.pi) - math.pi
            L = ylin / C[2] + C[1]
            E = ang / C[4] + C[3]
            j = _scan_index_nb(L, C[0]) - k0
            j = 0 if j < 0 else (K - 1 if j > K - 1 else j)               # 與正算相同：表外用表的兩端
            if j == i:
                break
            i = j
        return L, E

    @numba.njit(parallel=True, cache=True)
    def forward(line, ele, C, S, k0, lon, lat):
        for i in numba.prange(line.shape[0]):
            if not (math.isfinite(line[i]) and math.isfinite(ele[i])):
                lon[i], lat[i] = np.nan, np.nan
                continue
            lon[i], lat[i] = _fwd1(line[i], ele[i], C, S, k0)

    @numba.njit(parallel=True, cache=True)
    def inverse(lon, lat, C, S, k0, kstart, line, ele):
        for i in numba.prange(lon.shape[0]):
            if not (math.isfinite(lon[i]) and math.isfinite(lat[i])):
                line[i], ele[i] = np.nan, np.nan
                continue
            line[i], ele[i] = _inv1(lon[i], lat[i], C, S, k0, kstart)


# ---- numpy 版 (沒有 numba 時) -----------------------------------------------------------------
def _scan_index_np(L, numsen):
    ilin = np.floor(L + 0.5)
    return np.trunc((ilin - 1.) / numsen).astype(np.int64) + 1


def _forward_np(line, ele, C, S, k0):
    i = np.clip(_scan_index_np(line, C[0]) - k0, 0, S.shape[0] - 1)
    s = S[i]
    rot, ct, st, gd = s[:, 4], s[:, 5], s[:, 6], s[:, 7]
    ylin = (line - C[1]) * C[2]
    yele = (ele - C[3] + C[5] + gd) * C[4] - rot
    cl, sl, ce, se = np.cos(ylin), np.sin(ylin), np.cos(yele), np.sin(yele)
    e1, e2, eni = C[7] * cl - C[8] * sl, C[9] * cl - C[10] * sl, C[11] * cl - C[12] * sl
    eli, emi = ce * e1 + se * e2, -se * e1 + ce * e2
    elo = C[13] * eli + C[16] * emi + C[19] * eni
    emo = C[14] * eli + C[17] * emi + C[20] * eni
    eno = C[15] * eli + C[18] * emi + C[21] * eni
    xs, ys, zs = s[:, 1], s[:, 2], s[:, 3]
    basq = BSQ / ASQ
    aq = basq + (1. - basq) * eno * eno
    bq = 2. * ((elo * xs + emo * ys) * basq + eno * zs)
    cq = (xs * xs + ys * ys) * basq + zs * zs - BSQ
    rad = bq * bq - 4. * aq * cq
    with np.errstate(invalid='ignore'):
        sq = -(bq + np.sqrt(rad)) / (2. * aq)
        x, y, z = xs + elo * sq, ys + emo * sq, zs + eno * sq
        x1, y1 = ct * x + st * y, -st * x + ct * y
        a = np.arctan(z / np.sqrt(x1 * x1 + y1 * y1))
        lat = np.degrees(np.arctan2(ASQ * np.sin(a), BSQ * np.cos(a)))
        lon = np.degrees(np.arctan2(y1, x1))
    bad = (rad < 1.) | ~np.isfinite(line) | ~np.isfinite(ele)
    return np.where(bad, np.nan, lon), np.where(bad, np.nan, lat)


def _inverse_np(lon, lat, C, S, k0, kstart):
    ok = np.isfinite(lon) & np.isfinite(lat)
    la, lo = np.radians(np.where(ok, lat, 0.)), np.radians(np.where(ok, lon, 0.))
    ylat = np.arctan2(BSQ * np.sin(la), ASQ * np.cos(la))
    sn, cs = np.sin(ylat), np.cos(ylat)
    tn = (sn / cs) ** 2
    r = AB * np.sqrt((1. + tn) / (BSQ + ASQ * tn))
    xe, ye, ze = r * cs * np.cos(lo), r * cs * np.sin(lo), r * sn
    K = S.shape[0]
    i = np.full(lon.shape, kstart - k0, dtype=np.int64)
    done = ~ok
    L = np.full(lon.shape, np.nan)
    E = np.full(lon.shape, np.nan)
    for _ in range(6):
        i = np.clip(i, 0, K - 1)
        s = S[i]
        xs, ys, zs, rot, ct, st, gd = (s[:, j] for j in range(1, 8))
        x, y = ct * xe - st * ye, st * xe + ct * ye
        vx, vy, vz = x - xs, y - ys, ze - zs
        vis = (vx * x / ASQ + vy * y / ASQ + vz * ze / BSQ) < 0.
        n = np.sqrt(vx * vx + vy * vy + vz * vz)
        elo, emo, eno = vx / n, vy / n, vz / n
        eli = C[13] * elo + C[14] * emo + C[15] * eno
        emi = C[16] * elo + C[17] * emo + C[18] * eno
        eni = C[19] * elo + C[20] * emo + C[21] * eno
        q = eni / C[22]
        vis &= np.abs(q) <= 1.
        with np.errstate(invalid='ignore'):
            ylin = C[23] - np.arcsin(np.clip(q, -1., 1.))
        cl, sl = np.cos(ylin), np.sin(ylin)
        e1, e2 = C[7] * cl - C[8] * sl, C[9] * cl - C[10] * sl
        ang = np.arctan2(e2, e1) - np.arctan2(emi, eli) + rot - (C[5] + gd) * C[4]
        ang = (ang + np.pi) % (2. * np.pi) - np.pi
        Ln = ylin / C[2] + C[1]
        En = ang / C[4] + C[3]
        upd = ~done
        L = np.where(upd, np.where(vis, Ln, np.nan), L)
        E = np.where(upd, np.where(vis, En, np.nan), E)
        j = np.clip(_scan_index_np(np.where(np.isfinite(L), L, 0.), C[0]) - k0, 0, K - 1)
        done = done | ~vis | (j == i)
        i = np.where(done, i, j)
        if done.all():
            break
    return L, E


_NOTICE = set()


def _kernels(engine):
    if engine == 'cython':
        try:
            import mcgoes_nav_cy
            return mcgoes_nav_cy.forward, mcgoes_nav_cy.inverse
        except ImportError:
            if 'cy' not in _NOTICE:
                _NOTICE.add('cy')
                print(tr('  [導航] 沒有編譯好的 McIDAS GOES Cython 模組 (mcgoes_nav_cy)，改用 numba；'
                         '編譯方式：python sensors/geo/build_vissr_nav_cy.py',
                         '  [navigation] the McIDAS GOES Cython module (mcgoes_nav_cy) is not built; using numba '
                         'instead; to build: python sensors/geo/build_vissr_nav_cy.py'))
            engine = 'numba'
    if engine == 'numba' and numba is not None:
        return forward, inverse
    return None


class McGoesNav(Nav):
    """McIDAS 'GOES' 導航。陣列第 r 列、第 c 行 = 影像行號 first_line + r x line_step、
    像素號 first_ele + c x ele_step (影像座標以可見光線、像素為單位，由 1 起算，同 AREA 檔)。"""

    def __init__(self, shape, navblock, first_line, first_ele, line_step=1., ele_step=1., key=None, engine='numba'):
        self.shape = (int(shape[0]), int(shape[1]))
        self.p = params(navblock)
        self.C = constants(self.p)
        self.first_line, self.first_ele = float(first_line), float(first_ele)
        self.line_step, self.ele_step = float(line_step), float(ele_step)
        last = self.first_line + (self.shape[0] - 1) * self.line_step
        ns = self.C[0]
        # 影像涵蓋的掃描 (前後多留，反算疊代時不必外插)；整個地球最多 scans 圈
        self.k0 = max(-2, min(_scan_index(self.first_line, ns), _scan_index(last, ns)) - 4)
        k1 = max(_scan_index(self.first_line, ns), _scan_index(last, ns)) + 4
        self.S = scan_table(self.p, self.C, self.k0, k1)
        self.kmid = _scan_index(0.5 * (self.first_line + last), ns)
        self._key = key
        self._k = _kernels(engine)

    def key(self):
        return ('mcgoes', self._key) if self._key is not None else id(self)

    def subpoint(self):
        """影像中央時刻的子衛星點 (經度, 緯度) 與高度 (km)。"""
        s = self.S[self.kmid - self.k0]
        x, y = s[5] * s[1] + s[6] * s[2], -s[6] * s[1] + s[5] * s[2]
        z = s[3]
        a = math.atan(z / math.hypot(x, y))
        lat = math.degrees(math.atan2(ASQ * math.sin(a), BSQ * math.cos(a)))
        return math.degrees(math.atan2(y, x)), lat, math.sqrt(x * x + y * y + z * z) - math.sqrt(ASQ)

    def to_lonlat(self, row, col):
        row, col = np.broadcast_arrays(np.asarray(row, dtype=np.float64), np.asarray(col, dtype=np.float64))
        shape = row.shape
        line = np.ascontiguousarray(self.first_line + row.ravel() * self.line_step)
        ele = np.ascontiguousarray(self.first_ele + col.ravel() * self.ele_step)
        if self._k is not None:
            lon, lat = np.empty(line.size), np.empty(line.size)
            self._k[0](line, ele, self.C, self.S, self.k0, lon, lat)
        else:
            lon, lat = _forward_np(line, ele, self.C, self.S, self.k0)
        return lon.reshape(shape), lat.reshape(shape)

    def to_rowcol(self, lon, lat):
        lon, lat = np.broadcast_arrays(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        shape = lon.shape
        lon = np.ascontiguousarray(((lon.ravel() + 180.) % 360.) - 180.)
        lat = np.ascontiguousarray(lat.ravel())
        if self._k is not None:
            line, ele = np.empty(lon.size), np.empty(lon.size)
            self._k[1](lon, lat, self.C, self.S, self.k0, self.kmid, line, ele)
        else:
            line, ele = _inverse_np(lon, lat, self.C, self.S, self.k0, self.kmid)
        rows = (line - self.first_line) / self.line_step
        cols = (ele - self.first_ele) / self.ele_step
        return rows.reshape(shape), cols.reshape(shape)
