"""GOES-8~15 Imager (GVAR 時代) 共用：導航 (O&A 模型)、GVAR 數值格式、校正
— GOES-GVAR 與 GOES-AREA 讀取器共用 (檔名沒有 '-'，不會被當成讀取器載入)

導航：GOES I-M 的軌道與姿態 (O&A) 模型，同 McIDAS nvxgvar 與舊版 GOES_GVAR_NAV.py / GOES_NAV.py 的 le2ll：
  影像座標 (可見光行號、像素號，由 1 起算)
    -> 掃描角 (仰角 alpha = elvmax - (行號 - 4.5) x 28 µrad、掃描角 zeta = (像素號 - 1) x 16 µrad - scnmax)
    -> 儀器中心偏差 (doff) 與失準 (roll / pitch misalignment) 修正
    -> 儀器 -> 衛星本體 (roll、pitch、yaw 小角度矩陣) -> 地固座標 (由衛星經度、軌道傾角、軌道偏航求得)
    -> 與地球橢球 (6378.137 / 6356.7533 km) 的交點
  IMC (影像運動補償) 啟動時 (GOES-8~15 平常都是)，只用參考經度與參考姿態；
  IMC 關閉時才依影像時間計算軌道與姿態級數 (同 McIDAS，整張影像用同一時間)。
  反算 (經緯度 -> 行列號) 直接解上述旋轉與修正 (固定點疊代到收斂)，與正算互為精確反函數。
  計算方式依 [navigation] engine：numba (預設) / numpy (結果相同)；cython 時沒有編譯好的模組則改用 numba。
數值：GVAR Block 0 的浮點數為 Gould (IBM 類) 格式；AREA 導航區塊為放大 10^7 的整數。時間為 BCD。
校正：紅外依 conversion_table/GOES 的計數值 -> 亮溫表 (衛星、side、頻道、偵測器)；
      可見光 = 偵測器各自的 bias + gain x 計數值 (+ 二次項)，再乘以 radiance -> albedo 係數。
"""
import datetime
import math
import os

import numpy as np

from satima.core.i18n import tr
from satima.geo.nav import Nav

try:
    import numba
except ImportError:
    numba = None

# ---- McIDAS nvxgvar 常數 ------------------------------------------------------
NOMORB = 42164.365              # 衛星標稱地心距離 (km)
AEC = 6378.137                  # 地球赤道半徑 (km)
FERC = 1. - 6356.7533 / AEC     # 扁率
AEBE2C = 1. / (1. - FERC) ** 2
INCMAX = 6136                   # 每圈 (cycle) 的步數
ELVINC = 8.e-6                  # 每步仰角 (rad)
SCNINC = 16.e-6                 # 每步掃描角 (rad)
ELVLN = 28.e-6                  # 每條可見光偵測器線的仰角 (rad)
SCNPX = 16.e-6                  # 每個可見光像素的掃描角 (rad)
ELVMAX0 = 0.220896              # 預設的南北中心 (= 4.5 x incmax x elvinc)
SCNMAX0 = 0.24544               # 預設的東西中心 (= 2.5 x incmax x scninc)
EWNOM = 0.24544
_ITER = 6                       # 反算疊代次數 (失準與 doff 都是 1e-3 以下的小量，3 次即收斂到 1e-12)

ATTITUDES = ('roll', 'pitch', 'yaw', 'roll_mis', 'pitch_mis')
SATELLITE = {8: 'GOES-8', 9: 'GOES-9', 10: 'GOES-10', 11: 'GOES-11', 12: 'GOES-12', 13: 'GOES-13', 14: 'GOES-14',
             15: 'GOES-15'}

# 紅外每條線 / 每個像素的中心，相對於其涵蓋的第一條可見光線 / 第一個可見光像素 (可見光單位)
# GVAR 的紅外行號 (= 該掃描最北的可見光行 INSLN + 4 x 偵測器) 是 4 條可見光線的第一條，中心在 +1.5
IR_LINE_CENTER = 1.5
IR_PIXEL_CENTER = 0.

# 頻道：(中心波長 µm：GOES-8~11、12~13、14~15；種類；功能代稱)
_WAVELENGTH = {1: (0.65, 0.65, 0.63), 2: (3.9, 3.9, 3.9), 3: (6.75, 6.5, 6.5), 4: (10.7, 10.7, 10.7), 5: (12.0, 12.0, 12.0),
               6: (13.3, 13.3, 13.3)}
_KIND = {1: ('vis', 'Visible'), 2: ('swir', 'Shortwave Window'), 3: ('wv', 'Water Vapor'), 4: ('ir', 'LW Window'),
         5: ('ir', 'Dirty LW Window'), 6: ('ir', 'CO2 LW')}


def make_channel(number, satellite, lines_per_row=4):
    """GVAR 頻道 1~6 -> Channel (名稱 CH1~CH6、顯示 BAND-01~06；解析度：可見光 1 km、紅外 4 km、8 km 偵測器 8 km)。"""
    from satima.geo.reader import Channel
    n = int(satellite.split('-')[1])
    wl = _WAVELENGTH[number][0 if n <= 11 else (1 if n <= 13 else 2)]
    kind, alias = _KIND[number]
    res = 1.0 if number == 1 else float(lines_per_row)
    return Channel('CH%d' % number, wl, res, kind, band='BAND-%02d' % number, alias=alias)


def sector_label(georange):
    """區域掃描的經緯度範圍 -> (名稱, 是否自動帶入地圖範圍)。北 / 南半球等大範圍不自動帶入 (請自行輸入範圍)。"""
    s, n, w, e = georange
    if n - s < 30. and e - w < 45.:
        return '', True
    if s > -30. and n > 45.:
        return 'NH', False
    if n < 30. and s < -45.:
        return 'SH', False
    return '', False


def fit_line_time(lines, times):
    """各線 (可見光行號) 的時間做穩健直線擬合 -> (行號 lines[0] 的時間, 每條可見光線的分鐘數)；不足時 None。"""
    pairs = [(float(l), t) for l, t in zip(lines, times) if t is not None]
    if len(pairs) < 10:
        return None
    ref = pairs[len(pairs) // 2][1]
    ll = np.array([p[0] for p in pairs])
    sec = np.array([(p[1] - ref).total_seconds() for p in pairs])
    keep = np.abs(sec - np.median(sec)) < 3600.
    for _ in range(3):
        slope, icpt = np.polyfit(ll[keep], sec[keep], 1)
        res = sec - (slope * ll + icpt)
        keep = np.abs(res) < max(3. * np.median(np.abs(res[keep])) + 1e-3, 0.5)
    l0 = float(lines[0])
    return ref + datetime.timedelta(seconds=float(slope * l0 + icpt)), slope / 60.


# ==================================================================================================
#  數值格式
# ==================================================================================================
def gould(words):
    """Gould 32 位元浮點數 (GVAR Block 0) -> float64 (同舊版 gould2native，可傳陣列)。"""
    u = np.asarray(words, dtype=np.int64) & 0xFFFFFFFF
    neg = (u >> 31) & 1
    v = np.where(neg == 1, ((u ^ 0xFFFFFFFF) + 1) & 0xFFFFFFFF, u)
    exp = (v >> 24) & 0x7F
    exp = np.where(exp == 0, 64, exp)
    val = (v & 0xFFFFFF) / 16.0 ** (70 - exp)
    out = np.where(neg == 1, -val, val)
    return float(out) if out.ndim == 0 else out


def bcd_time(raw):
    """GVAR 8 bytes BCD 時間 -> datetime (日數百位的最高位元是時鐘飛輪旗標，先遮掉)；無效時 None。"""
    b = bytes(raw)
    if len(b) < 8:
        return None
    n = [(x >> 4, x & 15) for x in b[:8]]
    try:
        year = n[0][0] * 1000 + n[0][1] * 100 + n[1][0] * 10 + n[1][1]
        doy = (n[2][0] & 7) * 100 + n[2][1] * 10 + n[3][0]
        hour = n[3][1] * 10 + n[4][0]
        minute = n[4][1] * 10 + n[5][0]
        sec = n[5][1] * 10 + n[6][0]
        msec = n[6][1] * 100 + n[7][0] * 10 + n[7][1]
        if not (1970 <= year <= 2100 and 1 <= doy <= 366 and hour < 24 and minute < 60 and sec < 61 and msec < 1000):
            return None
        return (datetime.datetime(year, 1, 1) + datetime.timedelta(days=doy - 1, hours=hour, minutes=minute,
                                                                   seconds=sec, milliseconds=msec))
    except ValueError:
        return None


def unpack10(raw):
    """GVAR 10 位元字組 (每 5 bytes 4 個字，高位元在前) -> uint16 陣列。"""
    b = np.frombuffer(raw, np.uint8)
    n = len(b) // 5
    g = b[:n * 5].reshape(n, 5).astype(np.uint16)
    out = np.empty((n, 4), np.uint16)
    out[:, 0] = (g[:, 0] << 2) | (g[:, 1] >> 6)
    out[:, 1] = ((g[:, 1] & 0x3F) << 4) | (g[:, 2] >> 4)
    out[:, 2] = ((g[:, 2] & 0x0F) << 6) | (g[:, 3] >> 2)
    out[:, 3] = ((g[:, 3] & 0x03) << 8) | g[:, 4]
    return out.ravel()


# ==================================================================================================
#  導航參數 (自然單位：弧度、km、分鐘)
# ==================================================================================================
def _series(mean, exp_mag, exp_tc, sins, monos):
    return dict(mean=mean, exp_mag=exp_mag, exp_tc=exp_tc, sin=list(sins), mono=list(monos))


def params_from_block0(b0, b0b):
    """GVAR Block 0 (gvar_format._Block_0 與 _Block_0_I_N / _Block_0_O 的紀錄) -> 導航參數。"""
    g = lambda name: gould(int(b0[name]))
    iscan = int(b0['ISCAN'])
    p = dict(imc=bool(iscan >> (31 - 8) & 1), flip=-1. if iscan >> (31 - 16) & 1 else 1.,
             nad=tuple(int(b0b[k]) for k in ('IOFNC', 'IOFEC', 'IOFNI', 'IOFEI')),
             lam=g('REFLO'), dr=g('REFRA'), phi=g('REFLA'), psi=g('REFOR'),
             roll=g('REROL'), pitch=g('REPIT'), yaw=g('REYAW'),
             epoch=bcd_time(b0['EPODT']), time=bcd_time(b0['TCHED']),
             chlon=gould(np.frombuffer(bytes(b0['CHLON']), '>u4')), chrad=gould(np.frombuffer(bytes(b0['CHRAD']), '>u4')),
             sinla=gould(np.frombuffer(bytes(b0['SINLA']), '>u4')), siyaw=gould(np.frombuffer(bytes(b0['SIYAW']), '>u4')),
             dayso=g('DAYSO'), start=g('START'), imc_corr=(g('COROL'), g('COPIT'), g('COYAW')), att={})
    for key, s in zip(ATTITUDES, ('R', 'P', 'Y', 'RM', 'PM')):
        n_sin = min(int(b0['SINAN' + s]), 15)
        n_mono = min(int(b0['NUMMS' + s]), 4)
        p['att'][key] = _series(
            g('MATAN' + s), g('EXMAG' + s), g('EXTIC' + s),
            [(g('MAG%dS%s' % (j, s)), g('PA%dSD%s' % (j, s))) for j in range(1, n_sin + 1)],
            [(int(b0['OAPS%d%s' % (l, s)]), int(b0['O%dMSD%s' % (l, s)]), g('MMSD%d%s' % (l, s)), g('PAMS%d%s' % (l, s)),
              g('AEMZ%d%s' % (l, s))) for l in range(1, n_mono + 1)])
    return p


def params_from_area(nav):
    """McIDAS AREA 的 GVAR 導航區塊 (gvar_area_format._Navigation 紀錄) -> 導航參数 (數值 / 10^7)。"""
    i = lambda name: int(nav[name])
    s7 = lambda name: i(name) / 1e7
    arr = lambda name: np.frombuffer(bytes(nav[name]), '<i4') / 1e7
    epoch_words = np.frombuffer(bytes(nav['EPOCH']), '<i4')
    p = dict(imc=bool(i('Imager Scan Status') >> 7 & 1), flip=-1. if i('W4') >> 15 & 1 else 1.,
             nad=(i('Instrument nadir, north/south cycles'), i('Instrument nadir, east/west cycles'),
                  i('Instrument nadir, north/south increments'), i('Instrument nadir, east/west increments')),
             lam=s7('Reference longitude'), dr=s7('Reference distance from nominal'), phi=s7('Reference latitude'),
             psi=s7('Reference yaw'), roll=s7('Reference attitude roll'), pitch=s7('Reference attitude pitch'),
             yaw=s7('Reference attitude yaw'), epoch=bcd_time(epoch_words.astype('>i4').tobytes()),
             time=None, chlon=arr('Longitude delta from reference values'),
             chrad=arr('Radial distance delta from reference values'), sinla=arr('Geocentric latitude delta values'),
             siyaw=arr('Orbit yaw delta values'), dayso=s7('Daily solar rate'),
             start=i('Exponential start time from epoch') / 100.,
             imc_corr=(s7('Image motion compensation roll'), s7('Image motion compensation pitch'),
                       s7('Image motion compensation yaw')), att={})
    names = {'roll': 'roll attitude angle', 'pitch': 'pitch attitude angle', 'yaw': 'yaw attitude angle',
             'roll_mis': 'roll misalignment angle', 'pitch_mis': 'pitch misalignment angle'}
    for key, pat in names.items():
        if 'Mean attitude angle for ' + pat not in nav.dtype.names:
            p['att'][key] = _series(0., 0., 1., [], [])
            continue
        n_sin = min(i('Number of sinusoids/angles for ' + pat), 15)
        n_mono = min(i('Number of monomial sinusoid for ' + pat), 4)
        p['att'][key] = _series(
            s7('Mean attitude angle for ' + pat), s7('Exponential magnitude for ' + pat),
            i('Exponential time constant for ' + pat) / 100.,
            [(s7('Magnitude of %d-order sinusoid for %s' % (j, pat)), s7('Phase angle of %d-order sinusoid for %s' % (j, pat)))
             for j in range(1, n_sin + 1)],
            [(i('%d order of applicable sinusoid for %s' % (l, pat)), i('Order of %d monomial sinusoid for %s' % (l, pat)),
              s7('Magnitude of %d monomial sinusoid for %s' % (l, pat)),
              s7('Phase angle of %d monomial sinusoid for %s' % (l, pat)),
              s7('Angle from epoch at daily solar rate (%d order) for %s' % (l, pat))) for l in range(1, n_mono + 1)])
    return p


def _attitude(a, wa, te):
    """姿態角級數 (McIDAS GATT)：常數 + 指數 + 太陽軌道角的正弦級數 + 單項式正弦。"""
    att = a['mean']
    if te >= 0 and a['exp_tc']:
        att += a['exp_mag'] * math.exp(-te / a['exp_tc'])
    for j, (mag, ph) in enumerate(a['sin']):
        att += mag * math.cos(wa * (j + 1) + ph)
    for order_sin, order_mono, mag, ph, ang in a['mono']:
        att += mag * (wa - ang) ** order_mono * math.cos(order_sin * wa + ph)
    return att


def minutes_since_epoch(p, time):
    return (time - p['epoch']).total_seconds() / 60. if (time and p.get('epoch')) else 0.


def nav_state(p, ts=None):
    """導航參數 -> 計算用的量 (衛星位置、轉換矩陣、掃描角範圍、失準)。
    ts：IMC 關閉時的時間 (自 epoch 起的分鐘)；None 時用參數中的影像時間。"""
    nsc, ewc, nsi, ewi = p['nad']
    if nsc and ewc and nsi and ewi:
        elvmax = (nsc * INCMAX + nsi) * ELVINC
        scnmax = (ewc * INCMAX + ewi) * SCNINC
    else:
        elvmax, scnmax = ELVMAX0, SCNMAX0
    lam, dr, phi, psi = p['lam'], p['dr'], p['phi'], p['psi']
    roll, pitch, yaw = p['roll'], p['pitch'], p['yaw']
    rma = pma = 0.
    if not p['imc']:
        if ts is None:
            ts = minutes_since_epoch(p, p.get('time'))
        dr = phi = psi = 0.
        w = 0.729115e-4 * 60. * ts
        sw, cw, sw1, cw1 = math.sin(w), math.cos(w), math.sin(0.927 * w), math.cos(0.927 * w)
        s2w, c2w, sw3, cw3 = math.sin(2 * w), math.cos(2 * w), math.sin(1.9268 * w), math.cos(1.9268 * w)
        cl, cr, sl, sy = p['chlon'], p['chrad'], p['sinla'], p['siyaw']
        lam = lam + cl[0] + (cl[1] + cl[2] * w) * w + (cl[9] * sw1 + cl[10] * cw1 + cl[3] * sw + cl[4] * cw + cl[5] * s2w +
                                                       cl[6] * c2w + cl[7] * sw3 + cl[8] * cw3 +
                                                       w * (cl[11] * sw + cl[12] * cw)) * 2.
        dr = dr + cr[0] + cr[1] * cw + cr[2] * sw + cr[3] * c2w + cr[4] * s2w + cr[5] * cw3 + cr[6] * sw3 + \
            cr[7] * cw1 + cr[8] * sw1 + w * (cr[9] * cw + cr[10] * sw)
        dlat = sl[0] + sl[1] * cw + sl[2] * sw + sl[3] * c2w + sl[4] * s2w + w * (sl[5] * cw + sl[6] * sw) + \
            sl[7] * cw1 + sl[8] * sw1
        phi = phi + dlat * (1. + dlat * dlat / 6.)
        dyaw = sy[0] + sy[1] * sw + sy[2] * cw + sy[3] * s2w + sy[4] * c2w + w * (sy[5] * sw + sy[6] * cw) + \
            sy[7] * sw1 + sy[8] * cw1
        psi = psi + dyaw * (1. + dyaw * dyaw / 6.)
        wa = p['dayso'] * ts
        te = ts - p['start']
        roll += _attitude(p['att']['roll'], wa, te)
        pitch += _attitude(p['att']['pitch'], wa, te)
        yaw += _attitude(p['att']['yaw'], wa, te)
        rma = _attitude(p['att']['roll_mis'], wa, te)
        pma = _attitude(p['att']['pitch_mis'], wa, te)
        roll, pitch, yaw = roll + p['imc_corr'][0], pitch + p['imc_corr'][1], yaw + p['imc_corr'][2]
    # 衛星經度、軌道傾角 -> 衛星座標 -> 地固座標的矩陣 b (ref: GOES-PCC-TM-2473)
    slat, syaw = math.sin(phi), math.sin(psi)
    sinoi = slat * slat + syaw * syaw
    cosoi = math.sqrt(1. - sinoi)
    sinoi = math.sqrt(sinoi)
    u = 0. if (slat == 0. and syaw == 0.) else math.atan2(slat, syaw)
    sinu, cosu = math.sin(u), math.cos(u)
    asc = lam - u
    sinasc, cosasc = math.sin(asc), math.cos(asc)
    b = np.array([[-cosasc * sinu - sinasc * cosu * cosoi, -sinasc * sinoi, -cosasc * cosu + sinasc * sinu * cosoi],
                  [-sinasc * sinu + cosasc * cosu * cosoi, cosasc * sinoi, -sinasc * cosu - cosasc * sinu * cosoi],
                  [cosu * sinoi, -cosoi, -slat]])
    r = (NOMORB + dr) / AEC
    xs = -b[:, 2] * r
    # 儀器 -> 衛星本體 (roll、pitch、yaw 小角度近似，同 McIDAS INST2E)
    rs, ps, ys = roll * roll, pitch * pitch, yaw * yaw
    rpy = np.array([[1. - 0.5 * (ps + ys), -yaw, pitch],
                    [yaw + pitch * roll, 1. - 0.5 * (ys + rs), -roll],
                    [-pitch + roll * yaw, roll + pitch * yaw, 1. - 0.5 * (ps + rs)]])
    bt = b @ rpy
    return dict(elvmax=elvmax, scnmax=scnmax, doff=scnmax - EWNOM, ff=float(p['flip']), rma=rma, pma=pma,
                xs=xs, bt=bt, binv=np.linalg.inv(bt),
                sublon=math.degrees(asc + math.atan2(cosoi * sinu, cosu)),
                sublat=math.degrees(math.atan(AEBE2C * math.tan(phi))), height=NOMORB + dr - AEC)


# ==================================================================================================
#  逐點計算 (numba；沒有 numba 時用 numpy)
#  狀態表 S (K x 28)：[elvmax, scnmax, doff, ff, rma, pma, aebe2c, xs(3), bt(9), binv(9)]，依時間線性內插
#  (IMC 啟動時 K = 1，與時間無關)。每點的時間 ts = ta + tb x (行號 - L0) (自 epoch 起的分鐘)。
# ==================================================================================================
NSTATE = 28


def _pack(st):
    return np.ascontiguousarray(np.r_[st['elvmax'], st['scnmax'], st['doff'], st['ff'], st['rma'], st['pma'], AEBE2C,
                                      st['xs'], st['bt'].ravel(), st['binv'].ravel()]).astype(np.float64)


if numba is not None:
    @numba.njit(cache=True)
    def _state_at(S, g0, dg, ts, out):
        K = S.shape[0]
        if K == 1:
            for j in range(NSTATE):
                out[j] = S[0, j]
            return
        f = (ts - g0) / dg
        i = int(math.floor(f))
        i = 0 if i < 0 else (K - 2 if i > K - 2 else i)
        w = f - i
        for j in range(NSTATE):
            out[j] = S[i, j] + w * (S[i + 1, j] - S[i, j])

    @numba.njit(cache=True)
    def _fwd1(L, P, s):
        elvmax, scnmax, doff, ff, rma, pma, k = s[0], s[1], s[2], s[3], s[4], s[5], s[6]
        x0, x1, x2 = s[7], s[8], s[9]
        a0 = elvmax - (L - 4.5) * ELVLN
        z0 = (P - 1.) * SCNPX - scnmax
        al = a0 - a0 * z0 * doff
        ze = z0 + 0.5 * a0 * a0 * doff
        ca, sa, cz = math.cos(al), math.sin(al), math.cos(ze)
        da = al - pma * sa * (ff / cz + math.tan(ze)) - rma * (1. - ca / cz)
        dz = ze + ff * rma * sa
        g0, g1, g2 = math.sin(dz), -math.cos(dz) * math.sin(da), math.cos(dz) * math.cos(da)
        h0 = s[10] * g0 + s[11] * g1 + s[12] * g2
        h1 = s[13] * g0 + s[14] * g1 + s[15] * g2
        h2 = s[16] * g0 + s[17] * g1 + s[18] * g2
        q1 = h0 * h0 + h1 * h1 + k * h2 * h2
        q2 = x0 * h0 + x1 * h1 + k * x2 * h2
        q3 = x0 * x0 + x1 * x1 + k * x2 * x2 - 1.
        d = q2 * q2 - q1 * q3
        if abs(d) < 1e-9:
            d = 0.
        if d < 0.:
            return np.nan, np.nan
        hh = -(q2 + math.sqrt(d)) / q1
        u0, u1, u2 = x0 + hh * h0, x1 + hh * h1, x2 + hh * h2
        d1 = u2 / math.sqrt(u0 * u0 + u1 * u1 + u2 * u2)
        return math.degrees(math.atan2(u1, u0)), math.degrees(math.atan(k * d1 / math.sqrt(1. - d1 * d1)))

    @numba.njit(cache=True)
    def _inv1(u0, u1, u2, s, iters):
        elvmax, scnmax, doff, ff, rma, pma, k = s[0], s[1], s[2], s[3], s[4], s[5], s[6]
        x0, x1, x2 = s[7], s[8], s[9]
        if (x0 - u0) * u0 + (x1 - u1) * u1 + k * (x2 - u2) * u2 <= 0.:     # 地球背面
            return np.nan, np.nan
        v0, v1, v2 = u0 - x0, u1 - x1, u2 - x2
        g0 = s[19] * v0 + s[20] * v1 + s[21] * v2
        g1 = s[22] * v0 + s[23] * v1 + s[24] * v2
        g2 = s[25] * v0 + s[26] * v1 + s[27] * v2
        n = math.sqrt(g0 * g0 + g1 * g1 + g2 * g2)
        dz = math.asin(g0 / n)
        da = math.atan2(-g1, g2)
        al, ze = da, dz
        for _ in range(iters):                                         # 失準修正的反函數
            ze = dz - ff * rma * math.sin(al)
            al = da + pma * math.sin(al) * (ff / math.cos(ze) + math.tan(ze)) + rma * (1. - math.cos(al) / math.cos(ze))
        a0, z0 = al, ze
        for _ in range(iters):                                         # doff 修正的反函數
            a0 = al / (1. - z0 * doff)
            z0 = ze - 0.5 * a0 * a0 * doff
        return 4.5 + (elvmax - a0) / ELVLN, 1. + (z0 + scnmax) / SCNPX

    @numba.njit(parallel=True, cache=True)
    def forward(line, pixel, S, g0, dg, ta, tb, L0, lon, lat):
        for i in numba.prange(line.shape[0]):
            if not (math.isfinite(line[i]) and math.isfinite(pixel[i])):
                lon[i], lat[i] = np.nan, np.nan
                continue
            s = np.empty(NSTATE)
            _state_at(S, g0, dg, ta + tb * (line[i] - L0), s)
            lon[i], lat[i] = _fwd1(line[i], pixel[i], s)

    @numba.njit(parallel=True, cache=True)
    def inverse(lon, lat, S, g0, dg, ta, tb, L0, iters, line, pixel):
        k = AEBE2C
        timed = S.shape[0] > 1
        for i in numba.prange(lon.shape[0]):
            line[i], pixel[i] = np.nan, np.nan
            if not (math.isfinite(lon[i]) and math.isfinite(lat[i])):
                continue
            psi = math.atan(math.tan(math.radians(lat[i])) / k)          # 大地緯度 -> 地心緯度
            lam = math.radians(lon[i])
            r = 1. / math.sqrt(math.cos(psi) ** 2 + k * math.sin(psi) ** 2)
            u0, u1, u2 = r * math.cos(psi) * math.cos(lam), r * math.cos(psi) * math.sin(lam), r * math.sin(psi)
            s = np.empty(NSTATE)
            ts = ta
            L, P = np.nan, np.nan
            for _ in range(4 if timed else 1):                          # IMC 關閉：依行號的掃描時間反覆修正
                _state_at(S, g0, dg, ts, s)
                L, P = _inv1(u0, u1, u2, s, iters)
                if not math.isfinite(L):
                    break
                ts = ta + tb * (L - L0)
            line[i], pixel[i] = L, P


def _states_np(S, g0, dg, ts):
    if S.shape[0] == 1:
        return np.broadcast_to(S[0], (np.size(ts), NSTATE))
    f = (np.asarray(ts) - g0) / dg
    i = np.clip(np.floor(f).astype(np.int64), 0, S.shape[0] - 2)
    w = (f - i)[:, None]
    return S[i] + w * (S[i + 1] - S[i])


def _forward_np(line, pixel, S, g0, dg, ta, tb, L0):
    s = _states_np(S, g0, dg, ta + tb * (line - L0))
    elvmax, scnmax, doff, ff, rma, pma = (s[:, j] for j in range(6))
    xs, bt = s[:, 7:10].T, s[:, 10:19].reshape(-1, 3, 3)
    a0 = elvmax - (line - 4.5) * ELVLN
    z0 = (pixel - 1.) * SCNPX - scnmax
    al = a0 - a0 * z0 * doff
    ze = z0 + 0.5 * a0 * a0 * doff
    ca, sa, cz = np.cos(al), np.sin(al), np.cos(ze)
    da = al - pma * sa * (ff / cz + np.tan(ze)) - rma * (1. - ca / cz)
    dz = ze + ff * rma * sa
    g = np.stack([np.sin(dz), -np.cos(dz) * np.sin(da), np.cos(dz) * np.cos(da)], axis=-1)
    h = np.einsum('nij,nj->in', bt, g)
    q1 = h[0] ** 2 + h[1] ** 2 + AEBE2C * h[2] ** 2
    q2 = xs[0] * h[0] + xs[1] * h[1] + AEBE2C * xs[2] * h[2]
    q3 = xs[0] ** 2 + xs[1] ** 2 + AEBE2C * xs[2] ** 2 - 1.
    d = q2 * q2 - q1 * q3
    d = np.where(np.abs(d) < 1e-9, 0., d)
    with np.errstate(invalid='ignore'):
        hh = -(q2 + np.sqrt(d)) / q1
        u = xs + hh * h
        d1 = u[2] / np.sqrt((u ** 2).sum(0))
        lat = np.degrees(np.arctan(AEBE2C * d1 / np.sqrt(1. - d1 * d1)))
        lon = np.degrees(np.arctan2(u[1], u[0]))
    bad = (d < 0) | ~np.isfinite(line) | ~np.isfinite(pixel)
    return np.where(bad, np.nan, lon), np.where(bad, np.nan, lat)


def _inverse_np(lon, lat, S, g0, dg, ta, tb, L0, iters=_ITER):
    ok = np.isfinite(lon) & np.isfinite(lat)
    psi = np.arctan(np.tan(np.radians(np.where(ok, lat, 0.))) / AEBE2C)
    lam = np.radians(np.where(ok, lon, 0.))
    r = 1. / np.sqrt(np.cos(psi) ** 2 + AEBE2C * np.sin(psi) ** 2)
    u = np.stack([r * np.cos(psi) * np.cos(lam), r * np.cos(psi) * np.sin(lam), r * np.sin(psi)])
    ts = np.full(lon.shape, float(ta))
    line = pixel = None
    for _ in range(4 if S.shape[0] > 1 else 1):
        s = _states_np(S, g0, dg, ts)
        elvmax, scnmax, doff, ff, rma, pma = (s[:, j] for j in range(6))
        xs, binv = s[:, 7:10].T, s[:, 19:28].reshape(-1, 3, 3)
        visible = ok & (((xs[0] - u[0]) * u[0] + (xs[1] - u[1]) * u[1] + AEBE2C * (xs[2] - u[2]) * u[2]) > 0)
        g = np.einsum('nij,jn->in', binv, u - xs)
        g = g / np.sqrt((g ** 2).sum(0))
        dz = np.arcsin(np.clip(g[0], -1, 1))
        da = np.arctan2(-g[1], g[2])
        al, ze = da, dz
        for _ in range(iters):
            ze = dz - ff * rma * np.sin(al)
            al = da + pma * np.sin(al) * (ff / np.cos(ze) + np.tan(ze)) + rma * (1. - np.cos(al) / np.cos(ze))
        a0, z0 = al, ze
        for _ in range(iters):
            a0 = al / (1. - z0 * doff)
            z0 = ze - 0.5 * a0 * a0 * doff
        line = np.where(visible, 4.5 + (elvmax - a0) / ELVLN, np.nan)
        pixel = np.where(visible, 1. + (z0 + scnmax) / SCNPX, np.nan)
        ts = np.where(np.isfinite(line), ta + tb * (line - L0), ta)
    return line, pixel


_NOTICE = set()


def _kernels(engine):
    if engine == 'cython':
        try:
            import goes_nav_cy
            return goes_nav_cy.forward, goes_nav_cy.inverse
        except ImportError:
            if 'cy' not in _NOTICE:
                _NOTICE.add('cy')
                print(tr('  [導航] 沒有編譯好的 GOES Cython 模組 (goes_nav_cy)，改用 numba；'
                         '編譯方式：python sensors/geo/build_vissr_nav_cy.py',
                         '  [navigation] the GOES Cython module (goes_nav_cy) is not built; using numba instead; '
                         'to build: python sensors/geo/build_vissr_nav_cy.py'))
            engine = 'numba'
    if engine == 'numba' and numba is not None:
        return forward, inverse
    return None


class GoesNav(Nav):
    """GOES Imager 導航。陣列第 r 列、第 c 行 = 可見光行號 first_line + r x line_step、像素號 first_pixel + c x pixel_step
    (紅外等粗解析度波段以可見光為單位，line_step = 4 或 8、pixel_step = 4)。
    line_time = (行號 L0 的掃描時間 (datetime), 每條可見光線的分鐘數)：IMC 關閉時依每條線的掃描時間計算
    (IMC 啟動或沒有給時整張影像用同一組參數)。"""

    def __init__(self, shape, params, first_line, first_pixel, line_step=1., pixel_step=1., key=None, engine='numba',
                 line_time=None):
        self.shape = (int(shape[0]), int(shape[1]))
        self.params = params
        self.first_line, self.first_pixel = float(first_line), float(first_pixel)
        self.line_step, self.pixel_step = float(line_step), float(pixel_step)
        self._key = key
        self._k = _kernels(engine)
        if params['imc'] or line_time is None or not params.get('epoch'):
            self.state = nav_state(params)
            self._S = _pack(self.state)[None, :]
            self._t = (0., 1., 0., 0., 0.)
        else:
            t0, per_line = line_time
            L0 = self.first_line
            ta = minutes_since_epoch(params, t0)
            tb = float(per_line)
            span = [ta + tb * (L - L0) for L in (self.first_line, self.first_line + self.shape[0] * self.line_step)]
            g0, dg = min(span) - 1., 0.1                                  # 每 6 秒一組參數，前後多留 1 分鐘
            n = int((max(span) + 1. - g0) / dg) + 2
            self._S = np.ascontiguousarray(np.stack([_pack(nav_state(params, g0 + i * dg)) for i in range(n)]))
            self._t = (g0, dg, ta, tb, L0)
            self.state = nav_state(params, 0.5 * (span[0] + span[1]))

    def key(self):
        return ('goes', self._key) if self._key is not None else id(self)

    def to_lonlat(self, row, col):
        row, col = np.broadcast_arrays(np.asarray(row, dtype=np.float64), np.asarray(col, dtype=np.float64))
        shape = row.shape
        line = np.ascontiguousarray(self.first_line + row.ravel() * self.line_step)
        pixel = np.ascontiguousarray(self.first_pixel + col.ravel() * self.pixel_step)
        g0, dg, ta, tb, L0 = self._t
        if self._k is not None:
            lon, lat = np.empty(line.size), np.empty(line.size)
            self._k[0](line, pixel, self._S, g0, dg, ta, tb, L0, lon, lat)
        else:
            lon, lat = _forward_np(line, pixel, self._S, g0, dg, ta, tb, L0)
        return lon.reshape(shape), lat.reshape(shape)

    def to_rowcol(self, lon, lat):
        lon, lat = np.broadcast_arrays(np.asarray(lon, dtype=np.float64), np.asarray(lat, dtype=np.float64))
        shape = lon.shape
        lon, lat = np.ascontiguousarray(lon.ravel()), np.ascontiguousarray(lat.ravel())
        g0, dg, ta, tb, L0 = self._t
        if self._k is not None:
            line, pixel = np.empty(lon.size), np.empty(lon.size)
            self._k[1](lon, lat, self._S, g0, dg, ta, tb, L0, _ITER, line, pixel)
        else:
            line, pixel = _inverse_np(lon, lat, self._S, g0, dg, ta, tb, L0)
        rows = (line - self.first_line) / self.line_step
        cols = (pixel - self.first_pixel) / self.pixel_step
        return rows.reshape(shape), cols.reshape(shape)


# ==================================================================================================
#  校正
# ==================================================================================================
_IR_TABLES = {}
# 可見光發射前係數 (NOAA OSPO：Weinreb & Han「Calibration of the Visible Channels of the GOES Imagers and Sounders」，
# 2009/03 修訂版)：radiance = m x count + b (b = -m x 29)，albedo = k x radiance；偵測器 1~8。
# GOES-14/15 (GOES O-P 格式) 的 GVAR Block 0 沒有可見光係數，用此表；GOES-13 的 Block 0 有係數，此表只在缺時備用。
_VIS_PRELAUNCH = {
    'GOES-13': {'m': [6.120196e-01, 6.118504e-01, 6.096360e-01, 6.087055e-01, 6.132860e-01, 6.118208e-01,
                      6.122307e-01, 6.066968e-01],
                'b': [-17.749, -17.744, -17.769, -17.653, -17.785, -17.743, -17.755, -17.594], 'k': 1.89544e-03},
    'GOES-14': {'m': [5.874693e-01, 5.865367e-01, 5.862807e-01, 5.864086e-01, 5.857146e-01, 5.852004e-01,
                      5.860814e-01, 5.841697e-01],
                'b': [-17.037, -17.010, -17.002, -17.006, -16.986, -16.971, -16.996, -16.941], 'k': 1.88772e-03},
    'GOES-15': {'m': [5.851966e-01, 5.879772e-01, 5.856793e-01, 5.854250e-01, 5.866992e-01, 5.836241e-01,
                      5.846555e-01, 5.843753e-01],
                'b': [-16.9707, -17.0513, -16.9847, -16.9773, -17.0143, -16.9251, -16.9550, -16.9469],
                'k': 1.88852e-03},
}
# 發射後校正 (可見光感測器在軌道上持續衰減，發射前係數會使影像偏暗)：
#   GOES-11/12/13/15：NOAA STAR 每月公布的修正值 C (radiance = C x 發射前 radiance)，
#     conversion_table/GOES/GOES_Imager_Vis_OpCal.txt (STAR GOES Imager Visible Channel Operational Calibration，
#     https://www.star.nesdis.noaa.gov/data/GOES15CalData/GOES_Imager_Vis_OpCal/)；
#     某月產生的 C 自次月第二個星期二起使用；早於表中第一筆者不修正，晚於最後一筆者沿用最後一筆。
#   GOES-8、GOES-10：NOAA OSPO「Implementation of the post-launch vicarious calibration of the GOES imager visible
#     channel」的線性衰減式 (d = 發射後天數)：GOES-8 1.192 x (1 + 0.0001688 d)、GOES-10 1.049 x (1 + 0.0001022 d)。
#   GOES-14：NOAA 沒有公布修正值序列；用發射後測試期間 (2009/12) 以 MODIS 比對得到的 C = 1.0695 (GOES-14 PLT 週報摘要，
#     未見於 NOAA 正式文件)，所有時間都用此值。
#   GOES-9 及表外期間：沒有官方值，不修正。
_LAUNCH = {'GOES-8': (datetime.datetime(1994, 4, 13), 1.192, 0.0001688),
           'GOES-10': (datetime.datetime(1997, 4, 25), 1.049, 0.0001022)}
_FIXED = {'GOES-14': (1.0695, ('GOES-14 發射後測試 2009/12 的 MODIS 比對值',          # 說明：(中文, 英文)，使用時以 tr 選擇
                               'MODIS comparison during the GOES-14 post-launch test, 2009/12'))}
_OPCAL = {}
_MONTHS = {m: i for i, m in enumerate(('jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug', 'sep', 'oct', 'nov',
                                       'dec'), 1)}


def _second_tuesday(year, month):
    d = datetime.datetime(year, month, 1)
    d += datetime.timedelta(days=(1 - d.weekday()) % 7)        # 第一個星期二
    return d + datetime.timedelta(days=7)


def _opcal_table(table_dir):
    """STAR 修正值表 -> {衛星: [(生效時間, C), ...] 依時間排序}。"""
    path = os.path.join(table_dir, 'GOES', 'GOES_Imager_Vis_OpCal.txt')
    if path in _OPCAL:
        return _OPCAL[path]
    out = {}
    if os.path.exists(path):
        sats = None
        for line in open(path, encoding='utf-8', errors='replace'):
            if 'GOES-11' in line and 'GOES-15' in line:
                sats = ['GOES-%d' % int(s.split('-')[1]) for s in line.split() if s.startswith('GOES-')]
                continue
            parts = line.split()
            if sats is None or len(parts) < 2 + len(sats) or parts[0][:3].lower() not in _MONTHS:
                continue
            try:
                year = int(parts[1])
            except ValueError:
                continue
            month = _MONTHS[parts[0][:3].lower()]
            start = _second_tuesday(year + month // 12, month % 12 + 1)          # 次月第二個星期二起生效
            for sat, val in zip(sats, parts[2:2 + len(sats)]):
                val = val.rstrip('*')
                if val not in ('-', ''):
                    out.setdefault(sat, []).append((start, float(val)))
        for sat in out:
            out[sat].sort()
    _OPCAL[path] = out
    return out


def vis_scale(cfg, satellite, time):
    """依 [visible] post_launch_calibration 回傳可見光的發射後校正倍數 (關閉時為 1)，並印出來源。"""
    if not cfg.visible.get('post_launch_calibration', True):
        return 1.0
    c, why = vis_post_launch(cfg.general.conversion_table_dir, satellite, time)
    key = (satellite, time.strftime('%Y%m'))
    if key not in _NOTICE:
        _NOTICE.add(key)
        print(tr('  可見光發射後校正：x %.3f (%s)', '  Visible post-launch calibration: x %.3f (%s)') % (c, why) if c != 1.
              else tr('  可見光未做發射後校正：%s', '  No visible post-launch calibration applied: %s') % why)
    return c


def vis_post_launch(table_dir, satellite, time):
    """可見光發射後校正倍數 -> (C, 說明)；沒有官方值時 (1.0, 原因)。"""
    if satellite in _LAUNCH:
        launch, c0, rate = _LAUNCH[satellite]
        d = (time - launch).days
        return c0 * (1. + rate * d), tr('NOAA OSPO 線性衰減式 (發射後 %d 天)',
                                        'NOAA OSPO linear degradation formula (%d days after launch)') % d
    if satellite in _FIXED:
        c, why = _FIXED[satellite]
        return c, tr(*why)
    rows = _opcal_table(table_dir).get(satellite)
    if not rows:
        return 1.0, tr('沒有 %s 的官方發射後校正值', 'no official post-launch calibration for %s') % satellite
    valid = [c for start, c in rows if start <= time]
    if not valid:
        return 1.0, tr('%s 的 STAR 修正值自 %s 起才有', 'STAR corrections for %s are available only from %s') % (
            satellite, rows[0][0].strftime('%Y/%m'))
    return valid[-1], tr('NOAA STAR 每月修正值', 'NOAA STAR monthly correction')


def vis_coefficients(satellite, bias, gain1, gain2, albedo):
    """可見光 8 個偵測器的係數 -> (bias, gain1, gain2, albedo)；GVAR 中沒有 (全為 0) 時用 NOAA 發射前係數表。"""
    bias, gain1, gain2 = (np.asarray(x, dtype=np.float64) for x in (bias, gain1, gain2))
    if np.any(gain1) and albedo:
        return bias, gain1, gain2, float(albedo)
    t = _VIS_PRELAUNCH.get(satellite)
    if t is None:
        raise ValueError(tr('%s 的 GVAR 沒有可見光校正係數，也沒有 NOAA 發射前係數表',
                            'the GVAR of %s has no visible calibration coefficients, and there is no NOAA pre-launch '
                            'coefficient table for it') % satellite)
    return np.array(t['b']), np.array(t['m']), np.zeros(8), t['k']


C1, C2 = 1.191066e-5, 1.438833     # 輻射常數 (NOAA OSPO：mW/(m2 sr cm-4)、K cm)


def _teff(rad, n):
    with np.errstate(divide='ignore', invalid='ignore'):
        return C2 * n / np.log1p(C1 * n ** 3 / rad)


def _fit_planck(rad, tb, n0):
    """表中的 (輻射, 亮溫) -> (n, a, b)：T = a + b x Teff(n, R) (NOAA OSPO 的換算形式)，n 在標稱波數 ±5% 內以黃金分割搜尋。"""
    def sse(n):
        te = _teff(rad, n)
        A = np.c_[np.ones_like(te), te]
        coef = np.linalg.lstsq(A, tb, rcond=None)[0]
        return float(np.sum((A @ coef - tb) ** 2)), coef
    g = (math.sqrt(5.) - 1.) / 2.
    lo, hi = 0.95 * n0, 1.05 * n0
    x1, x2 = hi - g * (hi - lo), lo + g * (hi - lo)
    f1, f2 = sse(x1)[0], sse(x2)[0]
    for _ in range(60):
        if f1 < f2:
            hi, x2, f2 = x2, x1, f1
            x1 = hi - g * (hi - lo)
            f1 = sse(x1)[0]
        else:
            lo, x1, f1 = x1, x2, f2
            x2 = lo + g * (hi - lo)
            f2 = sse(x2)[0]
    n = 0.5 * (lo + hi)
    a, b = sse(n)[1]
    return n, a, b


def ir_lut(table_dir, satellite, side, channel, detector):
    """紅外計數值 (0~1023) -> 亮溫 (°C) 查表 (conversion_table/GOES/<衛星>_IMAGER_SIDE<s>_CH<c>_DET<d>.txt，
    NOAA OSPO 的表，與 OSPO 公布的最新係數相差 < 0.02 K)。
    該 side / 偵測器沒有表時依序改用另一個 side、偵測器 1；都沒有時回傳 None。
    NOAA 的表只列 180~340 K，計數值仍有效的高溫端 (3.9 µm 白天的沙漠、反光可到 343 K) 以表本身擬合的
    T = a + b x Teff(n, R) 延伸 (誤差 < 0.02 K)；低溫端不延伸 (低於 180 K 的計數值是太空)。"""
    sat = 'GOES-%02d' % int(satellite.split('-')[1])
    for s, d in ((side, detector), (3 - side, detector), (side, 1), (3 - side, 1)):
        path = os.path.join(table_dir, 'GOES', '%s_IMAGER_SIDE%d_CH%d_DET%d.txt' % (sat, s, channel, d))
        if path in _IR_TABLES:
            return _IR_TABLES[path]
        if os.path.exists(path):
            t = np.loadtxt(path)
            lut = np.full(1024, np.nan, dtype=np.float32)
            cnt = t[:, 0].astype(int)
            ok = (cnt >= 0) & (cnt < 1024) & (t[:, 2] > 0)
            lut[cnt[ok]] = t[ok, 2] - 273.15
            top = cnt[ok].max() if ok.any() else 1023
            hot = (cnt > top) & (cnt < 1024) & (t[:, 1] > 0)
            if hot.any():
                warm = ok & (t[:, 2] > 0.5 * (t[ok, 2].min() + t[ok, 2].max()))
                n, a, b = _fit_planck(t[warm, 1], t[warm, 2], 1e4 / _WAVELENGTH[channel][0 if int(sat[5:]) <= 11 else (
                    1 if int(sat[5:]) <= 13 else 2)])
                lut[cnt[hot]] = a + b * _teff(t[hot, 1], n) - 273.15
            _IR_TABLES[path] = lut
            return lut
    return None


