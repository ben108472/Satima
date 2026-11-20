"""真實色共用演算法

各衛星原本各自複製了一份相同的 CREFL 程式；數學上完全相同的部分集中在這裡，
各格式處理腳本保留原函式名稱 (例如 VIIRS 的 crefl_correct、MODIS 的 _chand) 並呼叫此處。
各衛星特有的部分 (波段係數、水氣吸收公式、地形修正、銳化方式、增強微調) 仍留在各自的腳本。

  csalbr / chand        NASA CREFL (Vermote, 6S) Rayleigh 球面反照率與路徑反射率
  crefl_core            CREFL 合成：(ρ/tO3 − ρray) / (T↑ T↓ tH2O)，再做多次散射 (球面反照率) 修正
  o3_transmittance      臭氧穿透率 exp(−m·UO3·aO3)
  h2o_transmittance_*   水氣穿透率 (VIIRS / MODIS 兩種 CREFL 公式)
  rayleigh_tau          Hansen & Travis (1974) 海平面 Rayleigh 光學厚度
  sunz_corr             太陽天頂角校正 1/cos(θs)，88° 起漸減、95° 歸零 (同 satpy sunz_corrected)
  cira_stretch          CIRA 對數增強 (satpy true_color 預設)
  worldview_stretch     NASA Worldview / MODIS Rapid Response 非線性增強
  srgb_encode           sRGB 顯示編碼 (IEC 61966-2-1)
"""
import numpy as np

SCALE_HEIGHT = 8000.      # m，CREFL 地形修正用的 Rayleigh 標高
MAXAIRMASS = 18.          # 空氣質量上限 (晨昏線附近避免爆量)
REFLMIN, REFLMAX = -0.01, 1.6

# NASA Worldview / MODIS Rapid Response 增強表 (Gumley, Descloitres & Schmaltz 2010, Table 3)
WV_REFL_MAX = 1.1
WV_XP = (0., 30., 60., 120., 190., 255.)
WV_FP = (0., 110., 160., 210., 240., 255.)


def rayleigh_tau(wl_um):
    """海平面 Rayleigh 光學厚度，Hansen & Travis (1974)；wl_um 單位 µm。"""
    return 0.008569 * wl_um ** -4 * (1. + 0.0113 * wl_um ** -2 + 0.00013 * wl_um ** -4)


def csalbr(tau):
    """Rayleigh 層球面反照率 (CREFL csalbr，指數積分近似)。"""
    a = [-.57721566, 0.99999193, -0.24991055, 0.05519968, -0.00976004, 0.00107857]
    fintexp1 = np.polyval(a[::-1], tau) - np.log(tau)
    fintexp3 = (np.exp(-tau) * (1.0 - tau) + tau ** 2 * fintexp1) / 2.0
    return (3.0 * tau - fintexp3 * (4.0 + 2.0 * tau) + 2.0 * np.exp(-tau)) / (4.0 + 3.0 * tau)


def chand(phi, muv, mus, taur):
    """Chandrasekhar/6S Rayleigh 路徑反射率與上下行直射穿透率 (CREFL chand)。
    phi: 太陽與衛星方位角差 (度)；mus/muv: 太陽/衛星天頂角餘弦；taur: Rayleigh 光學厚度。"""
    xfd = 0.958725775            # 去偏振因子 0.0279 所得
    xbeta2 = 0.5
    as0 = [0.33243832, 0.16285370, -0.30924818, -0.10324388, 0.11493334,
           -6.777104e-02, 1.577425e-03, -1.240906e-02, 3.241678e-02, -3.503695e-02]
    as1 = [0.19666292, -5.439061e-02]
    as2 = [0.14545937, -2.910845e-02]
    xph1 = 1.0 + (3.0 * mus * mus - 1.0) * (3.0 * muv * muv - 1.0) * xfd / 8.0
    xph2 = -xfd * xbeta2 * 1.5 * mus * muv * np.sqrt(1.0 - mus * mus) * np.sqrt(1.0 - muv * muv)
    xph3 = xfd * xbeta2 * 0.375 * (1.0 - mus * mus) * (1.0 - muv * muv)
    pl = (1., mus + muv, mus * muv, mus * mus + muv * muv, mus * mus * muv * muv)
    fs01 = sum(p * c for p, c in zip(pl, as0[:5]))
    fs02 = sum(p * c for p, c in zip(pl, as0[5:]))
    xlntaur = np.log(taur)
    fs0 = fs01 + fs02 * xlntaur
    fs1 = as1[0] + xlntaur * as1[1]
    fs2 = as2[0] + xlntaur * as2[1]
    trdown = np.exp(-taur / mus)
    trup = np.exp(-taur / muv)
    xitm1 = (1.0 - trdown * trup) / 4.0 / (mus + muv)
    xitm2 = (1.0 - trdown) * (1.0 - trup)
    phios = np.deg2rad(phi + 180.0)
    rhoray = (xph1 * (xitm1 + xitm2 * fs0)
              + 2.0 * xph2 * (xitm1 + xitm2 * fs1) * np.cos(phios)
              + 2.0 * xph3 * (xitm1 + xitm2 * fs2) * np.cos(2.0 * phios))
    return rhoray, trdown, trup


def o3_transmittance(air_mass, uo3, ao3):
    """臭氧穿透率；aO3 = 0 的波段回傳 1。"""
    if not ao3:
        return 1.0
    return np.exp(-air_mass * uo3 * ao3)


def h2o_transmittance_viirs(air_mass, uh2o, ah2o, bh2o):
    """VIIRS CREFL 水氣穿透率 exp(−aH2O·(m·UH2O)^bH2O)。"""
    if not ah2o or not bh2o:
        return 1.0
    return np.exp(-(ah2o * ((air_mass * uh2o) ** bh2o)))


def h2o_transmittance_modis(air_mass, uh2o, ah2o, bh2o):
    """MODIS CREFL 水氣穿透率 exp(−exp(aH2O + bH2O·ln(m·UH2O)))。"""
    if not bh2o:
        return 1.0
    return np.exp(-np.exp(ah2o + bh2o * np.log(air_mass * uh2o)))


def crefl_core(refl, mus, muv, phi, taur, t_o3=1.0, t_h2o=1.0):
    """CREFL 校正核心 (各衛星共用)。refl 需為已做太陽天頂角校正的 TOA 反射率。"""
    rhoray, trdown, trup = chand(phi, muv, mus, taur)
    t_up = ((2. / 3. + muv) + (2. / 3. - muv) * trup) / (4. / 3. + taur)
    t_down = ((2. / 3. + mus) + (2. / 3. - mus) * trdown) / (4. / 3. + taur)
    corr = (refl / t_o3 - rhoray) / (t_up * t_down * t_h2o)
    return corr / (1. + corr * csalbr(taur))


def sunz_corr(sza, limit=88., max_sza=95.):
    """太陽天頂角校正係數 1/cos(θs)；limit 以上逐漸減弱，max_sza 以上為 0 (同 satpy sunzen_corr_cos)。
    sza 為 NaN 的像元回傳 0。"""
    sza = np.asarray(sza, dtype=np.float64)
    limit_rad, max_rad = np.deg2rad(limit), np.deg2rad(max_sza)
    sza_rad = np.deg2rad(sza)
    cos_zen = np.cos(sza_rad)
    with np.errstate(invalid='ignore', divide='ignore'):
        grad = 1. - np.log((sza_rad - limit_rad) / (max_rad - limit_rad) + 1) / np.log(2)
        corr = np.where(cos_zen > np.cos(limit_rad), 1. / cos_zen,
                        np.clip(np.nan_to_num(grad), 0., None) / np.cos(limit_rad))
    corr = np.asarray(corr)
    corr[np.isnan(sza)] = 0.
    return corr


def srgb_encode(x):
    """標準 sRGB 顯示編碼 (IEC 61966-2-1)，線性反射率轉為螢幕亮度。"""
    x = np.nan_to_num(x)
    return np.where(x <= 0.0031308, 12.92 * x, 1.055 * np.power(np.clip(x, 0., None), 1. / 2.4) - 0.055)


def cira_stretch(refl, normalize_white=False):
    """CIRA 對數增強 (satpy cira_stretch)，輸入反射率 0–1，NaN 保留。
    normalize_white=True：標準 CIRA 把反射率 1.0 映射到約 0.83 (雲頂偏灰)，再除以該值使 1.0 = 滿白 (OLCI 微調)。"""
    log_root = np.log10(0.0223)
    denom = (1.0 - log_root) * 0.75
    with np.errstate(invalid='ignore', divide='ignore'):
        out = (np.log10(np.clip(refl, np.finfo(np.float32).eps, None)) - log_root) / denom
    if normalize_white:
        out = out / (-log_root / denom)
    return np.clip(out, 0., 1.)


def worldview_stretch(refl):
    """NASA Worldview 非線性增強：反射率 0–1.1 線性對應 0–255，再分段線性查表；NaN 保留。"""
    x = np.clip(np.nan_to_num(refl, nan=0.), 0., WV_REFL_MAX) * (255. / WV_REFL_MAX)
    out = np.interp(x, WV_XP, WV_FP) / 255.
    return np.where(np.isfinite(refl), out, np.nan)
