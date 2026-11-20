'''太陽與衛星觀測幾何 (以 pyorbital 計算，與 satpy 的作法相同)'''
import numpy as np
from pyorbital.astronomy import cos_zen, get_alt_az
from pyorbital.orbital import get_observer_look


def cos_sun_zenith(utc_time, lon, lat):
    '''太陽天頂角的餘弦值'''
    with np.errstate(invalid='ignore'):
        return cos_zen(utc_time, lon, lat)


def sun_angles(utc_time, lon, lat):
    '''回傳 (太陽天頂角, 太陽方位角)，單位為度，方位角範圍 0~360'''
    with np.errstate(invalid='ignore'):
        sunz = np.rad2deg(np.arccos(np.clip(cos_zen(utc_time, lon, lat), -1., 1.)))
        suna = np.rad2deg(get_alt_az(utc_time, lon, lat)[1]) % 360.
    return sunz, suna


def satellite_angles(sat_lon, sat_alt, utc_time, lon, lat, sat_lat=0.):
    '''回傳 (衛星天頂角, 衛星方位角)，單位為度。sat_alt 為衛星離地高度 (km)'''
    with np.errstate(invalid='ignore'):
        sata, satel = get_observer_look(sat_lon, sat_lat, sat_alt, utc_time, lon, lat, 0.)
    return 90. - satel, sata


def relative_azimuth(suna, sata):
    '''太陽與衛星的相對方位角 (0~180 度，0 為背向散射)'''
    diff = np.abs(suna - sata)
    return np.minimum(diff, 360. - diff)
