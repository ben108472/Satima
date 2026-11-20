"""JPSS ATMS 的 SDR 亮溫與 TDR 天線溫度 (NOAA CLASS / NOAA 開放資料的 SATMS_*、TATMS_*，HDF5)

檔名  <SATMS|TATMS>_<npp|j01|j02>_d<YYYYMMDD>_t<HHMMSSs>_e<HHMMSSs>_b<軌道>_c<產製時間>_<來源>_ops.h5
      npp = S-NPP、j01 = NOAA-20、j02 = NOAA-21；SATMS = SDR (亮溫)、TATMS = TDR (天線溫度)
      一個檔 (granule) 只有 32 秒、12 條掃描 x 96 點。同衛星、時間相接 (間隔 10 分鐘以內) 的檔自動接成一筆，
      資料列表的備註欄顯示檔數；中間缺檔照樣接、另外標出缺幾個 (接續的判定見 satima/core/segments.py)。
      同一個 granule 有多個檔 (重新下載、不同產製時間) 時只用產製時間最新的。
定位  **經緯度不在 SATMS / TATMS 裡**，在另外的 GATMO_* 檔 (ATMS-SDR-GEO)：須放在同一個資料夾，
      以「衛星_d日期_t起_e訖_b軌道」對應 (產製時間 c… 不必相同)；找不到定位檔的 granule 會略過。
      NOAA CLASS 訂資料時可一併勾選 GATMO；NOAA 的開放資料 (AWS S3，不需帳號) 也有：
        https://noaa-nesdis-n20-pds.s3.amazonaws.com/ATMS-SDR-GEO/年/月/日/   (S-NPP：noaa-nesdis-snpp-pds、
        NOAA-21：noaa-nesdis-n21-pds；同一個 bucket 的 ATMS-SDR/、ATMS-TDR/ 是 SATMS、TATMS)
      CLASS 把定位與資料包在同一個檔的合併檔 (GATMO-SATMS_*、GATMO-TATMS_*、GATMO-SATMS-TATMS_* 等，
      檔名前段列出所含的產品) 也可讀 (同樣的資料集；沒有樣本，未實測)。
      舊格式的 GATMO 沒有各頻段的波束經緯度 (BeamLatitude) 時，改用單一組的 Latitude / Longitude。
結構  All_Data/ATMS-SDR_All/BrightnessTemperature (掃描, 96, 22) uint16，x BrightnessTemperatureFactors[0] + [1] = K
      All_Data/ATMS-TDR_All/AntennaTemperature   同上，係數為 AntennaTemperatureFactors；65528 以上為無效值
      All_Data/ATMS-SDR-GEO_All/BeamLatitude、BeamLongitude (掃描, 96, 5)：K、Ka、V、W、G 五個頻段各一組
        (各頻段的波束指向略有不同，相差約 1 km)；StartTime (掃描,) 為 IET (1958 起算、含閏秒的微秒數)
頻道  22 個 (QV / QH = 準垂直 / 準水平極化，隨掃描角變化)：
        K  1 23.8QV                Ka 2 31.4QV
        V  3 50.3QH  4 51.76QH  5 52.8QH  6 53.596±0.115QH  7 54.4QH  8 54.94QH  9 55.5QH
           10 57.29QH (fo = 57.290344)  11 fo±0.217  12~15 fo±0.3222±0.048 / 0.022 / 0.010 / 0.0045   (溫度探空)
        W  16 88.2QV
        G  17 165.5QH  18~22 183.31±7 / 4.5 / 3 / 1.8 / 1 QH   (水氣探空)
      橫向掃描：星下點取樣約 16 km，往掃描帶邊緣變大；沒有雙極化，只有單頻道亮溫。
TDR 是天線溫度 (未做天線場型修正)：頻道附註 TA。
"""
import os
import re

import h5py
import numpy as np

from satima.core import segments
from satima.core.i18n import tr
from satima.mw.reader import Channel, Geo, Scene, SwathReader, clean_latlon, from_seconds, iet_to_unix

# <產品[-產品…]>_<衛星>_d<日期>_t<起>_e<訖>_b<軌道>_c<產製時間>_…；產品為 GATMO (定位)、SATMS (SDR)、TATMS (TDR)
_NAME = re.compile(r'^(?P<prods>(?:GATMO|SATMS|TATMS)(?:-(?:GATMO|SATMS|TATMS))*)_(?P<sat>[a-z0-9]+)_d(?P<date>\d{8})'
                   r'_t(?P<t>\d{7})_e(?P<e>\d{7})_b(?P<orbit>\d{5})_(?:c(?P<created>\d+))?', re.I)
_SATELLITES = {'npp': 'S-NPP', 'j01': 'NOAA-20', 'j02': 'NOAA-21', 'j03': 'NOAA-22', 'j04': 'NOAA-23'}
_KIND = {'S': ('SDR', 'All_Data/ATMS-SDR_All', 'BrightnessTemperature'),
         'T': ('TDR', 'All_Data/ATMS-TDR_All', 'AntennaTemperature')}
_GEO = 'All_Data/ATMS-SDR-GEO_All'
_FILL = 65528
GAP = 600.              # 兩個 granule 的間隔在幾秒內接成同一筆
# (名稱, 中心頻率 GHz, 極化, 頻段索引：K、Ka、V、W、G = 0~4)
CHANNELS = (
    ('23.8', 23.8, 'QV', 0), ('31.4', 31.4, 'QV', 1),
    ('50.3', 50.3, 'QH', 2), ('51.76', 51.76, 'QH', 2), ('52.8', 52.8, 'QH', 2), ('53.596±0.115', 53.596, 'QH', 2),
    ('54.4', 54.4, 'QH', 2), ('54.94', 54.94, 'QH', 2), ('55.5', 55.5, 'QH', 2), ('57.29', 57.29, 'QH', 2),
    ('57.29±0.217', 57.29, 'QH', 2), ('57.29±0.322±0.048', 57.29, 'QH', 2), ('57.29±0.322±0.022', 57.29, 'QH', 2),
    ('57.29±0.322±0.01', 57.29, 'QH', 2), ('57.29±0.322±0.0045', 57.29, 'QH', 2),
    ('88.2', 88.2, 'QV', 3),
    ('165.5', 165.5, 'QH', 4), ('183.31±7', 183.31, 'QH', 4), ('183.31±4.5', 183.31, 'QH', 4),
    ('183.31±3', 183.31, 'QH', 4), ('183.31±1.8', 183.31, 'QH', 4), ('183.31±1', 183.31, 'QH', 4),
)


def _granule_times(m):
    """檔名的 d、t、e (HHMMSS + 十分之一秒) -> (起, 訖) 1970 起算的秒數。"""
    import datetime
    day = datetime.datetime.strptime(m['date'], '%Y%m%d')

    def sec(s):
        return int(s[:2]) * 3600 + int(s[2:4]) * 60 + int(s[4:6]) + int(s[6]) / 10.

    a, b = sec(m['t']), sec(m['e'])
    if b < a:
        b += 86400.
    base = (day - datetime.datetime(1970, 1, 1)).total_seconds()
    return base + a, base + b


class ATMSReader(SwathReader):
    NAME = 'ATMS-SDR'
    SENSOR = 'ATMS'
    SOURCE = 'CLASS'
    DESCRIPTION = 'JPSS ATMS SDR 亮溫 / TDR 天線溫度 (NOAA CLASS SATMS、TATMS；需同資料夾的 GATMO 定位檔)'
    DESCRIPTION_EN = ('JPSS ATMS SDR brightness temperature / TDR antenna temperature (NOAA CLASS SATMS, TATMS; '
                      'needs the GATMO geolocation files in the same folder)')
    PATTERNS = ('SATMS_*.h5', 'TATMS_*.h5', 'SATMS-*.h5', 'TATMS-*.h5', 'GATMO-*.h5')

    def scenes(self, folder='.'):
        # 每個檔所含的產品：定位 (GATMO) 與資料 (SATMS / TATMS) 可能分開，也可能包在同一個檔
        geo, data = {}, {}          # granule -> (產製時間, 檔)；(種類, 衛星) -> {granule: (產製時間, 起, 訖, 軌道, 檔)}
        for name in sorted(os.listdir(folder)):
            m = _NAME.match(name)
            if not m or not name.lower().endswith('.h5'):
                continue
            path = os.path.join(folder, name)
            created = m['created'] or ''
            key = (m['sat'].lower(), m['date'], m['t'], m['e'], m['orbit'])
            start, end = _granule_times(m)
            for prod in m['prods'].upper().split('-'):
                if prod == 'GATMO':
                    if key not in geo or created > geo[key][0]:
                        geo[key] = (created, path)
                else:
                    items = data.setdefault((prod[0], m['sat'].lower()), {})
                    if key not in items or created > items[key][0]:
                        items[key] = (created, start, end, m['orbit'], path)
        out, nogeo = [], 0
        for (kind, sat), items in sorted(data.items()):
            grans = []
            for key, (_, start, end, orbit, path) in items.items():
                if key in geo:
                    grans.append((start, end, orbit, path, geo[key][1]))
                else:
                    nogeo += 1
            for chain, missing in segments.chain(grans, lambda g: g[0], lambda g: g[1], max_gap=GAP):
                out.append(Scene(_SATELLITES.get(sat, sat.upper()), 'ATMS', from_seconds(chain[0][0]),
                                 from_seconds(max(c[1] for c in chain)), file=chain[0][3], orbit=chain[0][2],
                                 format=_KIND[kind][0],
                                 extra={'kind': kind, 'granules': [(c[3], c[4]) for c in chain],
                                        'nfiles': len(chain), 'missing': missing}))
        if nogeo:
            print(tr('  [ATMS] %d 個 SATMS / TATMS 檔找不到對應的 GATMO 定位檔 (須在同一個資料夾)，已略過。',
                     '  [ATMS] %d SATMS / TATMS file(s) have no matching GATMO geolocation file (it must be in the '
                     'same folder); skipped.') % nogeo)
        return out

    def channels(self, scene):
        note = 'TA' if scene.extra['kind'] == 'T' else ''
        return [Channel(name + pol, freq, pol, str(band), key=k, note=note)
                for k, (name, freq, pol, band) in enumerate(CHANNELS)]

    def geo(self, scene, group):
        band = int(group)
        lats, lons, times = [], [], []
        for _, gpath in scene.extra['granules']:
            with h5py.File(gpath, 'r') as f:
                g = f[_GEO]
                if 'BeamLatitude' in g:
                    lats.append(g['BeamLatitude'][:, :, band])
                    lons.append(g['BeamLongitude'][:, :, band])
                else:
                    lats.append(g['Latitude'][:])
                    lons.append(g['Longitude'][:])
                times.append(iet_to_unix(g['StartTime'][:]))
        lat, lon = clean_latlon(np.concatenate(lats), np.concatenate(lons))
        return [Geo(lat, lon, np.concatenate(times), (scene.file, len(scene.extra['granules']), band))]

    def tb(self, scene, channel):
        _, group, var = _KIND[scene.extra['kind']]
        parts = []
        for path, _ in scene.extra['granules']:
            with h5py.File(path, 'r') as f:
                g = f[group]
                raw = g[var][:, :, channel.key]
                scale, offset = (float(v) for v in g[var + 'Factors'][:2])
            data = raw.astype(np.float32) * np.float32(scale) + np.float32(offset)
            data[(raw >= _FILL) | (raw == 0)] = np.nan
            parts.append(data)
        return [np.concatenate(parts)]


READER = ATMSReader
