"""極軌成像儀 (PO) 自己的設定：預設值與檢查

共同設定 (目錄、地圖、海岸線、圖片大小、標題、標註、檔名 …) 的預設值與讀取方式見 satima/core/config.py；
這裡的 DEFAULTS 只列 PO 才有的幾節，以及與共同預設值不同的項目。設定檔為 Satima.toml + Satima-PO.toml。
"""
from ..core.config import (HOME, MAX_SCALE, _is_num, check_colormaps, effective, fail, image_size,  # noqa: F401
                           unwrap)

KIND = 'po'
PROJECTIONS = ('merc', 'cyl', 'lcc')

_BW = {'name': 'BW', 'mpl': 'gray_r', 'vmin': -80.0, 'vmax': 50.0}

DEFAULTS = {
    # ---- 與共同預設值不同的項目 ----
    'general': {
        'range_warning': 20.0,               # 緯度 / 經度範圍 (中心往外的半寬) 超過此值時提醒記憶體用量；0 = 不提醒
    },
    'coastline': {'color': '#00BB00', 'linewidth': 0.5},
    'image': {'max_size': 10000},
    'title': {'ask': True},                  # 原本的 Picname：預設會詢問
    # ---- PO 才有的設定 ----
    'swath': {
        'merge': 'ask',                      # 一次選到時間相接的多個長檔 (整圈、半圈) 時合併成一張："ask" / true / false
                                             # (幾分鐘以內的短片段一律自動接成一筆，不受這裡影響)
        'merge_gap': 600.0,                  # 兩個檔的間隔在幾秒內視為同一段
    },
    'eyetemp': {
        'enable': True,                      # 紅外線頻道：印出中心點附近的最高亮溫
        'radius_deg': 0.5,
        'legacy_counts': False,
    },
    'colormap': {
        'ir': [dict(_BW),
               {'name': 'BD', 'file': 'bd'},
               {'name': 'CycloneCenter', 'file': 'cc'},
               {'name': 'RBTOP', 'file': 'rbtop'},
               {'name': 'CWB NPP_VIIRS', 'file': 'cwbm'},
               {'name': 'NRLCOLOR', 'file': 'color'},
               {'name': 'NHC', 'file': 'nhc'},
               {'name': 'RAMMB', 'file': 'rammb'},
               {'name': 'OTT', 'file': 'ott'},
               {'name': 'AVN', 'file': 'avn'},
               {'name': 'RB', 'file': 'rb'},
               {'name': 'MB', 'file': 'mb'},
               {'name': 'CA', 'file': 'ca'}],
        'wv': [dict(_BW),
               {'name': 'NRL', 'file': 'nrl'},
               {'name': 'SSD', 'file': 'ssd'},
               {'name': 'NOAA', 'file': 'noaa'}],
    },
}

# Satima-PO 舊設定檔的鍵 -> 整合版的寫法
_MOVED = (('output', 'facecolor', '[image] background'), ('output', 'ask_title', '[title] ask'),
          ('output', 'title', '[title] default'), ('coastline', 'resolution', '[map] resolution'))


def validate(cfg):
    """PO 各節的檢查 (共同的部分已由 core.config.validate 檢查)。"""
    for section, key, new in _MOVED:
        if key in cfg.get(section, {}):
            fail('[%s] %s 已改為 %s (Satima-PO 舊設定檔的寫法)', '[%s] %s is now %s (old Satima-PO config key)',
                 section, key, new)
    w = cfg.general.range_warning
    if not (_is_num(w) and w >= 0):
        fail('[general] range_warning 須為 0 以上的數值 (目前為 %r)', '[general] range_warning must be a number >= 0 (got %r)', w)
    sw = cfg.swath
    if not (isinstance(sw.merge, bool) or str(sw.merge).strip().lower() == 'ask'):
        fail('[swath] merge 須為 "ask"、true 或 false (目前為 %r)', '[swath] merge must be "ask", true or false (got %r)',
             sw.merge)
    if not (_is_num(sw.merge_gap) and sw.merge_gap >= 0):
        fail('[swath] merge_gap 須為 0 以上的秒數 (目前為 %r)', '[swath] merge_gap must be seconds >= 0 (got %r)',
             sw.merge_gap)
    e = cfg.eyetemp
    if not (_is_num(e.radius_deg) and e.radius_deg > 0):
        fail('[eyetemp] radius_deg 須為正數 (目前為 %r)', '[eyetemp] radius_deg must be a positive number (got %r)',
             e.radius_deg)
    check_colormaps(cfg, ('ir', 'wv'))
