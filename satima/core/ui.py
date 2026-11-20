"""互動式輸入 (GEO / MW / PO 共用)

輸入以 ast.literal_eval 解析 (不執行任意程式碼)；以空白或逗號分隔的多個數字也接受 (20 27 117 124)。
輸入錯誤時會重新詢問。
"""
import ast
import re

from .i18n import tr

_SPLIT = re.compile(r'[\s,]+')


def _read(prompt):
    """input() 並去掉前後空白與 BOM (管線或貼上文字時可能出現)。"""
    try:
        return input(prompt).replace('﻿', '').strip()
    except EOFError:
        raise SystemExit(tr('\n輸入結束。', '\nInput ended.'))


def parse(text):
    """把輸入文字轉成 Python 值；以空白或逗號分隔的多個值轉成 tuple。"""
    text = text.strip()
    try:
        return ast.literal_eval(text)
    except (ValueError, SyntaxError):
        pass
    items = [t for t in _SPLIT.split(text) if t]
    if len(items) > 1:
        try:
            return tuple(ast.literal_eval(t) for t in items)
        except (ValueError, SyntaxError):
            pass
    raise ValueError(text)


def ask(prompt, check=None, default=None):
    """詢問直到輸入合法；check(value) 回傳錯誤訊息字串 (合法時回傳 None)。"""
    while True:
        text = _read(prompt)
        if not text:
            if default is not None:
                return default
            continue
        try:
            value = parse(text)
        except ValueError:
            print(tr('  無法辨識的輸入：%s', '  Invalid input: %s') % text)
            continue
        if check is not None:
            err = check(value)
            if err:
                print('  ' + err)
                continue
        return value


def _is_int(v):
    return isinstance(v, int) and not isinstance(v, bool)


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _fmt_choices(choices):
    c = sorted(choices)
    if len(c) > 6 and c == list(range(c[0], c[-1] + 1)):
        return '%d-%d' % (c[0], c[-1])
    return ', '.join(str(i) for i in c)


def ask_int(prompt, choices, default=None):
    choices = list(choices)

    def check(v):
        if not _is_int(v) or v not in choices:
            return tr('請輸入：%s', 'Please input one of: %s') % _fmt_choices(choices)
        return None
    return ask(prompt, check, default)


def ask_choices(prompt, choices):
    """可複選 (以空白或逗號分隔，例如 13 17)；回傳不重複、依輸入順序的列表。"""
    choices = list(choices)

    def check(v):
        items = list(v) if isinstance(v, (tuple, list)) else [v]
        bad = [i for i in items if not _is_int(i) or i not in choices]
        if not items or bad:
            return tr('請輸入：%s (可複選，以空白分隔)',
                      'Please input: %s (several allowed, separated by spaces)') % _fmt_choices(choices)
        return None
    v = ask(prompt, check)
    items = list(v) if isinstance(v, (tuple, list)) else [v]
    return list(dict.fromkeys(items))


def ask_menu(header, options, prompt, allow_all=True):
    """印出單行選單 (例：'0.ALL 1.BW 2.BD ...') 並詢問，可複選；回傳所選的 key 列表 (0 = 全部)。
    只有一個選項時直接回傳，不詢問。"""
    keys = list(options)
    if len(keys) == 1:
        return keys
    items = (['0.ALL'] if allow_all else []) + ['%d.%s' % (k, options[k]) for k in keys]
    print(((header + ' ') if header else '') + ' '.join(items))
    picked = ask_choices(prompt, ([0] if allow_all else []) + keys)
    return keys if 0 in picked else picked


def ask_index(items, prompt=None):
    """選擇：3 = 單一；(3, 5) = 第 3 到第 5 個；[1, 3, 5] = 指定多個。回傳所選項目列表。"""
    n = len(items)
    if prompt is None:
        prompt = tr('輸入資料編號 (多筆可用 (3,5) 或 [1,3,5])：', 'Input the index (several: (3,5) or [1,3,5]): ')

    def check(v):
        idx = list(v) if isinstance(v, (tuple, list)) else [v]
        if not idx or not all(_is_int(k) and 1 <= k <= n for k in idx):
            return tr('編號須為 1~%d 的整數 (多筆可用 (3,5) 或 [1,3,5])',
                      'The index must be an integer between 1 and %d (several: (3,5) or [1,3,5])') % n
        return None

    v = ask(prompt, check)
    if isinstance(v, tuple) and len(v) == 2:
        return items[min(v) - 1:max(v)]
    if isinstance(v, (tuple, list)):
        return [items[k - 1] for k in v]
    return [items[v - 1]]


def ask_georange(mode='bounds'):
    """地圖範圍 (度)。回傳 ((latmin, latmax, lonmin, lonmax), center)。
    mode = 'bounds'：問「南 北 西 東」；'center'：先問中心 (緯度 經度)，再問範圍 (中心往外的緯度、經度半寬)。
    不論哪一種，輸入 4 個數字一律當成 南 北 西 東、2 個數字一律當成中心 (接著問範圍)。
    center = (緯度, 經度)：輸入中心時為該點，輸入四邊時為範圍中點 (極軌的中心時間、颱風眼溫度用)。"""
    def check(v):
        if not (isinstance(v, tuple) and len(v) in (2, 4) and all(_is_num(x) for x in v)):
            return tr('請輸入 4 個數字 (南 北 西 東，例如 20 27 117 124) 或 2 個數字 (中心緯度 經度，例如 23.5 121)',
                      'Please input 4 numbers (south north west east, e.g. 20 27 117 124) '
                      'or 2 numbers (center lat lon, e.g. 23.5 121)')
        if len(v) == 2:
            return None if -89.9 < v[0] < 89.9 else tr('中心緯度須在 -89.9 ~ 89.9 之間',
                                                      'The center latitude must be between -89.9 and 89.9')
        latmin, latmax, lonmin, lonmax = v
        if not (-89.9 < latmin < latmax < 89.9):
            return tr('緯度須滿足 -89.9 < 南 < 北 < 89.9', 'Latitudes must satisfy -89.9 < south < north < 89.9')
        if not lonmin < lonmax:
            return tr('經度須滿足 西 < 東 (跨換日線請用 0~360，例如 170 200)',
                      'Longitudes must satisfy west < east (across the date line use 0~360, e.g. 170 200)')
        return None

    prompt = (tr('地圖中心 (緯度 經度)：', 'Map center (lat lon): ') if mode == 'center'
              else tr('地圖範圍 (南 北 西 東)：', 'Map range (south north west east): '))
    v = ask(prompt, check)
    if len(v) == 4:
        r = tuple(float(x) for x in v)
        return r, (0.5 * (r[0] + r[1]), 0.5 * (r[2] + r[3]))
    clat, clon = (float(x) for x in v)

    def span(r):
        r = r if isinstance(r, tuple) else (r,)
        if not (len(r) in (1, 2) and all(_is_num(x) and x > 0 for x in r)):
            return tr('請輸入 1 或 2 個正數：緯度半寬 經度半寬 (只輸入一個 = 兩者相同)，例如 5 或 5 6',
                      'Please input 1 or 2 positive numbers: lat half-width, lon half-width '
                      '(one number = both), e.g. 5 or 5 6')
        if abs(clat) + r[0] >= 89.9:
            return tr('緯度範圍超出 ±89.9 度', 'The latitude range exceeds ±89.9 degrees')
        return None

    r = ask(tr('範圍 (中心往外的緯度、經度各幾度；只輸入一個 = 兩者相同)：',
               'Range (degrees of latitude / longitude from the center; one number = both): '), span)
    r = r if isinstance(r, tuple) else (r,)
    dlat, dlon = float(r[0]), float(r[-1])
    return (clat - dlat, clat + dlat, clon - dlon, clon + dlon), (clat, clon)


def ask_yes_no(prompt, default=False):
    """是 / 否；直接 Enter = default。"""
    while True:
        text = _read(prompt).lower()
        if not text:
            return default
        if text in ('y', 'yes', '1', '是'):
            return True
        if text in ('n', 'no', '0', '否'):
            return False
        print(tr('  請輸入 y 或 n', '  Please input y or n'))


def ask_projection(options, default):
    """地圖投影選單 (options = {代號: 顯示名稱})；直接 Enter = default。回傳代號。"""
    keys = list(options)
    print('Projection: ' + '  '.join('%d.%s' % (i, options[k]) for i, k in enumerate(keys, 1)))
    k = ask_int(tr('輸入投影編號 (Enter = %s)：', 'Input the projection (Enter = %s): ') % options[default],
                range(1, len(keys) + 1), default=keys.index(default) + 1)
    return keys[k - 1]


_SCALE = re.compile(r'^[xX×*]\s*(\d+)$|^(\d+)\s*[xX×*]$')


def ask_width(max_size, native_min, scale=1, max_scale=8, floor=200, native_zh='原始解析度', native_en='native resolution'):
    """圖片寬度 (原始解析度模式)。回傳 (寬度, 倍數)：
    直接按 Enter = 原始解析度 x scale ([image] size = "xN"；"native" 為 1)，回傳 (None, scale)，寬度至少 native_min；
    輸入 x3 (或 3x) = 這次用原始解析度 x 3，回傳 (None, 3)；
    輸入數字 = 指定寬度，回傳 (寬度, 1)，只受 floor~max_size 限制 (想要小圖時可以指定)。
    native_zh / native_en：「原始解析度」的說法 (微波為「原始取樣間距」)。"""
    times = ' x%d' % scale if scale > 1 else ''
    native = tr(native_zh, native_en) + times
    prompt = tr('圖片寬度 (%d~%d；Enter = %s，寬度至少 %d；x2~x%d = %s的倍數)：',
                'Image width (%d~%d; Enter = %s, at least %d wide; x2~x%d = multiples of the %s): ') % (
        floor, max_size, native, native_min, max_scale, tr(native_zh, native_en))
    while True:
        text = _read(prompt)
        if not text:
            return None, scale
        m = _SCALE.match(text)
        if m:
            k = int(m.group(1) or m.group(2))
            if 1 <= k <= max_scale:
                return None, k
            print(tr('  倍數須為 1~%d', '  The multiple must be 1~%d') % max_scale)
            continue
        try:
            v = parse(text)
        except ValueError:
            v = None
        if _is_int(v) and floor <= v <= max_size:
            return v, 1
        print(tr('  請輸入 %d~%d 的整數、x2~x%d，或直接按 Enter 使用%s',
                 '  Please input an integer %d~%d, x2~x%d, or just press Enter for the %s') % (
            floor, max_size, max_scale, native))


def ask_title(cfg_title):
    """圖片左下角的標題；設定 [title] ask = false 時使用 default。"""
    default = str(cfg_title.get('default', '') or '')
    if not cfg_title.get('ask', False):
        return default
    hint = ' [%s]' % default if default else ''
    return _read(tr('圖片標題%s：', 'Image title%s: ') % hint) or default


def print_columns(sections, indent=4, gap=3):
    """分組的多欄選單：sections = [(標題, [項目文字, ...])]；欄數依終端機寬度。"""
    import shutil
    texts = [t for _, items in sections for t in items]
    if not texts:
        return
    cell = max(len(t) for t in texts) + gap
    ncol = max(1, (shutil.get_terminal_size((120, 24)).columns - indent) // cell)
    for title, items in sections:
        print('  ' + title)
        for k in range(0, len(items), ncol):
            print(' ' * indent + ''.join(t.ljust(cell) for t in items[k:k + ncol]).rstrip())


def format_range(georange):
    """(南, 北, 西, 東) -> '13.4N-23.4N 123.9E-134.6E' (經度 0~360 也轉回東經/西經)。"""
    def lat(v):
        return '%.1f%s' % (abs(v), 'N' if v >= 0 else 'S')

    def lon(v):
        v = (v + 180.) % 360. - 180.
        return '%.1f%s' % (abs(v), 'E' if v >= 0 else 'W')
    latmin, latmax, lonmin, lonmax = georange
    return '%s-%s %s-%s' % (lat(latmin), lat(latmax), lon(lonmin), lon(lonmax))


def _width(text):
    """終端機上的顯示寬度 (中日韓全形字算 2 格)。"""
    import unicodedata
    return sum(2 if unicodedata.east_asian_width(c) in ('W', 'F') else 1 for c in text)


def _pad(text, width):
    return text + ' ' * max(0, width - _width(text))


def print_listing(rows, start=1, index_width=None):
    """對齊欄位印出列表：序號、各欄。start = 第一列的序號。"""
    rows = [[str(c) for c in r] for r in rows]
    if not rows:
        return
    ncol = max(len(r) for r in rows)
    widths = [max(_width(r[i]) for r in rows if i < len(r)) for i in range(ncol)]
    iw = index_width or len(str(start + len(rows) - 1))
    for i, r in enumerate(rows, start):
        cells = [_pad(c, widths[j]) for j, c in enumerate(r)]
        print('%*d   %s' % (iw, i, '   '.join(cells).rstrip()))


def print_groups(groups):
    """分段的列表：groups = [(標題, [列, ...])]；序號連續，各段欄寬各自對齊。只有一段時不印標題。"""
    groups = [(title, rows) for title, rows in groups if rows]
    total = sum(len(rows) for _, rows in groups)
    iw = len(str(total))
    start = 1
    for title, rows in groups:
        if len(groups) > 1:
            print('[%s]' % title)
        print_listing(rows, start, iw + (2 if len(groups) > 1 else 0))
        start += len(rows)
