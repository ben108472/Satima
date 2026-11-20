"""輸出檔名

樣板寫在設定檔 [output] filename，可用欄位 (某一類衛星沒有的欄位留空)：
  {title}      執行時輸入的標題
  {satellite}  衛星，例如 Himawari-9、GCOM-W1、NOAA-20
  {sensor}     儀器，例如 AHI、AMSR2、VIIRS
  {band}       波段 / 頻道 / 產品，例如 BAND-13、89V、I5、TrueColor
  {colorbar}   色表或樣式名稱，例如 BD、NRL (沒有時為空)
  {product}    {band} 加上 -{colorbar}，例如 BAND-13-OTT
  {time}       觀測時間 (極軌為掃描帶通過地圖範圍的時間；格式見 [output] time_format，也可寫 {time:%Y%m%d})
  {projection} 地圖投影：CYL、LCC、FD (麥卡托為空)
  {reader}     讀取器名稱，例如 AHI-HSD
  {dpi}        dpi 模式該張圖的 dpi (原始解析度模式為空)
  {suffix}     dpi 模式為 [[output.image]] 的 suffix；原始解析度模式為 [image] suffix
  {wavelength} GEO：中心波長，例如 10.4µm (檔名中為 10.4um)
  {area}       GEO：觀測範圍，例如 FLDK、R301、CONUS
  {frequency}  MW：頻率，例如 89GHz
  {orbit}      MW：軌道 / granule 編號
  {res}        PO：解析度 (MODIS、SGLI 有提供)
空白的欄位會連同它前面的分隔字元一起省略，例如沒有標題、沒有色表、麥卡托時
  {title}.{satellite}.{sensor}.{band}.{colorbar}.{projection}.{time}  ->  Himawari-9.AHI.BAND-03.202611080400Z
"""
import datetime
import os
import re
import string

from .i18n import tr

FIELDS = ('title', 'satellite', 'sensor', 'band', 'colorbar', 'product', 'time', 'projection', 'reader', 'dpi',
          'suffix', 'wavelength', 'area', 'frequency', 'orbit', 'res')
_ILLEGAL = re.compile(r'[<>:"/\\|?*\x00-\x1f]')


def clean(value, space='_'):
    """去掉檔名不允許的字元；空白改為 space。"""
    text = '' if value is None else str(value).strip()
    text = _ILLEGAL.sub('-', text)
    if space is not None:
        text = re.sub(r'\s+', space, text)
    return text


def fields_in(template):
    return [name for _, name, _, _ in string.Formatter().parse(template) if name]


def unknown_fields(template):
    return [name for name in fields_in(template) if name not in FIELDS]


def render(template, fields, time_format='%Y%m%d%H%MZ', space='_', clean_values=True):
    segments, trailing = [], ''
    for literal, name, spec, _ in string.Formatter().parse(template):
        if name is None:
            trailing += literal
            continue
        value = fields.get(name)
        if isinstance(value, datetime.datetime):
            value = value.strftime(spec or time_format)
        elif value not in (None, '') and spec:
            value = format(value, spec)
        segments.append((literal, clean(value, space) if clean_values else ('' if value is None else str(value))))
    prefix = segments[0][0] if segments else ''
    body = []
    for i, (literal, value) in enumerate(segments):
        if not value:
            continue
        body.append((literal if body and i > 0 else '') + value)
    return prefix + ''.join(body) + trailing


def output_path(out, fields, outdir):
    """依 [output] filename 樣板產生完整路徑；outdir 可含欄位 (例如 D:/png/{satellite}/{time:%Y%m%d})。"""
    space = out.get('space', '_')
    base = render(str(out.filename), fields, out.time_format, space) or 'satima'
    folder = render(str(outdir), fields, out.time_format, space) if '{' in str(outdir) else str(outdir)
    return os.path.join(folder, '%s.%s' % (base, str(out.format).lstrip('.')))


def output_paths(out, fields, outdir, images, annotation=True):
    """多張輸出 (dpi 模式的 [[output.image]]；原始解析度模式為單一張)：回傳 [(路徑, 該張的設定 dict, 是否加圖說)]。
    images：[{'dpi': 1000 或 None, 'suffix': '', 'annotation': 可省略}, ...]
    樣板中沒有 {suffix} / {dpi} 時，有 suffix 的圖自動在檔名後加上 '.suffix'；
    仍然同名 (多張都沒有 suffix) 時第 2 張起加上 dpi 或序號，避免互相覆蓋。"""
    template = str(out.filename)
    names = set(fields_in(template))
    ext = str(out.format).lstrip('.')
    space = out.get('space', '_')
    paths, seen = [], set()
    for i, img in enumerate(images):
        dpi = int(img['dpi']) if img.get('dpi') else None
        suffix = str(img.get('suffix', '') or '')
        f = dict(fields, dpi=dpi, suffix=suffix)
        base = render(template, f, out.time_format, space) or 'satima'
        if suffix and 'suffix' not in names:
            base += '.' + clean(suffix, space)
        folder = render(str(outdir), f, out.time_format, space) if '{' in str(outdir) else str(outdir)
        path = os.path.join(folder, '%s.%s' % (base, ext))
        if path in seen:
            base += '.' + ('%dDPI' % dpi if dpi else str(i + 1))
            path = os.path.join(folder, '%s.%s' % (base, ext))
        if path in seen:
            print(tr('注意：輸出檔名重複 %s (請檢查 [output] filename 與 suffix)',
                     'WARNING: duplicated output name %s (check [output] filename / suffix)') % path)
        seen.add(path)
        paths.append((path, img, bool(img.get('annotation', annotation))))
    return paths
