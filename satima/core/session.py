"""主流程：命令列、讀取器載入、分段的資料列表、選資料，再交給各類衛星 (GEO / MW / PO) 的繪圖流程

各類衛星在 satima/<類別>/family.py 提供：
  config                  該類別的設定模組 (DEFAULTS、PROJECTIONS、validate；見 core/config.py)
  collect(readers, cfg, name, data_dir, output_dir)   搜尋資料 -> [(讀取器物件, 場景)] (列表的順序)
  rows(found)             資料列表的各列 (最後一欄為備註)
  run(picked, cfg)        使用者選定資料之後的流程 (範圍、產品、繪圖、存檔)
"""
import importlib
import importlib.util
import os
import re
import sys
import traceback
from glob import glob

from . import config, i18n, maps, naming, ui
from .i18n import tr

FAMILIES = ('geo', 'mw', 'po')
TITLES = {'geo': 'GEO', 'mw': 'MW', 'po': 'PO'}


# ---- 讀取器載入 ------------------------------------------------------------
def discover(folder):
    """載入 folder 內所有檔名含 '-' 的 .py (例如 AHI-HSD.py)，取其 READER 類別。
    檔名不含 '-' 的 .py 視為輔助模組，不會被當成讀取器載入。載入失敗 (缺套件等) 只略過該格式。"""
    if folder not in sys.path:
        sys.path.insert(0, folder)
    readers = []
    for path in sorted(glob(os.path.join(folder, '*-*.py'))):
        stem = os.path.splitext(os.path.basename(path))[0]
        modname = 'satima_sensor_' + re.sub(r'\W', '_', stem)
        try:
            spec = importlib.util.spec_from_file_location(modname, path)
            module = importlib.util.module_from_spec(spec)
            sys.modules[modname] = module
            spec.loader.exec_module(module)
            cls = module.READER
        except Exception as err:
            print(tr('[略過] %s：%s: %s', '[skip] %s: %s: %s') % (stem, type(err).__name__, err))
            continue
        if not cls.NAME:
            cls.NAME = stem
        readers.append(cls)
    return readers


def collect(readers, cfg, name=None, data_dir=None, output_dir=None):
    """各讀取器在自己的資料目錄搜尋場景 (讀取器的 scenes())。回傳 [(讀取器物件, 場景)] (未排序)。"""
    found = []
    for cls in readers:
        if name and cls.NAME.lower() != name.lower():
            continue
        rcfg = config.effective(cfg, cls.NAME, data_dir, output_dir)
        folder = rcfg.general.data_dir
        if not os.path.isdir(folder) or not cls.find(folder):
            continue
        reader = cls(rcfg)
        try:
            scenes = reader.scenes(folder)
        except Exception as err:
            print(tr('[略過] %s：%s: %s', '[skip] %s: %s: %s') % (cls.NAME, type(err).__name__, err))
            continue
        found.extend((reader, s) for s in scenes)
    return found


def describe(cls):
    """讀取器的說明 (--list)：DESCRIPTION 為中文，DESCRIPTION_EN 為英文 (沒有時兩種語言都用 DESCRIPTION)。"""
    return tr(cls.DESCRIPTION, getattr(cls, 'DESCRIPTION_EN', '') or cls.DESCRIPTION)


# ---- 各類衛星共用的提問 ------------------------------------------------------
def ask_view(cfg, projections):
    """詢問投影與地圖範圍，回傳 (maps.View, 中心 (緯度, 經度))。projections：這一類可用的投影代號。"""
    m = cfg.map
    proj = str(m.projection).strip().lower()
    if m.ask_projection and len(projections) > 1:
        proj = ui.ask_projection({k: maps.PROJECTIONS[k] for k in projections}, proj)
    georange, center = ui.ask_georange(str(m.range_input).strip().lower())
    return make_view(proj, georange, cfg), center


def make_view(proj, georange, cfg):
    lcc = maps.lcc_params(georange, cfg.map) if proj == 'lcc' else None
    if lcc:
        print(tr('LCC：標準緯線 %.2f、%.2f，中央經線 %.2f',
                 'LCC: standard parallels %.2f, %.2f; central longitude %.2f') % (lcc[0], lcc[1], lcc[3]))
    return maps.View(proj, tuple(georange), lcc=lcc)


def ask_size(cfg, native_zh='原始解析度', native_en='native resolution'):
    """原始解析度模式的圖片大小：回傳 (固定寬度或 None, 倍數)。[image] size 為 "native" / "xN" 時執行時詢問
    (Enter = 設定檔的倍數)。dpi 模式不詢問 (大小由 [[output.image]] 的 dpi 決定)，回傳 (None, 1)。"""
    if cfg.image.mode == 'dpi':
        return None, 1
    fixed, scale = config.image_size(cfg.image.size)
    if fixed is None:
        fixed, scale = ui.ask_width(int(cfg.image.max_size), int(cfg.image.min_size), scale, config.MAX_SCALE,
                                    native_zh=native_zh, native_en=native_en)
    return fixed, scale


def images(cfg):
    """這次要輸出的圖：dpi 模式為 [[output.image]] 的各項；原始解析度模式只有一張 (dpi = None)。"""
    if cfg.image.mode == 'dpi':
        return [dict(item) for item in cfg.output.image]
    return [{'dpi': None, 'suffix': str(cfg.image.get('suffix', '') or '')}]


def size_for(canvas, img, native_res, fixed, scale):
    """某一張輸出的 (寬度像素, 圖寬英吋或 None)。"""
    if img.get('dpi'):
        return canvas.dpi_size(img['dpi'])
    return canvas.width_for(native_res, fixed, scale), None


# ---- 主流程 ----------------------------------------------------------------
def _family(kind):
    return importlib.import_module('satima.%s.family' % kind)


def main(argv=None):
    import argparse
    ap = argparse.ArgumentParser(prog='Satima', description='Satima: weather satellite imagery (GEO / MW / PO)')
    ap.add_argument('-c', '--config', help='common config file (default: Satima.toml)')
    ap.add_argument('-d', '--data-dir', help='data directory (overrides [general] data_dir)')
    ap.add_argument('-o', '--output-dir', help='output directory (overrides [general] output_dir)')
    ap.add_argument('-t', '--type', dest='kind', choices=FAMILIES, type=str.lower,
                    help='only this kind of satellite: geo / mw / po')
    ap.add_argument('-f', '--format', dest='reader', help='only this format, e.g. AHI-HSD')
    ap.add_argument('--lang', help='language: zh / en (overrides [general] language)')
    ap.add_argument('--mode', choices=config.MODES, type=str.lower,
                    help='image mode: native / dpi (overrides [image] mode)')
    ap.add_argument('--list', action='store_true', help='list available formats and exit')
    args = ap.parse_args(argv)

    lang = args.lang or config.peek_language(args.config) or 'zh'
    try:
        i18n.set_language(lang)
    except ValueError:
        raise SystemExit('--lang: zh / en (got %r)' % lang)
    overrides = {'image': {'mode': args.mode}} if args.mode else {}
    data_dir = os.path.abspath(args.data_dir) if args.data_dir else None
    output_dir = os.path.abspath(args.output_dir) if args.output_dir else None

    try:
        groups = []                     # (類別, family 模組, 設定, 讀取器類別)
        for kind in ([args.kind] if args.kind else FAMILIES):
            try:
                fam = _family(kind)
            except ImportError as err:      # 缺套件 (例如 PO 需要 scipy)：略過這一類，其他類照常使用
                print(tr('[略過] %s：%s', '[skip] %s: %s') % (TITLES[kind], err))
                continue
            cfg = config.load(kind, fam.config, args.config, overrides)
            if not args.lang:           # 類別檔也可以設定語言
                i18n.set_language(cfg.general.language)
            groups.append((kind, fam, cfg, discover(os.path.join(config.HOME, 'sensors', kind))))
        if args.list:
            ui.print_groups([(TITLES[kind], [(c.NAME, c.SENSOR, describe(c), '  '.join(c.PATTERNS))
                                             for c in readers]) for kind, _, _, readers in groups])
            return
        paths = list(dict.fromkeys(p for _, _, cfg, _ in groups for p in cfg['_paths']))
        print(tr('Satima  設定檔：%s', 'Satima  config: %s') % (', '.join(paths) or tr('(預設值)', '(defaults)')))
        for _, _, cfg, _ in groups:
            bad = naming.unknown_fields(cfg.output.filename) + naming.unknown_fields(cfg.general.output_dir)
            if bad:
                print(tr('注意：[output] filename / output_dir 有不認得的欄位 %s (會留空)',
                         'WARNING: unknown field(s) in [output] filename / output_dir: %s (they will be empty)')
                      % ', '.join(dict.fromkeys(bad)))
        names = [c.NAME for _, _, _, readers in groups for c in readers]
        if args.reader and args.reader.lower() not in (n.lower() for n in names):
            raise SystemExit(tr('沒有 %s 這個格式。可用：%s', 'Unknown format %s. Available: %s')
                             % (args.reader, ', '.join(names)))

        found = []                      # (類別, family, 設定, [(讀取器, 場景)])
        for kind, fam, cfg, readers in groups:
            items = fam.collect(readers, cfg, args.reader, data_dir, output_dir)
            if items:
                found.append((kind, fam, cfg, items))
        if not found:
            lines = [tr('找不到可讀取的衛星資料。各格式的搜尋樣式：',
                        'No supported satellite files found. Search patterns:')]
            for kind, _, cfg, readers in groups:
                for cls in readers:
                    if args.reader and cls.NAME.lower() != args.reader.lower():
                        continue
                    folder = config.effective(cfg, cls.NAME, data_dir, output_dir).general.data_dir
                    lines.append(tr('  %-4s %-12s %-50s 於 %s', '  %-4s %-12s %-50s in %s')
                                 % (TITLES[kind], cls.NAME, '  '.join(cls.PATTERNS), folder))
            raise SystemExit('\n'.join(lines))

        ui.print_groups([(TITLES[kind], fam.rows(items)) for kind, fam, _, items in found])
        flat = [(kind, fam, cfg, item) for kind, fam, cfg, items in found for item in items]
        while True:
            picked = ui.ask_index(flat)
            readers = list(dict.fromkeys(p[3][0].NAME for p in picked))
            if len({id(p[3][0]) for p in picked}) == 1:
                break
            print(tr('  一次只能選同一種格式的資料 (選到了 %s)，請重新選擇。',
                     '  Only one format can be selected at a time (got %s); please select again.')
                  % tr('、', ', ').join(readers))
        _, fam, cfg, _ = picked[0]
        fam.run([p[3] for p in picked], cfg)
    except KeyboardInterrupt:
        print(tr('\n已中斷。', '\nInterrupted.'))
    except SystemExit:
        raise
    except Exception:
        traceback.print_exc()
        raise SystemExit(1)
