"""極軌成像儀 (PO) 的流程：搜尋檔案並接成段、資料列表、範圍、頻道與色表、繪圖存檔"""
import contextlib
import dataclasses
import os
import time

from . import config                                # noqa: F401  (主流程讀設定檔用)
from ..core import naming, segments as seg, session, ui
from ..core.i18n import tr
from .canvas import Canvas


@contextlib.contextmanager
def _cwd(folder):
    """讀取器以相對於資料目錄的檔名找檔 (glob、對應的定位檔)：在資料目錄中執行，結束後還原。"""
    old = os.getcwd()
    os.chdir(folder)
    try:
        yield
    finally:
        os.chdir(old)


def collect(readers, cfg, name=None, data_dir=None, output_dir=None):
    """各讀取器在自己的資料目錄搜尋檔案並接成段。回傳 [(讀取器物件, 段)]，依時間排序。"""
    found = []
    for cls in readers:
        if name and cls.NAME.lower() != name.lower():
            continue
        rcfg = config.effective(cfg, cls.NAME, data_dir, output_dir)
        folder = rcfg.general.data_dir
        if not os.path.isdir(folder):
            continue
        with _cwd(folder):
            names = cls.find()
            if not names:
                continue
            reader = cls(rcfg)
            files = []
            for n in names:
                try:
                    files.append(reader.make_file(n))
                except Exception as err:
                    print(tr('[略過] %s：%s', '[skip] %s: %s') % (n, err))
            if not files:
                continue
            reader.all_files = files
            try:
                segs = reader.segments(files)
            except Exception as err:
                print(tr('[略過] %s：%s: %s', '[skip] %s: %s: %s') % (cls.NAME, type(err).__name__, err))
                continue
        found.extend((reader, s) for s in segs)
    found.sort(key=lambda rs: (rs[1].time, rs[0].NAME, rs[1].satellite, rs[1].product))
    return found


def rows(found):
    """資料列表：衛星、儀器、格式、時間、資料來源、頻道 / 產品、備註 (多檔接成一筆時的檔數與缺檔數)。"""
    return [(s.satellite, r.sensor_label(s), r.NAME.split('-', 1)[-1], _span(s), r.source_label(s), s.product,
             seg.note(s.slots, s.missing)) for r, s in found]


def _span(s):
    """列表的時間欄：2026/10/20 06:54-08:36Z (跨日時結束時間前加日期)；不知道結束時間時只列起始時間。"""
    a, b = s.time, s.end
    if b is None or b <= a:
        return a.strftime('%Y/%m/%d %H:%MZ')
    end = b.strftime('%H:%M') if a.date() == b.date() else b.strftime('%m/%d %H:%M')
    return '%s-%sZ' % (a.strftime('%Y/%m/%d %H:%M'), end)


def _merge(segs, reader, cfg):
    """一次選到多段時：時間相接的接成一段。短片段接成的段自動合併；整圈 / 半圈的長檔依 [swath] merge 詢問。"""
    if len(segs) < 2 or not reader.MULTI_FILE:
        return segs
    gap = float(cfg.swath.merge_gap)
    chains = []
    key = reader.merge_key
    for s in sorted(segs, key=lambda s: (key(s), s.time)):
        last = chains[-1][-1] if chains else None
        if (last is not None and key(last) == key(s)
                and (s.time - max(x.end for x in chains[-1])).total_seconds() <= gap):
            chains[-1].append(s)
        else:
            chains.append([s])
    if all(len(c) == 1 for c in chains):
        return segs
    mode = cfg.swath.merge
    long_chains = [c for c in chains if len(c) > 1 and not all(s.short for s in c)]
    merge_long = mode is True
    if long_chains and mode not in (True, False):
        for c in long_chains:
            print(tr('時間相接：%s %s  %s ~ %s (%d 個檔)', 'Contiguous: %s %s  %s ~ %s (%d files)') % (
                c[0].satellite, reader.sensor_label(c[0]), c[0].time.strftime('%Y/%m/%d %H:%MZ'),
                max(s.end for s in c).strftime('%H:%MZ'), sum(s.slots for s in c)))
        merge_long = ui.ask_yes_no(tr('時間相接的檔案合併成一張？(y/N)：',
                                      'Merge the contiguous files into one image? (y/N): '), False)
    out = []
    for c in chains:
        if len(c) == 1 or not (all(s.short for s in c) or merge_long):
            out.extend(c)
            continue
        files = [f for s in c for f in s.files]
        out.append(dataclasses.replace(c[0], end=max(s.end for s in c), files=files,
                                       slots=sum(s.slots for s in c), missing=sum(s.missing for s in c)))
    return sorted(out, key=lambda s: s.time)


def run(picked, cfg_all):
    reader = picked[0][0]
    cfg = reader.cfg
    segs = _merge([s for _, s in picked], reader, cfg)

    view, center = session.ask_view(cfg, config.PROJECTIONS)
    latmin, latmax, lonmin, lonmax = view.georange
    warn = float(cfg.general.range_warning)
    if warn and (latmax - latmin > 2 * warn or lonmax - lonmin > 2 * warn):
        print(tr('注意：範圍很大 (超過中心往外 %g 度)，高解析度的頻道可能用掉很多記憶體。',
                 'NOTE: large range (over %g degrees from the center); high-resolution channels may need a lot '
                 'of memory.') % warn)
    canvas = Canvas(view, cfg)
    reader.georange, reader.center = canvas.frame_range(), center

    os.chdir(cfg.general.data_dir)          # 讀取器以相對於資料目錄的檔名讀檔；輸出目錄已是絕對路徑
    reader.begin(segs[0])
    reader.select_band()
    variants = reader.select_variants() or [None]
    fixed, scale = session.ask_size(cfg, '原始取樣間距', 'native sampling interval')
    images = session.images(cfg)
    title = ui.ask_title(cfg.title)
    print()

    saved = []
    for n, segment in enumerate(segs):
        if n:
            reader.begin(segment)
        if len(segs) > 1:
            print(tr('==== %s %s %s ====', '==== %s %s %s ====') % (segment.satellite, reader.sensor_label(segment), _span(segment)))
        try:
            ready = reader.prepare()
        except (OSError, ValueError) as err:        # 缺檔 (定位檔、真實色的某個頻道)、檔案內容不符
            print(tr('[略過] %s', '[skip] %s') % err)
            continue
        if not ready:
            print(tr('地圖範圍內沒有資料 (緯度 %.2f~%.2f，經度 %.2f~%.2f)，略過。',
                     'No data inside the map range (lat %.2f~%.2f, lon %.2f~%.2f); skipped.') % reader.georange)
            continue
        for variant in variants:
            canvas.clear()
            try:
                info = reader.plot(canvas, variant)
            except (OSError, ValueError) as err:
                print(tr('[略過] %s', '[skip] %s') % err)
                continue
            if not info:
                continue
            band, colorbar = str(info.get('band', '') or ''), str(info.get('colorbar', '') or '')
            label = band + ('-' + colorbar if colorbar else '')
            if canvas.empty:
                print(tr('%s 沒有可畫的資料，略過。', 'Nothing to plot for %s; skipped.') % (label or tr('這個產品', 'this product')))
                continue
            fields = dict(title=title, satellite=reader.satellite, sensor=reader.sensor_label(segment), band=band,
                          colorbar=colorbar, product=label, time=reader.time, res=reader.res, reader=reader.NAME,
                          projection=view.tag)
            a = cfg.annotation
            text = naming.render(str(a.text), fields, a.time_format, space=None, clean_values=False)
            native = None
            for path, img, with_label in naming.output_paths(cfg.output, fields, cfg.general.output_dir, images,
                                                             a.draw):
                t0 = time.time()
                if native is None and not img.get('dpi'):
                    native = canvas.native_res()
                width, inch = session.size_for(canvas, img, native, fixed, scale)
                print(tr('  圖片大小 %dx%d，約 %.2f km/像素', '  Image size %dx%d, about %.2f km/pixel')
                      % (width, round(width * canvas.aspect), canvas.width_km / width))
                image = canvas.render(width)
                if image is None:
                    break
                fig = canvas.draw(image, width, title, text if with_label else '', inch=inch)
                del image
                os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
                canvas.save(fig, path)
                print(tr('  已存檔 %s (%.1f 秒)', '  Saved %s (%.1f s)') % (path, time.time() - t0))
                saved.append(path)
    return saved
