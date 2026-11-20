"""微波 (MW) 的流程：搜尋場景、資料列表、多檔合併、範圍、產品、繪圖存檔"""
import os
import time

from . import config, products                      # noqa: F401  (config：主流程讀設定檔用)
from ..core import naming, session, ui
from ..core.i18n import tr
from .canvas import Canvas
from .reader import from_seconds, merge_scenes


def collect(readers, cfg, name=None, data_dir=None, output_dir=None):
    """各讀取器在自己的資料目錄搜尋場景。回傳 [(讀取器物件, 場景)]，依時間排序。"""
    found = session.collect(readers, cfg, name, data_dir, output_dir)
    found.sort(key=lambda rs: (rs[1].time, rs[0].NAME, rs[1].satellite, rs[1].orbit))
    return found


def rows(found):
    """資料列表：衛星、儀器、格式、時間、資料來源、頻道數、備註 (多檔接成一筆時的檔數)。"""
    return [(s.satellite, r.sensor_label(s), r.format_label(s), _span(s), r.source_label(s),
             '%d ch' % len(r.channels(s)), r.note(s)) for r, s in found]


def _span(scene):
    """列表的時間欄：2026/10/20 06:54-08:36Z (跨日時結束時間前加日期)。"""
    a, b = scene.time, scene.end
    end = b.strftime('%H:%M') if a.date() == b.date() else b.strftime('%m/%d %H:%M')
    return '%s-%sZ' % (a.strftime('%Y/%m/%d %H:%M'), end)


def _merge(scenes, cfg):
    """同衛星、同儀器、時間相接的場景合併成一個 (依 [swath] merge)。"""
    mode = cfg.swath.merge
    if len(scenes) < 2 or mode is False:
        return scenes
    gap = float(cfg.swath.merge_gap)
    chains = []
    for s in sorted(scenes, key=lambda s: (s.satellite, s.sensor, s.format, s.time)):
        last = chains[-1][-1] if chains else None
        if (last is not None and (last.satellite, last.sensor, last.format) == (s.satellite, s.sensor, s.format)
                and (s.time - last.end).total_seconds() <= gap):
            chains[-1].append(s)
        else:
            chains.append([s])
    if all(len(c) == 1 for c in chains):
        return scenes
    if mode is not True:
        for c in chains:
            if len(c) > 1:
                print(tr('時間相接：%s %s  %s ~ %s (%d 個檔)', 'Contiguous: %s %s  %s ~ %s (%d files)') % (
                    c[0].satellite, c[0].sensor, c[0].time.strftime('%Y/%m/%d %H:%MZ'),
                    max(s.end for s in c).strftime('%H:%MZ'), len(c)))
        if not ui.ask_yes_no(tr('時間相接的檔案合併成一張？(y/N)：', 'Merge the contiguous files into one image? (y/N): '),
                             False):
            return scenes
    out = [c[0] if len(c) == 1 else merge_scenes(c) for c in chains]
    return sorted(out, key=lambda s: s.time)


def run(picked, cfg_all):
    reader = picked[0][0]
    cfg = reader.cfg
    scenes = _merge([s for _, s in picked], cfg)

    view, _ = session.ask_view(cfg, config.PROJECTIONS)
    prods = products.choose(reader, scenes[0], cfg)
    fixed, scale = session.ask_size(cfg, '原始取樣間距', 'native sampling interval')
    images = session.images(cfg)
    title = ui.ask_title(cfg.title)
    print()

    canvas = Canvas(view, cfg)
    saved = []
    for scene in scenes:
        start = scene.time.strftime('%Y/%m/%d %H:%MZ')
        for chosen in prods:
            # 產品依第一筆資料的選單選擇；其他資料依名稱 (89V、PCT89、NRL-37) 對應到自己的頻道
            prod = products.resolve(chosen, reader, scene, cfg)
            if prod is None:
                print(tr('[略過] %s %s %s 無法製作 %s', '[skip] %s %s %s cannot make %s')
                      % (scene.satellite, reader.sensor_label(scene), start, chosen.name))
                continue
            print(tr('正在製作 %s %s %s %s ...', 'Making %s %s %s %s ...')
                  % (scene.satellite, reader.sensor_label(scene), start, prod.label))
            native, when = None, None
            fields = info = paths = None
            for k, img in enumerate(images):
                t0 = time.time()
                try:
                    if native is None and not img.get('dpi'):
                        native = prod.native_res(reader, scene)
                    width, inch = session.size_for(canvas, img, native, fixed, scale)
                    if k == 0:
                        rs = products.resampler(canvas, width, cfg)
                        cover = 0.
                        for geo, _ in reader.swath(scene, prod.channels[0]):
                            c, t = rs.coverage(geo)
                            if c > cover:
                                cover, when = c, t
                        if cover < float(cfg.swath.min_coverage):
                            print(tr('  [略過] 掃描帶沒有通過這個地圖範圍',
                                     '  [skip] the swath does not pass over this map range'))
                            break
                        # 標註與檔名的時間：掃描帶通過地圖範圍的時間 (不是檔案的起始時間)
                        when = from_seconds(when) if when is not None else scene.time
                        print(tr('  通過時間 %s，涵蓋地圖 %.0f%%', '  Overpass time %s, covering %.0f%% of the map')
                              % (when.strftime('%Y/%m/%d %H:%M:%SZ'), 100. * cover))
                    print(tr('  圖片大小 %dx%d，約 %.2f km/像素', '  Image size %dx%d, about %.2f km/pixel')
                          % (width, round(width * canvas.aspect), canvas.width_km / width))
                    image = products.render(prod, reader, scene, canvas, width, cfg)
                except ValueError as err:
                    print(tr('  [略過] %s', '  [skip] %s') % err)
                    break
                if not image[..., 3].any():
                    # 掃描帶有通過，但這個頻道沒有任何有效的值 (故障的頻道：F18 SSMIS 的 150 GHz、AMSR-E 2004/11 後的 89A…)
                    print(tr('  [略過] 這個頻道在地圖範圍內沒有有效的資料',
                             '  [skip] this channel has no valid data inside the map range'))
                    break
                if paths is None:
                    fields = dict(title=title, satellite=scene.satellite, sensor=reader.sensor_label(scene),
                                  band=prod.display, colorbar=prod.colorbar, product=prod.label,
                                  frequency=prod.frequency, time=when, orbit=scene.orbit, reader=reader.NAME,
                                  projection=view.tag)
                    a = cfg.annotation
                    info = naming.render(str(a.text), fields, a.time_format, space=None, clean_values=False)
                    for key in ('band', 'product'):         # 檔名用 ASCII：183.31±7V -> 183.31+-7V
                        fields[key] = fields[key].replace('±', '+-')
                    paths = naming.output_paths(cfg.output, fields, cfg.general.output_dir, images, a.draw)
                path, _, with_label = paths[k]
                fig = canvas.draw(image, width, title, info if with_label else '', inch=inch)
                del image
                os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
                canvas.save(fig, path)
                print(tr('  已存檔 %s (%.1f 秒)', '  Saved %s (%.1f s)') % (path, time.time() - t0))
                saved.append(path)
    return saved
