"""靜止衛星 (GEO) 的流程：搜尋場景、資料列表、範圍、產品、繪圖存檔"""
import dataclasses
import os
import time

from . import config, products                      # noqa: F401  (config：主流程讀設定檔用)
from ..core import maps, naming, session, ui
from ..core.i18n import tr
from .canvas import Canvas


def collect(readers, cfg, name=None, data_dir=None, output_dir=None):
    """各讀取器在自己的資料目錄搜尋場景。回傳 [(讀取器物件, 場景)]，依時間排序。"""
    found = session.collect(readers, cfg, name, data_dir, output_dir)
    found.sort(key=lambda rs: (rs[1].time, rs[0].NAME, rs[1].area))
    return found


def rows(found):
    """資料列表：衛星、儀器、格式 (同一衛星、時間可能同時有多種格式，例如 GOES 的 GVAR / NSMC GVAR / AREA)、
    時間、觀測範圍、波段數、備註 (區域觀測的經緯度範圍)。"""
    return [(s.satellite, r.sensor_label(s), r.format_label(s), s.time.strftime('%Y/%m/%d %H:%MZ'), s.area_label,
             '%d bands' % len(r.channels(s)), _where(r, s)) for r, s in found]


def _where(reader, scene):
    """區域觀測 (目標區、mesoscale 等) 在資料列表顯示影像的經緯度範圍，才知道觀測的地點。"""
    if not scene.auto_range:
        return ''
    try:
        return ui.format_range(reader.auto_georange(scene))
    except (Exception, SystemExit):
        return '(range unknown)'


def _view(reader, scenes, cfg):
    """詢問投影與地圖範圍，回傳所有場景共用的 maps.View；全圓盤 (fd) 回傳 None。
    全圓盤不需要經緯度範圍：圖框依各場景、各產品的波段決定 (_fd_view)。
    其他投影：區域觀測直接用影像範圍，否則詢問。"""
    m = cfg.map
    proj = str(m.projection).strip().lower()
    if m.ask_projection:
        proj = ui.ask_projection({k: maps.PROJECTIONS[k] for k in config.PROJECTIONS}, proj)
    if proj == 'fd':
        return None
    if all(s.auto_range for s in scenes) and len({s.area for s in scenes}) == 1:
        georange = reader.auto_georange(scenes[0])
        print(tr('地圖範圍 (%s)：南 %.2f 北 %.2f 西 %.2f 東 %.2f',
                 'Map range (%s): south %.2f north %.2f west %.2f east %.2f') % ((scenes[0].area_label,) + georange))
    else:
        georange, _ = ui.ask_georange(str(m.range_input).strip().lower())
    return session.make_view(proj, georange, cfg)


def _fd_view(reader, scene, prod, cfg):
    """全圓盤：依產品的波段 (真實色用紅光) 決定圖框。"""
    ch = reader.CHANNELS[prod.recipe['red']] if prod.kind == 'truecolor' else prod.channel
    view = maps.fd_view(reader, scene, ch, str(cfg.map.backend).lower())
    lon0, h = view.geos[:2]
    print(tr('  全圓盤：星下點 %.2f%s，衛星高度 %.0f km', '  Full disk: sub-satellite point %.2f%s, altitude %.0f km')
          % (abs(lon0), 'E' if lon0 >= 0 else 'W', h / 1000.))
    return view


def run(picked, cfg_all):
    reader = picked[0][0]
    scenes = [s for _, s in picked]
    cfg = reader.cfg

    shared = _view(reader, scenes, cfg)
    prods = products.choose(reader, scenes[0], cfg)
    fixed, scale = session.ask_size(cfg)
    images = session.images(cfg)
    title = ui.ask_title(cfg.title)
    print()

    canvases = {}

    def canvas_for(view):
        if view not in canvases:
            canvases[view] = Canvas(view, cfg)
        return canvases[view]

    saved = []
    for scene in scenes:
        own = {ch.name: ch for ch in reader.channels(scene)}
        by_band = {ch.display: ch for ch in own.values()}
        tcs = reader.truecolor_products(scene)
        when = scene.time.strftime('%Y/%m/%d %H:%MZ')
        for prod in prods:
            # 產品依第一筆資料的選單選擇；其他資料依顯示名稱 (BAND-14、TrueColor) 對應到自己的波段 / 真實色配置
            # (例如同時選 MTG FDHSI 與 HRFI：BAND-14 在 HRFI 場景改用 ir_105_hr)
            if prod.kind == 'truecolor':
                if prod.name not in tcs:
                    print(tr('[略過] %s %s %s 無法製作 %s', '[skip] %s %s %s cannot make %s')
                          % (scene.satellite, when, scene.area, prod.name))
                    continue
                prod = dataclasses.replace(prod, recipe=tcs[prod.name])
            elif prod.name not in own:
                ch = by_band.get(prod.display)
                if ch is None:
                    print(tr('[略過] %s %s %s 沒有 %s', '[skip] %s %s %s has no %s')
                          % (scene.satellite, when, scene.area, prod.display))
                    continue
                prod = dataclasses.replace(prod, name=ch.name, channel=ch)
            band = ' '.join(t for t in (prod.display, scene.tag) if t)          # 例如 BAND-14 HRFI
            label = band + ('-' + prod.colorbar if prod.colorbar else '')
            print(tr('正在製作 %s %s %s ...', 'Making %s %s %s ...') % (scene.satellite, when, label))
            view = shared or _fd_view(reader, scene, prod, cfg)
            canvas = canvas_for(view)
            fields = dict(title=title, satellite=scene.satellite, sensor=reader.sensor_label(scene), band=band,
                          colorbar=prod.colorbar, product=label, wavelength=prod.wavelength,
                          time=scene.time, area=scene.area, reader=reader.NAME, projection=view.tag)
            a = cfg.annotation
            info = naming.render(str(a.text), fields, a.time_format, space=None, clean_values=False)
            fields['wavelength'] = prod.wavelength.replace('µ', 'u')        # 檔名用 ASCII：10.4µm -> 10.4um
            native = None
            for path, img, with_label in naming.output_paths(cfg.output, fields, cfg.general.output_dir, images,
                                                             a.draw):
                t0 = time.time()
                if native is None and not img.get('dpi'):
                    native = prod.native_res(reader, scene)
                width, inch = session.size_for(canvas, img, native, fixed, scale)
                print(tr('  圖片大小 %dx%d，約 %.2f km/像素', '  Image size %dx%d, about %.2f km/pixel')
                      % (width, round(width * canvas.aspect), canvas.width_km / width))
                try:
                    image = products.render(prod, reader, scene, canvas, width, cfg)
                except ValueError as err:
                    print(tr('  [略過] %s', '  [skip] %s') % err)
                    break
                fig = canvas.draw(image, width, title, info if with_label else '',
                                  lines=prod.kind != 'map' or bool(cfg[prod.basemap].lines), inch=inch)
                del image
                os.makedirs(os.path.dirname(path) or '.', exist_ok=True)
                canvas.save(fig, path)
                print(tr('  已存檔 %s (%.1f 秒)', '  Saved %s (%.1f s)') % (path, time.time() - t0))
                saved.append(path)
    return saved
