"""回歸測試：以相同的輸入分別執行原版 (Satima-GEO / Satima-MW / Satima-PO) 與整合版 (Satima)，逐像素比對輸出圖。

用法
  python tests/regress.py run orig [個案 ...]     執行原版，輸出到 <OUT>/orig/<個案>/
  python tests/regress.py run new  [個案 ...]     執行整合版，輸出到 <OUT>/new/<個案>/
  python tests/regress.py cmp orig new [個案 ...] 比對兩邊的圖 (檔名相同者逐像素比較)
  python tests/regress.py list                   列出個案

樣本資料與原版程式的位置寫在下方 (ROOT)；輸出目錄 OUT 預設為 D:/github/_baseline。
個案的 stdin 就是互動時依序輸入的內容；整合版多出的提問用 new_stdin 另外指定。
"""
import os
import re
import subprocess
import sys
import time

ROOT = r'D:\github'
OUT = os.path.join(ROOT, '_baseline')
GEO_DATA = os.path.join(ROOT, 'Satima-GEO', 'data')
L1B = os.path.join(ROOT, 'L1B-sample')
L1C = os.path.join(ROOT, 'L1C-sample')
PO_DATA = os.path.join(ROOT, 'Satima-PO', 'data')
VIIRS = r'D:\VIIRS'                         # I5 / I1 / DNB
VIIRS_TC = r'D:\VIIRS\Claude'               # M3 / M4 / M5 / I1 (真實色)

TARGETS = {
    'orig': {'geo': ('Satima-GEO', 'Satima-GEO.py', []),
             'mw': ('Satima-MW', 'Satima-MW.py', []),
             'po': ('Satima-PO', 'Satima-PO.py', [])},
    'new': {'geo': ('Satima', 'Satima.py', ['-t', 'geo']),
            'mw': ('Satima', 'Satima.py', ['-t', 'mw']),
            'po': ('Satima', 'Satima.py', ['-t', 'po'])},
}

# stdin：資料編號、投影、範圍、產品、色表、寬度 … (依各程式的提問順序)
CASES = {
    # ---- GEO ----
    'geo_ahi_merc': dict(family='geo', args=['-f', 'AHI-HSD', '-d', GEO_DATA],
                         stdin=['1', '1', '12 28 126 142', '6 7', '2', '']),
    'geo_ahi_fd': dict(family='geo', args=['-f', 'AHI-HSD', '-d', GEO_DATA],
                       stdin=['1', '4', '6', '1', '1 2', '1000']),
    'geo_mtsat_lcc': dict(family='geo', args=['-f', 'MTSAT-HRIT', '-d', GEO_DATA],
                          stdin=['1', '3', '5 35 115 150', '1 2', '1', '2', '1500']),
    # ---- MW ----
    'mw_amsr2': dict(family='mw', args=['-f', 'AMSR2-L1B', '-d', L1B],
                     stdin=['1', '1', '14 26 116 128', '13 18 19', '1', '1', '']),
    'mw_atms': dict(family='mw', args=['-f', 'ATMS-SDR', '-d', L1B],
                    stdin=['1', '1', '14 42 -82 -48', '16', '1', '']),
    'mw_ssmis_l1c': dict(family='mw', args=['-f', 'PPS-L1C', '-d', L1C],
                         stdin=['1', '1', '12 28 126 142', '1', '1', '']),
    # ---- PO ----
    'po_avhrr': dict(family='po', args=['-f', 'AVHRR-NAT', '-d', PO_DATA],
                     stdin=['1', '40 -5', '6 6', '5', '2', 'TEST'],
                     new_stdin=['1', '1', '40 -5', '6 6', '5', '2', '', 'TEST']),
    # VIIRS (真實檔案在 D:\VIIRS；原版一個檔一列、整合版一段一列，所以編號與提問順序不同)
    'po_viirs_i5': dict(family='po', args=['-f', 'VIIRS-HDF5', '-d', VIIRS],         # NOAA-20 2026/10/27 I5 色表 BD
                        stdin=['(21,23)', '15, 90', '5, 5', '0', '2', 'TEST'],
                        new_stdin=['3', '1', '15 90', '5 5', '2', '2', '', 'TEST']),
    'po_viirs_dnb': dict(family='po', args=['-f', 'VIIRS-HDF5', '-d', VIIRS],        # NOAA-21 2026/10/20 DNB ERF
                         stdin=['(1,7)', '16, 142', '4, 4', '0', '2', 'TEST'],
                         new_stdin=['1', '1', '16 142', '4 4', '2', '', 'TEST']),
    'po_viirs_tc': dict(family='po', args=['-f', 'VIIRS-HDF5', '-d', VIIRS_TC],      # NOAA-21 2026/10/28 TrueColor-I
                        stdin=['(18,21)', '12, 159', '4, 4', '2', '2', 'TEST'],
                        new_stdin=['1', '1', '12 159', '4 4', '6', '2', '', 'TEST']),
}


def run(target, names):
    for name in names:
        case = CASES[name]
        folder, script, extra = TARGETS[target][case['family']]
        out = os.path.join(OUT, target, name)
        os.makedirs(out, exist_ok=True)
        for f in os.listdir(out):
            if f.lower().endswith('.png'):
                os.remove(os.path.join(out, f))
        stdin = case.get(target + '_stdin', case['stdin'])
        cmd = [sys.executable, script] + extra + case.get(target + '_args', case['args']) + ['-o', out]
        env = dict(os.environ, PYTHONIOENCODING='utf-8')
        t0 = time.time()
        p = subprocess.run(cmd, cwd=os.path.join(ROOT, folder), input='\n'.join(stdin) + '\n', env=env,
                           capture_output=True, text=True, encoding='utf-8', errors='replace')
        with open(os.path.join(out, 'log.txt'), 'w', encoding='utf-8') as f:
            f.write('$ %s\n' % ' '.join(cmd))
            f.write(p.stdout)
            f.write('\n---- stderr ----\n' + p.stderr)
        pngs = sorted(f for f in os.listdir(out) if f.lower().endswith('.png'))
        print('%-16s %-5s exit=%d  %5.1fs  %d png' % (name, target, p.returncode, time.time() - t0, len(pngs)))
        if p.returncode or not pngs:
            print('    ' + '\n    '.join((p.stdout + p.stderr).strip().splitlines()[-8:]))


def _band_name(f):
    """檔名中的 PO 頻道名稱 .B04. / .B03A. 換成 .BAND-04. / .BAND-03A.。"""
    return re.sub(r'\.B(\d\d[AB]?)\.', r'.BAND-\1.', f)


def compare(a, b, names):
    import numpy as np
    from PIL import Image
    Image.MAX_IMAGE_PIXELS = None
    bad = 0
    for name in names:
        da, db = os.path.join(OUT, a, name), os.path.join(OUT, b, name)
        fa = sorted(f for f in os.listdir(da) if f.lower().endswith('.png')) if os.path.isdir(da) else []
        fb = sorted(f for f in os.listdir(db) if f.lower().endswith('.png')) if os.path.isdir(db) else []
        # 2026/11 起 PO 的頻道名稱由 B04 改為 BAND-04 (同 GEO)：比對時視為同名
        na, nb = {_band_name(f): f for f in fa}, {_band_name(f): f for f in fb}
        if not fa or sorted(na) != sorted(nb):
            print('%-16s 檔名不同：%s  vs  %s' % (name, fa, fb))
            bad += 1
            continue
        for key in sorted(na):
            f = na[key]
            xa = np.asarray(Image.open(os.path.join(da, f)).convert('RGBA')).astype(np.int16)
            xb = np.asarray(Image.open(os.path.join(db, nb[key])).convert('RGBA')).astype(np.int16)
            if xa.shape != xb.shape:
                print('%-16s %s 尺寸不同 %s vs %s' % (name, f, xa.shape[:2], xb.shape[:2]))
                bad += 1
                continue
            d = np.abs(xa - xb).max(axis=-1)
            n = int((d > 0).sum())
            if n:
                bad += 1
            print('%-16s %-70s %s' % (name, f, '相同' if not n else '不同：%d 像素 (%.4f%%)，最大差 %d' % (
                n, 100. * n / d.size, int(d.max()))))
    return bad


def main():
    for stream in (sys.stdout, sys.stderr):         # 主控台是 cp950 時 µ 等字元會出錯
        stream.reconfigure(encoding='utf-8', errors='replace')
    args = sys.argv[1:]
    if not args or args[0] == 'list':
        for name, case in CASES.items():
            print('%-16s %s  %s' % (name, case['family'], ' | '.join(case['stdin'])))
        return
    if args[0] == 'run':
        run(args[1], args[2:] or list(CASES))
    elif args[0] == 'cmp':
        sys.exit(1 if compare(args[1], args[2], args[3:] or list(CASES)) else 0)
    else:
        raise SystemExit(__doc__)


if __name__ == '__main__':
    main()
