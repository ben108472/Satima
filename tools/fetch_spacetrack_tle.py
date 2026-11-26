"""從 Space-Track 下載 NOAA 極軌衛星的歷史 TLE (AVHRR-L0 定位用)，存到 tle/spacetrack/<衛星>-<年>.txt。

1980~2004 年 tle/ 已有 CelesTrak 的歷史檔；這裡預設只下載 2005 年以後 (CelesTrak 不提供的部分)。

帳號 (Space-Track 免費註冊：https://www.space-track.org/auth/createAccount)：
  方法一：環境變數 SPACETRACK_IDENTITY (帳號 / email) 與 SPACETRACK_PASSWORD
  方法二：使用者家目錄的 .spacetrack 檔 (例如 C:\\Users\\<名稱>\\.spacetrack)，內容兩行：
            identity = 帳號
            password = 密碼
          (或只寫兩行：第一行帳號、第二行密碼)
  密碼不會印出，也不會寫進 Satima 的資料夾。

用法
  python tools/fetch_spacetrack_tle.py                     下載全部 (NOAA-12/14/15/16/17/18/19，2005 年起到停止運作)
  python tools/fetch_spacetrack_tle.py NOAA-18 NOAA-19     只下載指定衛星
  python tools/fetch_spacetrack_tle.py --from 2013 --to 2013
  python tools/fetch_spacetrack_tle.py --force             已下載的年份也重新下載
已存在的年份檔會略過 (中斷後重跑即可續傳)；今年的檔一律重新下載。
Space-Track 的流量限制是每分鐘 30 次、每小時 300 次查詢，這裡每次查詢間隔 3 秒 (全部約 85 次，約 5 分鐘)。
"""
import argparse
import datetime
import http.cookiejar
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = 'https://www.space-track.org'
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, 'tle', 'spacetrack')
PAUSE = 3.0
# 衛星 -> (NORAD 編號, 開始年, 停止運作的年)；只下載 2005 年以後
SATELLITES = {
    'NOAA-12': (21263, 1991, 2007), 'NOAA-14': (23455, 1994, 2007), 'NOAA-15': (25338, 1998, 2025),
    'NOAA-16': (26536, 2000, 2014), 'NOAA-17': (27453, 2002, 2013), 'NOAA-18': (28654, 2005, 2025),
    'NOAA-19': (33591, 2009, 2025),
}
FIRST_YEAR = 2005


def credentials():
    ident = os.environ.get('SPACETRACK_IDENTITY')
    pwd = os.environ.get('SPACETRACK_PASSWORD')
    if ident and pwd:
        return ident, pwd
    path = os.path.join(os.path.expanduser('~'), '.spacetrack')
    if os.path.isfile(path):
        kv = {}
        with open(path, encoding='utf-8-sig') as f:
            lines = [ln.strip() for ln in f if ln.strip() and not ln.lstrip().startswith('#')]
        for line in lines:
            if '=' in line:
                k, v = line.split('=', 1)
                kv[k.strip().lower()] = v.strip()
        if kv.get('identity') and kv.get('password'):
            return kv['identity'], kv['password']
        if not kv and len(lines) == 2:             # 只有兩行：第一行帳號、第二行密碼
            return lines[0], lines[1]
    sys.exit('找不到 Space-Track 帳號：請設定環境變數 SPACETRACK_IDENTITY / SPACETRACK_PASSWORD，或建立 %s' % path)


class Session:
    def __init__(self):
        self.opener = urllib.request.build_opener(urllib.request.HTTPCookieProcessor(http.cookiejar.CookieJar()))
        self.opener.addheaders = [('User-Agent', 'Satima-TLE-fetch/1.0')]

    def login(self, ident, pwd):
        data = urllib.parse.urlencode({'identity': ident, 'password': pwd}).encode()
        with self.opener.open(BASE + '/ajaxauth/login', data, timeout=60) as r:
            text = r.read().decode('utf-8', 'replace')
        if 'Failed' in text:
            sys.exit('Space-Track 登入失敗 (帳號或密碼錯誤，或帳號尚未核准)')

    def query(self, path):
        url = BASE + '/basicspacedata/query/' + urllib.parse.quote(path)
        with self.opener.open(url, timeout=300) as r:
            return r.read().decode('latin1')

    def logout(self):
        try:
            self.opener.open(BASE + '/ajaxauth/logout', timeout=30).read()
        except (urllib.error.URLError, OSError):
            pass


def main():
    ap = argparse.ArgumentParser(description='Download historical NOAA TLEs from Space-Track')
    ap.add_argument('satellites', nargs='*', help='e.g. NOAA-18 (default: all)')
    ap.add_argument('--from', dest='y0', type=int, default=FIRST_YEAR)
    ap.add_argument('--to', dest='y1', type=int, default=datetime.date.today().year)
    ap.add_argument('--force', action='store_true')
    args = ap.parse_args()
    names = [s.upper() for s in args.satellites] or list(SATELLITES)
    bad = [s for s in names if s not in SATELLITES]
    if bad:
        sys.exit('不認得的衛星：%s (可用：%s)' % (', '.join(bad), ', '.join(SATELLITES)))
    this_year = datetime.date.today().year
    jobs = []
    for name in names:
        norad, start, end = SATELLITES[name]
        for year in range(max(args.y0, start, FIRST_YEAR), min(args.y1, end) + 1):
            path = os.path.join(OUT, '%s-%d.txt' % (name.lower(), year))
            if not args.force and year < this_year and os.path.isfile(path) and os.path.getsize(path) > 0:
                continue
            jobs.append((name, norad, year, path))
    if not jobs:
        print('全部都已下載 (%s)' % OUT)
        return
    os.makedirs(OUT, exist_ok=True)
    s = Session()
    s.login(*credentials())
    print('登入成功，共 %d 次查詢 (每次間隔 %.0f 秒)' % (len(jobs), PAUSE))
    try:
        for k, (name, norad, year, path) in enumerate(jobs, 1):
            q = 'class/gp_history/NORAD_CAT_ID/%d/EPOCH/%d-01-01--%d-01-01/orderby/EPOCH asc/format/tle' % (
                norad, year, year + 1)
            text = s.query(q)
            lines = [ln for ln in text.splitlines() if ln.strip()]
            if lines and not lines[0].startswith('1 '):
                sys.exit('Space-Track 回傳的不是 TLE (%s %d)：%s' % (name, year, text[:200]))
            with open(path, 'w', encoding='ascii', newline='\n') as f:
                f.write('\n'.join(lines) + ('\n' if lines else ''))
            print('[%d/%d] %s %d：%d 組' % (k, len(jobs), name, year, len(lines) // 2))
            time.sleep(PAUSE)
    finally:
        s.logout()
    print('完成：%s' % OUT)


if __name__ == '__main__':
    main()
