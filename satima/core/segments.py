"""極軌資料的接續判定：把很短的片段 (granule) 接成「一段連續的觀測」

極軌的檔案長短差很多：ATMS 一檔 32 秒、METimage 約 1 分鐘、VIIRS 約 85 秒、MODIS 5 分鐘、JAXA L1B 半圈。
這裡把「相鄰」和「同一段」分成兩件事：
  相鄰    下一檔起始 - 前一檔結束 <= 相鄰門檻。門檻跟著片段長度走：該組檔案長度中位數的一半，限制在 2~120 秒
          (ATMS 16 秒、VIIRS 約 43 秒、5 分鐘以上的檔 120 秒)。缺一個檔的間隔等於一檔的長度，一定超過門檻，
          所以缺幾個檔都數得出來 (只用來在列表上標示，不影響能不能接)。
  同一段  間隔不超過 max_gap (預設 600 秒) 就接成同一段 (資料列表的一列)，中間缺檔照樣接：
          同一地點的兩次通過至少隔一圈 (約 100 分鐘)，10 分鐘內不可能是另一次通過；
          缺口處取樣時會留空 (掃描間距超過正常值即視為缺掃描)，不會被內插抹平。
軌道號不用來切段：軌道號在升交點過赤道時跳號，會剛好把低緯度颱風的通過切成兩半。
時間一律是 1970 起算的秒數或 datetime (同一次呼叫內一致)。
"""
import datetime
import statistics

from .i18n import tr

MAX_GAP = 600.              # 同一段：兩檔的間隔上限 (秒)
TOL_MIN, TOL_MAX = 2., 120.     # 相鄰門檻的下限 / 上限 (秒)
SHORT = 900.                # 一檔短於此 (秒) 視為短片段：自動接成一段，不詢問


def _seconds(dt):
    return dt.total_seconds() if isinstance(dt, datetime.timedelta) else float(dt)


def tolerance(lengths):
    """相鄰門檻 (秒)：片段長度中位數的一半，限制在 TOL_MIN~TOL_MAX。lengths：各檔的長度 (秒)。"""
    lengths = [v for v in lengths if v and v > 0]
    if not lengths:
        return TOL_MAX
    return min(TOL_MAX, max(TOL_MIN, 0.5 * statistics.median(lengths)))


def dedupe(items, start, newer):
    """同一片段重複的檔 (重新下載、不同產製時間) 只留最新的一份。
    start(item) = 起始時間；newer(item) = 可比較的新舊依據 (例如檔名中的產製時間)。回傳依起始時間排序的列表。"""
    best = {}
    for item in items:
        key = start(item)
        if key not in best or newer(item) > newer(best[key]):
            best[key] = item
    return [best[k] for k in sorted(best)]


def chain(items, start, end, max_gap=MAX_GAP):
    """把依時間排序的片段接成幾段。start(item)、end(item) 為起訖時間。
    回傳 [(該段的片段列表, 缺檔數)]；缺檔數 = 段內超過相鄰門檻的間隔各相當於幾個檔 (至少算 1 個)。"""
    items = sorted(items, key=start)
    if not items:
        return []
    lengths = [_seconds(end(g) - start(g)) for g in items]
    tol = tolerance(lengths)
    typical = statistics.median([v for v in lengths if v > 0] or [0.])
    out, members, missing = [], [items[0]], 0
    last = end(items[0])
    for g in items[1:]:
        gap = _seconds(start(g) - last)
        if gap > max_gap:
            out.append((members, missing))
            members, missing = [], 0
        elif gap > tol:
            missing += max(1, int(round(gap / typical))) if typical > 0 else 1
        members.append(g)
        if end(g) > last or not members[:-1]:
            last = end(g)
    out.append((members, missing))
    return out


def is_short(length):
    """一檔的長度 (秒) 是否算短片段 (自動接、不詢問)。"""
    return 0 < length < SHORT


def note(nfiles, missing=0):
    """資料列表的備註欄：多檔接成一筆時顯示檔數，有缺檔才加 (缺 N)；單一檔留空。"""
    if nfiles <= 1 and not missing:
        return ''
    text = tr('%d 檔', '%d files') % nfiles
    if missing:
        text += tr(' (缺 %d)', ' (%d missing)') % missing
    return text
