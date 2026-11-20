"""檢查 config-en/ 的英文註解設定檔是否與專案根目錄的設定檔同步 (鍵與值必須相同，只有註解的語言不同)。
Check that the English-commented templates in config-en/ hold the same keys and values as the root config files.

  python tests/check_config_sync.py

每個檔印一行：OK，或列出差異 (缺少 / 多出的鍵與值不同的鍵，以 . 連接的路徑表示；註解掉的範例行不同也會列出)。
四個檔都相同時結束碼為 0，否則為 1。只用標準函式庫。
"""
import os
import sys
import tomllib
from collections import Counter

NAMES = ("Satima.toml", "Satima-GEO.toml", "Satima-MW.toml", "Satima-PO.toml")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def diff(a, b, path=""):
    """逐層比較解析後的資料，回傳差異說明的清單 (a = 根目錄，b = config-en)。"""
    if isinstance(a, dict) and isinstance(b, dict):
        out = []
        for key in list(a) + [k for k in b if k not in a]:
            sub = f"{path}.{key}" if path else str(key)
            if key not in b:
                out.append(f"config-en 缺少 (missing) {sub}")
            elif key not in a:
                out.append(f"config-en 多出 (extra) {sub}")
            else:
                out += diff(a[key], b[key], sub)
        return out
    if isinstance(a, list) and isinstance(b, list) and len(a) == len(b):
        out = []
        for i, (x, y) in enumerate(zip(a, b)):
            out += diff(x, y, f"{path}[{i}]")
        return out
    if type(a) is not type(b) or a != b:        # 型別也要相同 (1、1.0、true 在 Python 中相等)
        return [f"值不同 (differs) {path}: {a!r} != {b!r}"]
    return []


def examples(text):
    """註解掉的 TOML 範例行 (# key = value、# [節]、# { ... },)。

    去掉 # 之後能被 tomllib 解析的才算，所以一般的說明文字不會被當成範例；行尾的註解由解析器略過，翻譯不影響比較。
    """
    out = []
    for line in text.splitlines():
        s = line.strip()
        if not s.startswith("#"):
            continue
        s = s.lstrip("#").strip()
        if not s:
            continue
        if s.startswith("{"):                   # 陣列中被註解掉的一個元素
            s = "_ = [" + s + "]"
        try:
            data = tomllib.loads(s)
        except tomllib.TOMLDecodeError:
            continue
        if data:
            out.append(repr(data))
    return out


def check(name):
    texts = []
    for folder in (ROOT, os.path.join(ROOT, "config-en")):
        with open(os.path.join(folder, name), encoding="utf-8-sig") as f:
            texts.append(f.read())
    out = diff(tomllib.loads(texts[0]), tomllib.loads(texts[1]))
    a, b = (Counter(examples(t)) for t in texts)
    out += [f"註解範例只在根目錄 (commented example only in root): {x}" for x in a - b]
    out += [f"註解範例只在 config-en (commented example only in config-en): {x}" for x in b - a]
    return out


def main():
    sys.stdout.reconfigure(errors="backslashreplace")   # 主控台編碼 (cp950) 印不出的字元不要讓程式中斷
    failed = 0
    for name in NAMES:
        try:
            out = check(name)
        except (OSError, tomllib.TOMLDecodeError) as e:
            out = [f"無法讀取或解析 (cannot read or parse): {e}"]
        failed += bool(out)
        print(f"{name}: " + ("; ".join(out) if out else "OK"))
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
