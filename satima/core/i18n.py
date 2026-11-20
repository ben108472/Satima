"""介面語言 (繁體中文 / English)

每則訊息的兩種語言寫在同一處：tr('找不到設定檔：%s', 'Config file not found: %s') % path
語言由 [general] language 或命令列 --lang 決定，一次只顯示一種。
產品與波段名稱、色表名稱、檔名、圖上的標註不翻譯，兩種語言的輸出檔完全相同。
"""
LANGUAGES = ('zh', 'en')
_ALIASES = {'zh': 'zh', 'zh-tw': 'zh', 'zh_tw': 'zh', 'tw': 'zh', 'cht': 'zh', 'chinese': 'zh', '中文': 'zh',
            '繁體中文': 'zh', 'en': 'en', 'en-us': 'en', 'en_us': 'en', 'eng': 'en', 'english': 'en'}
_lang = 'zh'


def normalize(code):
    """'EN'、'english'、'zh-TW' … -> 'en' / 'zh'；不認得時回傳 None。"""
    return _ALIASES.get(str(code).strip().lower())


def set_language(code):
    global _lang
    lang = normalize(code)
    if lang is None:
        raise ValueError(code)
    _lang = lang
    return lang


def language():
    return _lang


def tr(zh, en):
    """依目前語言回傳其中一個字串。"""
    return en if _lang == 'en' else zh
