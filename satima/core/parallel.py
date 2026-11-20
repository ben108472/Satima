"""多核心處理：以執行緒平行執行互相獨立的工作

numpy 的陣列運算、bz2 / zlib 解壓、scipy 濾波在計算時會釋放 GIL，所以用執行緒就能用滿多核心，
不必像多行程那樣複製大陣列。輸出影像分成多條 (resample.OutputGrid.strips)，各條同時計算再寫回各自的位置，
結果與逐條依序處理完全相同。

[memory] workers："auto" = CPU 邏輯核心數 (最多 AUTO_MAX)；1 = 不平行；或指定數字。
每個執行緒同時處理一條，記憶體用量約隨執行緒數增加 (每條約數百 MB)。
"""
import os
import threading
from concurrent.futures import ThreadPoolExecutor

AUTO_MAX = 16
LOCK = threading.RLock()     # 讀取器檔案存取等非執行緒安全的操作共用


def workers(cfg):
    w = cfg.memory.workers if 'workers' in cfg.memory else 'auto'
    if isinstance(w, str):
        return max(1, min(os.cpu_count() or 1, AUTO_MAX))
    return max(1, int(w))


def run(jobs, n):
    """執行 jobs (無參數函式的序列)，回傳結果串列 (順序同 jobs)；n <= 1 時依序執行。任一工作失敗即拋出例外。"""
    jobs = list(jobs)
    if n <= 1 or len(jobs) <= 1:
        return [job() for job in jobs]
    with ThreadPoolExecutor(min(n, len(jobs))) as ex:
        futures = [ex.submit(job) for job in jobs]
        return [f.result() for f in futures]
