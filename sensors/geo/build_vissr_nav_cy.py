"""編譯導航的 Cython 模組 ([navigation] engine = "cython" 用)：
  vissr_nav_cy.pyx -> GMS VISSR (GMS-VISSR、GMS-AREA)
  goes_nav_cy.pyx  -> GOES Imager (GOES-GVAR、GOES-AREA)
  mcgoes_nav_cy.pyx -> McIDAS 'GOES' 導航 (GOES-PREGVAR：SMS、GOES-1~7)

用法 (在任何目錄執行皆可)：
    python sensors/geo/build_vissr_nav_cy.py
需要 Cython 與 C 編譯器 (Windows：Visual Studio Build Tools 的 C++ 工具)。
編譯出的 .pyd 只適用於編譯時的 Python 版本，換版本後重新執行一次即可。
(檔名沒有 '-'，不會被當成讀取器載入)
"""
import importlib
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
MODULES = ('vissr_nav_cy', 'goes_nav_cy', 'mcgoes_nav_cy')


def main():
    from Cython.Build import cythonize
    from setuptools import Extension, setup

    win = sys.platform == 'win32'
    exts = [Extension(name, [name + '.pyx'],
                      extra_compile_args=['/openmp', '/O2'] if win else ['-fopenmp', '-O3'],
                      extra_link_args=[] if win else ['-fopenmp']) for name in MODULES]
    tmp = tempfile.mkdtemp(prefix='satima_build_')
    cwd = os.getcwd()
    os.chdir(HERE)
    try:
        setup(name='satima_nav_cy', ext_modules=cythonize(exts, language_level=3, quiet=True),
              script_args=['build_ext', '--inplace', '--build-temp', tmp, '--build-lib', tmp])
    finally:
        os.chdir(cwd)
        shutil.rmtree(tmp, ignore_errors=True)
        for name in MODULES:
            c_file = os.path.join(HERE, name + '.c')
            if os.path.exists(c_file):
                os.remove(c_file)
    sys.path.insert(0, HERE)
    for name in MODULES:
        print('完成：%s' % importlib.import_module(name).__file__)


if __name__ == '__main__':
    main()
