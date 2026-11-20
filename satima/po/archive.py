"""Sentinel-3 SAFE 產品 (.SEN3) 的 NetCDF 讀取，不解壓縮

SAFE 產品為一個 .SEN3 資料夾 (xfdumanifest.xml + 各頻道 NetCDF)，下載時通常再包成 zip 或 tar。
原本每次執行都把整包解壓到目前目錄；這裡改為直接從壓縮檔把需要的 NetCDF 讀進記憶體開啟
(netCDF4.Dataset(memory=...)，METimage 腳本原本即用此方式)，不會在磁碟上留下解壓檔。
已解壓的 .SEN3 資料夾也可以直接讀。
"""
import os
import tarfile
import zipfile

import netCDF4 as nc

from ..core.i18n import tr


def prefer_archives(names):
    """同時有 X.SEN3.tar(.zip) 與 X.SEN3 資料夾 (例如舊程式解壓出來的) 時只列壓縮檔，避免重複。"""
    names = list(names)
    have = set(names)
    return [n for n in names if not (n + '.tar' in have or n + '.zip' in have)]


class SafeProduct:
    """SAFE 產品 (zip / tar / 資料夾)。

        safe = SafeProduct('S3A_OL_1_EFR____....SEN3.zip')
        ds = safe.dataset('Oa08_radiance')     # 檔名包含此字串的第一個 NetCDF
        ...
        ds.close()
    """

    def __init__(self, fullname):
        self.fullname = fullname
        self._zip = self._tar = None
        self._members = {}
        if os.path.isdir(fullname):
            self.kind = 'dir'
            names = []
            for root, _, files in os.walk(fullname):
                for fn in files:
                    names.append(os.path.relpath(os.path.join(root, fn), fullname).replace('\\', '/'))
        elif zipfile.is_zipfile(fullname):
            self.kind = 'zip'
            self._zip = zipfile.ZipFile(fullname)
            names = [i.filename for i in self._zip.infolist() if not i.is_dir()]
        elif tarfile.is_tarfile(fullname):
            self.kind = 'tar'
            self._tar = tarfile.open(fullname, 'r:*')
            self._members = {m.name: m for m in self._tar.getmembers() if m.isfile()}
            names = list(self._members)
        else:
            raise ValueError(tr('不支援的 SAFE 產品：%s', 'Unsupported SAFE product: %s') % fullname)
        self.names = sorted(names)

    def find(self, key):
        """名稱包含 key 的第一個檔案 (NetCDF 優先)。"""
        found = [n for n in self.names if key in os.path.basename(n)]
        if not found:
            raise FileNotFoundError(tr('%s 不在 %s 裡面', '%s not found in %s') % (key, self.fullname))
        nc_files = [n for n in found if n.lower().endswith('.nc')]
        return (nc_files or found)[0]

    def dataset(self, key):
        """開啟名稱包含 key 的 NetCDF (zip / tar 直接讀進記憶體)。"""
        name = self.find(key)
        if self.kind == 'dir':
            return nc.Dataset(os.path.join(self.fullname, name), mode='r')
        if self.kind == 'zip':
            raw = self._zip.read(name)
        else:
            raw = self._tar.extractfile(self._members[name]).read()
        return nc.Dataset(os.path.basename(name), mode='r', memory=raw)

    def close(self):
        if self._zip is not None:
            self._zip.close()
        if self._tar is not None:
            self._tar.close()
