"""輸出格點 (同 Satima-GEO 的 resample.OutputGrid，去掉靜止衛星導航的部分)

輸出地圖 (satima/maps.py) 的像素中心格點；第一列為北邊。影像分條處理：每一條各自算經緯度、
取樣、上色後即釋放，多核心時各條由多個執行緒同時處理，結果與整張處理相同。
"""
import numpy as np

CHUNK = 1 << 21             # 分段計算的像素數
STRIP = 1 << 22             # 節省記憶體模式：每一條的輸出像素數 (約 400 萬)
MIN_STRIP = 1 << 18         # 平行處理時每一條的最少像素數


def strip_pixels(size, low_memory, workers):
    """每一條的像素數：單核心時節省記憶體為 STRIP，否則整張一條；
    多核心時切成約每個執行緒 4 條 (不超過 STRIP)，讓各執行緒的工作量平均。"""
    if workers <= 1:
        return STRIP if low_memory else None
    return min(STRIP, max(MIN_STRIP, -(-size // (4 * workers))))


class OutputGrid:
    """輸出地圖的像素中心格點；第一列為北邊 (imshow 請用 origin='upper')。"""

    def __init__(self, gmap, shape, pixel_km):
        ny, nx = int(shape[0]), int(shape[1])
        self.map = gmap
        self.shape = (ny, nx)
        self.size = ny * nx
        self.pixel_km = float(pixel_km)
        x0, x1, y0, y1 = gmap.extent
        dx, dy = (x1 - x0) / nx, (y1 - y0) / ny
        self.x = x0 + (np.arange(nx) + 0.5) * dx
        self.y = y1 - (np.arange(ny) + 0.5) * dy
        if gmap.separable:
            # 麥卡托、等距圓柱：經度只與 x、緯度只與 y 有關，反算一列與一行即可
            self._lon1, _ = gmap.inverse(self.x, np.full(nx, 0.5 * (y0 + y1)))
            _, self._lat1 = gmap.inverse(np.full(ny, 0.5 * (x0 + x1)), self.y)
        self._lonlat = None

    def strip_ranges(self, max_pixels=None):
        """各條的列範圍 [(i0, i1), ...]。"""
        if not max_pixels or self.size <= max_pixels:
            return [(0, self.shape[0])]
        step = max(1, int(max_pixels) // self.shape[1])
        return [(i0, min(self.shape[0], i0 + step)) for i0 in range(0, self.shape[0], step)]

    def lonlat_rows(self, i0, i1):
        """第 i0 ~ i1-1 列的經緯度 (看不到的位置為 NaN)。"""
        if self._lonlat is not None:
            return self._lonlat[0][i0:i1], self._lonlat[1][i0:i1]
        if self.map.separable:
            return np.meshgrid(self._lon1, self._lat1[i0:i1])
        x, y = np.meshgrid(self.x, self.y[i0:i1])
        return self.map.inverse(x, y)

    def lonlat_lattice(self, rows, cols):
        """指定的列 rows、行 cols 交叉處 (格子點) 的經緯度，形狀 (rows.size, cols.size)。"""
        if self.map.separable:
            return np.meshgrid(self._lon1[cols], self._lat1[rows])
        x, y = np.meshgrid(self.x[cols], self.y[rows])
        return self.map.inverse(x, y)

    def lonlat_points(self, rows, cols):
        """逐點 (rows[k], cols[k]) 的經緯度。"""
        if self.map.separable:
            return self._lon1[cols], self._lat1[rows]
        return self.map.inverse(self.x[cols], self.y[rows])

    def lonlat_sparse(self, step=8):
        """每隔 step 個像素取樣的經緯度 (估計地圖涵蓋範圍用)，含最後一列與最後一行。"""
        ny, nx = self.shape
        rows = np.unique(np.append(np.arange(0, ny, step), ny - 1))
        cols = np.unique(np.append(np.arange(0, nx, step), nx - 1))
        return self.lonlat_lattice(rows, cols)
