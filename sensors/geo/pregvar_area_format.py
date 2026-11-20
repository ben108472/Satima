"""McIDAS AREA (GVAR 以前：SMS、GOES-1~7 的 VISSR / VAS) 的結構 — 原樣複製自舊版 satellite/AREA_format.py
(D:/Satima/Satima/satellite，2024/06)，欄位名稱不變，讀取器以這些名稱保留完整檔頭。
結構為原生位元組順序；讀取器依 AREA 目錄第 2 個 word 判斷順序後以 dtype.newbyteorder() 轉換。
檔名沒有 '-'，不會被當成讀取器載入。最後的 _Calibration_AAA 為新增 (舊版直接讀陣列)。
"""

import numpy as np


_Control_Block = np.dtype([ ('W1', 'i4'),
                            ('Area Format', 'i4'),
                            ('SSS', 'i4'),
                            ('YYDDD', 'i4'),
                            ('HHMMSS', 'i4'),
                            ('UPPERLEFTLINE', 'i4'),
                            ('UPPERLEFTELE', 'i4'),
                            ('W8', 'i4'),
                            ('NLINES', 'i4'),
                            ('NELES', 'i4'),
                            ('ELESIZ', 'i4'),
                            ('LINERES', 'i4'),
                            ('ELERES', 'i4'),
                            ('NCHANS', 'i4'),
                            ('PRESIZ', 'i4'),
                            ('PROJ', 'i4'),
                            ('CREATION DATE', 'i4'),
                            ('CREATION TIME', 'i4'),
                            ('FILTER MAP', 'i4'),
                            ('W20', 'i4'),
                            ('W21', 'i4'),
                            ('W22', 'i4'),
                            ('W23', 'i4'),
                            ('W24', 'i4'), 
                            ('Memo', 'a32'),
                            ('Area file number', 'i4'),
                            ('W34', 'i4'),
                            ('W35', 'i4'),
                            ('Validity Code', 'i4'),
                            ('W37', 'i4'),
                            ('W38', 'i4'),
                            ('W39', 'i4'),
                            ('W40', 'i4'),
                            ('W41', 'i4'),
                            ('W42', 'i4'),
                            ('W43', 'i4'),
                            ('W44', 'i4'),
                            ('W45', 'i4'),
                            ('W46', 'i4'),
                            ('W47', 'i4'),
                            ('W48', 'i4'),
                            ('W49', 'i4'),
                            ('W50', 'i4'),
                            ('W51', 'i4'),
                            ('Image source type', 'a4'),
                            ('Calibration Type', 'a4'),
                            ('W54', 'i4'),
                            ('W55', 'i4'),
                            ('W56', 'i4'),
                            ('W57', 'i4'),
                            ('W58', 'i4'),
                            ('W59', 'i4'),
                            ('W60', 'i4'),
                            ('W61', 'i4'),
                            ('W62', 'i4'),
                            ('W63', 'i4'),
                            ('W64', 'i4'),
                            ])


                            
_Navigation = np.dtype([ ('Navigation type', 'a4'),
                            ('SSSYYDDD', 'i4'),
                            ('HHMMSS', 'i4'),
                            ('Orbit type', 'i4'),
                            ('ETIMY', 'i4'),
                            ('ETIMH', 'i4'),
                            ('SEMIMA', 'i4'),
                            ('ECCEN', 'i4'),
                            ('ORBINC', 'i4'),
                            ('MEANA', 'i4'),
                            ('PERIGEE', 'i4'),
                            ('ASNODE', 'i4'),
                            ('DECLIN', 'i4'),
                            ('RASCEN', 'i4'),
                            ('PICLIN', 'i4'),
                            ('SPINP', 'i4'),
                            ('DEGLIN', 'i4'),
                            ('LINTOT', 'i4'),
                            ('DEGELE', 'i4'),
                            ('ELETOT', 'i4'),
                            ('PITCH', 'i4'),
                            ('YAW', 'i4'),
                            ('ROLL', 'i4'),
                            ('Reserved', 'i4'),
                            ('IAJUST', 'i4'),
                            ('IAJTIM', 'i4'),
                            ('Reserved2', 'i4'),
                            ('ISEANG', 'i4'),
                            ('Reserved3', 'i4'),
                            ('Reserved4', 'i4'),
                            ('Scan line of the first beta', 'i4'),
                            ('Time of the above scan line (BEGINNING)', 'i4'),
                            ('Time of the above scan line (CONTINUED)', 'i4'),
                            ('Beta count 1', 'i4'),
                            ('Scan line of the second beta', 'i4'),
                            ('Time of the above scan line (BEGINNING)2', 'i4'),
                            ('Time of the above scan line (CONTINUED)2', 'i4'),
                            ('Beta count 2', 'i4'),
                            ('GAMMA', 'i4'),
                            ('GAMMA-DOT', 'i4'),
                            ('Rserved5', 'a320'),
                            ('Memo entry', 'a32'),
                            ])      
                            
_Calibration = np.dtype([ ('Satellite number', 'a4'),
                            ('YYDDD', 'i4'),
                            ('HHMMSS', 'i4'),])
                            

_Data_Block = np.dtype([('Validity code', 'i4'),
                        ('IR common documentation', 'i4', (128,)),
                        ('YYDDD', 'i4'),
                        ('HHMMSS', 'i4'),
                        ('Scan number', 'i4'),
                        ('Band data', 'i2', (13,4)),
                        ])                               
                        
                        
channel_to_band = {0: 0,
    1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9, 10: 10,
    11: 11, 12: 12, 13: 1, 14: 2, 15: 3, 16: 4, 17: 5, 18: 6, 19: 7,
    20: 8, 21: 9, 22: 10, 23: 11, 24: 12, 25: 3, 26: 4, 27: 5, 28: 7,
    29: 8, 30: 9, 31: 10, 32: 3, 33: 4, 34: 5, 35: 7, 36: 8, 37: 9, 38: 10
}


# ---- 新增：GOES-7 AAA 校正區塊 (128 words，McIDAS Programmer's Manual「AAA calibration block」) ----
# 輻射 R = (IAB[ch, 1] x P / 32 - IAB[ch, 0]) / 2^(15 - IFAB[ch])，ch = 頻道 (偵測器) 1~38 - 1
_Calibration_AAA = np.dtype([('Satellite number', 'i4'), ('YYDDD', 'i4'), ('HHMMSS', 'i4'),
                             ('IAB', 'i4', (38, 2)), ('IFAB', 'i4', (38,)), ('Spare', 'i4', (11,))])
