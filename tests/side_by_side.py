"""把兩張圖縮小後左右並排存成一張 (目視比對原版與整合版用)。

  python tests/side_by_side.py <左圖> <右圖> <輸出> [每張的寬度，預設 700]
"""
import sys

from PIL import Image

Image.MAX_IMAGE_PIXELS = None
a, b, out = sys.argv[1:4]
w = int(sys.argv[4]) if len(sys.argv) > 4 else 700
ims = []
for p in (a, b):
    im = Image.open(p).convert('RGB')
    ims.append(im.resize((w, round(im.height * w / im.width)), Image.LANCZOS))
canvas = Image.new('RGB', (2 * w + 10, max(i.height for i in ims)), 'white')
canvas.paste(ims[0], (0, 0))
canvas.paste(ims[1], (w + 10, 0))
canvas.save(out)
print([Image.open(p).size for p in (a, b)])
