import os
import sys
from PIL import Image, ImageDraw

d = sys.argv[1]
styles = sys.argv[2].split(",")
ligs = sys.argv[3].split(",")
s = 250
im = Image.new("RGB", (s * len(styles), s * len(ligs) + 20), "white")
dr = ImageDraw.Draw(im)
for j, st in enumerate(styles):
    dr.text((j * s + 5, 3), st, fill="black")
    for i, lig in enumerate(ligs):
        p = f"{d}/{lig}_{st}.png"
        if os.path.exists(p):
            im.paste(Image.open(p).convert("RGB").resize((s, s)), (j * s, 20 + i * s))
            dr.text((j * s + 5, 22 + i * s), lig, fill="grey")
im.save(f"{d}/sheet.png")
