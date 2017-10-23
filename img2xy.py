#!/usr/bin/python3
from PIL.Image import *
import sys

if len(sys.argv) < 2:
    print("%s image.png [treshold]" % sys.argv[0])
    exit(42)

treshold = 127
if len(sys.argv) > 2:
    treshold = (int)(sys.argv[2])

img=open(sys.argv[1])
(w, h) = img.size
for y in range(h):
    for x in range(w):
        r, g, b = img.getpixel((x, y))
        value = sum(img.getpixel((x, y)))
        if value < 3 * treshold:
            print(x, y)
