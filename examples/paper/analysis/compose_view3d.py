#!/usr/bin/env python3
"""Compose the perspective views of render_3d.py (transparent PNGs N1_3d.png ... N10_3d.png) into one
figure: crop all views to their common bounding box, find the smallest vertical pitch at which neighbouring
views keep a gap of GAP px, and write view3d.tex (N = 1...5 down the left column, 6...10 down the right one)
together with the cropped views c_N1_3d.png ... c_N10_3d.png into the same directory.
usage: compose_view3d.py figure_dir"""

import argparse
import os

import numpy as np
from PIL import Image

ap = argparse.ArgumentParser(description="Compose the perspective views of render_3d.py into view3d.tex.")
ap.add_argument("figure_dir", help="directory holding N1_3d.png ... N10_3d.png; view3d.tex is written there")
d = ap.parse_args().figure_dir
GAP = 14                       # px between the drawings of neighbouring views
WIDTH, COLGAP = 228, 10        # pt
N = range(1, 11)

img = {n: np.array(Image.open(os.path.join(d, "N%d_3d.png" % n))) for n in N}
alpha = {n: img[n][:, :, 3] > 8 for n in N}
union = np.logical_or.reduce([alpha[n] for n in N])
ys, xs = np.where(union)
y0, y1, x0, x1 = ys.min(), ys.max() + 1, xs.min(), xs.max() + 1
for n in N:
    Image.fromarray(img[n][y0:y1, x0:x1]).save(os.path.join(d, "c_N%d_3d.png" % n))
A = {n: alpha[n][y0:y1, x0:x1] for n in N}
H, W = A[1].shape


def dilate(a):
    b = a.copy()
    for k in range(1, GAP + 1):
        b[k:, :] |= a[:-k, :]
        b[:-k, :] |= a[k:, :]
    return b


def pitch(a, b):
    """smallest offset p (px) of view b below view a without contact"""
    db = dilate(A[b])
    for p in range(H + 1):
        if not np.any(A[a][p:, :] & db[:H - p, :]):
            return p
    return H


p = max(pitch(a, a + 1) for a in (1, 2, 3, 4, 6, 7, 8, 9))
dy = p * WIDTH / W
nodes = []
for k, n in enumerate(N):
    col, row = divmod(k, 5)
    nodes.append("  \\node[anchor=north west,inner sep=0] at (%d,%.2f) "
                 "{\\includegraphics[width=%dpt]{c_N%d_3d.png}};\n" % ((WIDTH + COLGAP) * col, -row * dy, WIDTH, n))
with open(os.path.join(d, "view3d.tex"), "w") as f:
    f.write("\\documentclass[tikz,border=0pt]{standalone}\n\\usepackage{graphicx}\n\\begin{document}\n"
            "\\begin{tikzpicture}[x=1pt,y=1pt]\n"
            "  %% transparent perspective views, N = 1...5 down the left column and 6...10 down the right one;\n"
            "  %% vertical pitch %d of %d px: the smallest one with a %d px gap between neighbouring views\n"
            % (p, H, GAP) + "".join(nodes) + "\\end{tikzpicture}\n\\end{document}\n")
print("common box %d x %d px, pitch %d px = %.2f pt" % (W, H, p, dy))
