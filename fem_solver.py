import getfem as gf
import numpy as np

def generate_mesh(rects, h, K):
    # rects is array or [[bottom corner coords], [top corner coords]]
    union = gf.MesherObject("rectangle", rects[0][0], rects[0][1])
    for i in range(1, len(rects) - 1):
        temp_rect = gf.MesherObject("rectangle", rects[i][0], rects[i][1])
        union = gf.MesherObject('union', union, temp_rect)

    mesh = gf.Mesh('generate', union, h, K)

    return mesh

