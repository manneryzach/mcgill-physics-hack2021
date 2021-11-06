import getfem as gf
import numpy as np

import pyvista as pv
from pyvirtualdisplay import Display

from fem_solver import *

# L shaped
rects = [ [ [0.,0.], [1.,1.] ], [ [0.,1.], [2.,2.] ] ]

mesh1 = generate_mesh(rects, 0.1, 1)

mfu = gf.MeshFem(mesh1, 1)

# Integration method
elements_degree = 2
mfu.set_classical_fem(elements_degree)

U = model(mesh1, 2.)
print(U)

sl = gf.Slice(("none",), mesh1, 1)
sl.export_to_vtk('sl.vtk', "ascii")



