import getfem as gf
import numpy as np

import pyvista as pv
from pyvirtualdisplay import Display

from fem_solver import *

rects = [ [ [0.,0.], [1.,1.] ], [ [0.,1.], [1.,2.] ] ]

mesh1 = generate_mesh(rects, 0.1, 1)
mfu = gf.MeshFem(mesh1, 1)

U = model(mesh1, 2.)
print(U)

sl = gf.Slice(("none",), mesh1, 1)
sl.export_to_vtk('sl.vtk', "ascii", mfu, U, "U")



