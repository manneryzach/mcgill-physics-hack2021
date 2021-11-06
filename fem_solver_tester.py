import getfem as gf
import numpy as np

import pyvista as pv
from pyvirtualdisplay import Display

from fem_solver import *

rects = [ [ [0.,0.], [1.,1.] ], [ [0.,1.], [1.,2.] ] ]

mesh1 = generate_mesh(rects, 0.1, 1)

# sl = gf.Slice(("none",), mesh1, 1)
# sl.export_to_vtk('sl.vtk', "ascii")

# # Display
# display = Display(visible=0, size=(1280, 1024))
# display.start()
# p = pv.Plotter()
# m = pv.read("./sl.vtk")
# p.add_mesh(m, show_edges=True)
# pts = m.points
# p.show(window_size=[512, 384], cpos="xy")
# display.stop()


