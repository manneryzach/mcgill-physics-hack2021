from dolfin import *
from dolfin.cpp.mesh import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from ufl import inner, grad, dx

from dolfin_solver import *

# Carre
# test_points = [ [0.,0.], [1.,0.], [1.,1.], [0.,1.], [0.,0.]]

# L
# test_points = [ [0.,0.], [1.,0.], [1.,1.], [2.,1.,], [2.,2.], [0.,2.], [0.,0.]]

# Star
# test_points = [ [0.,0.], [0.5, 0.5], [1.,0.], [0.75, 0.75], [1., 1.],
#                 [0.75, 1.], [0.5, 1.5], [0.25, 1.],
#                 [0.,1.],[0.25, 0.75], [0.,0.]]

# Canadian Flag

N = 60
lx = 0.2
ly = 0.1
lz = 1.0
# testmesh = UnitSquareMesh(100,100)
testmesh = generate_mesh_from_coords(test_points, 50)
plot(testmesh)
plt.show()

eigenvals, v = eigenpair_solver(testmesh, 25)

for i in range(0, 25):
    plt.subplot(5, 5, i+1)

    plot(v[i])

plt.show()
