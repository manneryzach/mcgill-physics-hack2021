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

# Test
test_points = [[358.0, 560.0], [563.0, 342.0], [485.0, 215.0], [348.0, 242.0], [256.0, 194.0], [190.0, 271.0], [338.0, 557.0]]


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
