from dolfin import *
from mshr import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from ufl import inner, grad, dx

from dolfin_solver import *

# Carre
# test_points = [ [0.,0.], [1.,0.], [1.,1.], [0.,1.], [0.,0.]]

# L
test_points = [ [0.,0.], [1.,0.], [1.,1.], [2.,1.,], [2.,2.], [0.,2.], [0.,0.]]

# Star
# test_points = [ [0.,0.], [0.5, 0.5], [1.,0.], [0.75, 0.75], [1., 1.],
#                 [0.75, 1.], [0.5, 1.5], [0.25, 1.],
#                 [0.,1.],[0.25, 0.75], [0.,0.]]

# Test
# test_points = [[358.0, 560.0], [563.0, 342.0], [485.0, 215.0], [348.0, 242.0], [256.0, 194.0], [190.0, 271.0], [338.0, 557.0]]


N = 60
lx = 0.2
ly = 0.1
lz = 1.0
# testmesh = UnitSquareMesh(100,100)
testmesh = generate_mesh_from_coords(test_points, 50)
plot(testmesh)
plt.show()

eigenvals, v = eigenpair_solver(testmesh, 25)


# Function space and mesh
V = FunctionSpace(testmesh, 'Lagrange', 1)

# Boundary conditions
boundary_mesh = BoundaryMesh(testmesh, "exterior")

def boundary(x, on_boundary):
    return on_boundary

bc = DirichletBC(V, 0, boundary)

w = Function(V)
w.vector()[:] = 0

xyz = V.tabulate_dof_coordinates()
x = xyz[:,0]
y = xyz[:,1]

mesh_points = []

It_mesh = vertices(testmesh)
for c in It_mesh:
    mesh_points.append(c.point().x())

dof_points = V.dofmap().dofs()

print(dof_points)
print(mesh_points)

# w.set_local()

l = w.vector()
l[[1500]] = 1
l.get_local()
w.vector()[:]=(l)

plot(w)
plt.show()

# for i in range(0, 25):
#     plt.subplot(5, 5, i+1)

#     plot(v[i])

# plt.show()
