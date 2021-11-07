from dolfin import *
from mshr import *
import numpy as np
from scipy.spatial import distance
from ufl import inner, grad, dx
import matplotlib.pyplot as plt

def strike_drum(coord, mesh, eig_funcs, eig_vals):

    # Function space and mesh
    V = FunctionSpace(testmesh, 'Lagrange', 1)

    # Boundary conditions
    boundary_mesh = BoundaryMesh(testmesh, "exterior")

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, 0, boundary)

    p = Function(V)
    p.vector()[:] = 0

    # get mesh and degree of freedom points
    mesh_points = []

    It_mesh = vertices(testmesh)
    for c in It_mesh:
        mesh_points.append([c.point().x(), c.point().y()])

    dof_points = V.dofmap().dofs()

    # find closest node
    def closest_node(coord, mesh_points):
        closest_index = distance.cdist([node], nodes).argmin()

    # Assign strike to map
    l = p.vector()
    l[[closest_index]] = 5
    l.get_local()
    p.vector()[:] = l

    # Calculate coeffs
    c_n = []
    for eig_func in eig_funcs:
        t_c_n = inner(eig_func, p)
        c_n.append(t_c_n)

    # animate
    dt = 0.1
    ntimesteps = 10
    for i in range(0, ntimesteps):
        t = i*dt

        # general solution
        gen_sol = Function(V)

        for j in range(len(c_n)):

            gen_sol = gen_sol + c_n[j]*eig_funcs[j]*math.cos(eig_vals[j]*t)

            plot(gen_sol)
            plt.show()



