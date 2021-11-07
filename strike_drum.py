from dolfin import *
from mshr import *
import numpy as np
from scipy.spatial import distance
from sklearn.neighbors import NearestNeighbors
from ufl import inner, grad, dx
import matplotlib.pyplot as plt
import math

def strike_drum(coord, mesh, eig_funcs, eig_vals, ntimesteps):

    # Function space and mesh
    V = FunctionSpace(mesh, 'Lagrange', 1)

    # Boundary conditions
    boundary_mesh = BoundaryMesh(mesh, "exterior")

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, 0, boundary)

    p = Function(V)
    p.vector()[:] = 0

    # get mesh and degree of freedom points
    mesh_points = []

    It_mesh = vertices(mesh)
    for c in It_mesh:
        mesh_points.append([c.point().x(), c.point().y()])

    # mesh_points.append(coord)

    dof_points = V.dofmap().dofs()

    # find closest node and nearest neighbors
    closest_index = [distance.cdist([coord], mesh_points).argmin()]

    # knn = NearestNeighbors(n_neighbors=5)
    # knn.fit(mesh_points)
    # closest_index = knn.kneighbors([mesh_points[-1]], return_distance=False)[0]

    # closest_index = np.delete(closest_index, 0)

    # Assign strike to map
    l = p.vector()
    l[closest_index] = 1
    l.get_local()
    p.vector()[:] = l

    # Calculate coeffs
    c_n = []
    for eig_func in eig_funcs:
        t_c_n = dot(eig_func, p)

        # t_c_n = eig_func * p

        c_n.append(t_c_n)

    # print(c_n)
    # animate
    dt = 1
    for i in range(0, ntimesteps):
        t = i*dt

        # general solution
        gen_sol = Function(V)

        for j in range(len(c_n)):

            gen_sol = gen_sol + c_n[j]*eig_funcs[j]*math.cos(eig_vals[j]*t)

        plot(gen_sol)
        plt.show()



