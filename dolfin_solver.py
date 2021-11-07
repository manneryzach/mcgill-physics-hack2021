from dolfin import *
from dolfin.cpp.mesh import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from ufl import inner, grad, dx

test_points = [ [0.,0.], [1.,0.], [1.,1.], [0.,1.], [0.,0.]]

# def generate_mesh(coords):

#     # convert array of coords to Dolfin point types
#     points = []
#     for coord in coords:
#         points.append(Point(coord[0], coord[1]))

#     PolygonalMeshGenerator.generate(mesh, domain_vertices, 0.1)

#     return mesh

N = 60
lx = 0.2
ly = 0.1
lz = 1.0
testmesh = UnitSquareMesh(100,100)

def eigenpair_solver(mesh):
    # Function space
    V = FunctionSpace(mesh, 'Lagrange', 1)

    # Boundary conditions
    boundary_mesh = BoundaryMesh(mesh, "exterior")

    def boundary(x, on_boundary):
        return on_boundary

    bc = DirichletBC(V, 0, boundary)

    # Test/Trial functions
    u_h = TrialFunction(V)
    v = TestFunction(V)

    # For variationnal problem
    a = inner(grad(u_h), grad(v)) * dx
    b = inner(u_h, v) *dx
    dummy = inner(1., v) * dx

    # Assemble system
    asm = SystemAssembler(b, dummy, bc)
    B = PETScMatrix()
    asm.assemble(B)

    diag_value = 1e6

    A = PETScMatrix()
    assemble(a, tensor=A)
    dummy_vec = assemble(dummy)

    bc.zero(A)
    bc.zero_columns(A, dummy_vec, diag_value)


    # Solve
    solver = SLEPcEigenSolver(A, B)

    solver.parameters['solver'] = 'krylov-schur'
    solver.parameters['spectrum'] = 'smallest magnitude'
    solver.parameters['problem_type'] = 'gen_hermitian'
    solver.parameters['tolerance'] = 1e-10

    n_eig = 20
    solver.solve(n_eig)

    w, v = [], []

    u = Function(V)
    for i in range(solver.get_number_converged()):
        r, _, rv, _ = solver.get_eigenpair(i)
        w.append(r)

        u.vector()[:] = rv

        plot(u)
        plt.show()

        v.append(u)
        print(u)


    w = np.array(w)
    v = np.array(u)

    return w, v

eigenvals, v = eigenpair_solver(testmesh)
print(eigenvals)
print(v)

# import matplotlib.pyplot as plt

# plt.plot(v[0])
# plt.show()

