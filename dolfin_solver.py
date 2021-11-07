from dolfin import *
from dolfin.cpp.mesh import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from ufl import inner, grad, dx

def generate_mesh(coords):

    # convert array of coords to Dolfin point types
    points = []
    for coord in coords:
        points.append(Point(coord[0], coord[1]))

    PolygonalMeshGenerator.generate(mesh, domain_vertices, 0.1)

    return mesh

def eigenpair_solver(mesh, n_eig):
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

    solver.solve(n_eig)

    w, v = [], []

    for i in range(solver.get_number_converged()):
        r, _, rv, _ = solver.get_eigenpair(i)
        w.append(r)

        u = Function(V)
        u.vector()[:] = rv

        v.append(u)

    w = np.array(w)

    return w, v

