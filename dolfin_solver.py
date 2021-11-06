from dolfin import *
from dolfin.cpp.mesh import *
import scipy
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

def eigenvalue_solver(mesh):
    # Function space
    V = FunctionSpace(mesh, 'Lagrange', 1)

    # Boundary conditions
    boundary = BoundaryMesh(mesh)
    bc = DirichletBC(V, 0, boundary)

    # Test/Trial functions
    u_h = TrialFunction(V)
    v = TestFunction(V)

    # For variationnal problem
    a = inner(grad(u_h), grad(v)) * dx
    b = inner(u_h, v) *dx
    dummy = inner(1., v) * dx

    # Assemble system
    asm = SystemAssembler(b, dummy, bcs)
    B = PETScMatrix()
    asm.assemble(B)

    diag_value = 1e6

    A = PETScMatrix()
    assemble(a, tensor=A)

    for bc in bcs:
        bc.zero(A)
        bc.zero_columns(A, dummy_vec, diag_values)


    # Scipy solve the eigenvalue equation
    A_array = scipy.sparse.csr_matrix(A.mat().getValuesCSR()[::-1])
    B_array = scipy.sparse.csr_matrix(B.mat().getValuesCSR()[::-1])

    k = 20
    which = 'SM'
    w, v = scipy.linalg.eigh(A_array, B_array, k=k, which=which)

    return w

eigeinvals = eigenvalue_solver(testmesh)
print(eigenvals)

