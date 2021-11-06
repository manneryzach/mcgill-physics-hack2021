import getfem as gf
import numpy as np

def generate_mesh(rects, h, K):
    # rects is array or [[bottom corner coords], [top corner coords]]
    union = gf.MesherObject("rectangle", rects[0][0], rects[0][1])
    for i in range(1, len(rects)):
        temp_rect = gf.MesherObject("rectangle", rects[i][0], rects[i][1])
        union = gf.MesherObject('union', union, temp_rect)

    # union of all the rectangles
    mesh = gf.Mesh('generate', union, h, K)

    # boundaries
    outer_faces = mesh.outer_faces()
    mesh.set_region(1, outer_faces)

    return mesh

def model(mesh, k):
    # Separation into finite elements
    mfu = gf.MeshFem(mesh, 1)

    # Integration method
    elements_degree = 2
    mfu.set_classical_fem(elements_degree)

    mim = gf.MeshIm(mesh, pow(elements_degree, 2))

    # Model
    md = gf.Model("complex")
    md.add_fem_variable("u", mfu)

    # Laplacian term of the equation
    md.add_Laplacian_brick(mim, "u")

    md.add_fem_data("k", mfu)

    # # Eigenvalue term of the equation
    md.add_Helmholtz_brick(mim, "u", "k")

    # Poisson equation (for testing)
    # F = 1
    # md.add_fem_data("F", mfu)
    # md.set_variable("F", np.repeat(F, mfu.nbdof()))
    # md.add_source_term_brick(mim, "u", "F")

    # Dirichlet boundary conditions
    md.add_Dirichlet_condition_with_multipliers(mim, "u", (elements_degree - 1), 1)

    md.display()
    # Solve the model
    md.solve()

    U = md.variable("u")

    return U
