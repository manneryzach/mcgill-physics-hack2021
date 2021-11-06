from dolfin import *

test_points = [ [0.,0.], [1.,0.], [1.,1.], [0.,1.], [0.,0.]]

def generate_mesh(coords):

    # convert array of coords to Dolfin point types
    points = []
    for coord in coords:
        points.append(Point(coord[0], coord[1]))

    PolygonalMeshGenerator.generate(mesh, domain_vertices, 0.1)

    return mesh

mesh = generate_mesh(test_points)
plot(mesh)
