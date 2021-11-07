from dolfin import *
from dolfin.cpp.mesh import *
import scipy
import numpy as np
import matplotlib.pyplot as plt
from ufl import inner, grad, dx

from dolfin_solver import *

test_points = [ [0.,0.], [1.,0.], [1.,1.], [0.,1.], [0.,0.]]

N = 60
lx = 0.2
ly = 0.1
lz = 1.0
testmesh = UnitSquareMesh(100,100)


eigenvals, v = eigenpair_solver(testmesh, 25)

for i in range(0, 25):
    plt.subplot(5, 5, i+1)

    plot(v[i])

plt.show()
