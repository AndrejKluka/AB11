from matplotlib import pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from skimage import measure

n=10
x = np.linspace(0, n, n+1)

X, Y, Z = np.meshgrid(x, x, x, indexing = 'ij')

def f(x, y, z):
    return (x**2+y**2+z**2)**0.5 - 3


verts, faces = measure.marching_cubes_classic(f(X, Y, Z), 5)

fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

mesh = Poly3DCollection(verts[faces])
ax.add_collection3d(mesh)

m=00

ax.set_xlim(-m, n+5)
ax.set_ylim(-m, n+5)
ax.set_zlim(-m, n+5)

plt.show()


