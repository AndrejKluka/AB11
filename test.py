'''from matplotlib import pyplot as plt
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

plt.show()'''
import time


def print_statusline(msg: str):
    last_msg_length = len(print_statusline.last_msg) if hasattr(print_statusline, 'last_msg') else 0
    print(' ' * last_msg_length, end='\r')
    print(msg, end='\r')
    print_statusline.last_msg = msg


    
for i in range(0, 20, 1):
    print_statusline("{}".format(i))
    time.sleep(0.15)
print('')
print(' nanh')    
for i in range(0, 20, 1):
    last_msg_length = len(print_statusline.last_msg) if hasattr(print_statusline, 'last_msg') else 0
    print(' ' * last_msg_length, end='\r')
    print("{}".format(i), end='\r')
    time.sleep(0.15)



'''
for x in range(10):
    for i in range(1000000):
        y=i**0.5
    print("Progress {:2.1%}".format(x / 10), end="\r")


for x in ['abc', 1]:
    for i in range(1000000):
        y=i**0.5
    print ('{}\r'.format(x)),
    
    
'''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    