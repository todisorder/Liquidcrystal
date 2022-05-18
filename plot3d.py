from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np


fig = plt.figure()
ax = fig.gca(projection='3d')


#with open('Sol3DUURe.txt', 'r') as f:
#    lines = f.readlines()

with open('Sol3DUUAbs.txt', 'r') as f:
    lines = f.readlines()


verts=[]
#print(len(lines))
#   data format: each line is a time, and is
#   x_1:x_2:x_3:...:x_n;y_1:y_2:y_3:...:y_n:
#   poly must be:
#   [[(x_1,0),(x_2,y_2),...(x_n,0)],[(x_1,0),(x_2,y_2),...(x_n,0)]]

for a in lines:
    a=a.strip()
    xs = [float(a) for a in a.split(';')[0][:-1].split(':')]
    ys = [float(a) for a in a.split(';')[1][:-1].split(':')]
    ys[0], ys[-1] = 0, 0
    layer = list(zip(xs,ys))
    verts.append(layer)

ts = list(range(len(verts)))

poly = PolyCollection(verts,facecolors='white',edgecolor='black',lw=1)
#poly.set_alpha(0.7)
ax.add_collection3d(poly, zs=ts, zdir='y')

ax.set_xlabel('X')
ax.set_xlim3d(-10, 10)
ax.set_ylabel('Y')
ax.set_ylim3d(-1, 10)
ax.set_zlabel('Z')
ax.set_zlim3d(-3, 3)






#plt.plot(x ,y,c='black')



plt.show()
