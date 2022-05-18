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
    trim = 300
    xs = [float(a) for a in a.split(';')[0][:-1].split(':')]
    ys = [float(a) for a in a.split(';')[1][:-1].split(':')]
    xs = xs[trim:-trim]
    ys = ys[trim:-trim]
    ys[0], ys[-1] = 0, 0
    layer = list(zip(xs,ys))
    verts.append(layer)


ts = list(range(len(verts)))

poly = PolyCollection(verts,facecolors='white',edgecolor='black',lw=1)
poly.set_alpha(0.7)
ax.add_collection3d(poly, zs=ts, zdir='y')

ax.xaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.yaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.zaxis.set_pane_color((1.0, 1.0, 1.0, 0.0))
ax.xaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.yaxis._axinfo["grid"]['color'] =  (1,1,1,0)
ax.zaxis._axinfo["grid"]['color'] =  (1,1,1,0)


ax.set_xlabel('x')
ax.set_xlim3d(verts[0][0][0], verts[0][-1][0])
ax.set_ylabel('t')
ax.set_ylim3d(-1, 10)
ax.set_zlabel('Abs u(x,t)')
ax.set_zlim3d(-3, 3)



fig.savefig("waterfall.png",dpi=500)


#plt.plot(x ,y,c='black')



plt.show()
