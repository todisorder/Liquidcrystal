from matplotlib import *
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
import matplotlib.pyplot as plt
from matplotlib import colors as mcolors
import numpy as np

fig=plt.figure()

with open('./UU-inire.txt', 'r') as f:
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
plt.plot(x ,y,':',lw=1.,c='black')

with open('./UU-iniabs.txt', 'r') as f:
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
plt.plot(x ,y,'-.',lw=1,c='black')

with open('./WW-ini.txt', 'r') as f:
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
plt.plot(x ,y,'--',lw=1,c='black')


with open('./RR-ini.txt', 'r') as f:
    lines = f.readlines()
    x = [float(line.split()[0]) for line in lines]
    y = [float(line.split()[1]) for line in lines]
plt.plot(x ,y,'-',lw=1,c='black')


plt.xlim(-3.1,2.1)
fig.savefig("asdf.png",dpi=500)


plt.show()
