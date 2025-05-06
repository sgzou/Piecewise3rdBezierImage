from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import os
from os.path import isfile, join

dirpath = os.path.dirname(os.path.abspath(__file__))
os.chdir(dirpath)
onlyfiles = [f for f in os.listdir("data") if isfile(join("data", f)) and "_coord.txt" in f]

for file in onlyfiles:

    fig = plt.figure()
    ax = plt.axes(projection='3d')
    ax.grid(True)

    lines = np.loadtxt("data/"+file)

    x = lines[:,0].tolist()
    y = lines[:,1].tolist()
    z = lines[:,2].tolist()
    ax.plot3D(x, y, z, 'green')

    ax.set_title(file)
    ax.set_xlabel("X")
    ax.set_ylabel("Y")

    plt.show()
