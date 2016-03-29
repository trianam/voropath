#!/bin/python

import sys
import pickle
from mpl_toolkits.mplot3d import  Axes3D
import matplotlib.pyplot as plt


if len(sys.argv) >= 2:

    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    fig = plt.figure()
    ax = Axes3D(fig)

    voronoi = record['voronoi']

    voronoi.plotSites(ax)
    voronoi.plotPolyhedrons(ax)
    #voronoi.plotGraph(ax, edges=False, labels=True)
    #voronoi.plotGraph(ax, pathExtremes=True)
    voronoi.plotGraph(ax)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.set_xlim3d(record['minX'], record['maxX'])
    ax.set_ylim3d(record['minY'], record['maxY'])
    ax.set_zlim3d(record['minZ'], record['maxZ'])

    plt.show()

else:
    print('use: {} sceneFile'.format(sys.argv[0]))
