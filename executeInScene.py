#!/bin/python

import sys
import numpy as np
import pickle
from mpl_toolkits.mplot3d import  Axes3D
import matplotlib.pyplot as plt

if len(sys.argv) >= 2:
    if len(sys.argv) == 4:
        i = 2
        startPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
        i += 1
        endPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
    else:
        startPoint = np.array(tuple(eval(input('Insert start point (x,y,z): '))),dtype=float)
        endPoint = np.array(tuple(eval(input('Insert end point (x,y,z): '))),dtype=float)

    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    fig = plt.figure()
    ax = Axes3D(fig)

    voronoi = record['voronoi']
    
    print('Calculate shortest path', flush=True)
    voronoi.calculateShortestPath(startPoint, endPoint, 'near', verbose=True, debug=False)

    print('Plot', flush=True)
    #voronoi.plotSites(ax)
    voronoi.plotPolyhedrons(ax)
    #voronoi.plotGraph(ax, edges=False, labels=True)
    #voronoi.plotGraph(ax, pathExtremes=True)
    #voronoi.plotGraph(ax)
    voronoi.plotShortestPath(ax)

    ax.set_xlabel('x')
    ax.set_ylabel('y')
    ax.set_zlabel('z')

    ax.set_xlim3d(record['minX'], record['maxX'])
    ax.set_ylim3d(record['minY'], record['maxY'])
    ax.set_zlim3d(record['minZ'], record['maxZ'])
    
    plt.show()

else:
    print('use: {} sceneFile [startPoint endPoint]'.format(sys.argv[0]))
