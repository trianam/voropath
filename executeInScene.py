#!/bin/python

import sys
import numpy as np
import pickle
from mpl_toolkits.mplot3d import  Axes3D
import matplotlib.pyplot as plt
import voronizator

if len(sys.argv) >= 2:
    if len(sys.argv) == 5 or len(sys.argv) == 6:
        i = 1
        startPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
        i += 1
        endPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
        i += 1
        distributePoints = bool(eval(sys.argv[i]))
        i += 1
        if distributePoints:
            maxEmptyArea = float(sys.argv[i])
            i += 1
    else:
        startPoint = np.array(tuple(eval(input('Insert start point (x,y,z): '))),dtype=float)
        endPoint = np.array(tuple(eval(input('Insert end point (x,y,z): '))),dtype=float)
        distributePoints = bool(eval(input('Do you want to (re)distribute points in obstacles surfaces? (True/False): ')))
        if distributePoints:
            maxEmptyArea = float(input('Insert max empty area: '))

    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    fig = plt.figure()
    ax = Axes3D(fig)

    voronoi = voronizator.Voronizator()
    
    for obstacle in record['obstacles']:
        if distributePoints:
            obstacle.distributePoints(maxEmptyArea)
        voronoi.addPolyhedron(obstacle)

    print('Set sites for Voronoi')
    voronoi.setPolyhedronsSites()
    print('Make pruned voronoi Graph')
    voronoi.makeVoroGraph()
    print('Calculate shortest path')
    voronoi.calculateShortestPath(startPoint, endPoint, 'near')

    print('Plot')
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
    print('use: {} sceneFile [startPoint endPoint (False|True maxEmptyArea)]'.format(sys.argv[0]))
