#!/bin/python

import sys
import numpy as np
import pickle
import plotter

if len(sys.argv) >= 2:
    if len(sys.argv) == 4:
        i = 2
        startPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
        i += 1
        endPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
    else:
        startPoint = np.array(tuple(eval(input('Insert start point (x,y,z): '))),dtype=float)
        endPoint = np.array(tuple(eval(input('Insert end point (x,y,z): '))),dtype=float)

    print('Load file', flush=True)
    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    voronoi = record['voronoi']
    voronoi.setBsplineDegree(2)

    print('Calculate shortest path', flush=True)
    voronoi.calculateShortestPath(startPoint, endPoint, 'near', postSimplify=False, verbose=True, debug=False)

    print('Build renderer, window and interactor', flush=True)
    plt = plotter.Plotter()
    
    #voronoi.plotSites(plt, verbose = True)
    voronoi.plotPolyhedrons(plt, verbose = True)
    #voronoi.plotGraph(plt, verbose = True)
    voronoi.plotShortestPath(plt, verbose = True)

    print('Render', flush=True)
    plt.draw()

else:
    print('use: {} sceneFile [startPoint endPoint]'.format(sys.argv[0]))
