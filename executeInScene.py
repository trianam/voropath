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
    plotter = plotter.Plotter()
    
    #voronoi.plotSites(plotter, verbose = True)
    voronoi.plotPolyhedrons(plotter, verbose = True)
    #voronoi.plotGraph(plotter, edges=False, labels=True)
    #voronoi.plotGraph(plotter, pathExtremes=True)
    #voronoi.plotGraph(plotter)
    voronoi.plotShortestPath(plotter, verbose = True)

    print('Render', flush=True)
    plotter.draw()

else:
    print('use: {} sceneFile [startPoint endPoint]'.format(sys.argv[0]))
