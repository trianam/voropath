#!/bin/python

import sys
import numpy as np
import pickle
import plotter

if len(sys.argv) >= 2:
    if len(sys.argv) == 8:
        i = 2
        startPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
        i += 1
        endPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
        i += 1
        bsplineDegree = int(sys.argv[i])
        i += 1
        useTrijkstra = bool(eval(sys.argv[i]))
        i += 1
        postSimplify = bool(eval(sys.argv[i]))
        i += 1
        adaptivePartition = bool(eval(sys.argv[i]))
    else:
        startPoint = np.array(tuple(eval(input('Insert start point (x,y,z): '))),dtype=float)
        endPoint = np.array(tuple(eval(input('Insert end point (x,y,z): '))),dtype=float)
        bsplineDegree = int(input('Insert B-spline degree (2 or 4): '))
        useTrijkstra = bool(eval(input('Do you want to use Trijkstra method? (True/False): ')))
        postSimplify = bool(eval(input('Do you want to simplify path? (True/False): ')))
        adaptivePartition = bool(eval(input('Do you want adaptive partition? (True/False): ')))

    print('Load file', flush=True)
    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    voronoi = record['voronoi']
    voronoi.setBsplineDegree(bsplineDegree)

    voronoi.calculateShortestPath(startPoint, endPoint, 'near', useTrijkstra=useTrijkstra, postSimplify=postSimplify, verbose=True, debug=False)

    print('Build renderer, window and interactor', flush=True)
    plt = plotter.Plotter()
    
    #voronoi.plotSites(plt, verbose = True)
    voronoi.plotPolyhedrons(plt, verbose = True)
    #voronoi.plotGraph(plt, verbose = True)
    voronoi.plotShortestPath(plt, adaptivePartition = adaptivePartition, verbose = True)

    print('Render', flush=True)
    plt.draw()

else:
    print('use: {} sceneFile [startPoint endPoint degree(2,4) useTrijkstra postSimplify adaptivePartition]'.format(sys.argv[0]))
