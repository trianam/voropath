#!/bin/python

import sys
import pickle
import plotter

if len(sys.argv) >= 2:
    if len(sys.argv) == 4:
        i = 2
        plotSites = bool(eval(sys.argv[i]))
        i += 1
        plotGraph = bool(eval(sys.argv[i]))
    else:
        plotSites = bool(eval(input('Do you want to plot Voronoi sites? (True/False): ')))
        plotGraph = bool(eval(input('Do you want to plot graph? (True/False): ')))


    print('Load file', flush=True)
    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    voronoi = record['voronoi']

    print('Build renderer, window and interactor', flush=True)
    plt = plotter.Plotter()
    
    voronoi.plotPolyhedrons(plt, verbose = True)
    if plotSites:
        voronoi.plotSites(plt, verbose = True)
    if plotGraph:
        voronoi.plotGraph(plt, verbose = True)

    print('Render', flush=True)
    plt.draw()
    
else:
    print('use: {} sceneFile [plotSites plotGraph]'.format(sys.argv[0]))
