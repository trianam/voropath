#!/bin/python

import sys
import numpy as np
import random
import math
import pickle
import voronizator
import convexHull

if len(sys.argv) == 2:
    fileName = sys.argv[1]

else:
    fileName = input('Insert file name: ')

minX = -1
minY = -1
minZ = -1
maxX = 3
maxY = 2
maxZ = 2
maxEmptyArea = 0.5

voronoi = voronizator.Voronizator()
points = np.array([[1.,0.,1.], [1.,1.,0], [0.2,0.8,0.2], [0.,1.,1.], [0.4,0.8,0.9], [1.5,0.5,0.5], [0.,0.,0.], [1.,1.,1.], [0.6,0.1,0.7], [2.,0.5,0.5], [0.,1.,0.], [0.2,0.5,0.3], [0.,0.,1.], [1.,0.,0.]])

newObstacle = convexHull.ConvexHull(points=points, distributePoints = True, maxEmptyArea = maxEmptyArea)
voronoi.addPolyhedron(newObstacle)

voronoi.addBoundingBox([minX, minY, minZ], [maxX, maxY, maxZ], maxEmptyArea, verbose=True)

voronoi.setPolyhedronsSites(verbose=True)
voronoi.makeVoroGraph(verbose=True)

print('Write file', flush=True)
record = {}
record['voronoi'] = voronoi
record['minX'] = minX
record['minY'] = minY
record['minZ'] = minZ
record['maxX'] = maxX
record['maxY'] = maxY
record['maxZ'] = maxZ
with open(fileName, 'wb') as f:
    pickle.dump(record, f)
