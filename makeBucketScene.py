#!/bin/python

import sys
import numpy as np
import random
import math
import pickle
import voronizator
import bucket

if len(sys.argv) == 9:
    i = 1
    minPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
    i += 1
    maxPoint = np.array(tuple(eval(sys.argv[i])),dtype=float)
    i += 1
    center = np.array(tuple(eval(sys.argv[i])),dtype=float)
    i += 1
    width = float(sys.argv[i])
    i += 1
    height = float(sys.argv[i])
    i += 1
    thickness = float(sys.argv[i])
    i += 1    
    maxEmptyArea = float(sys.argv[i])
    i += 1
    fileName = sys.argv[i]

else:
    minPoint = np.array(tuple(eval(input('Insert min point (x,y,z): '))),dtype=float)
    maxPoint = np.array(tuple(eval(input('Insert max point (x,y,z): '))),dtype=float)
    center = np.array(tuple(eval(input('Insert bucket center point (x,y,z): '))),dtype=float)
    width = float(input('Insert bucket width: '))    
    height = float(input('Insert bucket height: '))    
    thickness = float(input('Insert bucket thickness: '))    
    maxEmptyArea = float(input('Insert max empty area (for points distribution in obstacles): '))    
    fileName = input('Insert file name: ')

voronoi = voronizator.Voronizator()

print('Create bucket', flush=True)
voronoi.addPolyhedron(bucket.Bucket(center, width, height, thickness, distributePoints=True, maxEmptyArea=maxEmptyArea))
voronoi.addBoundingBox(minPoint, maxPoint, maxEmptyArea, verbose=True)
voronoi.setPolyhedronsSites(verbose=True)
voronoi.makeVoroGraph(verbose=True)

print('Write file', flush=True)
record = {}
record['voronoi'] = voronoi
# record['minX'] = minX
# record['minY'] = minY
# record['minZ'] = minZ
# record['maxX'] = maxX
# record['maxY'] = maxY
# record['maxZ'] = maxZ
with open(fileName, 'wb') as f:
    pickle.dump(record, f)
