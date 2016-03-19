#!/bin/python

import sys
import numpy as np
import random
import math
import pickle
import polyhedron

if len(sys.argv) >= 13 or len(sys.argv) <= 15:
    i = 1
    minX = float(sys.argv[i])
    i += 1
    minY = float(sys.argv[i])
    i += 1
    minZ = float(sys.argv[i])
    i += 1
    maxX = float(sys.argv[i])
    i += 1
    maxY = float(sys.argv[i])
    i += 1
    maxZ = float(sys.argv[i])
    i += 1
    fixedRadius = bool(eval(sys.argv[i]))
    i += 1
    if fixedRadius:
        radius = float(sys.argv[i])
        i += 1
    else:
        minRadius = float(sys.argv[i])
        i += 1
        maxRadius = float(sys.argv[i])
        i += 1
    avoidCollisions = bool(eval(sys.argv[i]))
    i += 1
    numObstacles = int(sys.argv[i])
    i += 1
    distributePoints = bool(eval(sys.argv[i]))
    i += 1
    if distributePoints:
        maxEmptyArea = float(sys.argv[i])
        i += 1
    fileName = sys.argv[i]
    
else:
    minX = float(input('Insert min scene X: '))
    minY = float(input('Insert min scene Y: '))
    minZ = float(input('Insert min scene Z: '))
    maxX = float(input('Insert max scene X: '))
    maxY = float(input('Insert max scene Y: '))
    maxZ = float(input('Insert max scene Z: '))
    fixedRadius = bool(eval(input('Do you want fixed obstacle radius? (True/False): ')))
    if fixedRadius:
        radius = float(input('Insert obstacle radius: '))
    else:
        minRadius = float(input('Insert min obstacle radius: '))
        maxRadius = float(input('Insert max obstacle radius: '))
    avoidCollisions = bool(eval(input('Do you want to avoid collisions between obstacles? (True/False): ')))
    numObstacles = int(input('Insert obstacles number: '))
    distributePoints = bool(eval(input('Do you want to distribute points in obstacles surfaces? (True/False): ')))
    if distributePoints:
        maxEmptyArea = float(input('Insert max empty area: '))
    fileName = input('Insert file name: ')

if not distributePoints:
    maxEmptyArea = 0.
    
obstacles = []

for ob in range(numObstacles):
    print('Creating obstacle {} '.format(ob+1), end='')
    ok = False
    while not ok:
        print('.', end='')
        if not fixedRadius:
            radius = random.uniform(minRadius,maxRadius)
        center = np.array([random.uniform(minX+radius,maxX-radius), random.uniform(minY+radius,maxY-radius), random.uniform(minZ+radius,maxZ-radius)])
        points = []
        for pt in range(4):
            elev = random.uniform(-math.pi/2., math.pi/2.)
            azim = random.uniform(0., 2.*math.pi)
            points[:0] = [center+np.array([
                radius*math.cos(elev)*math.cos(azim),
                radius*math.cos(elev)*math.sin(azim),
                radius*math.sin(elev)])]

        newObstacle = polyhedron.Polyhedron(a = points[0], b = points[1], c = points[2], d = points[3], distributePoints = distributePoints, maxEmptyArea = maxEmptyArea)

        ok = True
        if avoidCollisions:
            for obstacle in obstacles:
                if newObstacle.intersectPolyhedron(obstacle):
                    ok = False
                    break

        if ok:
            obstacles[:0] = [newObstacle]
            
    print(' done')
            
record = {}
record['obstacles'] = obstacles
record['minX'] = minX
record['minY'] = minY
record['minZ'] = minZ
record['maxX'] = maxX
record['maxY'] = maxY
record['maxZ'] = maxZ
with open(fileName, 'wb') as f:
    pickle.dump(record, f)
