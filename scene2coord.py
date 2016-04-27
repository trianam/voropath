#!/bin/python

import sys
import pickle
import xml.etree.cElementTree as ET

if len(sys.argv) == 3:

    print('Load file', flush=True)
    with open(sys.argv[1], 'rb') as f:
        record = pickle.load(f)

    voronoi = record['voronoi']

    print('Create XML', flush=True)
    xmlRoot = ET.Element('scene')
    voronoi.extractXmlTree(xmlRoot)
    xmlTree = ET.ElementTree(xmlRoot)

    print('Write file', flush=True)
    xmlTree.write(sys.argv[2])
        
else:
    print('use: {} sceneFile coordinateFile'.format(sys.argv[0]))
