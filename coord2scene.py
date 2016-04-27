#!/bin/python

import sys
import pickle
import xml.etree.cElementTree as ET
import voronizator

if len(sys.argv) == 4:
    xmlFileName = sys.argv[1]
    sceneFileName = sys.argv[2]
    maxEmptyArea = float(sys.argv[3])

    xmlRoot = ET.parse(xmlFileName).getroot()

    voronoi = voronizator.Voronizator()

    print('Import XML', flush=True)
    voronoi.importXmlTree(xmlRoot, maxEmptyArea)

    print('Set sites and make graph', flush=True)
    voronoi.setPolyhedronsSites(verbose=True)
    voronoi.makeVoroGraph(verbose=True)

    print('Write file', flush=True)
    record = {}
    record['voronoi'] = voronoi
    with open(sceneFileName, 'wb') as f:
        pickle.dump(record, f)

else:
    print('use: {} coordinateFile sceneFile maxEmptyArea'.format(sys.argv[0]))
    
