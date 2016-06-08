import numpy as np
import numpy.linalg
import scipy as sp
import scipy.spatial
import networkx as nx
import numpy.linalg
import polyhedron
import polyhedronsContainer
import path
import uuid
import xml.etree.cElementTree as ET

class Voronizator:
    def __init__(self, sites=np.array([]), bsplineDegree=4, adaptivePartition=False):
        self._shortestPath = path.Path(bsplineDegree, adaptivePartition)
        self._sites = sites
        self._graph = nx.Graph()
        self._tGraph = nx.DiGraph()
        self._startTriplet = None
        self._endTriplet = None
        self._polyhedronsContainer = polyhedronsContainer.PolyhedronsContainer()
        self._pathStart = np.array([])
        self._pathEnd = np.array([])
        self._startId = uuid.uuid4()
        self._endId = uuid.uuid4()
        self._bsplineDegree = bsplineDegree

    def setBsplineDegree(self, bsplineDegree):
        self._bsplineDegree = bsplineDegree
        self._shortestPath.setBsplineDegree(bsplineDegree)

    def setAdaptivePartition(self, adaptivePartition):
        self._shortestPath.setAdaptivePartition(adaptivePartition)

    def setCustomSites(self, sites):
        self._sites = sites

    def setRandomSites(self, number, seed=None):
        if seed != None:
            np.random.seed(0)
        self._sites = sp.rand(number,3)

    def addPolyhedron(self, polyhedron):
        self._polyhedronsContainer.addPolyhedron(polyhedron)

    def addBoundingBox(self, a, b, maxEmptyArea=1, invisible=True, verbose=False):
        if verbose:
            print('Add bounding box', flush=True)

        self._polyhedronsContainer.addBoundingBox(a,b,maxEmptyArea, invisible)

    def setPolyhedronsSites(self, verbose=False):
        if verbose:
            print('Set sites for Voronoi', flush=True)

        sites = []
        for polyhedron in self._polyhedronsContainer.polyhedrons:
            sites.extend(polyhedron.allPoints)

        self._sites = np.array(sites)

    def makeVoroGraph(self, prune=True, verbose=False, debug=False):
        if verbose:
            print('Calculate Voronoi cells', flush=True)
        ids = {}
        vor = sp.spatial.Voronoi(self._sites)

        if verbose:
            print('Make pruned Graph from cell edges ', end='', flush=True)
            printDotBunch = 0
        vorVer = vor.vertices
        for ridge in vor.ridge_vertices:
            if verbose:
                if printDotBunch == 0:
                    print('.', end='', flush=True)
                printDotBunch = (printDotBunch+1)%10

            for i in range(1, len(ridge)):
                for j in range(i):
                    if (ridge[i] != -1) and (ridge[j] != -1):
                        a = vorVer[ridge[i]]
                        b = vorVer[ridge[j]]
                        if (not prune) or (not self._polyhedronsContainer.segmentIntersectPolyhedrons(a,b)):
                            if tuple(a) in ids.keys():
                                idA = ids[tuple(a)]
                            else:
                                idA = uuid.uuid4()
                                self._graph.add_node(idA, coord=a)
                                ids[tuple(a)] = idA

                            if tuple(b) in ids.keys():
                                idB = ids[tuple(b)]
                            else:
                                idB = uuid.uuid4()
                                self._graph.add_node(idB, coord=b)
                                ids[tuple(b)] = idB

                            self._graph.add_edge(idA, idB, weight=np.linalg.norm(a-b))

        if verbose:
            print('', flush=True)

        self._createTripleGraph(verbose, debug)

    def calculateShortestPath(self, start, end, attachMode='near', prune=True, useMethod='cleanPath', postSimplify=True, verbose=False, debug=False):
        """
        useMethod: cleanPath; trijkstra; annealing
        """
        if verbose:
            print('Attach start and end points', flush=True)
        if attachMode=='near':
            self._attachToGraphNear(start, end, prune)
        elif attachMode=='all':
            self._attachToGraphAll(start, end, prune)
        else:
            self._attachToGraphNear(start, end, prune)

        self._attachSpecialStartEndTriples(verbose)

        self._pathStart = start
        self._pathEnd = end

        if useMethod == 'trijkstra':
            self._removeCollidingTriples(verbose, debug)

        triPath = self._dijkstra(verbose, debug)
        path = self._extractPath(triPath, verbose)
        self._shortestPath.assignValues(path, self._polyhedronsContainer)

        if useMethod == 'cleanPath':
            self._shortestPath.clean(verbose, debug)
        elif useMethod == 'annealing':
            self._shortestPath.anneal(verbose)
        if postSimplify:
            self._shortestPath.simplify(verbose, debug)

        #print(self._bsplineDegree)
        if useMethod != 'annealing':
            if self._bsplineDegree == 3:
                self._shortestPath.addNAlignedVertexes(1, verbose, debug)        
            if self._bsplineDegree == 4:
                self._shortestPath.addNAlignedVertexes(2, verbose, debug)
            
            

    def plotSites(self, plotter, verbose=False):
        if verbose:
            print('Plot Sites', end='', flush=True)
            
        if self._sites.size > 0:
            plotter.addPoints(self._sites, plotter.COLOR_SITES, thick=True)

    def plotPolyhedrons(self, plotter, verbose=False):
        if verbose:
            print('Plot Polyhedrons', end='', flush=True)
            
        for poly in self._polyhedronsContainer.polyhedrons:
            poly.plot(plotter)
            if verbose:
                print('.', end='', flush=True)

        if verbose:
            print('', flush=True)

    def plotShortestPath(self, plotter, verbose=False):
        if verbose:
            print('Plot shortest path', flush=True)
            
        if self._shortestPath.vertexes.size > 0:
            if self._polyhedronsContainer.hasBoundingBox:
                splineThickness = np.linalg.norm(np.array(self._polyhedronsContainer.boundingBoxB) - np.array(self._polyhedronsContainer.boundingBoxA)) / 1000.
                pointThickness = splineThickness * 2.
                lineThickness = splineThickness / 2.
                
                plotter.addPolyLine(self._shortestPath.vertexes, plotter.COLOR_CONTROL_POLIG, thick=True, thickness=lineThickness)
                plotter.addPoints(self._shortestPath.vertexes, plotter.COLOR_CONTROL_POINTS, thick=True, thickness=pointThickness)
                plotter.addBSpline(self._shortestPath, self._bsplineDegree, plotter.COLOR_PATH, thick=True, thickness=splineThickness)

            else:
                plotter.addPolyLine(self._shortestPath.vertexes, plotter.COLOR_CONTROL_POLIG, thick=True)
                plotter.addPoints(self._shortestPath.vertexes, plotter.COLOR_CONTROL_POINTS, thick=True)
                plotter.addBSpline(self._shortestPath, self._bsplineDegree, plotter.COLOR_PATH, thick=True)

    def plotGraph(self, plotter, verbose=False):
        if verbose:
            print('Plot shortest path', flush=True)

        plotter.addGraph(self._graph, plotter.COLOR_GRAPH)

    def extractXmlTree(self, root):
        if self._polyhedronsContainer.hasBoundingBox:
            xmlBoundingBox = ET.SubElement(root, 'boundingBox')
            ET.SubElement(xmlBoundingBox, 'a', x=str(self._polyhedronsContainer.boundingBoxA[0]), y=str(self._polyhedronsContainer.boundingBoxA[1]), z=str(self._polyhedronsContainer.boundingBoxA[2]))
            ET.SubElement(xmlBoundingBox, 'b', x=str(self._polyhedronsContainer.boundingBoxB[0]), y=str(self._polyhedronsContainer.boundingBoxB[1]), z=str(self._polyhedronsContainer.boundingBoxB[2]))

        xmlPolyhedrons = ET.SubElement(root, 'polyhedrons')
        for polyhedron in self._polyhedronsContainer.polyhedrons:
            xmlPolyhedron = polyhedron.extractXmlTree(xmlPolyhedrons)

    def importXmlTree(self, root, maxEmptyArea):
        xmlBoundingBox = root.find('boundingBox')
        if xmlBoundingBox:
            xmlA = xmlBoundingBox.find('a')
            xmlB = xmlBoundingBox.find('b')
            
            self._polyhedronsContainer.hasBoundingBox = True
            self._polyhedronsContainer.boundingBoxA = [float(xmlA.attrib['x']), float(xmlA.attrib['y']), float(xmlA.attrib['z'])]
            self._polyhedronsContainer.boundingBoxB = [float(xmlB.attrib['x']), float(xmlB.attrib['y']), float(xmlB.attrib['z'])]

        xmlPolyhedrons = root.find('polyhedrons')
        if xmlPolyhedrons:
            for xmlPolyhedron in xmlPolyhedrons.iter('polyhedron'):
                invisible = False
                if 'invisible' in xmlPolyhedron.attrib.keys():
                    invisible = bool(eval(xmlPolyhedron.attrib['invisible']))

                faces = []
                for xmlFace in xmlPolyhedron.iter('face'):
                    vertexes = []
                    for xmlVertex in xmlFace.iter('vertex'):
                        vertexes.append([float(xmlVertex.attrib['x']), float(xmlVertex.attrib['y']), float(xmlVertex.attrib['z'])])

                    faces.append(vertexes)

                newPolyhedron = polyhedron.Polyhedron(faces=np.array(faces), invisible=invisible, maxEmptyArea=maxEmptyArea)
                self_polyhedrons.addPolyhedron(newPolyhedron)

    def _attachToGraphNear(self, start, end, prune):
        firstS = True
        firstE = True
        minAttachS = None
        minAttachE = None
        minDistS = 0.
        minDistE = 0.
        for node,nodeAttr in self._graph.node.items():
            if (not prune) or (not self._polyhedronsContainer.segmentIntersectPolyhedrons(start,nodeAttr['coord'])):
                if firstS:
                    minAttachS = node
                    minDistS = np.linalg.norm(start - nodeAttr['coord'])
                    firstS = False
                else:
                    currDist = np.linalg.norm(start - nodeAttr['coord'])
                    if currDist < minDistS:
                        minAttachS = node
                        minDistS = currDist

            if (not prune) or (not self._polyhedronsContainer.segmentIntersectPolyhedrons(end, nodeAttr['coord'])):
                if firstE:
                    minAttachE = node
                    minDistE = np.linalg.norm(end - nodeAttr['coord'])
                    firstE = False
                else:
                    currDist = np.linalg.norm(end - nodeAttr['coord'])
                    if currDist < minDistE:
                        minAttachE = node
                        minDistE = currDist

        if minAttachS != None:
            self._addNodeToTGraph(self._startId, start, minAttachS, minDistS, rightDirection=True)
        if minAttachE != None:
            self._addNodeToTGraph(self._endId, end, minAttachE, minDistE, rightDirection=False)

    def _attachToGraphAll(self, start, end, prune):
        for node,nodeAttr in self._graph.node.items():
            if (not prune) or (not self._polyhedronsContainer.segmentIntersectPolyhedrons(start, nodeAttr['coord'])):
                self._addNodeToTGraph(self._startId, start, node, np.linalg.norm(start - nodeAttr['coord']), rightDirection=True)
            if (not prune) or (not self._polyhedronsContainer.segmentIntersectPolyhedrons(end, nodeAttr['coord'])):
                self._addNodeToTGraph(self._endId, end, node, np.linalg.norm(end - nodeAttr['coord']), rightDirection=False)

    def _addNodeToTGraph(self, newId, coord, attachId, dist, rightDirection):
        self._graph.add_node(newId, coord=coord)
        self._graph.add_edge(newId, attachId, weight=dist)
        for otherId in filter(lambda node: node != newId, self._graph.neighbors(attachId)):
            newTriplet = uuid.uuid4()
            if rightDirection:
                self._tGraph.add_node(newTriplet, triplet=[newId,attachId,otherId])
                self._tGraph.add_edges_from([(newTriplet, otherTriplet, {'weight':dist}) for otherTriplet in self._tGraph.nodes() if self._tGraph.node[otherTriplet]['triplet'][0] == attachId  and self._tGraph.node[otherTriplet]['triplet'][1] == otherId])

            else:
                self._tGraph.add_node(newTriplet, triplet=[otherId,attachId,newId])
                self._tGraph.add_edges_from([(otherTriplet, newTriplet, {'weight':dist}) for otherTriplet in self._tGraph.nodes() if self._tGraph.node[otherTriplet]['triplet'][1] == otherId  and self._tGraph.node[otherTriplet]['triplet'][2] == attachId])

    def _attachSpecialStartEndTriples(self, verbose):
        #attach special starting and ending triplet
        if verbose:
            print('Create starting and ending triplets', flush=True)

        self._startTriplet = uuid.uuid4()
        self._endTriplet = uuid.uuid4()
        self._tGraph.add_node(self._startTriplet, triplet = [self._startId,self._startId,self._startId], hit = False)
        self._tGraph.add_node(self._endTriplet, triplet = [self._endId,self._endId,self._endId], hit = False)
        self._tGraph.add_edges_from([(self._startTriplet, n, {'weight':0.}) for n in self._tGraph.nodes() if self._tGraph.node[n]['triplet'][0] == self._startId])
        self._tGraph.add_edges_from([(n, self._endTriplet, {'weight':0.}) for n in self._tGraph.nodes() if self._tGraph.node[n]['triplet'][2] == self._endId])

    def _createTripleGraph(self, verbose, debug):
        #create triplets

        if debug:
            triplets_file = open("triplets.txt","w")

        if verbose:
            print('Create triplets ', end='', flush=True)
            printDotBunch = 0

        tripletIdList = {}
        def getUniqueId(triplet):
            if tuple(triplet) in tripletIdList.keys():
                tripletId = tripletIdList[tuple(triplet)]
            else:
                tripletId = uuid.uuid4()
                tripletIdList[tuple(triplet)] = tripletId
                self._tGraph.add_node(tripletId, triplet = triplet)
            return tripletId

        for edge in self._graph.edges():
            if verbose:
                if printDotBunch == 0:
                    print('.', end='', flush=True)
                printDotBunch = (printDotBunch+1)%10


            tripletsSxOutgoing = []
            tripletsSxIngoing = []
            tripletsDxOutgoing = []
            tripletsDxIngoing = []

            for nodeSx in filter(lambda node: node != edge[1], self._graph.neighbors(edge[0])):
                tripletId = getUniqueId([nodeSx,edge[0],edge[1]])
                tripletsSxOutgoing.append(tripletId)
                if debug:
                    triplets_file.write('SxO: {}\n'.format(self._tGraph.node[tripletId]['triplet']))

                tripletId = getUniqueId([edge[1],edge[0],nodeSx])
                tripletsSxIngoing.append(tripletId)
                if debug:
                    triplets_file.write('SxI: {}\n'.format(self._tGraph.node[tripletId]['triplet']))

            for nodeDx in filter(lambda node: node != edge[0], self._graph.neighbors(edge[1])):
                tripletId = getUniqueId([nodeDx,edge[1],edge[0]])
                tripletsDxOutgoing.append(tripletId)
                if debug:
                    triplets_file.write('DxO: {}\n'.format(self._tGraph.node[tripletId]['triplet']))

                tripletId = getUniqueId([edge[0],edge[1],nodeDx])
                tripletsDxIngoing.append(tripletId)
                if debug:
                    triplets_file.write('DxI: {}\n'.format(self._tGraph.node[tripletId]['triplet']))

            for tripletSx in tripletsSxOutgoing:
                for tripletDx in tripletsDxIngoing:
                    self._tGraph.add_edge(tripletSx, tripletDx, {'weight':self._graph.edge[self._tGraph.node[tripletSx]['triplet'][0]][self._tGraph.node[tripletDx]['triplet'][0]]['weight']})

            for tripletDx in tripletsDxOutgoing:
                for tripletSx in tripletsSxIngoing:
                    self._tGraph.add_edge(tripletDx, tripletSx, {'weight':self._graph.edge[self._tGraph.node[tripletDx]['triplet'][0]][self._tGraph.node[tripletSx]['triplet'][0]]['weight']})

        if verbose:
            print('', flush=True)

        if debug:
            triplets_file.close()


    def _dijkstra(self, verbose, debug):
        try:
            if verbose:
                print('Dijkstra algorithm', flush=True)

            length,triPath=nx.bidirectional_dijkstra(self._tGraph, self._startTriplet, self._endTriplet)


        except (nx.NetworkXNoPath, nx.NetworkXError):
            print('ERROR: Impossible to find a path')
            triPath = []

        return triPath

    def _extractPath(self, triPath, verbose):
        if verbose:
            print('Extract path', flush=True)

        path = []
        for t in triPath:
            path.append(self._graph.node[self._tGraph.node[t]['triplet'][1]]['coord'])
        return np.array(path)


    def _removeCollidingTriples(self, verbose, debug):
        if verbose:
            print('Remove colliding triples', flush=True)
            printDotBunch = 0

        toRemove = []
        for triple in self._tGraph:
            if verbose:
                if printDotBunch == 0:
                    print('.', end='', flush=True)
                printDotBunch = (printDotBunch+1)%10

            a = self._graph.node[self._tGraph.node[triple]['triplet'][0]]['coord']
            b = self._graph.node[self._tGraph.node[triple]['triplet'][1]]['coord']
            c = self._graph.node[self._tGraph.node[triple]['triplet'][2]]['coord']
            intersect,intersectRes = self._polyhedronsContainer.triangleIntersectPolyhedrons(a, b, c)
            if intersect:
                toRemove.append(triple)
                
        if verbose:
            print ("", flush=True)

        for triple in toRemove:
            self._tGraph.remove_node(triple)
    
