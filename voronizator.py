import numpy as np
import numpy.linalg
import scipy as sp
import scipy.spatial
import networkx as nx
import numpy.linalg
import polyhedron
import uuid
import xml.etree.cElementTree as ET

class Voronizator:
    def __init__(self, sites=np.array([]), bsplineDegree=4):
        self._sites = sites
        self._shortestPath = np.array([])
        self._graph = nx.Graph()
        self._tGraph = nx.DiGraph()
        self._startTriplet = None
        self._endTriplet = None
        self._polyhedrons = []
        self._pathStart = np.array([])
        self._pathEnd = np.array([])
        self._startId = uuid.uuid4()
        self._endId = uuid.uuid4()
        self._bsplineDegree = bsplineDegree
        self._hasBoundingBox = False

    def setBsplineDegree(self, bsplineDegree):
        self._bsplineDegree = bsplineDegree

    def setCustomSites(self, sites):
        self._sites = sites

    def setRandomSites(self, number, seed=None):
        if seed != None:
            np.random.seed(0)
        self._sites = sp.rand(number,3)

    def addPolyhedron(self, polyhedron):
        self._polyhedrons.append(polyhedron)

    def addBoundingBox(self, a, b, maxEmptyArea=1, invisible=True, verbose=False):
        if verbose:
            print('Add bounding box', flush=True)

        self._hasBoundingBox = True
        self._boundingBoxA = a
        self._boundingBoxB = b

        c = [a[0], b[1], a[2]]
        d = [b[0], a[1], a[2]]
        e = [a[0], a[1], b[2]]
        f = [b[0], b[1], a[2]]
        g = [b[0], a[1], b[2]]
        h = [a[0], b[1], b[2]]

        self._polyhedrons.append(polyhedron.Polyhedron(faces=np.array([
            [a,g,e],[a,d,g],[d,f,g],[f,b,g],[f,b,h],[f,h,c],
            [h,a,e],[h,c,a],[e,h,g],[h,b,g],[a,d,f],[a,f,c]
            ]), invisible=invisible, maxEmptyArea=maxEmptyArea))


    def setPolyhedronsSites(self, verbose=False):
        if verbose:
            print('Set sites for Voronoi', flush=True)

        sites = []
        for polyhedron in self._polyhedrons:
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
                        if (not prune) or (not self._segmentIntersectPolyhedrons(a,b)):
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

    def calculateShortestPath(self, start, end, attachMode='near', prune=True, useTrijkstra=False, postSimplify=True, verbose=False, debug=False):
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

        if useTrijkstra:
            self._removeCollidingTriples(verbose, debug)

        triPath = self._dijkstra(verbose, debug)
        shortestPath = self._extractPath(triPath, verbose, debug)

        if not useTrijkstra:
            shortestPath = self._cleanPath(shortestPath, verbose, debug)
        if postSimplify:
            shortestPath = self._simplifyPath(shortestPath, verbose, debug)

        #print(self._bsplineDegree)
        if self._bsplineDegree == 3:
            shortestPath = self._addNAlignedVertexes(1, shortestPath, verbose, debug)        
        if self._bsplineDegree == 4:
            shortestPath = self._addNAlignedVertexes(2, shortestPath, verbose, debug)
            
        self._shortestPath = shortestPath
            

    def plotSites(self, plotter, verbose=False):
        if verbose:
            print('Plot Sites', end='', flush=True)
            
        if self._sites.size > 0:
            plotter.addPoints(self._sites, plotter.COLOR_SITES, thick=True)

    def plotPolyhedrons(self, plotter, verbose=False):
        if verbose:
            print('Plot Polyhedrons', end='', flush=True)
            
        for poly in self._polyhedrons:
            poly.plot(plotter)
            if verbose:
                print('.', end='', flush=True)

        if verbose:
            print('', flush=True)

    def plotShortestPath(self, plotter, adaptivePartition=False, verbose=False):
        if verbose:
            print('Plot shortest path', flush=True)
            
        if self._shortestPath.size > 0:
            if self._hasBoundingBox:
                splineThickness = np.linalg.norm(np.array(self._boundingBoxB) - np.array(self._boundingBoxA)) / 1000.
                pointThickness = splineThickness * 2.
                lineThickness = splineThickness / 2.
                
                plotter.addPolyLine(self._shortestPath, plotter.COLOR_CONTROL_POLIG, thick=True, thickness=lineThickness)
                plotter.addPoints(self._shortestPath, plotter.COLOR_CONTROL_POINTS, thick=True, thickness=pointThickness)
                plotter.addBSpline(self._shortestPath, self._bsplineDegree, adaptivePartition, plotter.COLOR_PATH, thick=True, thickness=splineThickness)

            else:
                plotter.addPolyLine(self._shortestPath, plotter.COLOR_CONTROL_POLIG, thick=True)
                plotter.addPoints(self._shortestPath, plotter.COLOR_CONTROL_POINTS, thick=True)
                plotter.addBSpline(self._shortestPath, self._bsplineDegree, adaptivePartition, plotter.COLOR_PATH, thick=True)

    def plotGraph(self, plotter, verbose=False):
        if verbose:
            print('Plot shortest path', flush=True)

        plotter.addGraph(self._graph, plotter.COLOR_GRAPH)

    def extractXmlTree(self, root):
        if self._hasBoundingBox:
            xmlBoundingBox = ET.SubElement(root, 'boundingBox')
            ET.SubElement(xmlBoundingBox, 'a', x=str(self._boundingBoxA[0]), y=str(self._boundingBoxA[1]), z=str(self._boundingBoxA[2]))
            ET.SubElement(xmlBoundingBox, 'b', x=str(self._boundingBoxB[0]), y=str(self._boundingBoxB[1]), z=str(self._boundingBoxB[2]))

        xmlPolyhedrons = ET.SubElement(root, 'polyhedrons')
        for polyhedron in self._polyhedrons:
            xmlPolyhedron = polyhedron.extractXmlTree(xmlPolyhedrons)

    def importXmlTree(self, root, maxEmptyArea):
        xmlBoundingBox = root.find('boundingBox')
        if xmlBoundingBox:
            xmlA = xmlBoundingBox.find('a')
            xmlB = xmlBoundingBox.find('b')
            
            self._hasBoundingBox = True
            self._boundingBoxA = [float(xmlA.attrib['x']), float(xmlA.attrib['y']), float(xmlA.attrib['z'])]
            self._boundingBoxB = [float(xmlB.attrib['x']), float(xmlB.attrib['y']), float(xmlB.attrib['z'])]

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
                self.addPolyhedron(newPolyhedron)

    def _segmentIntersectPolyhedrons(self, a, b):
        intersect = False
        if self._hasBoundingBox:
            if((a<self._boundingBoxA).any() or (a>self._boundingBoxB).any() or (b<self._boundingBoxA).any() or (b>self._boundingBoxB).any()):
                intersect = True

        if not intersect:
            for polyhedron in self._polyhedrons:
                if polyhedron.intersectSegment(a,b)[0]:
                    intersect = True
                    break
                    
        return intersect

    def _triangleIntersectPolyhedrons(self, a, b, c):
        triangle = polyhedron.Polyhedron(faces=np.array([[a,b,c]]), distributePoints = False)
        intersect = False
        result = np.array([])
        for currPolyhedron in self._polyhedrons:
            currIntersect,currResult = currPolyhedron.intersectPathTriple(triangle)
            if currIntersect and (not intersect or (currResult[1] > result[1])):
                intersect = True
                result = currResult

        return (intersect, result)

    def _attachToGraphNear(self, start, end, prune):
        firstS = True
        firstE = True
        minAttachS = None
        minAttachE = None
        minDistS = 0.
        minDistE = 0.
        for node,nodeAttr in self._graph.node.items():
            if (not prune) or (not self._segmentIntersectPolyhedrons(start,nodeAttr['coord'])):
                if firstS:
                    minAttachS = node
                    minDistS = np.linalg.norm(start - nodeAttr['coord'])
                    firstS = False
                else:
                    currDist = np.linalg.norm(start - nodeAttr['coord'])
                    if currDist < minDistS:
                        minAttachS = node
                        minDistS = currDist

            if (not prune) or (not self._segmentIntersectPolyhedrons(end, nodeAttr['coord'])):
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
            if (not prune) or (not self._segmentIntersectPolyhedrons(start, nodeAttr['coord'])):
                self._addNodeToTGraph(self._startId, start, node, np.linalg.norm(start - nodeAttr['coord']), rightDirection=True)
            if (not prune) or (not self._segmentIntersectPolyhedrons(end, nodeAttr['coord'])):
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
            intersect,intersectRes = self._triangleIntersectPolyhedrons(a, b, c)
            if intersect:
                toRemove.append(triple)
                
        if verbose:
            print ("", flush=True)

        for triple in toRemove:
            self._tGraph.remove_node(triple)
    
    def _extractPath(self, triPath, verbose, debug):
        if verbose:
            print('Extract path', flush=True)

        path = []
        for t in triPath:
            path.append(self._graph.node[self._tGraph.node[t]['triplet'][1]]['coord'])
        return np.array(path)

    def _cleanPath(self, path, verbose, debug):
        if verbose:
            print('Clean path (avoid obstacles)', flush=True)

        newPath = []
        if len(path) > 0:
            a = path[0]
            newPath.append(path[0])

        for i in range(1, len(path)-1):
            v = path[i]
            b = path[i+1]

            intersect,intersectRes = self._triangleIntersectPolyhedrons(a, v, b)
            if intersect:
                alpha = intersectRes[1]

                a1 = (1.-alpha)*a + alpha*v
                b1 = alpha*v + (1.-alpha)*b

                newPath.append(a1)
                newPath.append(v)
                newPath.append(b1)

                a = b1
            else:
                newPath.append(v)

                a = v

        if len(path) > 0:
            newPath.append(path[len(path)-1])

        return np.array(newPath)


    def _simplifyPath(self, path, verbose, debug):
        if verbose:
            print('Simplify path (remove useless triples)', flush=True)

        simplifiedPath = []
        if len(path) > 0:
            a = path[0]
            simplifiedPath.append(path[0])
        first = True
        for i in range(1,len(path)-1):
            v = path[i]
            b = path[i+1]
            keepV = False

            intersectCurr,nihil = self._triangleIntersectPolyhedrons(a, v, b)

            if not intersectCurr:
                if first:
                    intersectPrec = False
                else:
                    a1 = path[i-2]
                    intersectPrec,nihil = self._triangleIntersectPolyhedrons(a1, a, b)

                if i == len(path)-2:
                    intersectSucc = False
                else:
                    b1 = path[i+2]
                    intersectSucc,nihil = self._triangleIntersectPolyhedrons(a, b, b1)

                if intersectPrec or intersectSucc:
                    keepV = True

            else:
                keepV = True

            if keepV:
                first = False
                simplifiedPath.append(v)
                a = v

        if len(path) > 0:
            simplifiedPath.append(path[len(path)-1])

        return np.array(simplifiedPath)
            
    def _addNAlignedVertexes(self, numVertexes, path, verbose, debug):
        if verbose:
            print('Increase degree', flush=True)

        newPath = []
        for i in range(1, len(path)):
            a = path[i-1]
            b = path[i]
            newPath.append(a)

            if numVertexes == 1:
                n = 0.5 * a + 0.5 * b
                newPath.append(n)
                
            elif numVertexes == 2:
                n1 = 0.33 * a + 0.67 * b
                n2 = 0.33 * b + 0.67 * a
                newPath.append(n1)
                newPath.append(n2)
                
        if len(path) > 0:
            newPath.append(path[len(path)-1])

        return np.array(newPath)
    
