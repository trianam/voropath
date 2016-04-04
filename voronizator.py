import numpy as np
import numpy.linalg
import scipy as sp
import scipy.spatial
import scipy.interpolate
import networkx as nx
import numpy.linalg
import polyhedron
import uuid

class Voronizator:
    def __init__(self, sites=np.array([])):
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

    def setCustomSites(self, sites):
        self._sites = sites

    def setRandomSites(self, number, seed=None):
        if seed != None:
            np.random.seed(0)
        self._sites = sp.rand(number,3)

    def addPolyhedron(self, polyhedron):
        self._polyhedrons.append(polyhedron)

    def addBoundingBox(self, a, b, maxEmptyArea=1, invisible=True):
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

    def setPolyhedronsSites(self):
        sites = []
        for polyhedron in self._polyhedrons:
            sites.extend(polyhedron.allPoints)

        self._sites = np.array(sites)

    def makeVoroGraph(self, prune=True):
        ids = {}
        vor = sp.spatial.Voronoi(self._sites)
        vorVer = vor.vertices
        for ridge in vor.ridge_vertices:
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

    def calculateShortestPath(self, start, end, attachMode='near', prune=True, verbose=False, debug=False):
        if verbose:
            print('Attach start and end points', flush=True)
        if attachMode=='near':
            self._attachToGraphNear(start, end, prune)
        elif attachMode=='all':
            self._attachToGraphAll(start, end, prune)
        else:
            self._attachToGraphNear(start, end, prune)

        self._pathStart = start
        self._pathEnd = end

        self._createTripleGraph(verbose, debug)
        triPath = self._trijkstra(verbose, debug)
        self._shortestPath = self._extractPath(triPath, verbose,debug)

    def plotSites(self, plotter):
        if self._sites.size > 0:
            plotter.plot(self._sites[:,0], self._sites[:,1], self._sites[:,2], 'o')

    def plotPolyhedrons(self, plotter):
        for poly in self._polyhedrons:
            poly.plot(plotter)

    def plotShortestPath(self, plotter):
        if self._shortestPath.size > 0:
            x = self._shortestPath[:,0]
            y = self._shortestPath[:,1]
            z = self._shortestPath[:,2]

            t = range(len(self._shortestPath))
            ipl_t = np.linspace(0.0, len(self._shortestPath) - 1, 1000)
            #TODO: find a better way to substitute 100 above
            x_tup = sp.interpolate.splrep(t, x, k=5)
            y_tup = sp.interpolate.splrep(t, y, k=5)
            z_tup = sp.interpolate.splrep(t, z, k=5)

            x_list = list(x_tup)
            xl = x.tolist()
            x_list[1] = xl + [0.0, 0.0, 0.0, 0.0]

            y_list = list(y_tup)
            yl = y.tolist()
            y_list[1] = yl + [0.0, 0.0, 0.0, 0.0]

            z_list = list(z_tup)
            zl = z.tolist()
            z_list[1] = zl + [0.0, 0.0, 0.0, 0.0]

            x_i = sp.interpolate.splev(ipl_t, x_list)
            y_i = sp.interpolate.splev(ipl_t, y_list)
            z_i = sp.interpolate.splev(ipl_t, z_list)

            plotter.plot(x, y, z, 'r--')
            plotter.plot(x_i, y_i, z_i, 'r', lw=2)
            plotter.plot(x, y, z, 'ro')
        # if self._pathStart.size > 0:
        #     plotter.plot([self._pathStart[0]], [self._pathStart[1]], [self._pathStart[2]], 'ro')
        # if self._pathEnd.size > 0:
        #     plotter.plot([self._pathEnd[0]], [self._pathEnd[1]], [self._pathEnd[2]], 'ro')

    def plotGraph(self, plotter, vertexes=True, edges=True, labels=False, pathExtremes=False, showOnly=[]):
        if vertexes:
            for ver in self._graph.nodes():
                if not showOnly or ver in showOnly:
                    if (ver != self._startId and ver != self._endId):
                        plotter.plot([self._graph.node[ver]['coord'][0]], [self._graph.node[ver]['coord'][1]], [self._graph.node[ver]['coord'][2]], 'og')
                        if labels:
                                plotter.text(self._graph.node[ver]['coord'][0], self._graph.node[ver]['coord'][1], self._graph.node[ver]['coord'][2], ver, color='red')
                    elif pathExtremes==True:
                        plotter.plot([self._graph.node[ver]['coord'][0]], [self._graph.node[ver]['coord'][1]], [self._graph.node[ver]['coord'][2]], 'or')
                        if labels:
                                plotter.text(self._graph.node[ver]['coord'][0], self._graph.node[ver]['coord'][1], self._graph.node[ver]['coord'][2], ver, color='red')

        if edges:
            for edge in self._graph.edges():
                if not showOnly or (edge[0] in showOnly and edge[1] in showOnly):
                    if pathExtremes==True or (edge[0] != self._startId and edge[0] != self._endId and edge[1] != self._startId and edge[1] != self._endId):
                        plotter.plot([self._graph.node[edge[0]]['coord'][0], self._graph.node[edge[1]]['coord'][0]], [self._graph.node[edge[0]]['coord'][1], self._graph.node[edge[1]]['coord'][1]], [self._graph.node[edge[0]]['coord'][2], self._graph.node[edge[1]]['coord'][2]], 'k--')

    def _segmentIntersectPolyhedrons(self, a, b):
        for polyhedron in self._polyhedrons:
            if polyhedron.intersectSegment(a,b)[0]:
                return True
        return False

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
            self._graph.add_node(self._startId, coord=start)
            self._graph.add_edge(self._startId, minAttachS, weight=minDistS)
        if minAttachE != None:
            self._graph.add_node(self._endId, coord=end)
            self._graph.add_edge(self._endId, minAttachE, weight=minDistE)

    def _attachToGraphAll(self, start, end, prune):
        for node,nodeAttr in self._graph.node.items():
            if (not prune) or (not self._segmentIntersectPolyhedrons(start, nodeAttr['coord'])):
                self._graph.add_node(self._startId, coord=start)
                self._graph.add_edge(self._startId, node, weight=np.linalg.norm(start - nodeAttr['coord']))
            if (not prune) or (not self._segmentIntersectPolyhedrons(end, nodeAttr['coord'])):
                self._graph.add_node(self._endId, coord=end)
                self._graph.add_edge(self._endId, node, weight=np.linalg.norm(end - nodeAttr['coord']))

    def _createTripleGraph(self, verbose, debug):
        #create triplets
        if verbose:
            print('Create triplets ', end='', flush=True)

        if debug:
            hits_file = open("hits.txt","w")

        for node0 in self._graph.nodes():
            for node1 in self._graph.neighbors(node0):
                for node2 in filter(lambda node: node!=node0, self._graph.neighbors(node1)):
                    if verbose:
                        print('.', end='', flush=True)

                    tripletId = uuid.uuid4()
                    self._tGraph.add_node(tripletId, triplet = [node0,node1,node2], hit=False)
                    if debug:
                        hits_file.write(str(self._tGraph.node[tripletId]['triplet'])+"\n")


                    # if not self._tGraph.node[tripletId]['triplet'][::-1] in [n for n,attr in self._tGraph.node.items() if attr['hit'] == True]:
                    #     intersect,result = self._triangleIntersectPolyhedrons(self._graph.node[node0]['coord'], self._graph.node[node1]['coord'], self._graph.node[node2]['coord'])
                    #     if not intersect:
                    #         self._tGraph.node[tripletId]['hit'] = False
                    #         #Q.add(triplet, d)
                    #     else:
                    #         if debug:
                    #             hits_file.write(str(self._tGraph.node[tripletId]['triplet'])+"\n")
                    #         self._tGraph.node[tripletId]['hit'] = True
                    #         self._tGraph.node[tripletId]['hitRes'] = result[1]
                    # else:
                    #     if debug:
                    #         hits_file.write(str(self._tGraph.node[tripletId]['triplet'])+"\n")
                    #     self._tGraph.node[tripletId]['hit'] = True
                    #     self._tGraph.node[tripletId]['hitRes'] = [n for n,attr in self._tGraph.node.items() if attr['triplet'][::-1] == self._tGraph.node[tripletId]['triplet']][0]

        if verbose:
            print('', flush=True)

        if debug:
            hits_file.close()

        #create edges between triples
        if verbose:
            print('Create edges between triples', end='', flush=True)

        for triple in self._tGraph.nodes():
            #self._tGraph.add_edges_from([(triple,n,{'weight':self._graph.edge[self._tGraph.node[triple]['triplet'][0]][self._tGraph.node[n]['triplet'][0]]['weight']}) for n in self._tGraph.nodes() if self._tGraph.node[n]['triplet'][0] == self._tGraph.node[triple]['triplet'][1] and self._tGraph.node[n]['triplet'][1] == self._tGraph.node[triple]['triplet'][2]])
            for n in self._tGraph.nodes():
                if self._tGraph.node[n]['triplet'][0] == self._tGraph.node[triple]['triplet'][1] and self._tGraph.node[n]['triplet'][1] == self._tGraph.node[triple]['triplet'][2]:
                    if verbose:
                        print('.', end='', flush=True)
                    self._tGraph.add_edge(triple, n, {'weight':self._graph.edge[self._tGraph.node[triple]['triplet'][0]][self._tGraph.node[n]['triplet'][0]]['weight']})

        #attach special starting and ending triplet
        if verbose:
            print('', flush=True)
            print('Create starting and ending triplets', flush=True)

        self._startTriplet = uuid.uuid4()
        self._endTriplet = uuid.uuid4()
        self._tGraph.add_node(self._startTriplet, triplet = [self._startId,self._startId,self._startId], hit = False)
        self._tGraph.add_node(self._endTriplet, triplet = [self._endId,self._endId,self._endId], hit = False)
        self._tGraph.add_edges_from([(self._startTriplet, n, {'weight':0.}) for n in self._tGraph.nodes() if self._tGraph.node[n]['triplet'][0] == self._startId])
        self._tGraph.add_edges_from([(n, self._endTriplet, {'weight':0.}) for n in self._tGraph.nodes() if self._tGraph.node[n]['triplet'][2] == self._endId])

    def _trijkstra(self, verbose, debug):
        try:
            if verbose:
                print('Dijkstra algorithm', end='', flush=True)

            length,triPath=nx.bidirectional_dijkstra(self._tGraph, self._startTriplet, self._endTriplet)


        except (nx.NetworkXNoPath, nx.NetworkXError):
            triPath = []

        return triPath

    def _extractPath(self, triPath, verbose, debug):
        if verbose:
            print('', flush=True)
            print('Adjust hits and construct path', flush=True)

        path = []
        for i,t in [(i,t) for i,t in enumerate(triPath)]:
            if(t == self._startTriplet):
                path.append(self._graph.node[self._startId]['coord'])
            elif(t == self._endTriplet):
                path.append(self._graph.node[self._endId]['coord'])
            else:
                a = self._tGraph.node[t]['triplet'][0]
                v = self._tGraph.node[t]['triplet'][1]
                b = self._tGraph.node[t]['triplet'][2]
                intersect,alpha = self._triangleIntersectPolyhedrons(self._graph.node[a]['coord'], self._graph.node[v]['coord'], self._graph.node[b]['coord'])
                if intersect:
                    #adjust graph
                    #a1,b1 for avoiding obstacle
                    #a2,b2 for augmenting grade
                    a1 = uuid.uuid4()
                    b1 = uuid.uuid4()
                    a2 = uuid.uuid4()
                    b2 = uuid.uuid4()

                    self._graph.add_node(a1, coord = (1.-alpha)*self._graph.node[a]['coord'] + alpha*self._graph.node[v]['coord'])
                    self._graph.add_node(b1, coord = alpha*self._graph.node[v]['coord'] + (1.-alpha)*self._graph.node[b]['coord'])
                    self._graph.add_node(a2, coord = 0.25*self._graph.node[a1]['coord'] + 0.75*self._graph.node[v]['coord'])
                    self._graph.add_node(b2, coord = 0.25*self._graph.node[b1]['coord'] + 0.75*self._graph.node[v]['coord'])

                    self._graph.remove_edge(a,v)
                    self._graph.remove_edge(v,b)

                    self._graph.add_edge(a,a1, weight=np.linalg.norm(self._graph.node[a]['coord'] - self._graph.node[a1]['coord']))
                    self._graph.add_edge(a1,a2, weight=np.linalg.norm(self._graph.node[a1]['coord'] - self._graph.node[a2]['coord']))
                    self._graph.add_edge(a2,v, weight=np.linalg.norm(self._graph.node[a2]['coord'] - self._graph.node[v]['coord']))
                    self._graph.add_edge(v,b2, weight=np.linalg.norm(self._graph.node[v]['coord'] - self._graph.node[b2]['coord']))
                    self._graph.add_edge(b2,b1, weight=np.linalg.norm(self._graph.node[b2]['coord'] - self._graph.node[b1]['coord']))
                    self._graph.add_edge(b1,b, weight=np.linalg.norm(self._graph.node[b1]['coord'] - self._graph.node[b]['coord']))

                    #adjust next triple
                    if i < len(triPath)-1:
                        self._tGraph.node[triPath[i+1]]['triplet'][0] = b1

                    path.append(self._graph.node[a1]['coord'])
                    path.append(self._graph.node[a2]['coord'])
                    path.append(self._graph.node[v]['coord'])
                    path.append(self._graph.node[b2]['coord'])
                    path.append(self._graph.node[b1]['coord'])

                else:
                    # #delete triple
                    # if i < len(triPath)-1:
                    #     self._tGraph.node[triPath[i+1]]['triplet'][0] = a
                    #     self._graph.add_edge(a,self._tGraph.node[triPath[i+1]]['triplet'][1], weight=np.linalg.norm(self._graph.node[a]['coord'] - self._graph.node[self._tGraph.node[triPath[i+1]]['triplet'][1]]['coord']))

                    #a1,b1 for augmenting grade
                    a1 = uuid.uuid4()
                    b1 = uuid.uuid4()

                    self._graph.add_node(a1, coord = 0.25*self._graph.node[a]['coord'] + 0.75*self._graph.node[v]['coord'])
                    self._graph.add_node(b1, coord = 0.25*self._graph.node[b]['coord'] + 0.75*self._graph.node[v]['coord'])

                    self._graph.remove_edge(a,v)
                    self._graph.remove_edge(v,b)

                    self._graph.add_edge(a,a1, weight=np.linalg.norm(self._graph.node[a]['coord'] - self._graph.node[a1]['coord']))
                    self._graph.add_edge(a1,v, weight=np.linalg.norm(self._graph.node[a1]['coord'] - self._graph.node[v]['coord']))
                    self._graph.add_edge(v,b1, weight=np.linalg.norm(self._graph.node[v]['coord'] - self._graph.node[b1]['coord']))
                    self._graph.add_edge(b1,b, weight=np.linalg.norm(self._graph.node[b1]['coord'] - self._graph.node[b]['coord']))

                    #adjust next triple
                    if i < len(triPath)-1:
                        self._tGraph.node[triPath[i+1]]['triplet'][0] = b1

                    path.append(self._graph.node[a1]['coord'])
                    path.append(self._graph.node[v]['coord'])
                    path.append(self._graph.node[b1]['coord'])

        #path = [self._graph.node[self._tGraph.node[n]['triplet'][1]]['coord'] for n in triPath]
        return np.array(path)
