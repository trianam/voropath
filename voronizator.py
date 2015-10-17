import numpy as np
import numpy.linalg
import scipy as sp
import scipy.spatial
import networkx as nx
import numpy.linalg
import polyhedron
import smoothener
import priorityQueue

class Voronizator:
    def __init__(self, sites=np.array([])):
        self._smoothener = smoothener.Smoothener()
        self._sites = sites
        self._shortestPath = np.array([])
        self._graph = nx.Graph()
        self._polyhedrons = []
        self._pathStart = np.array([])
        self._pathEnd = np.array([])

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
        vor = sp.spatial.Voronoi(self._sites)
        vorVer = vor.vertices
        for ridge in vor.ridge_vertices:
            for i in range(len(ridge)):
                if (ridge[i] != -1) and (ridge[(i+1)%len(ridge)] != -1):
                    a = vorVer[ridge[i]]
                    b = vorVer[ridge[(i+1)%len(ridge)]]
                    if (not prune) or (not self._segmentIntersectPolyhedrons(a,b)):
                        self._graph.add_edge(tuple(a), tuple(b), weight=np.linalg.norm(a-b))
        i = 0
        for node in self._graph.nodes():
            self._graph.node[node]['index'] = i
            i = i + 1

    def calculateShortestPath(self, start, end, attachMode='near', prune=True):
        if attachMode=='near':
            self._attachToGraphNear(start, end, prune)
        elif attachMode=='all':
            self._attachToGraphAll(start, end, prune)
        else:
            self._attachToGraphNear(start, end, prune)
            
        self._pathStart = start
        self._pathEnd = end

        if tuple(start) in self._graph.nodes():
            self._graph.node[tuple(start)]['index'] = 's'
        if tuple(end) in self._graph.nodes():
            self._graph.node[tuple(end)]['index'] = 'e'

        self._shortestPath = self._trijkstra(start, end)
        
        # try:
        #     length,path=nx.bidirectional_dijkstra(self._graph, tuple(start), tuple(end))
        # except (nx.NetworkXNoPath, nx.NetworkXError):
        #     path = []
        #self._shortestPath = np.array(path)

    def plotSites(self, plotter):
        if self._sites.size > 0:
            plotter.plot(self._sites[:,0], self._sites[:,1], self._sites[:,2], 'o')

    def plotPolyhedrons(self, plotter):
        for poly in self._polyhedrons:
            poly.plot(plotter)
            
    def plotShortestPath(self, plotter):
        if self._shortestPath.size > 0:
            self._smoothener.plot(self._shortestPath, plotter)
        if self._pathStart.size > 0:
            plotter.plot([self._pathStart[0]], [self._pathStart[1]], [self._pathStart[2]], 'ro')
        if self._pathEnd.size > 0:
            plotter.plot([self._pathEnd[0]], [self._pathEnd[1]], [self._pathEnd[2]], 'ro')

    def plotGraph(self, plotter, vertexes=True, edges=True, labels=False, pathExtremes=False, showOnly=[]):
        if vertexes:
            for ver in self._graph.nodes():
                if not showOnly or self._graph.node[ver]['index'] in showOnly:
                    if pathExtremes==True or (ver!=tuple(self._pathStart) and ver!=tuple(self._pathEnd)):
                        plotter.plot([ver[0]], [ver[1]], [ver[2]], 'og')
                        if labels and ('index' in self._graph.node[ver]):
                                plotter.text(ver[0], ver[1], ver[2], self._graph.node[ver]['index'], color='red')

        if edges:
            for edge in self._graph.edges():
                if not showOnly or (self._graph.node[edge[0]]['index'] in showOnly and self._graph.node[edge[1]]['index'] in showOnly):
                    if pathExtremes==True or (edge[0]!=tuple(self._pathStart) and edge[0]!=tuple(self._pathEnd) and edge[1]!=tuple(self._pathStart) and edge[1]!=tuple(self._pathEnd)):
                        plotter.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], [edge[0][2], edge[1][2]], 'k--')

    def _segmentIntersectPolyhedrons(self, a, b):
        for polyhedron in self._polyhedrons:
            if polyhedron.intersectSegment(a,b):
                return True
        return False
                    
    def _triangleIntersectPolyhedrons(self, a, b, c):
        triangle = polyhedron.Polyhedron(faces=np.array([[a,b,c]]), distributePoints = False)
        for currPolyhedron in self._polyhedrons:
            if currPolyhedron.intersectPolyhedron(triangle):
                return True
        return False
                    
    def _attachToGraphNear(self, start, end, prune):
        firstS = True
        firstE = True
        minAttachS = None
        minAttachE = None
        minDistS = 0.
        minDistE = 0.
        for node in self._graph.nodes():
            if (not prune) or (not self._segmentIntersectPolyhedrons(start,np.array(node))):
                if firstS:
                    minAttachS = node
                    minDistS = np.linalg.norm(start-np.array(node))
                    firstS = False
                else:
                    currDist = np.linalg.norm(start-np.array(node))
                    if currDist < minDistS:
                        minAttachS = node
                        minDistS = currDist
                    
            if (not prune) or (not self._segmentIntersectPolyhedrons(end,np.array(node))):
                if firstE:
                    minAttachE = node
                    minDistE = np.linalg.norm(end-np.array(node))
                    firstE = False
                else:
                    currDist = np.linalg.norm(end-np.array(node))
                    if currDist < minDistE:
                        minAttachE = node
                        minDistE = currDist

        if minAttachS != None:
            self._graph.add_edge(tuple(start), minAttachS, weight=minDistS)
        if minAttachE != None:
            self._graph.add_edge(tuple(end), minAttachE, weight=minDistE)

    def _attachToGraphAll(self, start, end, prune):
        for node in self._graph.nodes():
            if (not prune) or (not self._segmentIntersectPolyhedrons(start,np.array(node))):
                self._graph.add_edge(tuple(start), node, weight=np.linalg.norm(start-np.array(node)))
            if (not prune) or (not self._segmentIntersectPolyhedrons(end,np.array(node))):
                self._graph.add_edge(tuple(end), node, weight=np.linalg.norm(end-np.array(node)))

    def _trijkstra(self, startA, endA):
        start = tuple(startA)
        end = tuple(endA)
        endTriplet = (end,end,end) #special triplet for termination
        inf = float("inf")
        path = []
        Q = priorityQueue.PQueue()
        dist = {}
        prev = {}
        hits = []

        #create triplets
        for node0 in self._graph.nodes():
            for node1 in self._graph.neighbors(node0):
                for node2 in filter(lambda node: node!=node0, self._graph.neighbors(node1)):
                    triplet = (node0,node1,node2)
                    if not triplet[::-1] in hits:
                        if not self._triangleIntersectPolyhedrons(np.array(node0), np.array(node1), np.array(node2)):
                            if node0 != start:
                                dist[triplet] = inf
                                Q.add(triplet, inf)
                            else:
                                dist[triplet] = 0
                                Q.add(triplet, 0)
                        else:
                            hits[:0] = [triplet]

        dist[endTriplet] = inf
        Q.add(endTriplet, inf)

        try:
            while True:
                u = Q.pop()

                if u == endTriplet or dist[u] == inf:
                    break

                for v in Q.filterGet(lambda tri: u[1] == tri[0] and u[2] == tri[1]):
                    tmpDist = dist[u] + self._graph[u[0]][u[1]]['weight']
                    if tmpDist < dist[v]:
                        dist[v] = tmpDist
                        prev[v] = u
                        Q.add(v, tmpDist)

                if u[2] == end:
                    tmpDist = dist[u] + self._graph[u[0]][u[1]]['weight'] + self._graph[u[1]][u[2]]['weight']
                    if tmpDist < dist[endTriplet]:
                        dist[endTriplet] = tmpDist
                        prev[endTriplet] = u
                        Q.add(v, tmpDist)
        except KeyError:
            pass
                        

        u = endTriplet
        while u in prev:
            u = prev[u]
            path[:0] = [u[1]]
            
        if path:
            path[len(path):] = [end]
            path[:0] = [start]

        return np.array(path)

    
