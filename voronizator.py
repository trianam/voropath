import numpy as np
import numpy.linalg
import scipy as sp
import scipy.spatial
import networkx as nx
import numpy.linalg
import polyhedron

class Voronizator:
    def __init__(self, sites=np.array([])):
        self._sites = sites
        self._shortestPath = np.array([])
        self._graph = nx.Graph()
        self._polyhedrons = []

    def setCustomSites(self, sites):
        self._sites = sites
        
    def setRandomSites(self, number, seed=None):
        if seed != None:
            np.random.seed(0)
        self._sites = sp.rand(number,3)

    def addPolyhedron(self, polyhedron):
        self._polyhedrons.append(polyhedron)

    def addBoundingBox(self, a, b):
        c = [a[0], b[1], a[2]]
        d = [b[0], a[1], a[2]]
        e = [a[0], a[1], b[2]]
        f = [b[0], b[1], a[2]]
        g = [b[0], a[1], b[2]]
        h = [a[0], b[1], b[2]]

        self._polyhedrons.append(polyhedron.Polyhedron(faces=np.array([
            [a,g,e],[a,d,g],[d,f,g],[f,b,g],[f,b,h],[f,h,c],
            [h,a,e],[h,c,a],[e,h,g],[h,b,g],[a,d,f],[a,f,c]
            ]), invisible=True))

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
                    if (not prune) or (not self._intersectPolyhedrons(a,b)):
                        self._graph.add_edge(tuple(a), tuple(b), weight=np.linalg.norm(a-b))

    def _intersectPolyhedrons(self, a, b):
        for polyhedron in self._polyhedrons:
            if polyhedron.intersectSegment(a,b):
                return True
        return False
                    
    def calculateShortestPath(self, start, end, prune=True):
        for node in self._graph.nodes():
            if (not prune) or (not self._intersectPolyhedrons(start,np.array(node))):
                self._graph.add_edge(tuple(start), node, weight=np.linalg.norm(start-node))
            if (not prune) or (not self._intersectPolyhedrons(end,np.array(node))):
                self._graph.add_edge(tuple(end), node, weight=np.linalg.norm(end-node))

        try:
            length,path=nx.bidirectional_dijkstra(self._graph, tuple(start), tuple(end))
        except nx.NetworkXNoPath:
            path = []
        self._shortestPath = np.array(path)

    def plotSites(self, plotter):
        plotter.plot(self._sites[:,0], self._sites[:,1], self._sites[:,2], 'o')

    def plotPolyhedrons(self, plotter):
        for poly in self._polyhedrons:
            poly.plot(plotter)
            
    def plotShortestPath(self, plotter):
        if self._shortestPath.size > 0:
            plotter.plot(self._shortestPath[:,0], self._shortestPath[:,1], self._shortestPath[:,2], 'r', lw=2)
            plotter.plot([self._shortestPath[0][0]], [self._shortestPath[0][1]], [self._shortestPath[0][2]], 'ro')
            plotter.plot([self._shortestPath[-1][0]], [self._shortestPath[-1][1]], [self._shortestPath[-1][2]], 'ro')

    def plotGraph(self, plotter, vertexes=True, edges=True, labels=True):
        if vertexes:
            i = 0
            for ver in self._graph.nodes():
                plotter.plot([ver[0]], [ver[1]], [ver[2]], 'og')
                if labels:
                    plotter.text(ver[0], ver[1], ver[2], i, color='red')
                    i = i+1

        if edges:
            for edge in self._graph.edges():
                plotter.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], [edge[0][2], edge[1][2]], 'k--')

