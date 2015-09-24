import numpy as np
import numpy.linalg
import scipy as sp
import scipy.spatial
import networkx as nx
import numpy.linalg

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

    def setPolyhedronsSites(self):
        sites = []
        for polyhedron in self._polyhedrons:
            sites.extend(polyhedron.allPoints)

        self._sites = np.array(sites)
        
    def makeVoroGraph(self):
        vor = sp.spatial.Voronoi(self._sites)
        vorVer = vor.vertices
        for ridge in vor.ridge_vertices:
            for i in range(len(ridge)):
                if (ridge[i] != -1) and (ridge[(i+1)%len(ridge)] != -1):
                    self._graph.add_edge(tuple(vorVer[ridge[i]]), tuple(vorVer[ridge[(i+1)%len(ridge)]]), weight=np.linalg.norm(vorVer[ridge[i]]-vorVer[ridge[(i+1)%len(ridge)]]))

    #def pruneVoroGraph(self):
                    
    def calculateShortestPath(self, start, end):
        startVertex, endVertex = self._calcStartEnd(start, end)

        length,path=nx.bidirectional_dijkstra(self._graph,startVertex,endVertex)
        path.insert(0,start)
        path.append(end)
        self._shortestPath = np.array(path)

    def plotSites(self, plotter):
        plotter.plot(self._sites[:,0], self._sites[:,1], self._sites[:,2], 'o')

    def plotShortestPath(self, plotter):
        plotter.plot(self._shortestPath[:,0], self._shortestPath[:,1], self._shortestPath[:,2], 'r', lw=2)

    def plotGraph(self, plotter):
        i = 0
        for ver in self._graph.nodes():
            plotter.plot([ver[0]], [ver[1]], [ver[2]], 'og')
            plotter.text(ver[0], ver[1], ver[2], i, color='red')
            i = i+1

        for edge in self._graph.edges():
            plotter.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], [edge[0][2], edge[1][2]], 'k--')

    def _calcStartEnd(self, start, end):
        minS = 0
        minE = 0
        first = True
        for ver in self._graph.nodes():
            if(first):
                startVertex = ver
                endVertex = ver
                minS = np.linalg.norm(start-ver)
                minE = np.linalg.norm(end-ver)
                first = False
            else:
                currS = np.linalg.norm(start-ver)
                currE = np.linalg.norm(end-ver)
                if(currS < minS):
                    startVertex = ver
                    minS = currS
                if(currE < minE):
                    endVertex = ver
                    minE = currE
        return (startVertex, endVertex)

