import numpy as np
import numpy.linalg
import scipy as sp
import scipy.spatial
import networkx as nx

def makeVoroGraph(sites):
    graph = nx.Graph()
    vor = sp.spatial.Voronoi(sites)
    vorVer = vor.vertices
    for ridge in vor.ridge_vertices:
        for i in range(len(ridge)):
            if (ridge[i] != -1) and (ridge[(i+1)%len(ridge)] != -1):
                graph.add_edge(tuple(vorVer[ridge[i]]), tuple(vorVer[ridge[(i+1)%len(ridge)]]), weight=np.linalg.norm(vorVer[ridge[i]]-vorVer[ridge[(i+1)%len(ridge)]]))

    return graph

