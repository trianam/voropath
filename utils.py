import numpy as np
import numpy.linalg


def calcStartEnd(graph, Vs, Ve):
    minS = 0
    minE = 0
    first = True
    for ver in graph.nodes():
        if(first):
            Vss = ver
            Vee = ver
            minS = np.linalg.norm(Vs-ver)
            minE = np.linalg.norm(Ve-ver)
            first = False
        else:
            currS = np.linalg.norm(Vs-ver)
            currE = np.linalg.norm(Ve-ver)
            if(currS < minS):
                Vss = ver
                minS = currS
            if(currE < minE):
                Vee = ver
                minE = currE
    return (Vss, Vee)
