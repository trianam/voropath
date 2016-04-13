import numpy as np
import polyhedron

class Tetrahedron(polyhedron.Polyhedron):
    def __init__(self, a, b, c, d, invisible=False, distributePoints=True, maxEmptyArea=0.1):
        super(Tetrahedron, self).__init__(np.array([[a,b,c],[a,b,d],[b,c,d],[c,a,d]]), invisible, distributePoints, maxEmptyArea)

        self._vertexes = [a,b,c,d]

    def plot(self, plotter):
        if self._invisible == False:
            plotter.addTetrahedron(self._vertexes, plotter.COLOR_OBSTACLE)

            
            


