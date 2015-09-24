import numpy as np
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class Polyhedron:
    _maxEmptyArea = 0.05
    
    def __init__(self, faces=np.array([])):
        self._faces = faces

        triangles = []
        allPoints = []

        for face in self._faces:
            triangles.append(face)
            
        while triangles:
            triangle = triangles.pop(0)
            a = triangle[0]
            b = triangle[1]
            c = triangle[2]
            if not any((a == x).all() for x in allPoints):
                allPoints.append(a)
            if not any((b == x).all() for x in allPoints):
                allPoints.append(b)
            if not any((c == x).all() for x in allPoints):
                allPoints.append(c)
            if (self._area(triangle) > Polyhedron._maxEmptyArea):
                ab = self._comb2(a,b)
                bc = self._comb2(b,c)
                ca = self._comb2(c,a)
                abc = self._comb3(a,b,c)
                
                triangles.append(np.array([a,ab,abc]))
                triangles.append(np.array([ab,b,abc]))
                triangles.append(np.array([b,bc,abc]))
                triangles.append(np.array([bc,c,abc]))
                triangles.append(np.array([c,ca,abc]))
                triangles.append(np.array([ca,a,abc]))

        self.allPoints = np.array(allPoints)
        
    def _area(self, triangle):
        a = np.linalg.norm(triangle[1]-triangle[0])
        b = np.linalg.norm(triangle[2]-triangle[1])
        c = np.linalg.norm(triangle[0]-triangle[2])
        s =  (a+b+c) / 2.
        return math.sqrt(s * (s-a) * (s-b) *(s-c))

    _comb2 = lambda self,a,b: 0.5*a + 0.5*b
    _comb3 = lambda self,a,b,c: 0.33*a + 0.33*b + 0.33*c
        
    def plotAllPoints(self, plotter):
        plotter.plot(self.allPoints[:,0], self.allPoints[:,1], self.allPoints[:,2], 'o')

    #def cross(self, pointA, pointB):
        #for face in self._faces:
        #    if
    #    return False

    def plot(self, plotter):
        col = Poly3DCollection(self._faces)
        col.set_color('y')
        col.set_edgecolor('k')
        plotter.add_collection3d(col)
