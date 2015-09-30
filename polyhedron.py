import numpy as np
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class Polyhedron:
    _maxEmptyArea = 0.1
    
    def __init__(self, **pars):
        """
        pars can be one of:
           faces -> an np.array of triangular faces
           a,b,c,d -> lists of the four vertexes of a triangular pyramid
        if invisible=True when plot will be called it will be useless
        """
        
        if 'faces' in pars:
            self._faces = pars['faces']
        elif 'a' in pars and 'b' in pars and 'c' in pars:
            a = pars['a']
            b = pars['b']
            c = pars['c']
            d = pars['d']
            self._faces = np.array([
                [a,b,c],
                [a,b,d],
                [b,c,d],
                [c,a,d]
            ])

        if 'invisible' in pars:
            self._invisible = pars['invisible']
        else:
            self._invisible = False
            
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
        plotter.plot(self.allPoints[:,0], self.allPoints[:,1], self.allPoints[:,2], 'ob')

    def intersectSegment(self, a, b):
        for triangle in self._faces:
            A = np.vstack(
                (np.vstack(
                    ((b-a), -1*triangle[0], -1*triangle[1], -1*triangle[2])
                ).T,np.array([0,1,1,1]))
            )
            B = np.append(-1*a,1)
            try:
                x = np.linalg.solve(A,B)
                if (x[0] >= 0.) and (x[0] <= 1.) and (x[1] >= 0.) and (x[2] >= 0.) and (x[3] >= 0.):
                    return True
            except np.linalg.linalg.LinAlgError:
                pass

        return False

    def plot(self, plotter):
        if self._invisible == False:
            col = Poly3DCollection(self._faces)
            col.set_color('y')
            col.set_edgecolor('k')
            plotter.add_collection3d(col)
