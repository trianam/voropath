import numpy as np
import math
from mpl_toolkits.mplot3d.art3d import Poly3DCollection

class Polyhedron:
    def __init__(self, **pars):
        """
        can be composed only by combined triangles
        pars can be one of:
           faces -> an np.array of triangular faces
           a,b,c,d -> lists of the four vertexes of a tetrahedron
        if invisible=True when plot will be called it will be useless
        """
        if 'maxEmptyArea' in pars:
            maxEmptyArea = pars['maxEmptyArea']
        else:
            maxEmptyArea = 0.1

        if 'faces' in pars:
            self._faces = pars['faces']
        elif 'a' in pars and 'b' in pars and 'c' in pars and 'd' in pars:
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
            
        if (not 'distributePoints' in pars) or (pars['distributePoints'] == True):
            self.distributePoints(maxEmptyArea)
        else:
            self.allPoints = np.array([])
        
    def _area(self, triangle):
        a = np.linalg.norm(triangle[1]-triangle[0])
        b = np.linalg.norm(triangle[2]-triangle[1])
        c = np.linalg.norm(triangle[0]-triangle[2])
        s =  (a+b+c) / 2.
        return math.sqrt(s * (s-a) * (s-b) *(s-c))

    _comb2 = lambda self,a,b: 0.5*a + 0.5*b
    _comb3 = lambda self,a,b,c: 0.33*a + 0.33*b + 0.33*c

    def distributePoints(self, maxEmptyArea):
        allPoints = []
        triangles = []

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
            if (self._area(triangle) > maxEmptyArea):
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

    
    def plotAllPoints(self, plotter):
        if self.allPoints.size > 0:
            plotter.plot(self.allPoints[:,0], self.allPoints[:,1], self.allPoints[:,2], 'ob')

    def intersectSegment(self, a, b):
        for triangle in self._faces:
            #solve {
            #        a+k(b-a) = v*triangle[0] + w*triangle[1] + s*triangle[2]
            #        v+w+s = 1
            #      }
            # for variables k, v, w, s

            #simplified in
            #      a+k(b-a) = (1-w-s)*triangle[0] + w*triangle[1] + s*triangle[2]
            # for variables k, w, s
            
            diffba = b-a
            difft0t1 = triangle[0] - triangle[1]
            difft0t2 = triangle[0] - triangle[2]
            difft0a = triangle[0] - a

            A = np.array([
                [diffba[0], difft0t1[0], difft0t2[0]],
                [diffba[1], difft0t1[1], difft0t2[1]],
                [diffba[2], difft0t1[2], difft0t2[2]]])
            B = np.array([difft0a[0], difft0a[1], difft0a[2]])
            
            try:
                x = np.linalg.solve(A,B)
                # check if
                #          0 <= k <= 1,
                #          w >= 0
                #          s >= 0
                #          w+s <= 1
                if (x[0] >= 0.) and (x[0] <= 1.) and (x[1] >= 0.) and (x[2] >= 0.) and (x[1]+x[2] <= 1.):
                    return (True, x)
            except np.linalg.linalg.LinAlgError:
                pass

        return (False,np.array([]))

    def intersectPolyhedron(self, polyhedron):
        """alert, not case of one polyhedron inside other"""
        for otherFace in polyhedron._faces:
            for myFace in self._faces:
                if (
                        self.intersectSegment(otherFace[0],otherFace[1])[0] or
                        self.intersectSegment(otherFace[1],otherFace[2])[0] or
                        self.intersectSegment(otherFace[2],otherFace[0])[0] or
                        polyhedron.intersectSegment(myFace[0], myFace[1])[0] or
                        polyhedron.intersectSegment(myFace[1], myFace[2])[0] or
                        polyhedron.intersectSegment(myFace[2], myFace[0])[0]):
                    return True
        return False
                
    def intersectPathTriple(self, triple):
        """alert, not case of one polyhedron inside other, and only
        check if the segments of self intersect the triple."""
        intersect = False
        result = np.array([])
        for myFace in self._faces:
            intersect1, result1 = triple.intersectSegment(myFace[0], myFace[1])
            intersect2, result2 = triple.intersectSegment(myFace[1], myFace[2])
            intersect3, result3 = triple.intersectSegment(myFace[2], myFace[0])
            if intersect1:
                intersect = True
                result = result1
            if intersect2 and (not intersect or (result2[1] > result[1])):
                intersect = True
                result = result2
            if intersect3 and (not intersect or (result3[1] > result[1])):
                intersect = True
                result = result3

        return intersect,result
                
                    
    def plot(self, plotter):
        if self._invisible == False:
            col = Poly3DCollection(self._faces)
            col.set_color('y')
            col.set_edgecolor('k')
            plotter.add_collection3d(col)
