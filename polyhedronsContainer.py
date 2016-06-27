import numpy as np
import scipy as sp
import scipy.spatial
import polyhedron

class PolyhedronsContainer:
    def __init__(self):
        self._polyhedrons = []
        self._hasBoundingBox = False
        self._boundingBoxA = None
        self._boundingBoxB = None

    @property
    def polyhedrons(self):
        return self._polyhedrons

    @property
    def hasBoundingBox(self):
        return self._hasBoundingBox

    @hasBoundingBox.setter
    def hasBoundingBox(self, value):
        self._hasBoundingBox = value

    @property
    def boundingBoxA(self):
        return self._boundingBoxA
        
    @boundingBoxA.setter
    def boundingBoxA(self, value):
        self._boundingBoxA = value

    @property
    def boundingBoxB(self):
        return self._boundingBoxB
        
    @boundingBoxB.setter
    def boundingBoxB(self, value):
        self._boundingBoxB = value

    def addPolyhedron(self, polyhedron):
        self._polyhedrons.append(polyhedron)

    def addBoundingBox(self, a, b, maxEmptyArea, invisible):
        self._hasBoundingBox = True
        self._boundingBoxA = a
        self._boundingBoxB = b

        c = [a[0], b[1], a[2]]
        d = [b[0], a[1], a[2]]
        e = [a[0], a[1], b[2]]
        f = [b[0], b[1], a[2]]
        g = [b[0], a[1], b[2]]
        h = [a[0], b[1], b[2]]

        self.addPolyhedron(polyhedron.Polyhedron(faces=np.array([
            [a,g,e],[a,d,g],[d,f,g],[f,b,g],[f,b,h],[f,h,c],
            [h,a,e],[h,c,a],[e,h,g],[h,b,g],[a,d,f],[a,f,c]
            ]), invisible=invisible, maxEmptyArea=maxEmptyArea, boundingBox=True))

    def pointInsidePolyhedron(self, p):
        inside = False
        if self._hasBoundingBox:
            if (p<self._boundingBoxA).any() or (p>self._boundingBoxB).any():
                inside = True

        if not inside:
            for polyhedron in self._polyhedrons:
                if (not polyhedron.isBoundingBox()) and polyhedron.hasPointInside(p):
                    inside = True
                    break

        return inside
        

    def segmentIntersectPolyhedrons(self, a, b, intersectionMargin = 0.):
        intersect = False
        if self._hasBoundingBox:
            if((a<self._boundingBoxA).any() or (a>self._boundingBoxB).any() or (b<self._boundingBoxA).any() or (b>self._boundingBoxB).any()):
                intersect = True

        if not intersect:
            minS = np.array([min(a[0],b[0]),min(a[1],b[1]),min(a[2],b[2])])
            maxS = np.array([max(a[0],b[0]),max(a[1],b[1]),max(a[2],b[2])])

            for polyhedron in self._polyhedrons:
                if polyhedron.intersectSegment(a,b,minS,maxS, intersectionMargin=intersectionMargin)[0]:
                    intersect = True
                    break
                    
        return intersect

    def triangleIntersectPolyhedrons(self, a, b, c):
        triangle = polyhedron.Polyhedron(faces=np.array([[a,b,c]]), distributePoints = False)
        intersect = False
        result = np.array([])
        for currPolyhedron in self._polyhedrons:
            currIntersect,currResult = currPolyhedron.intersectPathTriple(triangle)
            if currIntersect and (not intersect or (currResult[1] > result[1])):
                intersect = True
                result = currResult

        return (intersect, result)

    def convexHullIntersectsPolyhedrons(self, vertexes):
        convHull = sp.spatial.ConvexHull(vertexes, qhull_options="QJ Pp")
        for simplex in convHull.simplices:
            if self.triangleIntersectPolyhedrons(convHull.points[simplex[0]], convHull.points[simplex[1]], convHull.points[simplex[2]])[0]:
                return True

        return False
