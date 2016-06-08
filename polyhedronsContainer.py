import numpy as np
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

    @property
    def boundingBoxA(self):
        return self._boundingBoxA
        
    @property
    def boundingBoxB(self):
        return self._boundingBoxB
        
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
            ]), invisible=invisible, maxEmptyArea=maxEmptyArea))

    def pointInsidePolyhedron(self, p):
        inside = False
        if self._hasBoundingBox:
            if (p<self._boundingBoxA).any() or (p>self._boundingBoxB).any():
                inside = True

        if not inside:
            for polyhedron in self._polyhedrons:
                if polyhedron.hasPointInside(p):
                    inside = True
                    break
                    
        return inside
        

    def segmentIntersectPolyhedrons(self, a, b):
        intersect = False
        if self._hasBoundingBox:
            if((a<self._boundingBoxA).any() or (a>self._boundingBoxB).any() or (b<self._boundingBoxA).any() or (b>self._boundingBoxB).any()):
                intersect = True

        if not intersect:
            for polyhedron in self._polyhedrons:
                if polyhedron.intersectSegment(a,b)[0]:
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

