import numpy as np
import scipy as sp
import scipy.spatial
import math
import xml.etree.cElementTree as ET

class Polyhedron:
    def __init__(self, faces, invisible=False, distributePoints=True, maxEmptyArea=0.1, boundingBox=False):
        """
        can be composed only by combined triangles
        faces -> an np.array of triangular faces
        if invisible=True when plot will be called it will be useless
        """
        self._faces = faces
        self._invisible = invisible
        self._boundingBox = boundingBox

        self._minV = np.array([float('inf'),float('inf'),float('inf')])
        self._maxV = np.array([float('-inf'),float('-inf'),float('-inf')])
        for face in self._faces:
            for vertex in face:
                for i in range(len(vertex)):
                    if vertex[i] < self._minV[i]:
                        self._minV[i] = vertex[i]

                for i in range(len(vertex)):
                    if vertex[i] > self._maxV[i]:
                        self._maxV[i] = vertex[i]
                    
        if distributePoints:
            self.distributePoints(maxEmptyArea)
        else:
            self.allPoints = np.array([])

    @property
    def minV(self):
        return self._minV

    @property
    def maxV(self):
        return self._maxV
    
    def isBoundingBox(self):
        return self._boundingBox
        
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

    def hasPointInside(self, p):
        """
        check if a point is inside the convex hull of obstacle vertexes
        """
        outside = True
        if (p>self._minV).all() and (p<self._maxV).all():
            vertexes = [p]
            for triangle in self._faces:
                vertexes.append(triangle[0])
                vertexes.append(triangle[1])
                vertexes.append(triangle[2])

            chull = sp.spatial.ConvexHull(np.array(vertexes))
            outside = False
            for vertex in chull.vertices:
                if (p == chull.points[vertex]).all():
                    outside = True
                    break
            # for simplex in chull.simplices:
            #     if (p == chull.points[simplex[0]]).all() or (p == chull.points[simplex[1]]).all() or (p == chull.points[simplex[2]]).all():
            #         outside = True
            #         break

        return not outside
        
    def intersectSegment(self, a, b, minS=None, maxS=None):
        if minS is None or maxS is None:
            minS = np.array([min(a[0],b[0]),min(a[1],b[1]),min(a[2],b[2])])
            maxS = np.array([max(a[0],b[0]),max(a[1],b[1]),max(a[2],b[2])])
            
        if not ((self._minV > maxS).any() or (self._maxV < minS).any()):
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
                    #          0 < k < 1,
                    #          w > 0
                    #          s > 0
                    #          w+s < 1
                    if (x[0] > 0.) and (x[0] < 1.) and (x[1] > 0.) and (x[2] > 0.) and (x[1]+x[2] < 1.):
                        return (True, x)
                except np.linalg.linalg.LinAlgError:
                    pass

        return (False,np.array([]))

    def intersectPolyhedron(self, polyhedron):
        """alert, not case of one polyhedron inside other"""
        if not ((self._minV > polyhedron.maxV).any() or (self._maxV < polyhedron.minV).any()):
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

        if not ((self._minV > triple.maxV).any() or (self._maxV < triple.minV).any()):
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
                
                    
    def plotAllPoints(self, plotter):
        if self.allPoints.size > 0:
            plotter.addPoints(self.allPoints, plotter.COLOR_SITES)

    def plot(self, plotter):
        if self._invisible == False:
            plotter.addTriangles(self._faces, plotter.COLOR_OBSTACLE)

    def extractXmlTree(self, root):
        xmlPolyhedron = ET.SubElement(root, 'polyhedron', invisible=str(self._invisible))
        for face in self._faces:
            xmlFace = ET.SubElement(xmlPolyhedron, 'face')
            for vertex in face:
                xmlVertex = ET.SubElement(xmlFace, 'vertex', x=str(vertex[0]), y=str(vertex[1]), z=str(vertex[2]))

