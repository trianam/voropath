import numpy as np
import polyhedron

class CompositePolyhedron(polyhedron.Polyhedron):
    def __init__(self, components):
        self._components = components

        self._boundingBox = False

        self._minV = np.array([float('inf'),float('inf'),float('inf')])
        self._maxV = np.array([float('-inf'),float('-inf'),float('-inf')])

        for component in self._components:
            for i in range(3):
                if component.minV[i] < self._minV[i]:
                    self._minV[i] = component.minV[i]

                if component.maxV[i] > self._maxV[i]:
                    self._maxV[i] = component.maxV[i]

    @property
    def allPoints(self):
        allPoints = []
        for component in self._components:
            allPoints.extend(list(component.allPoints))
            
        return np.array(allPoints)
            
    @property
    def minV(self):
        return self._minV

    @property
    def maxV(self):
        return self._maxV

                    
    def distributePoints(self, maxEmptyArea):
        for component in self._components:
            component.distributePoints(maxEmptyArea)

    def hasPointInside(self, p):
        hasPI = False
        for component in self._components:
            if component.hasPointInside(p):
                hasPI = True
                break

        return hasPI

    def intersectSegment(self, a, b, minS=None, maxS=None, intersectionMargin=0.):
        intersect = (False,np.array([]))
        for component in self._components:
            current = component.intersectSegment(a,b,minS,maxS,intersectionMargin)
            if current[0]:
                intersect = current
                break

        return intersect

    def intersectPolyhedron(self, polyhedron):
        intersect = False
        for component in self._components:
            if component.intersectPolyhedron(polyhedron):
                intersect = True
                break

        return intersect

    def intersectPathTriple(self, triple):
        intersect = (False, np.array([]))
        for component in self._components:
            current = component.intersectPathTriple(triple)
            if current[0]:
                intersect = current
                break

        return intersect

    def plotAllPoints(self, plotter):
        for component in self._components:
            component.plotAllPoints(plotter)

    def plot(self, plotter):
        for component in self._components:
            component.plot(plotter)

    def extractXmlTree(self, root):
        for component in self._components:
            component.extractXmlTree(root)
        

