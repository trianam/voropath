import random
import math
import numpy as np
import scipy as sp
import scipy.interpolate

class Path:
    _initialTemperature = 10#1000
    _trials = 10#100
    _warmingRatio = 0.9#0.9
    _minTemperature=0.00001#0.00000001
    _minDeltaEnergy=0.000001
    _maxVlambdaPert = 1000.
    _maxVertexPertFactor = 100.
    _initialVlambda = 0.
    _changeVlambdaProbability = 0.05
    #====1
    #_useArcLen = True
    #_ratioCurvTorsLen = [0.1, 0.1, 0.8]
    #====2
    _useArcLen = False
    _ratioCurvTorsLen = [0.1, 0.1, 0.8]
    #====3
    #_useArcLen = True
    #_ratioCurvTorsLen = [0.3, 0.3, 0.4]

    def __init__(self, bsplineDegree, adaptivePartition):
        self._bsplineDegree = bsplineDegree
        self._adaptivePartition = adaptivePartition
        self._vertexes = np.array([])
        self._dimC = 0
        self._polyhedronsContainer = []
        self._vlambda = self._initialVlambda
        

    @property
    def vertexes(self):
        return self._vertexes

    def assignValues(self, path, polyhedronsContainer):
        self._vertexes = path
        self._dimC = self._vertexes.shape[1]

        self._polyhedronsContainer = polyhedronsContainer
        tau, u, spline, splineD1, splineD2, splineD3, curv, tors, arcLength, polLength = self._splinePoints(self._vertexes)
        self._maxVertexPert = polLength / self._maxVertexPertFactor


    def setBsplineDegree(self, bsplineDegree):
        self._bsplineDegree = bsplineDegree

    def setAdaptivePartition(self, adaptivePartition):
        self._adaptivePartition = adaptivePartition

    def clean(self, verbose, debug):
        if verbose:
            print('Clean path (avoid obstacles)', flush=True)

        newPath = []
        if len(self._vertexes) > 0:
            a = self._vertexes[0]
            newPath.append(self._vertexes[0])

        for i in range(1, len(self._vertexes)-1):
            v = self._vertexes[i]
            b = self._vertexes[i+1]

            intersect,intersectRes = self._polyhedronsContainer.triangleIntersectPolyhedrons(a, v, b)
            if intersect:
                alpha = intersectRes[1]

                a1 = (1.-alpha)*a + alpha*v
                b1 = alpha*v + (1.-alpha)*b

                newPath.append(a1)
                newPath.append(v)
                newPath.append(b1)

                a = b1
            else:
                newPath.append(v)

                a = v

        if len(self._vertexes) > 0:
            newPath.append(self._vertexes[len(self._vertexes)-1])

        self._vertexes = np.array(newPath)


    def anneal(self, verbose):
        if verbose:
            print('Anneal path', flush=True)

        tau, u, self._spline, splineD1, splineD2, splineD3, curv, tors, arcLength, polLength = self._splinePoints(self._vertexes)
        self._currentEnergy, self._maxCurvatureLength, self._currentConstraints = self._initializePathEnergy(self._vertexes, self._spline, splineD1, splineD2, self._vlambda)


        temperature = self._initialTemperature
        while True:
            initialEnergy = self._currentEnergy
            numMovedLambda = 0
            numMovedVertex = 0
            for i in range(self._trials):
                movedLambda,movedVertex = self._tryMove(temperature)
                if movedLambda:
                    numMovedLambda += 1
                if movedVertex:
                    numMovedVertex += 1
            deltaEnergy = abs(initialEnergy - self._currentEnergy)
            temperature = temperature * self._warmingRatio
            if verbose:
                print("T:{}; E:{}; DE:{}; L:{}; C:{}; ML:{}; MV:{}".format(temperature, self._currentEnergy, deltaEnergy, self._vlambda, self._currentConstraints, numMovedLambda, numMovedVertex), flush=True)
                #print(self._vertexes)

            if (temperature < self._minTemperature) or (numMovedVertex > 0 and (deltaEnergy < self._minDeltaEnergy) and self._currentConstraints == 0.):
                break


    def simplify(self, verbose, debug):
        if verbose:
            print('Simplify path (remove useless triples)', flush=True)
        if self._bsplineDegree == 2:
            self._simplify2()
        elif self._bsplineDegree == 3:
            self._simplify3()
        elif self._bsplineDegree == 4:
            self._simplify4()
            
    def _simplify2(self):
        simplifiedPath = []
        if len(self._vertexes) > 0:
            a = self._vertexes[0]
            simplifiedPath.append(self._vertexes[0])
        first = True
        for i in range(1,len(self._vertexes)-1):
            v = self._vertexes[i]
            b = self._vertexes[i+1]
            keepV = False

            intersectCurr,nihil = self._polyhedronsContainer.triangleIntersectPolyhedrons(a, v, b)

            if not intersectCurr:
                if first:
                    intersectPrec = False
                else:
                    #a1 = self._vertexes[i-2]
                    intersectPrec,nihil = self._polyhedronsContainer.triangleIntersectPolyhedrons(a1, a, b)

                if i == len(self._vertexes)-2:
                    intersectSucc = False
                else:
                    b1 = self._vertexes[i+2]
                    intersectSucc,nihil = self._polyhedronsContainer.triangleIntersectPolyhedrons(a, b, b1)

                if intersectPrec or intersectSucc:
                    keepV = True

            else:
                keepV = True

            if keepV:
                first = False
                simplifiedPath.append(v)
                a1 = a
                a = v

        if len(self._vertexes) > 0:
            simplifiedPath.append(self._vertexes[len(self._vertexes)-1])

        self._vertexes = np.array(simplifiedPath)

    def _simplify3(self):
        simp = list(self._vertexes)
        toEval = 0
        while toEval < len(simp)-2:
            toEval += 1
            if toEval >= 3:
                if toEval < len(simp)-1:
                    if self._polyhedronsContainer.convexHullIntersectsPolyhedrons([simp[toEval-3], simp[toEval-2], simp[toEval-1], simp[toEval+1]]):
                        continue

            if toEval >= 2:
                if toEval < len(simp)-2:
                    if self._polyhedronsContainer.convexHullIntersectsPolyhedrons([simp[toEval-2], simp[toEval-1], simp[toEval+1], simp[toEval+2]]):
                        continue

            if toEval < len(simp)-3:
                if self._polyhedronsContainer.convexHullIntersectsPolyhedrons([simp[toEval-1], simp[toEval+1], simp[toEval+2], simp[toEval+3]]):
                        continue

            del simp[toEval]
            toEval -= 1
            
        self._vertexes = np.array(simp)

    def _simplify4(self):
        simp = list(self._vertexes)
        toEval = 0
        while toEval < len(simp)-2:
            toEval += 1
            if toEval >= 4:
                if toEval < len(simp)-1:
                    if self._polyhedronsContainer.convexHullIntersectsPolyhedrons([simp[toEval-4], simp[toEval-3], simp[toEval-2], simp[toEval-1], simp[toEval+1]]):
                        continue

            if toEval >= 3:
                if toEval < len(simp)-2:
                    if self._polyhedronsContainer.convexHullIntersectsPolyhedrons([simp[toEval-3], simp[toEval-2], simp[toEval-1], simp[toEval+1], simp[toEval+2]]):
                        continue

            if toEval >= 2:
                if toEval < len(simp)-3:
                    if self._polyhedronsContainer.convexHullIntersectsPolyhedrons([simp[toEval-2], simp[toEval-1], simp[toEval+1], simp[toEval+2], simp[toEval+3]]):
                        continue

            if toEval < len(simp)-4:
                if self._polyhedronsContainer.convexHullIntersectsPolyhedrons([simp[toEval-1], simp[toEval+1], simp[toEval+2], simp[toEval+3], simp[toEval+4]]):
                        continue

            del(simp[toEval])
            toEval -= 1

        self._vertexes = np.array(simp)
        
    def addNAlignedVertexes(self, numVertexes, verbose, debug):
        if verbose:
            print('Increase degree', flush=True)

        newPath = []
        for i in range(1, len(self._vertexes)):
            a = self._vertexes[i-1]
            b = self._vertexes[i]
            newPath.append(a)

            if numVertexes == 1:
                n = 0.5 * a + 0.5 * b
                newPath.append(n)
                
            elif numVertexes == 2:
                n1 = 0.33 * b + 0.67 * a
                n2 = 0.33 * a + 0.67 * b
                newPath.append(n1)
                newPath.append(n2)
                
        if len(self._vertexes) > 0:
            newPath.append(self._vertexes[len(self._vertexes)-1])

        self._vertexes = np.array(newPath)

    def splinePoints(self):
        return self._splinePoints(self._vertexes)
    
    def _splinePoints(self, vertexes):
        
        x = vertexes[:,0]
        y = vertexes[:,1]
        z = vertexes[:,2]

        polLen = self._calculatePolyLength(vertexes)

        tau,t = self._createKnotPartition(vertexes)

        #[knots, coeff, degree]
        tck = [t,[x,y,z], self._bsplineDegree]

        u=np.linspace(0,1,(max(polLen*5,1000)),endpoint=True)

        out = sp.interpolate.splev(u, tck)
        outD1 = sp.interpolate.splev(u, tck, 1)
        outD2 = sp.interpolate.splev(u, tck, 2)

        spline = np.stack(out).T
        splineD1 = np.stack(outD1).T
        splineD2 = np.stack(outD2).T

        if self._bsplineDegree >= 3:
            outD3 = sp.interpolate.splev(u, tck, 3)
            splineD3 = np.stack(outD3).T
        else:
            splineD3 = None

        curv = []
        tors = []
        arcLength = 0.
        for i in range(len(u)):
            d1Xd2 = np.cross(splineD1[i], splineD2[i])
            Nd1Xd2 = np.linalg.norm(d1Xd2)
            Nd1 = np.linalg.norm(splineD1[i])

            currCurv = 0.
            if Nd1 >0.: #>= 1.:
                currCurv = Nd1Xd2 / math.pow(Nd1,3)

            currTors = 0.
            if self._bsplineDegree >= 3 and Nd1Xd2 > 0.: #>= 1.:
                try:
                    currTors = np.dot(d1Xd2, splineD3[i]) / math.pow(Nd1Xd2, 2)
                except RuntimeWarning:
                    currTors = 0.

            curv.append(currCurv)
            tors.append(currTors)

            if i >= 1:
                dMin = min(prevNd1, Nd1)
                dMax = max(prevNd1, Nd1)
                arcLength += (u[i]-u[i-1]) * (dMin + ((dMax-dMin) / 2.))
                
            prevNd1 = Nd1
                
        return (tau, u, spline, splineD1, splineD2, splineD3, curv, tors, arcLength, polLen)

    def _createKnotPartition(self, controlPolygon):
        nv = len(controlPolygon)
        nn = nv - self._bsplineDegree + 1

        if not self._adaptivePartition:
            T = np.linspace(0,1,nv-self._bsplineDegree+1,endpoint=True)
        else:
            d = [0]
            for j in range(1, nv):
                d.append(d[j-1] + np.linalg.norm(controlPolygon[j] - controlPolygon[j-1]))
            t = []
            for i in range(nn-1):
                a = i * (nv-1) / (nn-1)
                ai = math.floor(a)
                ad = a - ai
                p = ad * controlPolygon[ai+1] + (1-ad) * controlPolygon[ai]
                l = d[ai] + np.linalg.norm(p - controlPolygon[ai])
                t.append(l / d[nv-1])

            t.append(1.)

            T = np.array(t)

        Text = np.append([0]*self._bsplineDegree, T)
        Text = np.append(Text, [1]*self._bsplineDegree)

        return (T,Text)

    def _initializePathEnergy(self, vertexes, spline, splineD1, splineD2, vlambda):
        tau, u, spline, splineD1, splineD2, splineD3, curv, tors, arcLength, polLength = self._splinePoints(vertexes)
        if self._useArcLen:
            length = arcLength
        else:
            length = polLength
            
        self._initialLength = length
        maxCurvatureLength = self._calculateMaxCurvatureLength(length, curv, tors)
        constraints = self._calculateConstraints(spline)
        energy = maxCurvatureLength + vlambda * constraints

        return (energy, maxCurvatureLength, constraints)

    def _tryMove(self, temperature):
        """
        Move the path or lambda multipiers in a neighbouring state,
        with a certain acceptance probability.
        Pick a random vertex (except extremes), and move
        it in a random direction (with a maximum perturbance).
        Use a lagrangian relaxation because we need to evaluate
        min(measure(path)) given the constraint that all quadrilaters
        formed by 4 consecutive points in the path must be collision
        free; where measure(path) is, depending of the choose method,
        the length of the path or the mean
        of the supplementary angles of each pair of edges of the path.
        If neighbourMode=0 then move the node uniformly, if
        neighbourMode=1 then move the node with gaussian probabilities
        with mean in the perpendicular direction respect to the
        previous-next nodes axis.
        """

        movedLambda = False
        movedVertex = False
        moveVlambda = random.random() < self._changeVlambdaProbability
        if moveVlambda:
            newVlambda = self._vlambda
            newVlambda = newVlambda + (random.uniform(-1.,1.) * self._maxVlambdaPert)

            newEnergy = self._calculatePathEnergyLambda(newVlambda)

            #attention, different formula from below
            if (newEnergy > self._currentEnergy) or (math.exp(-(self._currentEnergy-newEnergy)/temperature) >= random.random()):
                self._vlambda = newVlambda
                self._currentEnergy = newEnergy
                movedLambda = True
        
        else:
            newVertexes = np.copy(self._vertexes)
            movedV = random.randint(1,len(self._vertexes) - 2) #don't change extremes

            moveC = random.randint(0,self._dimC - 1)
            newVertexes[movedV][moveC] = newVertexes[movedV][moveC] + (random.uniform(-1.,1.) * self._maxVertexPert)
                
            newEnergy,newMaxCurvatureLength,newConstraints = self._calculatePathEnergyVertex(newVertexes)

            #attention, different formula from above
            if (newEnergy < self._currentEnergy) or (math.exp(-(newEnergy-self._currentEnergy)/temperature) >= random.random()):
                self._vertexes = newVertexes
                self._currentEnergy = newEnergy
                self._currentMaxCurvatureLength = newMaxCurvatureLength
                self._currentConstraints = newConstraints
                movedVertex = True

        return (movedLambda, movedVertex)

    def _calculatePathEnergyLambda(self, vlambda):
        """
        calculate the energy when lambda is moved.
        """
        return (self._currentEnergy - (self._vlambda * self._currentConstraints) + (vlambda * self._currentConstraints))
    
    def _calculatePathEnergyVertex(self, vertexes):
        """
        calculate the energy when a vertex is moved and returns it.
        """
        tau, u, spline, splineD1, splineD2, splineD3, curv, tors, arcLength, polLength = self._splinePoints(vertexes)
        if self._useArcLen:
            length = arcLength
        else:
            length = polLength
            
        constraints = self._calculateConstraints(spline)#this is bottleneck
        maxCurvatureLength = self._calculateMaxCurvatureLength(length, curv, tors)

        energy = maxCurvatureLength + self._vlambda * constraints
            
        return (energy, maxCurvatureLength, constraints)

    def _calculatePolyLength(self, vertexes):
        length = 0.
        for i in range(1, len(vertexes)):
            length += sp.spatial.distance.euclidean(vertexes[i-1], vertexes[i])
            #length += np.linalg.norm(np.subtract(vertexes[i], vertexes[i-1]))
        return length

    def _calculateMaxCurvatureLength(self, length, curv, tors):
        normLength = length/self._initialLength * 100 #for making the ratio indipendent of the initial length

        maxCurvature = 0.
        maxTorsion = 0.
        for i in range(0, len(curv)):
            currCurv = curv[i]
            currTors = abs(tors[i])
            if currCurv > maxCurvature:
                maxCurvature = currCurv
            if currTors > maxTorsion:
                maxTorsion = currTors

        return self._ratioCurvTorsLen[0]*maxCurvature + self._ratioCurvTorsLen[1]*maxTorsion + self._ratioCurvTorsLen[2]*normLength

    def _calculateConstraints(self, spline):
        """
        calculate the constraints function. Is the ratio of the points
        of the calculated spline that are inside obstacles respect the
        total number of points of the spline.
        """
        pointsInside = 0
        for p in spline:
            if self._polyhedronsContainer.pointInsidePolyhedron(p):
                pointsInside = pointsInside + 1

        constraints = pointsInside / len(spline)

        return constraints
