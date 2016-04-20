import numpy as np
import scipy as sp
import scipy.spatial
import polyhedron

class ConvexHull(polyhedron.Polyhedron):
    def __init__(self, points, invisible=False, distributePoints=True, maxEmptyArea=0.1):
        convHull = sp.spatial.ConvexHull(points)
        faces = []
        for simplex in convHull.simplices:
            faces.append([convHull.points[simplex[0]], convHull.points[simplex[1]], convHull.points[simplex[2]]])

        super(ConvexHull, self).__init__(np.array(faces), invisible, distributePoints, maxEmptyArea)

