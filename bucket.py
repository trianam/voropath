import numpy as np
import compositePolyhedron
import parallelepiped

class Bucket(compositePolyhedron.CompositePolyhedron):
    def __init__(self, center, width, height, thickness, invisible=False, distributePoints=True, maxEmptyArea=0.1, boundingBox=False):
        c = center
        l = width
        h = height
        d = thickness
        parallelepipeds = []
        
        parallelepipeds.append(parallelepiped.Parallelepiped(
            np.array([c[0]-(l/2), c[1]-(l/2), c[2]-(h/2)]),\
            np.array([c[0]+(l/2), c[1]+(l/2), c[2]-(h/2)+d]), invisible, distributePoints, maxEmptyArea, boundingBox))
        
        parallelepipeds.append(parallelepiped.Parallelepiped(
            np.array([c[0]-(l/2), c[1]+(l/2)-d, c[2]-(h/2)+d]),\
            np.array([c[0]+(l/2), c[1]+(l/2), c[2]+(h/2)]), invisible, distributePoints, maxEmptyArea, boundingBox))

        parallelepipeds.append(parallelepiped.Parallelepiped(
            np.array([c[0]-(l/2), c[1]-(l/2), c[2]-(h/2)+d]),\
            np.array([c[0]+(l/2), c[1]-(l/2)+d, c[2]+(h/2)]), invisible, distributePoints, maxEmptyArea, boundingBox))

        parallelepipeds.append(parallelepiped.Parallelepiped(
            np.array([c[0]-(l/2), c[1]-(l/2)+d, c[2]-(h/2)+d]),\
            np.array([c[0]-(l/2)+d, c[1]+(l/2)-d, c[2]+(h/2)]), invisible, distributePoints, maxEmptyArea, boundingBox))

        parallelepipeds.append(parallelepiped.Parallelepiped(
            np.array([c[0]+(l/2)-d, c[1]-(l/2)+d, c[2]-(h/2)+d]),\
            np.array([c[0]+(l/2), c[1]+(l/2)-d, c[2]+(h/2)]), invisible, distributePoints, maxEmptyArea, boundingBox))

        super(Bucket, self).__init__(parallelepipeds)
        
