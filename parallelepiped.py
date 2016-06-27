import numpy as np
import polyhedron

class Parallelepiped(polyhedron.Polyhedron):
    def __init__(self, a, b, invisible=False, distributePoints=True, maxEmptyArea=0.1, boundingBox=False):
        
        c = [a[0], b[1], a[2]]
        d = [b[0], a[1], a[2]]
        e = [a[0], a[1], b[2]]
        f = [b[0], b[1], a[2]]
        g = [b[0], a[1], b[2]]
        h = [a[0], b[1], b[2]]

        super(Parallelepiped, self).__init__(faces=np.array([
            [a,g,e],[a,d,g],[d,f,g],[f,b,g],[f,b,h],[f,h,c],
            [h,a,e],[h,c,a],[e,h,g],[h,b,g],[a,d,f],[a,f,c]
            ]), invisible=invisible, distributePoints=distributePoints, maxEmptyArea=maxEmptyArea, boundingBox=boundingBox)

