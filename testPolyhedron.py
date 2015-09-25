#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import polyhedron

poly1 = polyhedron.Polyhedron(faces = np.array(
    [
        [[0.,0.,0.],[0.,1.,0.],[1.,0.,0.]],
        [[0.,0.,0.],[0.5,0.5,1.],[0.,1.,0.]],
        [[0.,0.,0.],[0.5,0.5,1.],[1.,0.,0.]]
    ]))

poly2 = polyhedron.Polyhedron(a=[2.,2.,2.], b=[3.,2.,2.], c=[2.5,3.,2.], d=[2.5,2.5,3.])

fig = plt.figure()
ax = fig.gca(projection='3d')

poly1.plot(ax)
poly1.plotAllPoints(ax)

poly2.plot(ax)
poly2.plotAllPoints(ax)
            
plt.show()
