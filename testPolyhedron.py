#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import polyhedron

poly1 = polyhedron.Polyhedron(np.array(
    [
        [[0.,0.,0.],[0.,1.,0.],[1.,0.,0.]],
        [[0.,0.,0.],[0.5,0.5,1.],[0.,1.,0.]],
        [[0.,1.,0.],[0.5,0.5,1.],[1.,0.,0.]],
        [[0.,0.,0.],[0.5,0.5,1.],[1.,0.,0.]]
    ]))

# poly1 = polyhedron.Polyhedron([
#         [(0.,0.,0.),(0.,1.,0.),(1.,0.,0.)],
#         [(0.,0.,0.),(0.5,0.5,1.),(0.,1.,0.)],
#         [(0.,1.,0.),(0.5,0.5,1.),(1.,0.,0.)],
#         [(0.,0.,0.),(0.5,0.5,1.),(1.,0.,0.)]
#     ]
#)

fig = plt.figure()
ax = fig.gca(projection='3d')

poly1.plot(ax)
poly1.plotAllPoints(ax)

# ax.set_xlim3d(0., 1.)
# ax.set_ylim3d(0., 1.)
# ax.set_zlim3d(0., 1.)
            
plt.show()
