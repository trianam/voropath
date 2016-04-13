#!/usr/bin/python

import numpy as np
import polyhedron
import tetrahedron
import plotter

poly1 = polyhedron.Polyhedron(faces = np.array(
    [
        [[0.,0.,0.],[0.,1.,0.],[1.,0.,0.]],
        [[0.,0.,0.],[0.5,0.5,1.],[0.,1.,0.]],
        [[0.,0.,0.],[0.5,0.5,1.],[1.,0.,0.]],
        [[1.,0.,0.],[0.,1.,0.],[0.5,0.5,1.]]
    ]))

poly2 = tetrahedron.Tetrahedron(a=[2.,2.,2.], b=[3.,2.,2.], c=[2.5,3.,2.], d=[2.5,2.5,3.])

plt = plotter.Plotter()

poly1.plot(plt)
poly1.plotAllPoints(plt)

poly2.plot(plt)
poly2.plotAllPoints(plt)
            
plt.draw()

