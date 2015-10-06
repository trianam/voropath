#!/usr/bin/python

from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import matplotlib.pyplot as plt
import voronizator
import polyhedron

voronoi = voronizator.Voronizator()

fig = plt.figure()
ax = fig.gca(projection='3d')

voronoi.addPolyhedron(polyhedron.Polyhedron(a = [0.1,0.1,0.1], b = [0.3,0.3,0.1], c = [0.4,0.2,0.2], d = [0.1,0.2,0.2]))

#make test graph
a = (0., 0. ,0.)
b = (0.25, 0.15, 0.18)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (0.25, 0.2, 0.22)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (0.25, 0.25, 0.18)
#b = (0.25, 0.25, 0.4)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (1., 1., 1.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = (0., 0. ,0.)
b = (0.25, 0.15, 0.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (0.4, 0.4, 0.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (0.8, 0., 0.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))
        
a = b
b = (0.8, 1., 0.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (0.8, 1., 1.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (0., 1., 1.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (1., 0., 1.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))

a = b
b = (1., 1., 1.)
voronoi._graph.add_edge(a, b, weight=np.linalg.norm(np.array(a)-np.array(b)))


#calculate test shortest path
voronoi._pathStart = np.array([0.,0.,0.])
voronoi._pathEnd = np.array([1.,1.,1.])

voronoi._shortestPath = voronoi._dijkstraPlus(voronoi._pathStart, voronoi._pathEnd)


#plot
voronoi.plotPolyhedrons(ax)
voronoi.plotGraph(ax, pathExtremes=True)
voronoi.plotShortestPath(ax)

# ax.set_xlim3d(0., 1.)
# ax.set_ylim3d(0., 1.)
# ax.set_zlim3d(0., 1.)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
