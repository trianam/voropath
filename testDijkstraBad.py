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

ns = (0., 0., 0.)
n1 = (0.25, 0.15, 0.18)
n2 = (0.25, 0.2, 0.3)
n3 = (0.25, 0.25, 0.18)
n4 = (0.15, 0.1, 0.3)
ne = (0.4, 0.4, 0.4)

voronoi._graph.add_edge(ns, n1, weight=np.linalg.norm(np.array(ns)-np.array(n1)))
voronoi._graph.add_edge(n1, n2, weight=np.linalg.norm(np.array(n1)-np.array(n2)))
voronoi._graph.add_edge(n2, n3, weight=np.linalg.norm(np.array(n2)-np.array(n3)))
voronoi._graph.add_edge(n3, ne, weight=np.linalg.norm(np.array(n3)-np.array(ne)))
voronoi._graph.add_edge(ns, n4, weight=np.linalg.norm(np.array(ns)-np.array(n4)))
voronoi._graph.add_edge(n4, n2, weight=np.linalg.norm(np.array(n4)-np.array(n2)))

voronoi._graph.node[ns]['index'] = 's'
voronoi._graph.node[n1]['index'] = 1
voronoi._graph.node[n2]['index'] = 2
voronoi._graph.node[n3]['index'] = 3
voronoi._graph.node[n4]['index'] = 4
voronoi._graph.node[ne]['index'] = 'e'


#calculate test shortest path
voronoi._pathStart = np.array(ns)
voronoi._pathEnd = np.array(ne)

#voronoi._shortestPath = voronoi._dijkstraPlus(voronoi._pathStart, voronoi._pathEnd)
voronoi._shortestPath = voronoi._trijkstra(voronoi._pathStart, voronoi._pathEnd)


#plot
voronoi.plotPolyhedrons(ax)
voronoi.plotGraph(ax, pathExtremes=True, labels=True)
voronoi.plotShortestPath(ax)

ax.set_xlim3d(0., .4)
ax.set_ylim3d(0., .4)
ax.set_zlim3d(0., .4)

ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')

plt.show()
