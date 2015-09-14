from mpl_toolkits.mplot3d import  axes3d,Axes3D
import numpy as np
import scipy as sp
#import scipy.misc
#import math
import matplotlib.pyplot as plt
import scipy.spatial
import networkx as nx
import voronoi
import utils

np.random.seed(0)

fig = plt.figure()
ax = fig.gca(projection='3d')

#V = np.array([[0.,0.,0.],[0.,1.,0.],[1.,1.,1.],[2.,1.,2.],[3.,2.,1.],[4.,3.,6.],[2.,3.,1.],[5.,3.,1.],[6.,5.,3.]])
V = sp.rand(10,3)
Vs = np.array([0.,0.,0.])
Ve = np.array([1.,1.,1.])

ax.plot(V[:,0], V[:,1], V[:,2], 'o')

graph = voronoi.makeVoroGraph(V)

Vss, Vee = utils.calcStartEnd(graph, Vs, Ve)

length,path=nx.bidirectional_dijkstra(graph,Vss,Vee)
pathArr = np.array(path)

ax.plot(pathArr[:,0], pathArr[:,1], pathArr[:,2], 'r', lw=2)

i = 0
for ver in graph.nodes():
    ax.plot([ver[0]], [ver[1]], [ver[2]], 'og')
    ax.text(ver[0], ver[1], ver[2], i, color='red')
    i = i+1

for edge in graph.edges():
    ax.plot([edge[0][0], edge[1][0]], [edge[0][1], edge[1][1]], [edge[0][2], edge[1][2]], 'k--')

ax.plot([Vs[0], Vss[0]], [Vs[1], Vss[1]], [Vs[2], Vss[2]], 'r', lw=2)
ax.plot([Ve[0], Vee[0]], [Ve[1], Vee[1]], [Ve[2], Vee[2]], 'r', lw=2)

ax.set_xlim3d(0., 1.)
ax.set_ylim3d(0., 1.)
ax.set_zlim3d(0., 1.)
            
plt.show()
