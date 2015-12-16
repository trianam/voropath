#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import scipy.spatial

numSites = 10

fig = plt.figure()

np.random.seed(0)
sites = sp.rand(numSites,2)

plt.plot(sites[:,0], sites[:,1], 'b.')

vor = sp.spatial.Voronoi(sites)
sp.spatial.voronoi_plot_2d(vor)

plt.show()
