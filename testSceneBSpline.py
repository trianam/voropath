#!/bin/python

import numpy as np
import plotter

points = np.array([[3 , 1, 0], [2.5, 4, 0], [0, 1, 0], [-2.5, 4, 0],[-3, 0, 0], [-2.5, -4, 0], [0, -1, 0], [2.5, -4, 0], [3, -1, 0]])

plt = plotter.Plotter()

plt.addPoints(points, plt.COLOR_CONTROL_POINTS, thick=True, thickness=0.02)
plt.addPolyLine(points, plt.COLOR_CONTROL_POLIG, thick=True, thickness=0.01)
plt.addBSpline(points, 3, plt.COLOR_PATH, thick=True, thickness=0.01)
#plt.addBSplineDEPRECATED(points, 3, plt.COLOR_PATH, thick=True, thickness=0.01)

plt.draw()
