#!/bin/python

import sys
import numpy as np
import plotter

degree = 3
#degree = 4
#controlPolygon = np.array([[0.,0.,0.], [0.,6.,1.], [-1.,5.,2.], [-3.,8.,3.], [-1.,14.,4.],[2.,14.,5.],[4.,8.,6.], [2.,5.,7.], [1.,6.,8.], [1.,0.,9.]])
controlPolygon = np.array([[0.,0.,0.], [0.,6.,10.], [-1.,5.,2.], [-3.,8.,13.], [-1.,14.,4.],[2.,14.,15.],[4.,8.,6.], [2.,5.,17.], [1.,6.,8.], [1.,0.,19.]])


plt = plotter.Plotter()
plt.addBSpline(controlPolygon, degree, plt.COLOR_PATH, thick=True)    
plt.draw()

