#!/usr/bin/env ipython -wthread
from enthought.mayavi.mlab import *

#===============================================================================
# Create a mayavi figure
#===============================================================================
fig = figure(bgcolor=(0.1, 0.1, 0.1))

##==============================================================================
# Show multiple spheres
##==============================================================================
# x,y and z have to be the same length
# x,y and z are python lists with 8 Entries
x = [0,0,0,0,1,1,1,1]
y = [0,0,1,1,0,0,1,1]
z = [0,1,0,1,0,1,0,1]
pts = points3d(x, y, z,
               color=(0.3, 1.0, 0.3))

