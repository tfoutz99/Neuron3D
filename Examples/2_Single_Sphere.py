#!/usr/bin/env ipython -wthread
# Run from neuron3d directory
from enthought.mayavi.mlab import *

#===============================================================================
# Create a mayavi figure
#===============================================================================
fig = figure(bgcolor=(0.1, 0.1, 0.1))

#===============================================================================
# Show a single sphere
#===============================================================================

# Define XYZ coordinates
x = [0] # Python List with 1 Entry
y = [0]
z = [0]

# Mayavi's points3d function
pts = points3d(x, y, z,
               color=(0.3, 1.0, 0.3))
pts.glyph.glyph_source.glyph_source.theta_resolution = 50
pts.glyph.glyph_source.glyph_source.phi_resolution = 50