#!/usr/bin/env ipython -wthread
from enthought.mayavi.mlab import *
from neuron3d import *

#===============================================================================
# Create a mayavi figure
#===============================================================================
fig = figure(bgcolor=(0.1, 0.1, 0.1))

##==============================================================================
# Simple Electrode
##==============================================================================
electrode = Electrode()
electrode.display()