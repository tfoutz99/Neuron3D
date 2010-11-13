#!/usr/bin/env ipython -wthread
from enthought.mayavi.mlab import *
from neuron3d import *

#===============================================================================
# Create a mayavi figure
#===============================================================================
fig = figure(bgcolor=(0.1, 0.1, 0.1))

##==============================================================================
# Deep Brain Stimulation
##==============================================================================
dbs = DBS()
dbs.display()
#dbs?
#dbs.<tab>
#dbs??
dbs.edit_traits()
