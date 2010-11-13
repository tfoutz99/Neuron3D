#!/usr/bin/env ipython -wthread
from enthought.mayavi.mlab import *
from neuron import h as nrn
from neuron import gui
from neuron3d import *

##==============================================================================
# All together
##==============================================================================
all_together = All()
all_together.dbs.edit_traits()
all_together.edit_traits()