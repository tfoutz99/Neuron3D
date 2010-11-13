#!/usr/bin/env ipython -wthread
from enthought.mayavi.mlab import *
from neuron import h as nrn
from neuron import gui
from neuron3d import *

##==============================================================================
# Network
##==============================================================================
# olfactory bulb model - Migliore
network = Network('data/somas-large-203-10s-20pc-g3.spk')
network.edit_traits()