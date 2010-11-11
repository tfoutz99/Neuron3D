#!/usr/bin/env ipython -wthread
from enthought.mayavi.mlab import *
from neuron import h as nrn
from neuron import gui
from neuron3d import *

##==============================================================================
# Cell
##==============================================================================
nrn.load_file("stn.hoc")
cell = Cell(name="Subthalamic Neuron")
cell.edit_traits()
