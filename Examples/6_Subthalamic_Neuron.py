from mayavi.mlab import *
from neuron import h as nrn
from neuron import gui
from neuron3d import *

##==============================================================================
# Cell
##==============================================================================
load_model()
cell = Cell(name="Subthalamic Neuron")
cell.edit_traits()

show()
