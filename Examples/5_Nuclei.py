from mayavi.mlab import *
from neuron3d import *

#===============================================================================
# Create a mayavi figure
#===============================================================================
fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1))

##==============================================================================
# Nuclei
##==============================================================================
gpi = Nucleus(ptsFile='data/gpi.pts',trisFile='data/gpi.tri',
              color=(0.7, 0.4, 0.4),
              name='Globus Pallidus Pars Interna')
gpi.display()

stn = Nucleus(ptsFile='data/stn.pts',trisFile='data/stn.tri',
              color=(0.4, 0.4, 0.7),
              name='Subthalamic Nucleus')
stn.display()

show()
