#!/usr/bin/env ipython -wthread
from enthought.mayavi.mlab import *
from neuron import h as nrn
from neuron import gui
from neuron3d import *

#===============================================================================
# Mayavi2 and Enthought
#===============================================================================
#   - Installation issues on mac/windows/linux:
#       - Make a guide
#   - Licenses:
#       - Acadamic is free (Now 64-bit)
#   - IPython:
#       - Syntax highlighting
#       - Access to help with ? and ??
#       - -wthread switch
#           - IPython starts running a separate thread for the graphical
#               toolkit's operation
#           - Allows opening and control of graphical elements
#               from within IPython without blocking
#       - -pylab switch
#           - Like "from pylab import *"
#           - Special support for matplotlib
#===============================================================================

output = module.submodule(parameter = (1.0, 1.0, 1.0))

(1.0, 1.0, 1.0) is a tuple, RGB

parameter names

modules

#from neuron import h as nrn
# nrn is an interface to NEURON
#nrn.tstop
#nrn.dt
#nrn.<tab>

from neuron import h as nrn
from neuron import gui # What is this?

I want to build a nice movie about a cell firing, with nuclei, a DBS electrode
 - Tools: 
        - Cpaste
 - Simple objects
 - DBS Electrode
 - 
I want to build Network firing movie

Describe tools

Sacrifice details for conceptual clarity
- let's not worry about the complexity of inner workings of code

Talk about how to keep Mayavi and NEURON working together
#===============================================================================
# How to create a mayavi figure
#===============================================================================
fig = figure(bgcolor = (1,1,1))

#===============================================================================
# Show a single sphere
#===============================================================================
x = [0] # Lists
y = [0]
z = [0]
pts = points3d(x, y, z,
               color=(0.3, 1.0, 0.3))

##==============================================================================
# Show multiple spheres
##==============================================================================
fig = figure(bgcolor = (0.1, 0.1, 0.1))
x = [0,0,0,0,1,1,1,1] # Lists
y = [0,0,1,1,0,0,1,1] # They have to be the same length
z = [0,1,0,1,0,1,0,1] # (XYZ in 3d)
pts = points3d(x, y, z,
               color=(0.3, 1.0, 0.3))

##==============================================================================
# Demonstrate Mayavi tools
##==============================================================================
#  - View the mayavi pipeline
#      - Change opacity:
#          - Glyph: Actor: Actor: Opacity: 1.0 --> 0.5
#          - Show overlap
#          - Change back to 1.0
#  - Record changes in order to add to script
#      - Change resolution
#          - Glyph: Glyph: Glyph Source:
#              - Phi Resolution 8   --> 20
#              - Theta resolution 8 --> 20
fig = figure(bgcolor = (0.1, 0.1, 0.1))
x = [0,0,0,0,1,1,1,1]
y = [0,0,1,1,0,0,1,1]
z = [0,1,0,1,0,1,0,1]
pts = points3d(x, y, z,
               color=(0.3, 1.0, 0.3))
pts.glypnrn.glyph_source.glyph_source.phi_resolution = 20
pts.glypnrn.glyph_source.glyph_source.theta_resolution = 20

##==============================================================================
# Simple Electrode
##==============================================================================
fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1),
                      size = (1920, 1080))
electrode = Electrode()
electrode.display()
# View pipeline, demonstrate that it is just a cylinder and a sphere
# Electrode?
# Electrode.<tab>
# Electrode.display?
# Electrode.display??

##==============================================================================
# Deep Brain Stimulation
##==============================================================================
fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1),
                      size = (1920, 1080))
dbs = DBS()
dbs.display()
#dbs?
dbs.update_contacts([-1,0,0,1])


##==============================================================================
# Nuclei
##==============================================================================
# Triangular mesh, give 3 points in 3D, get a plane
# pts, tris
# Add a wiremesh
fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1),
                      size = (1920, 1080))
# STN
stn = Nucleus(ptsFile='data/stn.pts',trisFile='data/stn.tri',
              name='Subthalamic Nucleus')
stn.display()

# GPi
gpi = Nucleus(ptsFile='data/gpi.pts',trisFile='data/gpi.tri',
              name='Globus Pallidus Pars Externa')
gpi.display()

# GPe
gpe = Nucleus(ptsFile='data/gpe.pts',trisFile='data/gpe.tri',
              name='Globus Pallidus Pars Externa')
gpe.display()

##==============================================================================
# Cell
##==============================================================================
fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1),
                      size = (1920, 1080))
nrn.load_file("j4a_mrg.hoc")

# Put thi in j4a_mrg.hoc
h('objref ic[3]')
h('soma ic[0] = new IClamp(0.5)')
h('soma ic[1] = new IClamp(0.5)')
h('soma ic[2] = new IClamp(0.5)')

# Look into this
stim = nrn.IClamp(sec = nrn.soma)

stim = nrn.IClamp(sec = nrn.soma)

# Work
soma = nrn.Section(name='soma')
stim = nrn.IClamp(soma(0.5))

# Doesn't work?
stim = nrn.IClamp(0.5, sec=soma)
stim = nrn.IClamp(seg=0.5, sec=soma)
stim = nrn.IClamp(loc=0.5, sec=soma)

# Write neuron in python
for ic in nrn.ic:
    ic.dur = 3
    ic.amp = 10
nrn.ic[0].delay = 0.5
nrn.ic[1].delay = 3.5
nrn.ic[2].delay = 6.5

# Describe CPaste?

# Have a slide that gives a top-level view of things
#  - Neuron + Mayavi (Python)
# Make a slide for each project
# ~ 6 main points

cell = Cell()
cell.display()
mlab.view(azimuth=88,
          elevation=58,
          distance=2000,
          focalpoint=array([-320,  75,  -60]),
          reset_roll=False)
# Show thresholding with mayavi tools (Mask), conductance

##==============================================================================
# Animation
##==============================================================================
# During a run (Using continuerun)
# Can't interact
#cell.play()
post-processing
##==============================================================================
# All together
##==============================================================================
fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1),
                      size = (1920, 1080))
dbs = DBS()
dbs.display()
nrn.load_file("j4a_mrg.hoc")
cell = Cell()
cell.display()
stn = Nucleus(ptsFile='data/stn.pts',trisFile='data/stn.tri',
              name='Subthalamic Nucleus')
stn.display()

##==============================================================================
# Network
##==============================================================================
fig = mlab.figure(bgcolor = (0.1, 0.1, 0.1),
                      size = (1920, 1080))
# olfactory bulb model - Migliore
network = Network('data/somas-large-203-10s-20pc-g3.spk')
network.display()
10000 green granule
500 mitral
just somas
flattened representation
##==============================================================================
# Animation
##==============================================================================
#network.play('Network')

# Summary

# Resources
