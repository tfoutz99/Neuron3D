""" Display NEURON models, DBS electrodes, and anatomical nuclei in 3D with Mayavi """
from pathlib import Path

import numpy
import neuron
from neuron import h as nrn
from mayavi import mlab
from mayavi.mlab import pipeline
from numpy import array, vstack, sqrt, dot
from functions import *
from traits.api import Float, HasTraits, Button, Bool
from traitsui.api import View, Item

MODULE_DIR = Path(__file__).resolve().parent
DATA_DIR = MODULE_DIR / 'data'
IMG_DIR = MODULE_DIR / 'img'

nrn.load_file('stdrun.hoc')

def load_model(hoc_file=MODULE_DIR / 'stn.hoc'):
    """ Load the compiled mechanisms in mod/ (if needed) and a hoc model """
    # NEURON loads mechanisms compiled in the working directory at import;
    # trainIClamp is defined in mod/train.mod
    if not hasattr(nrn, 'trainIClamp'):
        if not neuron.load_mechanisms(str(MODULE_DIR)):
            raise RuntimeError("Compiled mechanisms not found. "
                               "Run 'nrnivmodl mod' in %s" % MODULE_DIR)
    nrn.load_file(str(hoc_file))

#----------------------------------------------------------------------
class Electrode(HasTraits):
    """ Create an intracellular electrode
    """
    def __init__(self):
        super().__init__()
        # Electrode tip
        #self.sphere_center = [0, 0, 0]
        self.sphere_center = array([-7.34,
                              7.02,
                              -4.25])*1000 - array([0, 0, 2000])
        self.diameter = 1270

        # Electrode Shaft
        #self.direction = [0, 0, 1]
        self.top_electrode = [0, 0, 10000]
        self.shaft = array([self.sphere_center,
                            self.top_electrode])
        shaft_traj = self.top_electrode-self.sphere_center
        self.direction = shaft_traj/sqrt(dot(shaft_traj, shaft_traj))

    def display(self):
        """ Display Electrode in Mayavi """
        f = mlab.figure(figure=mlab.gcf(),
                        bgcolor=(0,0,0),
                        size=(800,600))

        # Build Electrode tip
        x, y, z = self.shaft[0,:].T
        self.mlab_tip = mlab.points3d(array(x),
                                      array(y),
                                      array(z),
                                      scale_factor=self.diameter,
                                      color=(0.5, 0.5, 0.5),
                                      name='Electrode Tip')
        self.mlab_tip.glyph.glyph_source.glyph_source.theta_resolution = 30
        self.mlab_tip.glyph.glyph_source.glyph_source.phi_resolution   = 16

        # Build Electrode Shaft
        x, y, z = self.shaft.T
        self.mlab_shaft = mlab.plot3d(x.T, y.T, z.T,
                               color=(0.5, 0.5, 0.5),
                               tube_radius=self.diameter / 2,
                               tube_sides=40,
                               name='Electrode')
        self.mlab_shaft.parent.parent.filter.capping=True

#----------------------------------------------------------------------
class DBS(Electrode):
    """ DBS Electrode Model """
    contact0 = Float(1, desc="Contact 0 Voltage")
    contact1 = Float(-0.5, desc="Contact 1 Voltage")
    contact2  = Float(0, desc="Contact 2 Voltage")
    contact3  = Float(-1, desc="Contact 3 Voltage")
    view = View( Item("contact0"),
                 Item("contact1"),
                 Item("contact2"),
                 Item("contact3"))

    # Initialize
    def __init__(self):
        """ Creates a DBS electrode with 4 contacts """

        # Call Super Class
        super().__init__()

        # Create Contacts in DBS electrode
        self.contact_height = 1000.0
        self.mlab_contacts = None

    def display(self):
        """ Display Electrode in Mayavi """
        f = mlab.figure(figure=mlab.gcf(),
                        bgcolor=(0,0,0),
                        size=(800,600))

        # Call Super Class
        Electrode.display(self)

        # Build electrode contacts
        c0 = find_new_point(self.sphere_center, self.direction, self.contact_height * 3 / 2)
        c1 = find_new_point(self.sphere_center, self.direction, self.contact_height * 7 / 2)
        c2 = find_new_point(self.sphere_center, self.direction, self.contact_height * 11 / 2)
        c3 = find_new_point(self.sphere_center, self.direction, self.contact_height * 15 / 2)

        contactCenters = array([c0, c1, c2, c3])
        contact = []
        for ii in range(4):
            topContact = find_new_point(contactCenters[ii,:], self.direction, self.contact_height / 2)
            bottomContact = find_new_point(contactCenters[ii,:], -self.direction, self.contact_height / 2)
            contact.append([bottomContact, topContact])

        contact_array = array(contact)
        x, y, z = contact_array.T

        activity_array = array([
            [self.contact0, self.contact0],
            [self.contact1, self.contact1],
            [self.contact2, self.contact2],
            [self.contact3, self.contact3]])
        src = pipeline.scalar_scatter(x.T, y.T, z.T, activity_array)
        src.mlab_source.dataset.lines = array([[0, 1], [2, 3], [4, 5], [6, 7]])
        src.name = 'Contacts'

        self.mlab_contacts = pipeline.tube(src,
                                                tube_radius=self.diameter * 1.03 / 2,
                                                tube_sides=40,
                                                name='DBS Contacts',)
        stripper = pipeline.stripper(self.mlab_contacts)
        surface = pipeline.surface(stripper)
        self.mlab_contacts.filter.capping = True
        surface.module_manager.scalar_lut_manager.lut_mode = 'RdBu'
        surface.module_manager.scalar_lut_manager.reverse_lut = True
        #surface.module_manager.scalar_lut_manager.use_default_range = False
        surface.module_manager.scalar_lut_manager.use_default_range = True
        surface.module_manager.scalar_lut_manager.data_range = array([-1., 1.])

    def _contact0_changed(self, value):
        self.set_contact(0, value)
    def _contact1_changed(self, value):
        self.set_contact(1, value)
    def _contact2_changed(self, value):
        self.set_contact(2, value)
    def _contact3_changed(self, value):
        self.set_contact(3, value)

    def set_contact(self, index, value):
        """ Recolor a contact once the electrode is displayed """
        if self.mlab_contacts is None:
            return
        self.mlab_contacts.mlab_source.scalars[2*index:2*index+2] = value
        self.mlab_contacts.mlab_source.update()

#----------------------------------------------------------------------
class Cell(HasTraits):
    """ multi-compartment neuron in hoc """
    play_dt = Float(0.5, desc="Time step for Mayavi simulations")
    play_tstop = Float(20, desc="Stopping time to update Mayavi")

    display_NEURON = Button("Display in NEURON", desc = "Display Neuron in NEURON")
    run_NEURON = Button("Run simulation (NEURON)", desc = "Run a neuron simulation")

    display_MayaVi = Button("Display in MayaVi", desc = "Display Neuron in Mayavi")
    save_img = Bool(0, desc="Save images to a file")
    run_mayavi = Button("Run simulation (NEURON+Mayavi)", desc = "Run a neuron simulation and update in mayavi")

    view = View(Item(name = 'play_dt'),
                Item(name = 'play_tstop'),
                Item(name = 'display_NEURON'),
                Item(name = 'run_NEURON'),
                Item(name = 'display_MayaVi'),
                Item(name = 'save_img'),
                Item(name = 'run_mayavi'))
    def _display_NEURON_fired(self):
        from neuron import gui
        nrn.xopen(str(MODULE_DIR / 'ses.ses'))
    def _display_MayaVi_fired(self):
        self.display()
    def _run_NEURON_fired(self):
        nrn.run()
    def _run_mayavi_fired(self):
        self.play()

    def __init__(self, name='Cell', var='v'):
        super().__init__()
        self.var = var
        self.name = name
        self.offset = array([0., 0., 0.])
        self.mlab_pts = None
        self.mlab_cell = None
        self.build_tree()

    def root_section(self):
        return nrn.SectionRef().root

    def sections(self):
        """ Sections of the tree, each parent before its children """
        def walk(sec):
            yield sec
            for child in nrn.SectionRef(sec=sec).child:
                yield from walk(child)
        return walk(self.root_section())

    def calculate_voltage(self):
        """ Calculate the voltage at this time """
        v = array([sec.v for sec in self.point_sections])
        self.xyzdv[:,4] = v
        return v

    def build_tree(self):
        """ Find [x, y, z, diam, v] at each 3d point and the edges between them """
        if self.var != 'v':
            raise AttributeError('Variable %s not implemented' % self.var)

        xyzdv = []
        connections = []
        self.point_sections = []  # section that each point belongs to
        point_ids = {}            # section -> xyzdv index of each of its 3d points
        for sec in self.sections():
            parent_seg = sec.parentseg()
            if parent_seg is None:
                xyzdv.append([sec.x3d(0), sec.y3d(0), sec.z3d(0), sec.diam3d(0), sec.v])
                self.point_sections.append(sec)
                ids = [0]
            else:
                # A child's first 3d point is where it attaches to its parent
                parent = parent_seg.sec
                arc = array([parent.arc3d(ii) for ii in range(int(parent.n3d()))])
                ids = [point_ids[parent][numpy.abs(arc - parent_seg.x * parent.L).argmin()]]
            for ii in range(1, int(sec.n3d())):
                xyzdv.append([sec.x3d(ii), sec.y3d(ii), sec.z3d(ii), sec.diam3d(ii), sec.v])
                self.point_sections.append(sec)
                ids.append(len(xyzdv)-1)
                connections.append([ids[-1], ids[-2]])
            point_ids[sec] = ids
        self.xyzdv = array(xyzdv)
        self.connections = array(connections)

    def display(self, var='v', scaling=1):
        ''' Display current cell in mayavi'''
        if self.mlab_pts is not None:
            try:
                self.mlab_pts.remove()
            except Exception:
                pass  # the figure was closed
        f = mlab.figure(figure=mlab.gcf(),
                        bgcolor=(0,0,0),
                        size=(800,600))

        ### Turn off vtk warnings # # # # # # # # # # # # # # # # # # # # # # #
        from vtk import vtkObject
        o = vtkObject
        o.GetGlobalWarningDisplay()
        o.SetGlobalWarningDisplay(0) # Turn it off.

        self.var = var
        self.build_tree()
        xs = self.xyzdv[:,0] + self.offset[0]
        ys = self.xyzdv[:,1] + self.offset[1]
        zs = self.xyzdv[:,2] + self.offset[2]
        diams = self.xyzdv[:,3] * scaling # larger scaling makes neurons more visible
        data = self.xyzdv[:,4]
        edges = self.connections

        # Display in mayavi
        pts = pipeline.scalar_scatter(xs, ys, zs, diams/2.0,
                                      name=self.name)
        dataset = pts.mlab_source.dataset
        dataset.point_data.get_array(0).name = 'diameter'
        dataset.lines = vstack(edges)

        array_id = dataset.point_data.add_array(data.T.ravel())
        dataset.point_data.get_array(array_id).name = 'data'
        dataset.point_data.update()

        #### Create tube with diameter data
        src = pipeline.set_active_attribute(pts,
                                            point_scalars='diameter')
        stripper = pipeline.stripper(src)
        tube = pipeline.tube(stripper,
                             tube_sides = 8,
                             tube_radius = 1)
        tube.filter.capping = True
        tube.filter.use_default_normal = False
        tube.filter.vary_radius = 'vary_radius_by_absolute_scalar'
        #tube.filter.radius_factor = 90.0 # just for making movies
        src2 = pipeline.set_active_attribute(tube, point_scalars='data')

        lines = pipeline.surface(src2)
        self.mlab_pts = pts
        self.mlab_cell = lines

        if self.var == 'v':
            self.mlab_cell.module_manager.scalar_lut_manager.use_default_range=False
            self.mlab_cell.module_manager.scalar_lut_manager.data_range = array([-70, 0])
        # Look at the root section (soma) from far enough to see the dendrites
        mlab.view(azimuth = 133.48567814586244,
                  elevation = 72.824281412546014,
                  distance = 2500,
                  focalpoint = self.xyzdv[0,:3] + self.offset)


    def move(self, xyz):
        """ Place the cell at an offset of xyz from its hoc coordinates """
        delta = array(xyz, dtype=float) - self.offset
        self.offset = self.offset + delta
        if self.mlab_cell is not None:
            source = self.mlab_cell.mlab_source
            source.x = source.x + delta[0]
            source.y = source.y + delta[1]
            source.z = source.z + delta[2]


    def retrieve_coordinates(self, sec):
        return [[sec.x3d(ii), sec.y3d(ii), sec.z3d(ii), sec.diam3d(ii)]
                for ii in range(int(sec.n3d()))]

    def play(self, fileroot='cell', show_colorbar=True, show_title=False,
             dt=None, tstop=None, save_img=None):
        ''' Step through cell response over time

        dt, tstop, and save_img default to the play_dt, play_tstop, and
        save_img traits. '''
        dt = self.play_dt if dt is None else dt
        tstop = self.play_tstop if tstop is None else tstop
        save_img = self.save_img if save_img is None else save_img
        nrn.cvode_active(0)
        nrn.stdinit()
        img_counter=0
        f = mlab.gcf()
        ftitle = None
        if show_colorbar: mlab.colorbar(self.mlab_cell)
        for x in range(int(round(tstop/dt)) + 1):
            timestamp = "TIME: %.1f" % (x*dt)
            print(timestamp)
            if show_title:
                if ftitle is None:
                    ftitle = mlab.title(timestamp)
                else:
                    ftitle.text = timestamp
            nrn.continuerun(x*dt)

            dataset = self.mlab_cell.mlab_source.dataset
            v = array(self.calculate_voltage())
            dataset.point_data.remove_array('data')
            array_id = dataset.point_data.add_array(v.T.ravel())
            dataset.point_data.get_array(array_id).name = 'data'
            dataset.point_data.update()
            self.mlab_cell.update_data()
            self.mlab_cell.update_pipeline()

            if save_img:
                f.scene.save_png(str(IMG_DIR / ('%s_%03d.png' % (fileroot, img_counter))))

            img_counter += 1

#----------------------------------------------------------------------
class Nucleus:
    def __init__(self, ptsFile=None, trisFile=None, color=(0.3, 0.3, 0.3), name=''):
        """ An anatomical nucleus """
        self.pts = None
        self.tris = None
        self.mlab_mesh = None
        self.color = color
        self.name = name
        self.ptsFile=ptsFile
        self.trisFile=trisFile
        if ptsFile and trisFile:
            self.create_pts_tris_from_file(ptsFile, trisFile)

    def create_pts_tris_from_file(self, ptsFile, trisFile):
        """ Load points (mm) and triangles (1-based point indices) from file """
        self.ptsFile = ptsFile
        self.trisFile = trisFile
        self.pts = numpy.loadtxt(ptsFile) * 1e3
        self.tris = numpy.loadtxt(trisFile, dtype=int) - 1

    def display(self):
        """ Display Nucleus in Mayavi """
        # Create and display a triangular mesh
        f = mlab.figure(figure=mlab.gcf(),
                        bgcolor=(0,0,0),
                        size=(800,600))
        self.mlab_mesh = mlab.triangular_mesh(self.pts[:, 0],
                                         self.pts[:, 1],
                                         self.pts[:, 2],
                                         self.tris,
                                         color=self.color,
                                         name=self.name)

        # Visual Tweaks
        #self.mlab_mesh.actor.property.backface_culling = True
        self.mlab_mesh.parent.parent.filter.splitting = False
    def set_opacity(self, opacity):
        if self.mlab_mesh is None:
            return
        self.mlab_mesh.actor.property.opacity = opacity
        if opacity:
            # I think transparent meshes look better with backface culling
            self.mlab_mesh.actor.property.backface_culling = True

#----------------------------------------------------------------------
class All(HasTraits):
    play_dt = Float(0.5, desc="Time step for Mayavi simulations")
    play_tstop = Float(20, desc="Stopping time to update Mayavi")

    display_dbs = Button("Display DBS", desc="Display DBS Electrode")
    display_nucleus = Button("Display STN Nucleus", desc="Display Subthalamic Nucleus")
    nucleus_opacity = Float(1.0, desc="Nucleus Opacity")
    display_neuron = Button("Display STN Neuron", desc="Display Subthalamic Nucleus Projection Neuron")
    x_neuron = Float(0, desc="Soma X-coordinate")
    y_neuron = Float(0, desc="Soma Y-coordinate")
    z_neuron = Float(-2000, desc="Soma Z-coordinate")
    diam_neuron = Float(1, desc="Relative diameter")

    run_NEURON = Button("Run simulation (NEURON)", desc = "Run a neuron simulation")

    save_img = Bool(0, desc="Save images to a file")
    run_mayavi = Button("Run simulation (NEURON+Mayavi)", desc = "Run a neuron simulation and update in mayavi")

    view = View(Item(name = 'play_dt'),
                Item(name = 'play_tstop'),
                Item(name = 'display_dbs'),
                Item(name = 'display_nucleus'),
                Item(name = 'nucleus_opacity'),
                Item(name = 'display_neuron'),
                Item(name = 'x_neuron'),
                Item(name = 'y_neuron'),
                Item(name = 'z_neuron'),
                Item(name = 'diam_neuron'),
                Item(name = 'run_NEURON'),
                Item(name = 'save_img'),
                Item(name = 'run_mayavi'))

    def _display_dbs_fired(self):
        self.dbs.display()
        self.reset_view()
    def _display_nucleus_fired(self):
        self.nucleus.display()
        self.nucleus_opacity=0.5
        self.reset_view()
    def _nucleus_opacity_changed(self, value):
        self.nucleus.set_opacity(value)

    def _display_neuron_fired(self):
        self.neuron.display(scaling=self.diam_neuron)
        self.reset_view()

    def _x_neuron_changed(self):
        self.move_neuron()
    def _y_neuron_changed(self):
        self.move_neuron()
    def _z_neuron_changed(self):
        self.move_neuron()
    def _diam_neuron_changed(self, value):
        self.neuron.display(scaling=value)
        self.reset_view()

    def _run_NEURON_fired(self):
        nrn.run()

    def _run_mayavi_fired(self):
        self.neuron.play(dt=self.play_dt,
                         tstop=self.play_tstop,
                         save_img=self.save_img)

    def __init__(self):
        super().__init__()
        self.dbs = DBS()
        self.nucleus = Nucleus(ptsFile=DATA_DIR / 'stn.pts',
                               trisFile=DATA_DIR / 'stn.tri',
                               name='Subthalamic Nucleus',
                               color= (0.75, 0.75, 0.75))
        load_model()
        self.neuron = Cell(name='stn')
        self.move_neuron()
        self.azimuth = 36
        self.elevation = 57
        self.distance = 6103
        self.focalpoint = array([-4480, 5135,  -4373])
        self.reset_roll = False

    def move_neuron(self):
        self.neuron.move([self.x_neuron,
                          self.y_neuron,
                          self.z_neuron])

    def reset_view(self):
        mlab.view(azimuth=self.azimuth,
                  elevation=self.elevation,
                  distance=self.distance,
                  focalpoint=self.focalpoint,
                  reset_roll=self.reset_roll)

#----------------------------------------------------------------------
class Network(HasTraits):
    t = Float(0, desc="Time")
    sim_start = Float(900, desc="Start Time")
    sim_end = Float(1000, desc="End Time")
    display_mayavi=Button("Display Network", desc="Display Network in Mayavi")
    save_img = Bool(0, desc="Save images to a file")
    run_mayavi = Button("Run simulation (NEURON+Mayavi)", desc = "Run a neuron simulation and update in mayavi")
    view = View(Item(name = 'display_mayavi'),
                Item(name = 'sim_start'),
                Item(name = 'sim_end'),
                Item(name = 't'),
                Item(name = 'save_img'),
                Item(name = 'run_mayavi'))
    def _display_mayavi_fired(self):
        self.display()
    def _run_mayavi_fired(self):
        self.play()
    #---------------------------------------------------------------------------
    def __init__(self, spikefile=DATA_DIR / 'somas-large-203-10s-20pc-g3.spk'):
    #--------------------------------------------------------------
        super().__init__()
        self.spikefile=spikefile
        self.nummit = 500
        self.numgcspermit = 20
        self.numgran = self.nummit * self.numgcspermit

        self.num_arcs = self.nummit
        self.mcsize = 1./75.
        self.gcsize = 1./float(self.numgcspermit*2*3)

        self.mcsize *= 2
        self.gcsize *= 2
        self.slice_height = .1
        self.mspikes, self.gspikes = self.read_file_to_vector(self.spikefile)
        self.create_cell_coordinates()

    #--------------------------------------------------------------
    def read_file_to_vector(self, file_name, spike_ids=None):
        """
        Read the spike data formatted from NEURON where each line is a timestamp
        of spike followed by a cell (or other spike generator) id that gave the
        spike.

        :param file_name: Name of the spike file to read.
        :param spike_ids: If specified, a subset of ids to load.

        :return: Mitral and granule cell spikes, each an array of (time, id)
            rows. Granule cell ids are renumbered to start from 0.
        """
        spikes = numpy.loadtxt(file_name, ndmin=2)
        if spike_ids is not None:
            spikes = spikes[numpy.isin(spikes[:,1], spike_ids)]
        mspikes = spikes[spikes[:,1] < self.nummit] # I assume this is how you ordered things?
        gspikes = spikes[spikes[:,1] >= self.nummit] # I assume this is how you ordered things?
        gspikes[:,1] -= self.nummit

        # Mitral and Granule Cell Spike times
        self.mspikes=mspikes
        self.gspikes=gspikes
        return mspikes, gspikes

    #---------------------------------------------------------------------------
    def create_cell_coordinates(self):
        def create_granule_cells(a=1, b=1, h=1, num_arcs=6, numgcspermit=2, \
                gcsize=None, color=None):
            if gcsize is None:
                gcsize = 1./float(numgcspermit*2*5)
            if color is None:
                color = (0, 1./3., 0)
            num_rings = numgcspermit + 1
            incr = 2*numpy.pi / self.num_arcs
            u = numpy.linspace(incr/2.,2*numpy.pi-incr/2., self.num_arcs)
            x = []; y = []; z = []
            for uu in u:
                for i in numpy.arange(3,(num_rings*2),2):
                    ratio = numpy.sqrt(float(i)/float(num_rings*2))
                    angle_noise = numpy.random.uniform(-incr/2.1, incr/2.1)
                    mag_noise = numpy.random.uniform(
                            numpy.sqrt(float(i-1)/float(num_rings*2)),
                            numpy.sqrt(float(i+1)/float(num_rings*2)))
                    x.append(numpy.sin(uu+angle_noise)*mag_noise*a)
                    y.append(numpy.cos(uu+angle_noise)*mag_noise*b)
                    z.append(numpy.random.uniform(0, h))
            #gcs = mlab.points3d(x,y,z, color=color, scale_factor = gcsize)
            gx=array(x)
            gy=array(y)
            gz=array(z)
            return gx, gy, gz
        #---------------------------------------------------------------------------
        def create_mitral_cells(a=1, b=1, h=1, mcsize=None, color=None):
            if mcsize is None:
                mcsize = 1./75.
            if color is None:
                color = (0, 0., 0.5)
            incr = 2 * numpy.pi / self.num_arcs
            u = numpy.linspace(incr/2., 2 * numpy.pi - incr/2., self.num_arcs)
            x = []; y = []; z = []
            angle_noise = numpy.random.uniform(-incr/2.1, incr/2.1, len(u))
            x.extend(numpy.sin(u + angle_noise) * a * 1.02) # Add 2% to separate from GCs
            y.extend(numpy.cos(u + angle_noise) * b * 1.02)
            z.extend(numpy.random.uniform(0, h, len(u)))

            #mcs = mlab.points3d(x,y,z, color=color, scale_factor = mcsize)
            mx=array(x)
            my=array(y)
            mz=array(z)
            return [mx, my, mz]
        self.gx, self.gy, self.gz = create_granule_cells(.5,
                                                         h = self.slice_height,
                                                         num_arcs = self.num_arcs,
                                                         numgcspermit = self.numgcspermit,
                                                         gcsize = self.gcsize)
        self.mx, self.my, self.mz = create_mitral_cells(.5,
                                                        h = self.slice_height,
                                                        mcsize = self.mcsize)

    #---------------------------------------------------------------------------
    def plot_points(self, x, y, z, spikes=[], t=0, color=(0.5, 0.5, 0.5), csize=1./75.):
        """ Plot spheres representing each node in network """
        if t==0:
            pts = mlab.points3d(x,y,z,
                                color=color,
                                scale_factor = csize)
        else:
            ind = spikes[spikes[:,0]==t, 1].astype(int)
            #ind = ind[ind<=numgran] # looks like there are some indices greater than numgran
            pts = mlab.points3d(x[ind],
                                y[ind],
                                z[ind],
                                color=color,
                                scale_factor = csize)
        pts.glyph.scale_mode = 'data_scaling_off'
        return pts

    #---------------------------------------------------------------------------
    def display(self):
        """ Plot background field """
        f = mlab.figure(figure=mlab.gcf(),
                bgcolor=(0,0,0),
                size=(800,600))

        self.mpts0 = self.plot_points(self.mx,
                                      self.my,
                                      self.mz,
                                      color=(0.1, 0.1, 0.3),
                                      csize=self.mcsize)
        self.gpts0 = self.plot_points(self.gx,
                                      self.gy,
                                      self.gz,
                                      color=(0.1, 0.3, 0.1),
                                      csize=self.gcsize)
        self.mpts0.parent.parent.scene.z_plus_view()

    #---------------------------------------------------------------------------
    def play(self, fileroot='Network', mspikes=None, gspikes=None, save_png=None,
             sim_step=10, windowsize=10):
        view = mlab.view()
        f = mlab.gcf()
        if mspikes is None:
            mspikes = self.mspikes
        if gspikes is None:
            gspikes = self.gspikes
        if save_png is None:
            save_png = self.save_img
        img_counter = 0
        ftitle = None
        ts = numpy.unique(gspikes[:,0]) # can use either spike set
        ts = ts[(ts > self.sim_start) & (ts < self.sim_end)]
        mqueue = []; gqueue = [];
        for t in ts[::sim_step]:
            self.t = t
            mlab.gcf().scene.disable_render=True
            timestamp = u"Time: %.1f" % (t)
            print(timestamp)
            if save_png:
                # Display time stamp
                f.name = timestamp
                if ftitle is None:
                    ftitle = mlab.title(timestamp)
                else:
                    ftitle.text = timestamp

            # Delete old spheres
            if len(mqueue) >= windowsize:
                mpts=mqueue.pop(0)
                mpts.parent.parent.remove()
                gpts=gqueue.pop(0)
                gpts.parent.parent.remove()

            # It would be great to make prevoius arrays dimmer

            # Plot activate spheres
            mqueue.append(self.plot_points(self.mx,
                                           self.my,
                                           self.mz,
                                           mspikes,
                                           t=t,
                                           color=(1., 1., 1.),
                                           csize = self.mcsize))
            gqueue.append(self.plot_points(self.gx,
                                           self.gy,
                                           self.gz,
                                           gspikes,
                                           t=t,
                                           color=(1., 1., 1.),
                                           csize = self.gcsize))
            mlab.view(view [0], view [1], view [2], view [3])
            mlab.gcf().scene.disable_render=False
            if save_png:
                f.scene.save_png(str(IMG_DIR / ('%s_%03d.png' % (fileroot, img_counter))))
                img_counter += 1
        return mqueue, gqueue

#----------------------------------------------------------------------
if __name__=='__main__':
    pass
    ## Electrode
    #electrode = Electrode()
    #electrode.display()

    ##Deep Brain Stimulation Electrode
    #dbs = DBS()
    #dbs.display()
    #dbs.edit_traits()

    ## Nucleus
    #stn = Nucleus(ptsFile='data/stn.pts',trisFile='data/stn.tri',
    #              name='Subthalamic Nucleus')
    #stn.display()
    #gpi = Nucleus(ptsFile='data/gpi.pts',trisFile='data/gpi.tri',
    #              name='Globus Pallidus Pars Interna')
    #gpi.display()
    #gpe = Nucleus(ptsFile='data/gpe.pts',trisFile='data/gpe.tri',
    #              name='Globus Pallidus Pars Externa')
    #gpe.display()

    ## Cell
    #nrn.load_file("stn.hoc")
    #cell = Cell()
    #cell.display()
    #cell.edit_traits()
    #mlab.view(azimuth=88,
              #elevation=58,
              #distance=2000,
              #focalpoint=array([-320,  75,  -60]),
              #reset_roll=False)
    #h('objref ic[3]')
    #h('soma ic[0] = new IClamp(0.5)')
    #h('soma ic[1] = new IClamp(0.5)')
    #h('soma ic[2] = new IClamp(0.5)')
    #for ic in nrn.ic:
        #ic.dur = 3
        #ic.amp = 10
    #nrn.ic[0].delay = 0.5
    #nrn.ic[1].delay = 3.5
    #nrn.ic[2].delay = 6.5
    #cell.play()


    ## Network
    #network = Network('data/somas-large-203-10s-20pc-g3.spk')
    #network.display()
    #network.play('Network')