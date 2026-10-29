import numpy as np
import pytest

import neuron3d
from neuron3d import DBS, All, Cell, Network, Nucleus
from functions import find_new_point

SPIKEFILE = neuron3d.DATA_DIR / "somas-large-203-10s-20pc-g3.spk"


def test_find_new_point():
    p = find_new_point([1, 2, 3], [0, 0, 2], 5)
    np.testing.assert_allclose(p, [1, 2, 8])


def test_nucleus_mesh():
    stn = Nucleus(ptsFile=neuron3d.DATA_DIR / "stn.pts", trisFile=neuron3d.DATA_DIR / "stn.tri")
    assert stn.pts.shape == (662, 3)
    assert stn.tris.shape == (1320, 3)
    assert stn.tris.min() == 0 and stn.tris.max() == len(stn.pts) - 1
    # mm in the file, um in the model
    np.testing.assert_allclose(stn.pts[0], [-7576.105, 7019.770, -3506.849])


def test_dbs_contacts_update_after_display():
    dbs = DBS()
    dbs.contact1 = 0.25   # before display: nothing to update yet
    dbs.display()
    dbs.contact3 = 0.75
    scalars = dbs.mlab_contacts.mlab_source.scalars.ravel()
    np.testing.assert_allclose(scalars[2:4], 0.25)
    np.testing.assert_allclose(scalars[6:8], 0.75)


#----------------------------------------------------------------------
# Network
@pytest.fixture(scope="module")
def network():
    return Network(SPIKEFILE)


def test_network_spike_times_match_file(network):
    raw = np.loadtxt(SPIKEFILE)
    granule = raw[raw[:, 1] >= network.nummit]
    np.testing.assert_array_equal(network.gspikes[:, 0], granule[:, 0])
    np.testing.assert_array_equal(network.gspikes[:, 1], granule[:, 1] - network.nummit)
    assert network.gspikes[:, 1].max() < network.numgran


def test_network_spike_ids_subset(network):
    mspikes, gspikes = network.read_file_to_vector(SPIKEFILE, spike_ids=[247, 600])
    assert set(mspikes[:, 1]) == {247}
    assert set(gspikes[:, 1]) == {600 - network.nummit}
    network.read_file_to_vector(SPIKEFILE)   # restore the full data set


def test_network_default_spikefile():
    assert len(Network().mspikes) > 0


def test_network_play(network, img_dir):
    network.display()
    network.sim_start, network.sim_end = 900, 905
    network.save_img = True
    mqueue, gqueue = network.play(mspikes=network.mspikes, gspikes=network.gspikes,
                                  sim_step=1, windowsize=3)
    assert len(mqueue) == len(gqueue) == 3
    frames = sorted(img_dir.glob("Network_*.png"))
    assert len(frames) == len(np.unique(network.gspikes[(network.gspikes[:, 0] > 900) &
                                                         (network.gspikes[:, 0] < 905), 0]))


#----------------------------------------------------------------------
# Cell
@pytest.fixture
def cell(stn_model):
    return Cell(name="stn")


def test_cell_tree(cell, stn_model):
    h = stn_model
    tree = [sec for sec in h.allsec() if sec.n3d() > 0]
    assert len(cell.xyzdv) == sum(int(sec.n3d()) for sec in tree) - (len(tree) - 1)
    assert len(cell.connections) == len(cell.xyzdv) - 1


def test_cell_children_connect_where_they_start(cell, stn_model):
    xyz = cell.xyzdv[:, :3]
    parent_of = dict(cell.connections)
    for sec in stn_model.allsec():
        if sec.parentseg() is None or sec.n3d() < 2:
            continue
        first = np.array([sec.x3d(1), sec.y3d(1), sec.z3d(1)])
        start = np.array([sec.x3d(0), sec.y3d(0), sec.z3d(0)])
        dist = np.linalg.norm(xyz - first, axis=1)
        assert dist.min() < 1e-6
        np.testing.assert_allclose(xyz[parent_of[dist.argmin()]], start, atol=1e-3,
                                   err_msg=sec.name())


def test_calculate_voltage_matches_tree_order(cell, stn_model):
    h = stn_model
    stim = h.IClamp(h.soma[1](0.5))
    stim.delay, stim.dur, stim.amp = 0, 5, 0.5
    h.finitialize(h.v_init)
    h.continuerun(3)
    v = cell.calculate_voltage().copy()
    cell.build_tree()
    np.testing.assert_array_equal(v, cell.xyzdv[:, 4])
    assert v.std() > 0


def test_cell_move_is_absolute(cell):
    cell.display()
    x0 = cell.mlab_cell.mlab_source.x.copy()
    cell.move([100, 0, 0])
    cell.move([100, 0, 0])
    np.testing.assert_allclose(cell.mlab_cell.mlab_source.x, x0 + 100)


def test_cell_move_before_display(cell):
    cell.move([0, 0, -2000])
    cell.display()
    np.testing.assert_allclose(cell.mlab_cell.mlab_source.z, cell.xyzdv[:, 2] - 2000)


def test_cell_play_reaches_tstop(cell, stn_model, img_dir):
    cell.display()
    cell.play_dt, cell.play_tstop, cell.save_img = 0.5, 2.0, True
    cell.play()
    assert stn_model.t == pytest.approx(2.0)
    assert len(list(img_dir.glob("cell_*.png"))) == 5


#----------------------------------------------------------------------
# All
def test_all_neuron_position_follows_panel(stn_model):
    scene = All()
    scene._display_neuron_fired()
    z0 = scene.neuron.xyzdv[:, 2]
    np.testing.assert_allclose(scene.neuron.mlab_cell.mlab_source.z, z0 + scene.z_neuron)
    scene.x_neuron = 50
    scene.z_neuron = -1000
    scene.diam_neuron = 2   # redraws the neuron
    src = scene.neuron.mlab_cell.mlab_source
    np.testing.assert_allclose(src.x, scene.neuron.xyzdv[:, 0] + 50)
    np.testing.assert_allclose(src.z, z0 - 1000)


def test_all_run_uses_panel_settings(stn_model, img_dir):
    scene = All()
    scene._display_neuron_fired()
    scene.play_dt, scene.play_tstop, scene.save_img = 1.0, 3.0, True
    scene._run_mayavi_fired()
    assert stn_model.t == pytest.approx(3.0)
    assert len(list(img_dir.glob("cell_*.png"))) == 4
