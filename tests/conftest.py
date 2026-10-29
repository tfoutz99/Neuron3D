import os

os.environ.setdefault("ETS_TOOLKIT", "null")

import pytest
from mayavi import mlab

mlab.options.offscreen = True


@pytest.fixture(autouse=True)
def fresh_figure():
    """ Each test draws into its own offscreen figure """
    mlab.figure()
    yield
    mlab.close(all=True)


@pytest.fixture(scope="session")
def stn_model():
    import neuron3d
    neuron3d.load_model()
    return neuron3d.nrn


@pytest.fixture
def img_dir(tmp_path, monkeypatch):
    import neuron3d
    monkeypatch.setattr(neuron3d, "IMG_DIR", tmp_path)
    return tmp_path
