# Neuron3D

3D visualization of [NEURON](https://nrn.readthedocs.io/) simulations with
[Mayavi](https://docs.enthought.com/mayavi/mayavi/): a multi-compartment
subthalamic nucleus (STN) projection neuron, a deep brain stimulation (DBS)
electrode, basal ganglia nuclei, and an olfactory bulb network. Simulated
membrane voltage is drawn on the 3D morphology as the simulation runs, and
frames can be saved to make movies. Presented at SfN 2010
([talk slides](docs/SFN_2010_Talk.pdf)).

## Setup

Requires Python 3.10+ and a C compiler for NEURON's mechanisms (on macOS,
the Xcode Command Line Tools).

```sh
uv sync                  # or: pip install neuron mayavi PyQt5 numpy
uv run nrnivmodl mod     # compile the ion channel mechanisms (creates arm64/ or x86_64/)
```

If `nrnivmodl` fails with *unable to locate the Python shared library*
(seen with uv-managed Pythons), point it at libpython and rerun:

```sh
export NMODL_PYLIB=$(uv run python -c "import sys, sysconfig; print(sys.base_prefix + '/lib/' + sysconfig.get_config_var('LDLIBRARY'))")
```

## Examples

```sh
uv run python main.py       # interactive menu
uv run python main.py 7     # run example 7, then show the menu
```

| # | Example | Shows |
|---|---------|-------|
| 1 | Spherical harmonics | Mayavi sanity check |
| 2 | Single sphere | Basic `points3d` |
| 3 | Electrode | Electrode shaft and tip |
| 4 | DBS electrode | Four-contact DBS lead with a panel for contact voltages |
| 5 | Nuclei | Globus pallidus pars interna and STN surface meshes |
| 6 | Subthalamic neuron | STN neuron with buttons to display it and run a simulation |
| 7 | All together | DBS lead, STN, and STN neuron in one scene |
| 8 | Network | Olfactory bulb network activity from recorded spike times |

Close an example's windows to return to the menu. To keep a Python prompt
while the windows are open, use IPython with Qt event loop integration from
the project directory:

```sh
uv run --with ipython ipython --gui=qt
In [1]: %run Examples/7_All_Together.py
```

## Making movies

Check *Save images* in a panel before running a simulation. Frames are
written to `img/`; turn them into a movie with ffmpeg:

```sh
cd img
./makecellmovie cell cell     # cell_000.png, ... -> cell.mov
./makenetworkmovie            # Network_000.png, ... -> network.mov
```

## Tests

```sh
uv run pytest
```

The tests need the compiled mechanisms and render offscreen, so no windows
open.

## Credits

- STN neuron model (`stn.hoc`): Svjetlana Miocinovic, 2004. Sources for each
  ion channel mechanism are cited in its file in `mod/`.
- Olfactory bulb network spike data: Migliore olfactory bulb model.

Licensed under the GNU GPL v3; see [LICENSE.txt](LICENSE.txt).
