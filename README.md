# espeon

Photonic lantern simulation code and linking with `hcipy`. 

This package allows you to simulate new photonic lanterns using lightbeam (setting the design parameters, and running multiple wavelengths in one file) via the function `save_lantern_design`. This will generate an HDF5 file with the entrance modes and design parameters under `pl_designs/`. This file can then be loaded using the `PhotonicLanternOptics` class, which has methods to return per-port intensities and full images of the single-mode end in response to focal-plane images generated with `hcipy`.

For practical reasons I've gotten rid of the git history in this repo - if you need anything in there it's under `umbreon`. This is mostly the same as the earlier version of this package, `plsim`, but with multiple wavelengths in a single file and with much smaller files, as I'm no longer saving full-frame images of the single-mode end.

Instructions, assuming you have python and pip:

1. Clone this repository from GitHub by running `git clone https://github.com/ucolabforadaptiveoptics/espeon`, and open a terminal window at the folder that you get.
2. (optional, but recommended) Create a virtual environment with `conda` or by running `python3 -m venv .venv` + `source .venv/bin/activate`; I currently use `uv` (https://pypi.org/project/uv/) but go with whatever's easiest.
3. Install all the requirements with `pip install -r requirements.txt`.
4. Optionally, install the _lightbeam_ requirements with `pip install -r requirements_lb.txt`. This is necessary to make your own PL designs and to visualize the full images, but if you don't have this, you can work with design files someone else has run. (I'm using my own version of lightbeam for this, but Jon's should still work.)
4. Install the local package with `pip install -e .`
5. Confirm that you can run `scripts/spie24_image.py` correctly.
