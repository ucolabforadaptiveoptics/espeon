# espeon

Photonic lantern simulation code and linking with `hcipy`. 

This package allows you to simulate new photonic lanterns using lightbeam (setting the design parameters, and running multiple wavelengths in one file) via the function `save_lantern_design`. This will generate an HDF5 file with the entrance modes and design parameters under `pl_designs/`. This file can then be loaded using the `PhotonicLanternOptics` class, which has methods to return per-port intensities and full images of the single-mode end in response to focal-plane images generated with `hcipy`.

For practical reasons I've gotten rid of the git history in this repo - if you need anything in there it's under `umbreon`. This is mostly the same as the earlier version of this package, `plsim`, but with multiple wavelengths in a single file and with much smaller files, as I'm no longer saving full-frame images of the single-mode end.

This package is pip installable with `pip install git+https://github.com/ucolabforadaptiveoptics/espeon`. Requirements for only getting intensities per port as a function of the PSF are in `requirements.txt`. Full-frame images additionally require `lightbeam`, for which the requirements are in `requirements_lb.txt`.

