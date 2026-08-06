# espeon

Photonic lantern simulation code and linking with `hcipy`. 

This package allows you to simulate new photonic lanterns using lightbeam (setting the design parameters, and running multiple wavelengths in one file) via the function `save_lantern_design`. This will generate an HDF5 file with the entrance modes and design parameters under `pl_designs/`. This file can then be loaded using the `PhotonicLanternOptics` class, which has methods to return per-port intensities and full images of the single-mode end in response to focal-plane images generated with `hcipy`.

For practical reasons I've gotten rid of the git history in this repo - if you need anything in there it's under `umbreon`. This is mostly the same as the earlier version of this package, `plsim`, but with multiple wavelengths in a single file and with much smaller files, as I'm no longer saving full-frame images of the single-mode end.

Install the core package with:

```sh
pip install git+https://github.com/ucolabforadaptiveoptics/espeon
```

Install only the simulation backend you need:

```sh
pip install "espeon[lightbeam] @ git+https://github.com/ucolabforadaptiveoptics/espeon"
pip install "espeon[cbeam] @ git+https://github.com/ucolabforadaptiveoptics/espeon"
```

For a source checkout, the equivalent requirement files are `requirements_lb.txt` and `requirements_cbeam.txt`. Backend imports are lazy, so Lightbeam and CBeam do not need to be installed together.

This repo includes code generated with GPT-5.6 Sol via Codex. My procedure for validating the outputs of AI-generated code is available upon request and the prompts used are under `prompts/`.