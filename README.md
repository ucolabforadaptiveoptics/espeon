# espeon

Photonic lantern simulation code and linking with _hcipy_. 

For practical reasons I've gotten rid of the git history in this repo - if you need anything in there it's under `umbreon`.

Instructions, assuming you have python and pip:

1. Clone this repository from GitHub by running `git clone https://github.com/ucolabforadaptiveoptics/espeon`, and open a terminal window at the folder that you get.
2. (optional, but recommended) Create a virtual environment with `conda` or by running `python3 -m venv .venv` + `source .venv/bin/activate`; I currently use `uv` (https://pypi.org/project/uv/) but go with whatever's easiest.
3. Install all the requirements with `pip install -r requirements.txt`.
4. Optionally, install the _lightbeam_ requirements with `pip install -r requirements_lb.txt`. This is necessary to make your own PL designs, but if you don't have this, you can work with design files someone else has run. (I'm using my own version of lightbeam for this, but Jon's should still work.)
4. Install the local package with `pip install -e .`
5. Confirm that you can run `scripts/spie24_image.py` correctly.
