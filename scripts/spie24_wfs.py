# %%
import numpy as np
import hcipy

from plsim.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("spie24")

pupil_grid = hcipy.make_pupil_grid(256, 1.0)
aperture = hcipy.make_circular_aperture(1.0)(pupil_grid)

# %%
