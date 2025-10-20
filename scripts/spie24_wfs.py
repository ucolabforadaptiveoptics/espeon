# %%
import numpy as np
import hcipy

from plsim.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("spie24")

pupil_grid = hcipy.make_pupil_grid(256, 1.0)
aperture = hcipy.make_circular_aperture(1.0)(pupil_grid)
flat_wavefront_pupil = hcipy.Wavefront(aperture, wavelength=1.55e-6)
prop = hcipy.FraunhoferPropagator(pupil_grid, pl.focal_grid)
focal_wf = prop(flat_wavefront_pupil)
pl.coeffs(focal_wf)
pl.show_lantern_modes()
# %%
