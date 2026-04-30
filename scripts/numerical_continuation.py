# %%
import numpy as np
import hcipy
from espeon.simulation import PhotonicLanternOptics
from scipy.optimize import root

pl = PhotonicLanternOptics("wavedriverlike_3port")
D = 1.0
wl_cen = 0.545e-6
pupil_grid = hcipy.make_pupil_grid(256, D)
aperture = hcipy.make_circular_aperture(D)(pupil_grid)
flat_wavefront_pupil = hcipy.Wavefront(aperture, wavelength=wl_cen)
prop = hcipy.FraunhoferPropagator(pupil_grid, pl.focal_grid, focal_length=10.0)
focal_wf = prop(flat_wavefront_pupil)
zernikes = hcipy.mode_basis.make_zernike_basis(2, D, pupil_grid, starting_mode=2)

def f_numcont(aberration):
    phase = sum([a * z for (a, z) in zip(aberration, zernikes)])
    wf = hcipy.Wavefront(aperture * np.exp(1j * phase), wavelength=wl_cen)
    focal_wf = prop(wf)
    pl_readout = pl.readout(focal_wf)
    return pl_readout / np.sum(pl_readout)

f_numcont([0.0, 0.2])
# %%
root(lambda x: f_numcont(x) - np.array([0.25, 0.25, 0.5]), x0=[-1,-1], method='lm')
# %%
