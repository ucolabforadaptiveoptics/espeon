# %%
import numpy as np
import hcipy

from espeon.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("spie24")
D = 1.0
wl_cen = 1.55e-6
pupil_grid = hcipy.make_pupil_grid(256, D)
aperture = hcipy.make_circular_aperture(D)(pupil_grid)
flat_wavefront_pupil = hcipy.Wavefront(aperture, wavelength=wl_cen)
prop = hcipy.FraunhoferPropagator(pupil_grid, pl.focal_grid, focal_length=7.5)
focal_wf = prop(flat_wavefront_pupil)
# %%
pl.coeffs(focal_wf) # length-19 complex valued coefficients on each SMF port
# %%
pl.readout(focal_wf) # as above, but norm-squared to make intensities
# %%
pl.show_modes(1)
# %%
# if you have lightbeam:
pl.generate_launch_fields(scaleup=3) 
# "scaleup" makes the ports artificially larger for better visualization
pl.image(focal_wf) # "readout" but embedded in the full image

# %%
pl.show_image(focal_wf)
# %%
# looking at aberrated outputs
zernikes = hcipy.make_zernike_basis(3, D, pupil_grid, starting_mode=2)
for z in zernikes:
    pl.show_image(
        prop(hcipy.Wavefront(aperture * np.exp(1j * z), wavelength=wl_cen))
    )

# %%
