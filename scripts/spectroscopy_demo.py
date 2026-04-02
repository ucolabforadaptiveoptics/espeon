# %%
import numpy as np
import hcipy

from espeon.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("spie24")
D_pupil = 1.0
um = 1e-6
D_focus = 52.5 * um
pupil_grid = hcipy.make_pupil_grid(256, D_pupil)
aperture = hcipy.make_circular_aperture(D_pupil)(pupil_grid)
prop = hcipy.FraunhoferPropagator(pupil_grid, pl.focal_grid, focal_length=7.5)
singlemode_fiber_core_radius = 3.1 * um
fiber_NA = 0.13
fiber_length = 10

import matplotlib.pyplot as plt

plt.figure(figsize=(7, 4.5))
for j in range(2):
    for (i, wl) in enumerate([1.54, 1.55, 1.56]):
        wl_cen = wl * um
        flat_wavefront_pupil = hcipy.Wavefront(aperture, wavelength=wl_cen)
        focal_wf = prop(flat_wavefront_pupil)

        single_mode_fiber = hcipy.StepIndexFiber(singlemode_fiber_core_radius, fiber_NA, fiber_length)
        if j == 1:
            focal_wf = single_mode_fiber.forward(focal_wf)
            t = "smf"
        else:
            t = "psf"

        plt.subplot(2, 3, j*3 + i + 1)
        pl.generate_launch_fields(scaleup=5) 
        plt.imshow(pl.image(focal_wf))
        plt.xticks([])
        plt.yticks([])
        plt.title(f"{(wl_cen / um):.2f} microns, {t}")

plt.subplots_adjust(hspace=0.6)
# %%
