# %%
import numpy as np
from espeon.design_cbeam import save_lantern_design_cbeam
from espeon.utils import ngon_pattern
from cbeam.waveguide import PhotonicLantern
from cbeam.propagator import Propagator
import hcipy
# %%

port_positions = ngon_pattern(5, 2, 10*2/3)
nports = 6
core_radius_um = np.repeat(2.2, 6)
cladding_radius_um = 10
jacket_radius_um = 30
z_extent_um=10_000
scale=8
n_clad=1.444
n_core=1.444+8.8e-3
n_jacket=1.444-5.5e-3

lant = PhotonicLantern(port_positions, core_radius_um/scale, cladding_radius_um, jacket_radius_um, [n_core]*nports, n_clad, n_jacket, z_extent_um, scale)
prop = Propagator(1.55, lant, nports)
prop.characterize()

# %%

save_lantern_design_cbeam(
    design_name="kim-6port-260730",
    port_positions=ngon_pattern(5, 2, 10*2/3), 
    core_radius_um=np.repeat(2.2, 6), cladding_radius_um=10, jacket_radius_um=30, z_extent_um=10_000, scale=8,
    n_clad=1.444, n_core=1.444+8.8e-3, n_jacket=1.444-5.5e-3, 
    wavelengths_um=[1.55],
    do_plot=False,
    force_overwrite=True
)
# %%
from espeon.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("kim-6port-260729")
D = 1.0
wl_cen = 1.55e-6
pupil_grid = hcipy.make_pupil_grid(256, D)
aperture = hcipy.make_circular_aperture(D)(pupil_grid)
flat_wavefront_pupil = hcipy.Wavefront(aperture, wavelength=wl_cen)
prop = hcipy.FraunhoferPropagator(pupil_grid, pl.focal_grid, focal_length=6.0)
focal_wf = prop(flat_wavefront_pupil)
pl.readout(focal_wf)
#print(pl._backend, pl.readout(focal_wf).sum() / focal_wf.total_power)
# %%
