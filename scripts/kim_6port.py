# %%
import numpy as np
from espeon.design import save_lantern_design, ngon_pattern
from cbeam.waveguide import PhotonicLantern

save_lantern_design(
    design_name="kim-6port-260729",
    port_positions=ngon_pattern(5, 2, 10*2/3), 
    core_radius_um=np.repeat(2.2, 6), cladding_radius_um=10, z_extent_um=10_000, scale=8,
    n_clad=1.444, n_core=1.444+8.8e-3, n_jacket=1.444-5.5e-3, 
    wavelengths_um=[1.55],
    do_plot=False,
    simulation_params = {
		"name" : "lightbeam",
		"mesh_extent_um" : 350,
		"mesh_spacing_um" : 0.1,
		"dz_um" : 10,
		"PML" : 8
	},
    force_overwrite=True
)
# %%
from espeon.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("kim-6port")
# %%
