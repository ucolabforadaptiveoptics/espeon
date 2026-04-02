# %%
import numpy as np
from espeon.design import setup_lantern, save_lantern_design, ngon_pattern
# %%
save_lantern_design(
    design_name="ames",
    port_positions=ngon_pattern(6, 2, 18.554/8, theta_init=np.pi/12), 
    core_radius_um=np.repeat(4.0, 7), cladding_radius_um=55.6/8, z_extent_um=10_000, scale=8,
    n_clad=1.44342, n_core=1.4512, n_jacket=1.4406, 
    wavelengths_um=[0.633, 1.550],
    simulation_params = {
		"name" : "lightbeam",
		"mesh_extent_um" : 150,
		"mesh_spacing_um" : 0.2,
		"dz_um" : 25,
		"PML" : 8
	}
)
# %%
