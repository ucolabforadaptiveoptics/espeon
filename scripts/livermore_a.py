# %%
import numpy as np
from espeon.design import save_lantern_design, ngon_pattern

n_clad = 1.4584 # fused silica
numerical_aperture = 0.117
n_core = np.sqrt(n_clad ** 2 + numerical_aperture ** 2)
n_jacket = np.sqrt(n_clad ** 2 - 0.125**2)

save_lantern_design(
    design_name="livermore_a_t",
    port_positions=ngon_pattern(6, 3, 27.4 / 7.14, theta_init=np.pi/12), 
    core_radius_um=np.repeat(3.11/2, 19), cladding_radius_um=62.7/7.14, z_extent_um=10_000, scale=7.14,
    n_clad=n_clad, n_core=n_core, n_jacket=n_jacket, 
    wavelengths_um=[0.633],
    simulation_params = {
		"name" : "lightbeam",
		"mesh_extent_um" : 150,
		"mesh_spacing_um" : 0.1,
		"dz_um" : 25,
		"PML" : 8
	}
)
# %%
from espeon.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("livermore_a_rotated")
# %%
