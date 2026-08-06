# %%
import numpy as np
from matplotlib import pyplot as plt
from espeon.design import save_lantern_design, ngon_pattern, setup_lantern

n_clad = 1.4584 
numerical_aperture = 0.117
n_core = np.sqrt(n_clad ** 2 + numerical_aperture ** 2)
n_jacket = np.sqrt(n_clad ** 2 - 0.125**2)

j = jr_radius_um = 4.04 / 2
m = ma_radius_um = 4.60 / 2
p = pa_radius_um = 5.22 / 2
core_radii_livermore_b = [j, p, m, p, m, p, m, m, j, p, j, m, j, p, j, m, j, p, j]

save_lantern_design(
    design_name="livermore_b_5over4",
    port_positions=ngon_pattern(6, 3, 28.4 / 7.14, theta_init=np.pi/12), 
    core_radius_um=core_radii_livermore_b, cladding_radius_um=65.0*5/4/7.14,z_extent_um=10_000, scale=7.14,
    n_clad=n_clad, n_core=n_core, n_jacket=n_jacket, 
    wavelengths_um=[0.633],
    simulation_params = {
		"name" : "lightbeam",
		"mesh_extent_um" : 175,
		"mesh_spacing_um" : 0.2,
		"dz_um" : 25,
		"PML" : 8
	},
)
# %%
