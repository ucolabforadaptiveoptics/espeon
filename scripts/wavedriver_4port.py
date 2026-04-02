# %%
import numpy as np
from espeon.design import setup_lantern, save_lantern_design, ngon_pattern

n_clad = 1.4584 # fused silica
numerical_aperture = 0.117
n_core = np.sqrt(n_clad ** 2 + numerical_aperture ** 2)
n_jacket = np.sqrt(n_clad ** 2 - 0.125**2)

save_lantern_design(
    "wavedriver_4port_secondlook_multiwl",
    port_positions=ngon_pattern(3, 2, 27.4/7.14, theta_init=np.pi/12), 
    core_radius_um=np.repeat(1.77, 4), cladding_radius_um=45.0/7.14, z_extent_um=10_000, scale=7.14,
    n_clad=n_clad, n_core=n_core, n_jacket=n_jacket, 
    wavelengths_um=np.array([0.5025, 0.5238, 0.5450, 0.5663, 0.5875]),
    do_plot=True,
        simulation_params = {
        "name" : "lightbeam",
        "mesh_extent_um" : 150,
        "mesh_spacing_um" : 0.1,
        "dz_um" : 10,
        "PML" : 8
    }
)
# %%
