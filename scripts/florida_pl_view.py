# %%
# n_core = 1.45213 and n_clad = 1.44692 for SMF28 https://iopscience.iop.org/article/10.1088/1742-6596/1655/1/012160/pdf
# sqrt(1.44692^2 - n_jacket^2) = 0.13 (design NA from Caleb) -> n_jacket = 1.4410
# other numbers set by Caleb 
# input SMFs radii are about 9um
# pretend all uniform cladding, but I should be close to that anyway
# 24um 

import numpy as np
from espeon.design import save_lantern_design, ngon_pattern
from itertools import product

pos = ngon_pattern(6, 3, 125 / (633.6/67.22), theta_init=np.pi/12)

save_lantern_design(
    "florida_test_view", 
    port_positions=pos, 
    core_radius_um=9/2, cladding_radius_um=67.22/2, z_extent_um=50000, scale=633.6/67.22,
    n_clad=1.44692, n_core=1.45213, n_jacket=1.441, 
    wavelengths_um=[1.55],
    do_plot=True,
    force_overwrite=True,
    simulation_params = {
        "name" : "lightbeam",
        "mesh_extent_um" : 700,
        "mesh_spacing_um" : 1.0,
        "dz_um" : 50,
        "PML" : 32
    }
)

# %%
