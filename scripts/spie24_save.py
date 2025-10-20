import numpy as np
from plsim.design import save_lantern_design, ngon_pattern

# Make the PL design from Aditya's 2024 SPIE paper (arxiv 2406.07771)

save_lantern_design(
    "spie24", 
    port_positions=ngon_pattern(6, 3, 7.4), 
    core_radius_um=2.2, cladding_radius_um=18.5, z_extent_um=60_000, scale=8,
    n_clad=1.4504, n_core=1.46076, n_jacket=np.sqrt(1.4504**2 - 0.125**2),
    wavelengths_um=[1.54, 1.55, 1.56],
    do_plot=True,
    force_overwrite=True
)

