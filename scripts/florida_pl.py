# %%
# 67.22 micron MMF cladding
# 15.64 micron MMF core (eyeballed off image)
# 6.5 micron SMF core
# 125 micron SMF core-to-core
# 10x taper factor (67.22 micron MMF cladding to 633.6 micron SMF jacket)
# 10cm taper distance (arbitrary)
# Slightly uneven packing, which I'll ignore for now but might 

import numpy as np
from espeon.design import save_lantern_design, ngon_pattern

save_lantern_design(
    "florida", 
    port_positions=ngon_pattern(6, 3, 2.8), 
    core_radius_um=6.5/2, cladding_radius_um=15.64, z_extent_um=10_000, scale=10,
    n_clad=1.444, n_core=1.444+8.8e-3, n_jacket=1.444-5.5e-3, 
    wavelengths_um=[1.55],
    do_plot=True,
    force_overwrite=True
)


# %%
