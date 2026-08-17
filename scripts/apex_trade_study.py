# hmm I need to first run a high res PL sim
# R~300 at lambda = 1550nm
# R = lambda / dlambda
# dlambda = 1550 / 300 ~ 5nm
# over like, 300nm
# so that's 1300-1600?

# %%
import numpy as np
from espeon.design_cbeam import save_lantern_design_cbeam
from espeon.utils import ngon_pattern
from cbeam.waveguide import PhotonicLantern
from cbeam.propagator import Propagator
import hcipy
# %%

save_lantern_design_cbeam(
    design_name="kim-6port-260730",
    port_positions=ngon_pattern(5, 2, 10*2/3), 
    core_radius_um=np.repeat(2.2, 6), cladding_radius_um=10, jacket_radius_um=30, z_extent_um=10_000, scale=8,
    n_clad=1.444, n_core=1.444+8.8e-3, n_jacket=1.444-5.5e-3, 
    wavelengths_um=np.arange(1.3, 1.6, 0.0025),
    do_plot=False,
    force_overwrite=True
)
# %%
