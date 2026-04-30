# %%
from espeon.simulation import PhotonicLanternOptics

pl = PhotonicLanternOptics("wavedriver_4port_secondlook_multiwl")
# %%
for k in range(5):
    pl.show_modes(k, nrows=2)
# %%
