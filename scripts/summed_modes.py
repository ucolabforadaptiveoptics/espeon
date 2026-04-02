# %%
import numpy as np
from espeon.simulation import PhotonicLanternOptics, input_to_2d

pl = PhotonicLanternOptics("livermore_a_rotated")
summed_mode = 0
for o in pl.output[0,:,:]:
   summed_mode += np.abs(input_to_2d(o, pl.input_footprint, pl.extent))

plt.imshow(summed_mode)
# %%
