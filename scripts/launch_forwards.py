# %%
import numpy as np
from espeon import PhotonicLanternOptics, input_to_2d
import lightbeam as lb

pl = PhotonicLanternOptics("wavedriver_4port_secondlook")
wavelength_um = pl.wavelengths_um[0] or 1.55
scale = pl.attributes["scale"]
n_core, n_clad, n_jacket = pl.attributes["n_core"], pl.attributes["n_clad"], pl.attributes["n_jacket"]
z_extent_um = pl.attributes["z_extent_um"]
PML = pl.attributes["sim_PML"]

mesh = lb.RectMesh3D(
    xw = pl.attributes["sim_mesh_extent_um"],
    yw = pl.attributes["sim_mesh_extent_um"],
    zw = z_extent_um,
    ds = pl.attributes["sim_mesh_spacing_um"],
    dz = pl.attributes["sim_dz_um"],
    PML = pl.attributes["sim_PML"]
)
lant = lb.optics.Lantern(
    pl.attributes["port_positions"], pl.attributes["core_radius_um"]/scale, pl.attributes["cladding_radius_um"], 0, z_extent_um, (n_core, n_clad, n_jacket), final_scale = scale
)

# %%
lant.set_sampling(mesh.xy)
xg, yg = mesh.grids_without_pml()
lbprop = lb.Prop3D(wavelength_um, mesh, lant, n_clad)
out_unpml = np.zeros_like(xg)
lant.set_IORsq(out_unpml, 0)
plt.imshow(out_unpml, vmin=n_jacket**2, vmax=n_core**2)
# %%

lant.set_IORsq(out_unpml, 0)
out = out_unpml[PML:-PML,PML:-PML]
input_footprint = np.where(out_unpml >= n_clad ** 2)
complement_mask = np.ones_like(out, dtype=bool)
for (i, j) in zip(input_footprint[0], input_footprint[1]):
    complement_mask[i,j] = False
extent = [
    [np.min(x), np.max(x)]
    for x in input_footprint
]
# %%
# pull the principal modes as the thing we're launching here
for k in range(pl.nports):
    launch_mode = np.zeros_like(out_unpml, dtype=np.complex64)
    launch_mode[input_footprint[0]-PML, input_footprint[1]-PML] = pl.output[0,k,:]
    plt.imshow(np.abs(lbprop.prop2end(launch_mode)) ** 2)
    plt.show()

# This works great for the 4-port PL
# Which makes me think that one is adiabatic/well designed
# and the Kim 6-port one isn't 
# (Well, it probably was, but we have an unconstrained parameter that's messing with this)
# %%
