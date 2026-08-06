"""Calculate and validate the LP-to-port transfer matrix of a saved lantern."""

from pathlib import Path

import hcipy as hc
import lightbeam as lb
import matplotlib.pyplot as plt
import numpy as np
from cbeam.propagator import Propagator, construct_B
from cbeam.waveguide import PhotonicLantern
from scipy.optimize import linear_sum_assignment

from espeon import PhotonicLanternOptics


tag = "kim-6port-260729"
pl = PhotonicLanternOptics(tag)

specs = ((0, 1, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0), (2, 1, 1), (0, 2, 0))

def make_lp_modes(x, y):
	"""Return the real LP01, LP11, LP21, and LP02 basis on arbitrary points."""
	phi, modes = np.arctan2(y, x), []
	for l, m, sine in specs:
		field = lb.lpfield(x, y, l, m, pl.attributes["cladding_radius_um"],
		                   pl.wavelengths_um[0], pl.attributes["n_clad"], pl.attributes["n_jacket"])
		field *= ([np.cos(l * phi), np.sin(l * phi)][sine]) if l else 1
		modes.append(field)
	return np.asarray(modes)


# Make the LP modes on the pixels saved by save_lantern_design.
i, j = pl.input_footprint
dx, width = pl.attributes["sim_mesh_spacing_um"], pl.attributes["sim_mesh_extent_um"]
x, y = (i - width / dx / 2) * dx, (j - width / dx / 2) * dx
lp_modes = make_lp_modes(x, y)
lp_modes /= np.linalg.norm(lp_modes, axis=1)[:, None]

# Both matrices below act on column vectors; transfer @ lp_coeffs = port_coeffs.
principal_projector = pl.projectors[0]
lp_projector = np.linalg.pinv(lp_modes.T)
transfer = principal_projector @ lp_modes.T

# Use a real focal-plane field, and compare with PhotonicLanternOptics.coeffs().
pupil_grid = hc.make_pupil_grid(256, 1)
pupil = hc.make_circular_aperture(1)(pupil_grid)
focal = hc.FraunhoferPropagator(pupil_grid, pl.focal_grid, focal_length=7.5)(
	hc.Wavefront(pupil, wavelength=pl.wavelengths_um[0] * 1e-6)
)
direct = pl.coeffs(focal)
via_transfer = transfer @ lp_projector @ focal.electric_field.shaped[pl.input_footprint]
projection_error = np.linalg.norm(via_transfer - direct) / np.linalg.norm(direct)

# Run CBeam directly, then express its instantaneous input eigenmodes in the
# same analytic LP basis used above.
n = pl.nports
a, scale = pl.attributes, pl.attributes["scale"]
core_radii = np.broadcast_to(a["core_radius_um"], n) / scale
cbeam_lantern = PhotonicLantern(
	a["port_positions"], core_radii, a["cladding_radius_um"], 3 * a["cladding_radius_um"],
	np.repeat(a["n_core"], n), a["n_clad"], a["n_jacket"], a["z_extent_um"], scale,
	core_res=20, clad_res=64, jack_res=48, core_mesh_size=.09, clad_mesh_size=.8,
)
cbeam = Propagator(pl.wavelengths_um[0], cbeam_lantern, Nmax=n,
                   save_dir="outputs/cbeam_cache")
cbeam.fixed_zstep = 20
cbeam.characterize(save=False)
cbeam_transfer = cbeam.compute_transfer_matrix(channel_basis=True)
cbeam_mesh = cbeam.make_mesh_at_z(0)
mass = construct_B(cbeam_mesh, sparse=True)
cbeam_lp_modes = make_lp_modes(*cbeam_mesh.points[:, :2].T)
cbeam_lp_modes /= np.sqrt([np.vdot(mode, mass @ mode) for mode in cbeam_lp_modes])[:, None]
cbeam_transfer = cbeam_transfer @ cbeam.inner_product(cbeam.vs[0], cbeam_lp_modes, mass)

# CBeam's numerical eigenmodes and isolated ports have arbitrary ordering and
# phase. Put its rows in Lightbeam port order and align those unobservable phases.
cost = np.linalg.norm(
	np.abs(cbeam_transfer)[:, None, :] - np.abs(transfer)[None, :, :], axis=2
)
cbeam_rows, lightbeam_rows = linear_sum_assignment(cost)
ordered = np.empty_like(cbeam_transfer)
ordered[lightbeam_rows] = cbeam_transfer[cbeam_rows]
ordered *= np.exp(1j * np.angle(np.sum(ordered.conj() * transfer.conj(), axis=1)))[:, None]
cbeam_transfer = ordered.conj()
cbeam_error = np.linalg.norm(cbeam_transfer - transfer) / np.linalg.norm(transfer)

# Compare predicted port intensities and the magnitude and phase of both matrices.
lp_coeffs = lp_projector @ focal.electric_field.shaped[pl.input_footprint]
powers = [np.abs(coeffs) ** 2 for coeffs in (direct, via_transfer, cbeam_transfer @ lp_coeffs)]
powers = [power / power.sum() for power in powers]
matrices = (transfer, cbeam_transfer)
fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
for row, (name, power, matrix) in enumerate(zip(("Lightbeam", "CBeam"), powers[1:], matrices)):
	ports = np.arange(1, n + 1)
	axes[row, 0].bar(ports - .2, powers[0], .4, label="Direct projection")
	axes[row, 0].bar(ports + .2, power, .4, label=f"{name} transfer")
	axes[row, 0].set(title=f"{name} predicted output", xlabel="Output port",
	                 ylabel="Normalized intensity", ylim=(0, max(powers[0]) * 1.15))
	axes[row, 0].legend()
	for col, (values, title, limits, cmap) in enumerate((
		(np.abs(matrix), "|Transfer|", (0, max(np.abs(m).max() for m in matrices)), "viridis"),
		(np.angle(matrix), "Transfer phase", (-np.pi, np.pi), "twilight"),
	), 1):
		image = axes[row, col].imshow(values, vmin=limits[0], vmax=limits[1], cmap=cmap)
		axes[row, col].set(
			title=f"{name} {title}",
			xlabel="LP input mode", ylabel="Output port",
			xticks=range(6), yticks=range(6),
			xticklabels=("01", "11c", "11s", "21c", "21s", "02"),
			yticklabels=range(1, 7),
		)
		fig.colorbar(image, ax=axes[row, col], shrink=.8)
output_path = Path("outputs/transfer_matrix_comparison.png")
output_path.parent.mkdir(exist_ok=True)
fig.savefig(output_path, dpi=180)

print(f"focal-field projection error: {projection_error:.3%}")
print(f"saved comparison plot to {output_path}")
