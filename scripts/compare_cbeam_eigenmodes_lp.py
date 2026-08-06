"""Compare CBeam entrance eigenmodes with the analytic LP-mode basis."""

from pathlib import Path

import lightbeam as lb
import matplotlib.pyplot as plt
import numpy as np
from cbeam import FEval
from cbeam.propagator import Propagator, construct_B
from cbeam.waveguide import PhotonicLantern
from scipy.optimize import linear_sum_assignment

from espeon import PhotonicLanternOptics


pl = PhotonicLanternOptics("kim-6port-260729")
labels = ("LP01", "LP11c", "LP11s", "LP21c", "LP21s", "LP02")
specs = ((0, 1, 0), (1, 1, 0), (1, 1, 1), (2, 1, 0), (2, 1, 1), (0, 2, 0))


def lp_modes(x, y):
	phi, fields = np.arctan2(y, x), []
	for l, m, sine in specs:
		field = lb.lpfield(x, y, l, m, pl.attributes["cladding_radius_um"],
		                   pl.wavelengths_um[0], pl.attributes["n_clad"], pl.attributes["n_jacket"])
		field *= ([np.cos(l * phi), np.sin(l * phi)][sine] * np.exp(-1j * l * phi)) if l else 1
		fields.append(np.real_if_close(field))
	return np.asarray(fields)


# Solve only the multimode entrance; no longitudinal characterization is needed.
n, a, scale = pl.nports, pl.attributes, pl.attributes["scale"]
lantern = PhotonicLantern(
	a["port_positions"], np.broadcast_to(a["core_radius_um"], n) / scale,
	a["cladding_radius_um"], 3 * a["cladding_radius_um"],
	np.repeat(a["n_core"], n), a["n_clad"], a["n_jacket"],
	a["z_extent_um"], scale, core_res=20, clad_res=64, jack_res=48,
	core_mesh_size=.09, clad_mesh_size=.8,
)
solver = Propagator(pl.wavelengths_um[0], lantern, Nmax=n,
                    save_dir="outputs/cbeam_eigenmode_cache")
neff, eigenmodes = solver.solve_at(0)
mesh, mass = solver.mesh, construct_B(solver.mesh, sparse=True)
node_lp = lp_modes(*mesh.points[:, :2].T)
node_lp /= np.sqrt([np.vdot(field, mass @ field) for field in node_lp])[:, None]
eigenmodes /= np.sqrt([np.vdot(field, mass @ field) for field in eigenmodes])[:, None]

# Test the user's reordering hypothesis.
overlap = eigenmodes.conj() @ mass @ node_lp.T
eigen_index, lp_index = linear_sum_assignment(-np.abs(overlap) ** 2)
reordered = np.empty_like(eigenmodes, dtype=complex)
reorder_fidelity = np.empty(n)
for eig, lp in zip(eigen_index, lp_index):
	phase = np.exp(1j * np.angle(overlap[eig, lp]))
	reordered[lp] = eigenmodes[eig] * phase
	reorder_fidelity[lp] = np.abs(overlap[eig, lp]) ** 2

# Degenerate eigenspaces are only defined up to unitary rotations. Align the
# LP01, LP11, and LP21/LP02 groups separately without mixing distinct neff.
aligned = np.empty_like(eigenmodes, dtype=complex)
subspace_fidelity = np.empty(n)
for group in ((0,), (1, 2), (3, 4, 5)):
	block = overlap[np.ix_(group, group)]
	u, singular_values, vh = np.linalg.svd(block)
	aligned[list(group)] = (u @ vh).conj().T @ eigenmodes[list(group)]
	subspace_fidelity[list(group)] = singular_values ** 2

# Resample all fields onto one regular grid for a direct visual comparison.
axis = np.linspace(-12, 12, 121)
xg, yg = np.meshgrid(axis, axis)
FEval.sort_mesh(mesh)
evaluate = lambda fields: np.asarray([
	np.asarray(FEval.evaluate_grid(axis, axis, np.real(field), mesh.tree)).T
	for field in fields
])
images = (lp_modes(xg, yg).real, evaluate(reordered), evaluate(aligned))
images = tuple(fields / np.sqrt(np.sum(fields ** 2, axis=(1, 2)))[:, None, None] for fields in images)
residual = images[2] - images[0]

fig, axes = plt.subplots(n, 4, figsize=(10, 15), constrained_layout=True)
headers = ("Analytic LP mode", "Best reordering", "Degenerate-group rotation", "Residual")
for column, title in enumerate(headers):
	axes[0, column].set_title(title)
for mode in range(n):
	limit = max(np.abs(fields[mode]).max() for fields in images)
	for column, field in enumerate((*images, residual)):
		axes[mode, column].imshow(field[mode], origin="lower", extent=(-12, 12, -12, 12),
		                          cmap="RdBu_r", vmin=-limit, vmax=limit)
		axes[mode, column].set(xticks=[], yticks=[])
	axes[mode, 0].set_ylabel(labels[mode])
	axes[mode, 1].text(.5, .04, f"overlap²={reorder_fidelity[mode]:.3f}",
	                   transform=axes[mode, 1].transAxes, ha="center")
	axes[mode, 2].text(.5, .04, f"subspace={subspace_fidelity[mode]:.3f}",
	                   transform=axes[mode, 2].transAxes, ha="center")
	axes[mode, 3].text(.5, .04, f"NRMSE={np.linalg.norm(residual[mode]):.3f}",
	                   transform=axes[mode, 3].transAxes, ha="center")

output = Path("outputs/cbeam_eigenmodes_vs_lp.png")
output.parent.mkdir(exist_ok=True)
fig.savefig(output, dpi=180)

fig, (ax_overlap, ax_fidelity) = plt.subplots(1, 2, figsize=(9, 4), constrained_layout=True)
image = ax_overlap.imshow(np.abs(overlap) ** 2, cmap="viridis", vmin=0, vmax=1)
ax_overlap.set(title="Raw pairwise overlap²", xlabel="Analytic LP mode",
               ylabel="CBeam eigenmode", xticks=range(n), yticks=range(n),
               xticklabels=labels, yticklabels=range(1, n + 1))
fig.colorbar(image, ax=ax_overlap)
x = np.arange(n)
ax_fidelity.bar(x - .2, reorder_fidelity, .4, label="Reordering only")
ax_fidelity.bar(x + .2, subspace_fidelity, .4, label="Degenerate rotation")
ax_fidelity.set(title="Mode/subspace fidelity", ylabel="Fidelity", ylim=(0, 1.05),
                xticks=x, xticklabels=labels)
ax_fidelity.legend()
metrics_output = Path("outputs/cbeam_eigenmode_lp_metrics.png")
fig.savefig(metrics_output, dpi=180)

print("CBeam entrance neff:", neff)
print("Best-reordering overlap²:", reorder_fidelity)
print("Degenerate-subspace fidelity:", subspace_fidelity)
print("Saved", output, "and", metrics_output)
