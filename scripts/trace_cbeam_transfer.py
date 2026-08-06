"""Explain why CBeam's raw and sandbox-reconstructed matrices differ."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


sandbox = Path.home() / "projects/plwfs_sim_sandbox/outputs/kim_6port_comparison"
raw = np.load(sandbox / "cbeam.npz")["transfer_to_ports"]
metrics = np.load(sandbox / "basis_metrics.npz")
coefficients = metrics["lp_cbeam_coefficients"]
labels = metrics["lp_labels"]

# The saved `raw` matrix was already calculated with channel_basis=True. That
# flag changes only its output coordinates:
#
#   channel_basis=False: eigenmodes(z=0) -> eigenmodes(z=L)
#   channel_basis=True:  eigenmodes(z=0) -> isolated output ports
#
# It does not change the input columns from CBeam eigenmodes to analytic LP
# modes, so the input basis conversion below is still required.
#
# The sandbox constructs normalized back-propagated port fields and projects
# them onto LP modes, so (up to its sampling/normalization approximations)
#
#   coefficients = effective_eigen_to_lp @ inv(raw).
#
# Therefore inv(coefficients) is raw expressed in the LP input basis. Its
# conjugate additionally changes from CBeam's phasor convention to Lightbeam's.
# "Effective" matters: the sandbox also resamples, masks, least-squares fits,
# and independently normalizes each port field, all of which are folded in here.
effective_eigen_to_lp = coefficients @ raw
lp_transfer = raw @ np.linalg.inv(effective_eigen_to_lp)
lightbeam_convention = lp_transfer.conj()
sandbox_reconstruction = np.linalg.inv(coefficients).conj()

relative = lambda a, b: np.linalg.norm(a - b) / np.linalg.norm(b)
errors = (
	relative(sandbox_reconstruction, raw),
	relative(np.linalg.inv(coefficients), raw),
	relative(lp_transfer, np.linalg.inv(coefficients)),
	relative(lightbeam_convention, sandbox_reconstruction),
)
residuals = metrics["residual_fractions"][1, 1]

fig, axes = plt.subplots(2, 3, figsize=(12, 7), constrained_layout=True)
plots = (
	(raw, "Raw CBeam transfer", "CBeam eigenmode", False),
	(sandbox_reconstruction, r"$\mathrm{inv}(C_\mathrm{LP})^*$", "LP mode", False),
	(effective_eigen_to_lp, r"Effective input map $A=C_\mathrm{LP}T$", "CBeam eigenmode", False),
	(raw, "Raw CBeam phase", "CBeam eigenmode", True),
	(sandbox_reconstruction, "Reconstructed phase", "LP mode", True),
)
vmax = max(np.abs(raw).max(), np.abs(sandbox_reconstruction).max())
for ax, (matrix, title, xlabel, phase) in zip(axes.flat, plots):
	is_basis_map = "input map" in title
	magnitude_limit = np.abs(matrix).max() if is_basis_map else vmax
	image = ax.imshow(np.angle(matrix) if phase else np.abs(matrix),
	                  cmap="twilight" if phase else "viridis",
	                  vmin=-np.pi if phase else 0,
	                  vmax=np.pi if phase else magnitude_limit)
	ax.set(title=title, xlabel=xlabel, ylabel="LP mode" if is_basis_map else "Output port",
	       xticks=range(6), yticks=range(6))
	if not is_basis_map:
		ax.set_yticklabels(range(1, 7))
	if xlabel == "LP mode":
		ax.set_xticklabels(labels)
	if is_basis_map:
		ax.set_yticklabels(labels)
	if xlabel == "CBeam eigenmode":
		ax.set_xticklabels(range(1, 7))
	fig.colorbar(image, ax=ax, shrink=.8)

ax = axes[1, 2]
names = ("naive\n+ conjugate", "naive", "basis\nchanged", "basis +\nconjugate")
ax.bar(names, np.asarray(errors) * 100)
ax.set(title="Relative matrix error", ylabel="Error [%]",
       ylim=(0, max(errors) * 115))
for index, value in enumerate(errors):
	ax.text(index, value * 100, f"{value:.1%}", ha="center", va="bottom")

output = Path("outputs/cbeam_transfer_trace.png")
output.parent.mkdir(exist_ok=True)
fig.savefig(output, dpi=180)

print("The matrices differ because their columns use different input bases.")
print("channel_basis=True is already present; it changes only the output basis.")
print("Off-diagonal fraction of eigenmode-to-LP map:",
      np.linalg.norm(effective_eigen_to_lp - np.diag(np.diag(effective_eigen_to_lp)))
      / np.linalg.norm(effective_eigen_to_lp))
print("Eigenmode-to-LP singular values:",
      np.linalg.svd(effective_eigen_to_lp, compute_uv=False))
print("LP fit residuals of sandbox port fields:", residuals)
print("Errors [naive conjugated, naive, basis-corrected, fully corrected]:", errors)
print("Saved", output)

assert errors[2] < 1e-12 and errors[3] < 1e-12
