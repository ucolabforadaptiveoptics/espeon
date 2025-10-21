# Common linearity tools so I don't keep rewriting them.
# Technically, these aren't PL simulation tools per se because they're WFS-agnostic
# but I end up using them a lot

import numpy as np
from matplotlib import pyplot as plt
from .utils import zernike_names

def interpolate_weights(arr: np.ndarray, x: float):
	"""
	Given a sorted array arr and some arr[0] < x < arr[-1], finds the left index idx and the left weight w such that x = w * arr[idx] + (1 - w) * arr[idx+1].
	
	Arguments:
	arr - np.ndarray, (n,)
		The sorted array to interpolate.
	x - float
		The value whose interpolation we want.
		
	Returns:
	idx - int
		The left-hand index.
	w - float
		The weight to apply to the left-hand index.
	"""
	idx = np.where(arr > x)[0][0] - 1
	w = 1 - (x - arr[idx]) / (arr[idx+1] - arr[idx])
	return idx, w

def make_postprocessed_interaction_matrix(amplitudes, mode_sweep, poke_amplitude=0.1):
	"""
	Makes a modal interaction matrix given a linearity sweep dataset.
	
	Arguments:
	amplitudes - np.ndarray, (namp,)
		The amplitudes poked for each Zernike mode.
	mode_sweep - np.ndarray, (nzern, namp, nwfs,)
		The WFS response for each Zernike number and amplitude.
	poke_amplitude - float
		The poke amplitude used for the push/pull interaction. If this isn't in "amplitudes", the given readings will be interpolated.
		
	Returns:
	interaction_matrix - np.ndarray, (nwfs, nzern,)
		The interaction matrix.
	"""
	nzern, namp, nwfs = mode_sweep.shape
	assert len(amplitudes.shape) == 1 and len(amplitudes) == namp, "Malformed input; make sure the amplitude arrays match."
	interaction_matrix = np.zeros((nwfs, nzern))
	idx_pos, w_pos = interpolate_weights(amplitudes, poke_amplitude)
	idx_neg, w_neg = interpolate_weights(amplitudes, -poke_amplitude)
	for i in range(nzern):
		s_push = w_pos * mode_sweep[i,idx_pos,:] + (1 - w_pos) * mode_sweep[i,idx_pos+1,:]
		s_pull = w_neg * mode_sweep[i,idx_neg,:] + (1 - w_neg) * mode_sweep[i,idx_neg+1,:]
		s = (s_push - s_pull) / (2 * poke_amplitude)
		interaction_matrix[:,i] = s.ravel()
	
	return interaction_matrix

def plot_linearity(amplitudes, responses, title_mod="", savepath="", showplot=True, zernike_modes=False, zernike_offset=0):
	nzern = responses.shape[0]
	nrow = min(nzern, 3)
	ncol = int(np.ceil(nzern / 3))
	fig, axs = plt.subplots(ncol, nrow, sharex=True, sharey=True, figsize=(3 * nrow, 3 * ncol))
	if nrow == 1 and ncol < 3:
		axs = [axs]
	title = "Linearity curves"
	if title_mod != "":
		title += f", {title_mod}"
	plt.suptitle(title, y=1.12)
	for i in range(nzern):
		r, c = i // 3, i % 3
		if ncol > 1:
			ax = axs[r][c]
		else:
			ax = axs[c]
		ax.set_ylim([min(amplitudes), max(amplitudes)])
		if not zernike_modes:
			ax.title.set_text(f"Mode {i + 1}")
		else:
			ax.title.set_text(zernike_names[i+zernike_offset])
		ax.plot(amplitudes, amplitudes, '--k')
		for j in range(nzern):
			alpha = 1 if i == j else 0.25
			ax.plot(amplitudes, responses[i,:,j], alpha=alpha)
	for k in range(i + 1, nrow * ncol):
		r, c = k // 3, k % 3
		fig.delaxes(axs[r][c])
	if len(savepath) > 0:
		plt.savefig(savepath, bbox_inches='tight')
	if showplot:
		plt.show()
