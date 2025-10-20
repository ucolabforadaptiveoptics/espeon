# Implementation of the PL in the style of hcipy's pyramid WFS.

import h5py
import numpy as np
import hcipy as hc
from os import path
from hcipy import imshow_field
from matplotlib import pyplot as plt

from .utils import input_to_2d, PL_DESIGNS_PATH

class PhotonicLanternOptics(hc.wavefront_sensing.WavefrontSensorOptics):
	"""
	The optical elements for a simulated photonic lantern. This mostly consists of accurately defining grids that are compatible with the simulation that's been run.
	"""
	def __init__(self, tag):
		path_to_pl = path.join(PL_DESIGNS_PATH, f"{tag}.hdf5")
		with h5py.File(path_to_pl) as f:
			self.output = np.array(f["pl_output"])
			self.wavelengths_um = np.array(f["wavelengths_um"])
			self.design_name = f["pl_output"].attrs["design_name"]
			mesh_spacing_um = f["pl_output"].attrs["sim_mesh_spacing_um"]
			mesh_extent_um = f["pl_output"].attrs["sim_mesh_extent_um"]
			self.input_footprint = (np.array(f["input_footprint_x"]), np.array(f["input_footprint_y"]))
			self.extent = [[f[f"input_footprint_{l}"].attrs[f"{l}min"], f[f"input_footprint_{l}"].attrs[f"{l}max"]] for l in ["x", "y"]]
		self.nports = self.output.shape[1]
		self.projectors = [np.linalg.pinv(self.output[i,:,:].T) for i in range(self.output.shape[0])]
		delta = mesh_spacing_um * 1e-6 * np.ones(2)
		dims = (mesh_extent_um // mesh_spacing_um + 1) * np.ones(2)
		zero = delta * (-dims / 2 + np.mod(dims, 2) * 0.5)
		self.focal_grid = hc.CartesianGrid(hc.RegularCoords(delta, dims, zero))
		# generate launch fields here
		
	def coeffs(self, focal_wavefront):
		"""
  		Takes in a PSF at the lantern entrance and returns the coefficients of a projection onto the lantern basis.
		If you want the lantern reading from here, do np.abs(_) ** 2 on the output of this function.
		"""
		wf_wl_um = focal_wavefront.wavelength * 1e6
		pl_run_index = np.argmin(np.abs(wf_wl_um - self.wavelengths_um))
		pl_wl_um = self.wavelengths_um[pl_run_index]
		assert np.abs(wf_wl_um - pl_wl_um) < 1e-3, f"PL simulation was run at a different wavelength than the input wavefront: closest PL simulation was at {pl_wl_um} microns vs. input at {wf_wl_um} microns."
		profile_to_project = focal_wavefront.electric_field.shaped[self.input_footprint]
		return self.projectors[pl_run_index] @ profile_to_project

	def lantern_output(self, focal_wavefront):
		"""
		Pairs the response of readout with the field of what the output looks like for plotting.
		"""
		coeff_vals = self.coeffs(focal_field)
		lantern_reading = sum(c * lf for (c, lf) in zip(coeff_vals, self.launch_fields))
		return coeff_vals, np.abs(lantern_reading) ** 2
		
	def show_lantern_output(self, focal_wavefront):
		coeffs, lantern_reading = self.plotting_lantern_output(focal_wavefront)
		fig, axs = plt.subplots(1, 2)
		fig.subplots_adjust(top=1.4, bottom=0.0)
		for ax in axs:
			ax.set_xticks([])
			ax.set_yticks([])
		imshow_field(np.log10(focal_field.intensity), ax=axs[0])
		axs[0].set_title("Lantern input")
		axs[1].imshow(np.abs(lantern_reading))
		axs[1].set_title("Lantern output")
		plt.show()

	def show_lantern_modes(self, wl_index=0, nrows=4, crop=1):
		rm, cm = nrows, int(np.ceil(self.nports / nrows))
		fig, axs = plt.subplots(rm, cm)
		plt.suptitle(f"Photonic lantern entrance modes, {self.design_name}, wavelength = {self.wavelengths_um[wl_index]} microns")
		plt.subplots_adjust(wspace=0.05, hspace=0.05)
		for (i, o) in enumerate(self.output[wl_index,:,:]):
			r, c = i // cm, i % cm
			axs[r][c].imshow(np.abs(input_to_2d(o, self.input_footprint, self.extent))[crop:-crop,crop:-crop])
			axs[r][c].set_xticks([])
			axs[r][c].set_yticks([])
		for i in range(self.nports, rm * cm):
			r, c = i // cm, i % cm
			fig.delaxes(axs[r][c])
		# plt.show()