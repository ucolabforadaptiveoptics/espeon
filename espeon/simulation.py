# Implementation of the PL in the style of hcipy's pyramid WFS.

import h5py
import hcipy
import numpy as np
import hcipy as hc
from os import path
from copy import copy
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
			att = f["pl_output"].attrs
			self.attributes = {k: att[k] for k in att}
			self.design_name = self.attributes["design_name"]
			self.wavelengths_um = np.array(f["wavelengths_um"])
			self.input_footprint = (np.array(f["input_footprint_x"]), np.array(f["input_footprint_y"]))
			self.extent = [[f[f"input_footprint_{l}"].attrs[f"{l}min"], f[f"input_footprint_{l}"].attrs[f"{l}max"]] for l in ["x", "y"]]
		self.nports = self.output.shape[1]
		self.projectors = [np.linalg.pinv(self.output[i,:,:].T) for i in range(self.output.shape[0])]
		delta = self.attributes["sim_mesh_spacing_um"] * 1e-6 * np.ones(2)
		dims = (self.attributes["sim_mesh_extent_um"] // self.attributes["sim_mesh_spacing_um"] + 1) * np.ones(2)
		zero = delta * (-dims / 2 + np.mod(dims, 2) * 0.5)
		self.focal_grid = hc.CartesianGrid(hc.RegularCoords(delta, dims, zero))
		# generate launch fields here

	def generate_launch_fields(self, scaleup=3):
		# this is if you want the full-frame images
		# requires lightbeam integration for the LP mode calculator
		import lightbeam as lb
		dx = self.attributes["sim_mesh_spacing_um"]
		extent = self.attributes["sim_mesh_extent_um"]
		n = int(extent // dx) + 1
		coords = np.linspace(-extent / 2.0, extent / 2.0, n)
		self.xg, self.yg = np.meshgrid(coords, coords, indexing='xy')

		self.launch_fields = [
			lb.normalize(lb.lpfield(self.xg-pos[0], self.yg-pos[1], 0, 1, corerad * scaleup, np.median(self.wavelengths_um), self.attributes["n_core"], self.attributes["n_clad"]))
			for (pos, corerad) in zip(self.attributes["port_positions"] * self.attributes["scale"], self.attributes["core_radius_um"])
		]
		
	def coeffs(self, focal_wavefront: hcipy.Wavefront) -> np.ndarray:
		"""
  		Takes in a PSF at the lantern entrance and returns the coefficients of a projection onto the lantern basis.

		Input: focal_wavefront (hcipy.Wavefront)
		Output: coeffs (np.ndarray, length = self.nports)
		"""
		wf_wl_um = focal_wavefront.wavelength * 1e6
		pl_run_index = np.argmin(np.abs(wf_wl_um - self.wavelengths_um))
		pl_wl_um = self.wavelengths_um[pl_run_index]
		assert np.abs(wf_wl_um - pl_wl_um) < 1e-3, f"PL simulation was run at a different wavelength than the input wavefront: closest PL simulation was at {pl_wl_um} microns vs. input at {wf_wl_um} microns."
		profile_to_project = focal_wavefront.electric_field.shaped[self.input_footprint]
		return np.array(self.projectors[pl_run_index] @ profile_to_project)

	def readout(self, focal_wavefront: hcipy.Wavefront) -> np.ndarray:
		"""
		Returns the intensity per port of a PL in response to a PSF at the lantern entrance.

		Input - focal_wavefront: hcipy.Wavefront
		Output - readout: np.ndarray, shape = (self.nports,)
		"""
		return np.abs(self.coeffs(focal_wavefront) ** 2)

	def image(self, focal_wavefront):
		"""
		Converts a PL readout into an image of the SMF end.

		Input - focal_wavefront: hcipy.Wavefront
		Output - image: np.ndarray, shape = self.xg.shape
		"""
		coeff_vals = self.coeffs(focal_wavefront)
		lantern_reading = sum(c * lf for (c, lf) in zip(coeff_vals, self.launch_fields))
		return np.abs(lantern_reading) ** 2
		
	def show_image(self, focal_wavefront):
		img = self.image(focal_wavefront)
		fig, axs = plt.subplots(1, 2)
		fig.subplots_adjust(top=1.4, bottom=0.0)
		for ax in axs:
			ax.set_xticks([])
			ax.set_yticks([])
		hcipy.imshow_field(np.log10(focal_wavefront.intensity), ax=axs[0])
		axs[0].set_title("Lantern input")
		axs[1].imshow(np.abs(img))
		axs[1].set_title("Lantern output")
		plt.show()

	def show_modes(self, wl_index=0, nrows=4, crop=1):
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
		plt.show()
