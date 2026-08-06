# Implementation of the PL in the style of hcipy's pyramid WFS.

import h5py
import hcipy
import hcipy as hc
import numpy as np
from matplotlib import pyplot as plt
from os import path

from .utils import input_to_2d, PL_DESIGNS_PATH, normalize


def _make_lightbeam_input_grid(attributes):
	delta = attributes["sim_mesh_spacing_um"] * 1e-6 * np.ones(2)
	dim = int(attributes["sim_mesh_extent_um"] // attributes["sim_mesh_spacing_um"] + 1)
	dims = dim * np.ones(2)
	zero = delta * (-dims / 2 + np.mod(dims, 2) * 0.5)
	return hc.CartesianGrid(hc.RegularCoords(delta, dims, zero))


def _make_cbeam_input_grid(attributes):
	# The focal grid is also used to measure the incident wavefront power. A
	# grid limited to the multimode cladding diameter crops most of the first
	# Airy ring and makes coupling efficiencies appear artificially high. Use
	# the lantern's largest transverse simulation boundary for the focal grid;
	# the LP projection itself is still restricted to the cladding below.
	# This padding captures the PSF power that misses the multimode entrance,
	# so ``readout.sum() / focal_wavefront.total_power`` is a coupling
	# efficiency rather than a conditional efficiency inside the core crop.
	radius = (
		attributes["jacket_radius_um"] * attributes["scale"] * 1e-6
	)
	return hc.make_uniform_grid(256, np.full(2, 2 * radius))


class _LightbeamBackend:
	def __init__(
		self, file, output, attributes, wavelengths_um, input_grid=None
	):
		self.output = output
		self.nports = output.shape[1]
		self.attributes = attributes
		self.input_footprint = (
			np.array(file["input_footprint_x"]),
			np.array(file["input_footprint_y"])
		)
		self.extent = [
			[
				file[f"input_footprint_{axis}"].attrs[f"{axis}min"],
				file[f"input_footprint_{axis}"].attrs[f"{axis}max"]
			]
			for axis in ("x", "y")
		]
		self.projectors = [
			np.linalg.pinv(self.output[index].T)
			for index in range(len(wavelengths_um))
		]
		self.input_grid = _make_lightbeam_input_grid(attributes)
		if input_grid is not None and input_grid != self.input_grid:
			raise ValueError(
				"A Lightbeam design must use the input grid saved with it."
			)
		if "total_power" in file:
			self.total_power = np.array(file["total_power"])

	def coeffs(self, focal_wavefront, wavelength_index):
		# Lightbeam arrays are indexed [x, y], while HCIPy shaped fields use
		# image order [y, x]. Sample the same physical points in both layouts.
		x_indices, y_indices = self.input_footprint
		profile = focal_wavefront.electric_field.shaped[
			y_indices, x_indices
		]
		return np.asarray(self.projectors[wavelength_index] @ profile)

	def readout(self, focal_wavefront, wavelength_index):
		coefficients = self.coeffs(focal_wavefront, wavelength_index)
		return np.abs(coefficients) ** 2 * self.input_grid.weights

	def mode_images(self, wavelength_index):
		return [
			input_to_2d(mode, self.input_footprint, self.extent)
			for mode in self.output[wavelength_index]
		]

	def make_launch_fields(self, scaleup, wavelength_um):
		import lightbeam as lb

		dx = self.attributes["sim_mesh_spacing_um"]
		extent = self.attributes["sim_mesh_extent_um"]
		n = int(extent // dx) + 1
		coords = np.linspace(-extent / 2.0, extent / 2.0, n)
		x, y = np.meshgrid(coords, coords, indexing="xy")
		core_radii = np.broadcast_to(
			self.attributes["core_radius_um"], self.nports
		)
		fields = [
			np.real(normalize(lb.lpfield(
				x - position[0], y - position[1], 0, 1,
				core_radius * scaleup, wavelength_um,
				self.attributes["n_core"], self.attributes["n_clad"]
			)))
			for position, core_radius in zip(
				self.attributes["port_positions"] * self.attributes["scale"],
				core_radii
			)
		]
		return x, y, fields


class _CBeamBackend:
	def __init__(
		self, file, output, attributes, wavelengths_um, input_grid=None
	):
		self.output = output
		self.attributes = attributes
		self.wavelengths_um = wavelengths_um
		self.nports = output[0].shape[0]
		if input_grid is None:
			input_grid = _make_cbeam_input_grid(attributes)
		self.input_grid = input_grid
		radius = attributes["cladding_radius_um"] * 1e-6
		self.input_mask = np.asarray(
			input_grid.x**2 + input_grid.y**2 <= radius**2
		)
		self.projection_grid = input_grid.subset(self.input_mask)
		self.bases = [
			self._make_basis(self.projection_grid, wavelength_index)
			for wavelength_index in range(len(wavelengths_um))
		]

	def _make_basis(self, grid, wavelength_index):
		radius = self.attributes["cladding_radius_um"] * 1e-6
		wavelength = self.wavelengths_um[wavelength_index] * 1e-6
		V_number = (
			2 * np.pi / wavelength * radius
			* np.sqrt(
				self.attributes["n_clad"] ** 2
				- self.attributes["n_jacket"] ** 2
			)
		)
		return hc.make_lp_modes(grid, V_number, radius)

	def _matrix(self, wavelength_index):
		return self.output[wavelength_index]

	def coeffs(self, focal_wavefront, wavelength_index):
		field = focal_wavefront.electric_field
		basis = self.bases[wavelength_index]
		if field.grid != self.input_grid:
			radius = self.attributes["cladding_radius_um"] * 1e-6
			mask = np.asarray(field.grid.x**2 + field.grid.y**2 <= radius**2)
			basis = self._make_basis(field.grid.subset(mask), wavelength_index)
		else:
			mask = self.input_mask
		lp_coeffs = basis.coefficients_for(np.asarray(field)[mask])
		matrix = self._matrix(wavelength_index)
		return np.asarray(matrix @ lp_coeffs)

	def readout(self, focal_wavefront, wavelength_index):
		return np.abs(self.coeffs(focal_wavefront, wavelength_index)) ** 2

	def mode_images(self, wavelength_index):
		basis = self.bases[wavelength_index]
		matrix = self._matrix(wavelength_index)
		projector = np.linalg.pinv(basis.transformation_matrix)
		port_modes = np.conj(matrix @ projector)
		images = []
		for mode in port_modes:
			field = np.zeros(self.input_grid.size, dtype=complex)
			field[self.input_mask] = mode
			images.append(hc.Field(field, self.input_grid).shaped)
		return images

	def make_launch_fields(self, scaleup, wavelength_um):
		diameter = (
			2 * self.attributes["jacket_radius_um"]
			* self.attributes["scale"] * 1e-6
		)
		grid = hc.make_pupil_grid(256, diameter)
		core_radii = np.broadcast_to(
			self.attributes["core_radius_um"], self.nports
		)
		fields = []
		for position, core_radius in zip(
			self.attributes["port_positions"] * self.attributes["scale"] * 1e-6,
			core_radii * scaleup * 1e-6
		):
			V_number = (
				2 * np.pi / (wavelength_um * 1e-6) * core_radius
				* np.sqrt(
					self.attributes["n_core"] ** 2
					- self.attributes["n_clad"] ** 2
				)
			)
			basis = hc.make_lp_modes(
				grid.shifted(-position), V_number, core_radius
			)
			fields.append(np.real(normalize(basis[0].shaped)))
		return grid.x.shaped * 1e6, grid.y.shaped * 1e6, fields


class PhotonicLanternOptics(hc.wavefront_sensing.WavefrontSensorOptics):
	"""Optical model for a saved Lightbeam or CBeam photonic lantern.

	For CBeam designs, ``input_grid`` sets the entrance-plane sampling; LP
	projection uses the Lightbeam convention of retaining the cladding disk.
	Lightbeam designs use the grid stored with the simulation.
	"""
	def __init__(
		self, tag, despath=None, wavelength_um=None, input_grid=None
	):
		if despath is None:
			despath = PL_DESIGNS_PATH
		path_to_pl = path.join(despath, f"{tag}.hdf5")
		with h5py.File(path_to_pl) as file:
			output_node = file["pl_output"]
			attributes = dict(output_node.attrs)
			try:
				wavelengths_um = np.array(file["wavelengths_um"])
			except KeyError:
				wavelengths_um = np.array([wavelength_um])
			if isinstance(output_node, h5py.Group):
				output = [
					np.array(output_node[str(index)])
					for index in range(len(wavelengths_um))
				]
			else:
				output = np.array(output_node)
			if "wavelengths_um" not in file:
				output = output[None]
			backend_name = attributes.get("sim_name", "lightbeam")
			if isinstance(backend_name, bytes):
				backend_name = backend_name.decode()
			backends = {
				"lightbeam": _LightbeamBackend,
				"cbeam": _CBeamBackend,
			}
			try:
				backend_class = backends[backend_name]
			except KeyError:
				raise ValueError(f"Unsupported photonic-lantern backend: {backend_name}")
			self._backend = backend_class(
				file, output, attributes, wavelengths_um, input_grid
			)

		self.output = self._backend.output
		self.attributes = attributes
		self.design_name = attributes["design_name"]
		self.wavelengths_um = wavelengths_um
		self.nports = self._backend.nports
		self.input_grid = self._backend.input_grid
		self.focal_grid = self.input_grid
		for name in ("input_footprint", "extent", "projectors", "total_power"):
			if hasattr(self._backend, name):
				setattr(self, name, getattr(self._backend, name))

	def generate_launch_fields(self, scaleup=1):
		self.xg, self.yg, self.launch_fields = self._backend.make_launch_fields(
			scaleup, np.median(self.wavelengths_um)
		)

	def _wavelength_index(self, focal_wavefront):
		wavelength_um = focal_wavefront.wavelength * 1e6
		index = np.argmin(np.abs(wavelength_um - self.wavelengths_um))
		closest_um = self.wavelengths_um[index]
		assert np.abs(wavelength_um - closest_um) < 1e-3, (
			"PL simulation was run at a different wavelength than the input "
			f"wavefront: closest PL simulation was at {closest_um} microns vs. "
			f"input at {wavelength_um} microns."
		)
		return index

	def coeffs(self, focal_wavefront: hcipy.Wavefront) -> np.ndarray:
		"""Return the complex output-port coefficients."""
		index = self._wavelength_index(focal_wavefront)
		return self._backend.coeffs(focal_wavefront, index)

	def readout(self, focal_wavefront: hcipy.Wavefront) -> np.ndarray:
		"""Return the intensity in each output port."""
		index = self._wavelength_index(focal_wavefront)
		return self._backend.readout(focal_wavefront, index)

	def image(self, focal_wavefront):
		readout_vals = self.readout(focal_wavefront)
		return sum([
			coefficient * launch_field
			for coefficient, launch_field in zip(readout_vals, self.launch_fields)
		])

	def show_image(self, focal_wavefront):
		image = self.image(focal_wavefront)
		fig, axs = plt.subplots(1, 2)
		fig.subplots_adjust(top=1.4, bottom=0.0)
		for ax in axs:
			ax.set_xticks([])
			ax.set_yticks([])
		hcipy.imshow_field(np.log10(focal_wavefront.intensity), ax=axs[0])
		axs[0].set_title("Lantern input")
		axs[1].imshow(np.abs(image))
		axs[1].set_title("Lantern output")
		plt.show()

	def show_modes(self, wl_index=0, nrows=4, crop=1, fn=np.abs):
		images = self._backend.mode_images(wl_index)
		nrows = min(nrows, len(images))
		ncols = int(np.ceil(len(images) / nrows))
		fig, axs = plt.subplots(nrows, ncols, squeeze=False)
		plt.subplots_adjust(wspace=0.05, hspace=0.05)
		for index, image in enumerate(images):
			row, column = divmod(index, ncols)
			image = fn(image)
			if crop:
				image = image[crop:-crop, crop:-crop]
			axs[row, column].imshow(image)
			axs[row, column].set_xticks([])
			axs[row, column].set_yticks([])
		for index in range(len(images), nrows * ncols):
			row, column = divmod(index, ncols)
			fig.delaxes(axs[row, column])

	def V(self, wl_um):
		return (
			2 * np.pi / wl_um * self.attributes["cladding_radius_um"]
			* np.sqrt(
				self.attributes["n_clad"] ** 2
				- self.attributes["n_jacket"] ** 2
			)
		)
