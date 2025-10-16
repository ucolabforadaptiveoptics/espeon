# Makes photonic lantern data files

import h5py
import numpy as np
from matplotlib import pyplot as plt
from os import path

import lightbeam as lb
from .utils import PROJECT_ROOT, input_to_2d

def triangle_pattern(core_offset):
	theta = np.array([0, 2 * np.pi / 3, 4 * np.pi / 3])
	return core_offset * np.vstack((np.cos(theta), np.sin(theta))).T

def ngon_pattern(n, nrings, core_offset, theta_init=0):
	"""
	Generate an n-gon pattern of PL port positions. The pattern starts with a single port at the center and expands outward in concentric 
 rings.

	Args:
		n (int): the number of sides of the n-gon.
		nrings (int): The number of concentric hexagonal rings to generate. The total number of ports increases with the number of rings.
		core_offset (float): The distance between adjacent rings in the hexagonal pattern.

	Returns:
		numpy.ndarray: A 2D array of shape (nports, 2), where `nports` is the total number of ports. Each row represents the (x, y) coordinates of a port in the hexan-gongonal pattern.
	"""
	theta_step = 2 * np.pi / n
	nports = int(1 + (n / 2) * nrings * (nrings - 1))
	port_positions = np.zeros((nports, 2))
	nports_so_far = 0
	for i in range(nrings):
		nports_per_ring = max(1, n*i)
		theta = theta_init
		current_position = i * core_offset * np.array([np.cos(theta), np.sin(theta)])
		next_position = i * core_offset * np.array([np.cos(theta + theta_step), np.sin(theta + theta_step)])
		for j in range(nports_per_ring):
			if i > 0 and j % i == 0:
				theta += theta_step
				current_position = next_position
				next_position = i * core_offset * np.array([np.cos(theta + theta_step), np.sin(theta + theta_step)])
			cvx_coeff = 0 if i == 0 else (j % i) / i
			port_positions[nports_so_far,:] = (1 - cvx_coeff) * current_position + cvx_coeff * next_position
			nports_so_far += 1
	return port_positions

def setup_lantern(
	port_positions, 
	core_radius_um, cladding_radius_um, z_extent_um, scale,
	n_clad, n_core, n_jacket, 
	wavelength_um,
	simulation_params = {
		"name" : "lightbeam",
		"mesh_extent_um" : 512,
		"mesh_spacing_um" : 1,
		"dz_um" : 50,
		"PML" : 8
	}
):
	PML = simulation_params["PML"]
	mesh = lb.RectMesh3D(
		xw = simulation_params["mesh_extent_um"],
		yw = simulation_params["mesh_extent_um"],
		zw = z_extent_um,
		ds = simulation_params["mesh_spacing_um"],
		dz = simulation_params["dz_um"],
		PML = PML
	)
	# set up beam propagation IOR profile
	
	lant = lb.optics.Lantern(
		port_positions * scale, core_radius_um, cladding_radius_um * scale, 0, z_extent_um, (n_core, n_clad, n_jacket), final_scale = 1/scale
	)
	lant.set_sampling(mesh.xy)
	xg, yg = mesh.grids_without_pml()
	launch_fields = [
		lb.normalize(lb.lpfield(xg-pos[0], yg-pos[1], 0, 1, corerad, wavelength_um, n_core, n_clad))
		for (pos, corerad) in zip(lant.init_core_locs, core_radius_um)
	]
	lbprop = lb.Prop3D(wavelength_um, mesh, lant, n_clad)
 	# set up grid definition so we can unravel a 2D PL output into a 1D column of our transformation matrix
	out = np.zeros_like(xg)
	lant.set_IORsq(out, z_extent_um)
	out = out[PML:-PML,PML:-PML]
	input_footprint = np.where(out >= n_clad ** 2)
	complement_mask = np.ones_like(out, dtype=bool)
	for (i, j) in zip(input_footprint[0], input_footprint[1]):
		complement_mask[i,j] = False
	extent = [
		[np.min(x), np.max(x)]
		for x in input_footprint
	]
	return lbprop, input_footprint, extent, launch_fields
	
def save_lantern_design(
	design_name, 
	port_positions, 
	core_radius_um, cladding_radius_um, z_extent_um, scale,
	n_clad, n_core, n_jacket, 
	wavelength_um,
	design_core_radius_um=None,
	simulation_params = {
		"name" : "lightbeam",
		"mesh_extent_um" : 512,
		"mesh_spacing_um" : 1,
		"dz_um" : 50,
		"PML" : 8
	},
	force_overwrite=False,
	do_plot=True
):
	"""
	Runs a PL simulation and saves the results with design metadata to an HDF5 file.
	
	Parameters
	----------
	design_name : str
		The user-provided identifier for this PL design.
		
	port_positions : np.array(nports, 2)
		The locations of the PL ports.
		
	core_radius_um, cladding_radius_um, z_extent_um : float
		Dimensions on the multi-mode end, all in microns: radii for the core and cladding, and the overall taper length.
		core_radius_um can also be a list of length nports.
		
	scale : float
		The (dimensionless) ratio of the single-mode end diameter to the multi-mode end diameter.
		
	n_clad, n_core, n_jacket : float
		Indices of refraction for the cladding, core, and jacket.
		
	wavelength_um : float
		The wavelength, in microns, at which the simulation should be run.
		
	simulator : str
		Which simulation to run. Either "lightbeam" or "cbeam"; currently only implementing "lightbeam".
		
	force_overwrite : bool
		Whether or not to overwrite an existing PL design with the same design_name.
		
	do_plot : bool
		Whether or not to plot simulation output live.
	"""
	if design_core_radius_um is None:
		design_core_radius_um = core_radius_um
  
	if isinstance(core_radius_um, float):
		core_radius_um = [core_radius_um for _ in range(port_positions.shape[0])]

	filepath = path.join(PROJECT_ROOT, "pl_designs", design_name + ".hdf5")
	if not force_overwrite and path.exists(filepath):
		raise OSError(f"File already exists at {filepath}, change design_name or force overwrite by passing in force_overwrite=True")
	# setup for lightbeam
	PML = simulation_params["PML"]
	lbprop, input_footprint, extent, launch_fields = setup_lantern(
		port_positions, core_radius_um, cladding_radius_um, z_extent_um, scale, n_clad, n_core, n_jacket, wavelength_um, simulation_params
	)
	pl_output = []
	# the backwards beam-propagation runs
	for (i, lf) in enumerate(launch_fields):
		print(f"Illuminating core {i}")
		if do_plot:
			plt.imshow(np.abs(lf) ** 2)
			plt.show()
		u = lbprop.prop2end(lf)[PML:-PML,PML:-PML] # this step takes ~minutes
		u = u[input_footprint]
		u /= np.linalg.norm(u)
		pl_output.append(u)
		if do_plot:
			output_intensity = np.abs(input_to_2d(u, input_footprint, extent)) ** 2
			plt.imshow(output_intensity)
			plt.show()

	# write output to hdf5
	with h5py.File(filepath, "w") as f:
		pl_output_dset = f.create_dataset("pl_output", data=pl_output)
		for k in ["design_name", "port_positions", "core_radius_um", "cladding_radius_um", "z_extent_um", "scale", "n_core", "n_clad", "n_jacket", "wavelength_um"]:
			pl_output_dset.attrs[k] = eval(k)
		for k in simulation_params:
			pl_output_dset.attrs["sim_" + k] = simulation_params[k]
		input_footprint_x_dset = f.create_dataset("input_footprint_x", data=input_footprint[0])
		input_footprint_x_dset.attrs["xmin"] = extent[0][0]
		input_footprint_x_dset.attrs["xmax"] = extent[0][1]
		input_footprint_y_dset = f.create_dataset("input_footprint_y", data=input_footprint[1])
		input_footprint_y_dset.attrs["ymin"] = extent[1][0]
		input_footprint_y_dset.attrs["ymax"] = extent[1][1]
		