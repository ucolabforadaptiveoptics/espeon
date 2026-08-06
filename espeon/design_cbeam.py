# Makes photonic lantern data files with `cbeam` as the simulation backend

from os import path

import h5py
import hcipy as hc
import numpy as np
from cbeam.propagator import Propagator, construct_B
from cbeam.waveguide import PhotonicLantern
from matplotlib import pyplot as plt

from .utils import PROJECT_ROOT


def _make_lp_basis(mesh, mass, wavelength_um, radius_um, n_core, n_clad):
    weights = np.asarray(mass @ np.ones(len(mesh.points))).ravel()
    grid = hc.CartesianGrid(
        hc.UnstructuredCoords(mesh.points[:, :2].T),
        weights=weights
    )
    V_number = (
        2 * np.pi / wavelength_um * radius_um
        * np.sqrt(n_core**2 - n_clad**2)
    )
    basis = hc.make_lp_modes(grid, V_number, radius_um)
    modes = basis.transformation_matrix.T
    norms = np.sqrt([
        np.real(np.vdot(mode, mass @ mode)) for mode in modes
    ])
    return modes / norms[:, None]


def _set_output_channel_order(prop, z_extent_um):
    mesh = prop.make_mesh_at_z(z_extent_um)
    isolated = prop.compute_isolated_basis(z_extent_um)
    channel_indices = {
        primitive.label: index
        for index, primitive in enumerate(prop.wvg.prim3Dgroups[-1])
    }
    order = [
        channel_indices[label]
        for label in mesh.cell_sets
        if label in channel_indices
    ]
    order.extend(index for index in range(len(isolated)) if index not in order)
    prop.compute_change_of_basis(isolated[order], z_extent_um)


def _compute_transfer_matrix(prop):
    columns = []
    for mode in range(prop.Nmax):
        input_coefficients = np.zeros(prop.Nmax)
        input_coefficients[mode] = 1
        _, _, output_coefficients = prop.propagate(input_coefficients)
        columns.append(prop.to_channel_basis(output_coefficients))
    return np.column_stack(columns)


def _write_output(file, transfer_matrices):
    shapes = {matrix.shape for matrix in transfer_matrices}
    if len(shapes) == 1:
        return file.create_dataset("pl_output", data=np.asarray(transfer_matrices))

    output = file.create_group("pl_output")
    for index, matrix in enumerate(transfer_matrices):
        output.create_dataset(str(index), data=matrix)
    return output


def save_lantern_design_cbeam(
    design_name,
    port_positions,
    core_radius_um, cladding_radius_um, jacket_radius_um, z_extent_um, scale,
    n_clad, n_core, n_jacket,
    wavelengths_um,
    force_overwrite=False,
    do_plot=True,
    **cbeam_kwargs
):
    """Run a CBeam simulation and save its LP-to-port transfer matrices."""
    port_positions = np.asarray(port_positions)
    nports = len(port_positions)
    core_radius_um = np.broadcast_to(np.asarray(core_radius_um, dtype=float), nports)
    wavelengths_um = np.atleast_1d(wavelengths_um)
    filepath = path.join(PROJECT_ROOT, "pl_designs", design_name + ".hdf5")
    if not force_overwrite and path.exists(filepath):
        raise OSError(
            f"File already exists at {filepath}, change design_name or force "
            "overwrite by passing in force_overwrite=True"
        )

    lantern = PhotonicLantern(
        port_positions, core_radius_um / scale, cladding_radius_um,
        jacket_radius_um, np.repeat(n_core, nports), n_clad, n_jacket,
        z_extent_um, scale, **cbeam_kwargs
    )
    transfer_matrices = []
    for wavelength_um in wavelengths_um:
        prop = Propagator(wavelength_um, lantern)
        mesh = prop.make_mesh_at_z(0)
        mass = construct_B(mesh, sparse=True)
        lp_modes = _make_lp_basis(
            mesh, mass, wavelength_um, cladding_radius_um, n_clad, n_jacket
        )
        prop.Nmax = len(lp_modes)
        prop.characterize(mesh=mesh)
        _set_output_channel_order(prop, z_extent_um)
        eigen_to_ports = _compute_transfer_matrix(prop)
        eigen_from_lp = prop.inner_product(prop.vs[0], lp_modes, mass)
        transfer_matrices.append((eigen_to_ports @ eigen_from_lp).conj())

    # Persist the completed simulation before opening any interactive plots.
    # ``plt.show()`` may block indefinitely, which previously meant the HDF5
    # file was never created until every plot window had been closed.
    with h5py.File(filepath, "w") as file:
        output = _write_output(file, transfer_matrices)
        values = {
            "design_name": design_name,
            "port_positions": port_positions,
            "core_radius_um": core_radius_um,
            "cladding_radius_um": cladding_radius_um,
            "jacket_radius_um": jacket_radius_um,
            "z_extent_um": z_extent_um,
            "scale": scale,
            "n_core": n_core,
            "n_clad": n_clad,
            "n_jacket": n_jacket,
            "sim_name": "cbeam",
        }
        for key, value in values.items():
            output.attrs[key] = value
        for key, value in cbeam_kwargs.items():
            if value is not None:
                output.attrs["sim_" + key] = value
        file.create_dataset("wavelengths_um", data=wavelengths_um)

    if do_plot:
        for wavelength_um, transfer in zip(wavelengths_um, transfer_matrices):
            plt.figure()
            plt.imshow(np.abs(transfer))
            plt.title(f"{design_name}, {wavelength_um} µm")
            plt.colorbar()
            plt.show()

    if len({matrix.shape for matrix in transfer_matrices}) == 1:
        return np.asarray(transfer_matrices)
    return transfer_matrices
