from datetime import datetime
import numpy as np
import os
from os import path
from copy import copy
from hcipy import Field, imshow_field
from matplotlib import pyplot as plt

PROJECT_ROOT = path.dirname(path.dirname(path.abspath(__file__)))
PL_DESIGNS_PATH = path.join(PROJECT_ROOT, "pl_designs")
if not os.path.isdir(PL_DESIGNS_PATH):
    os.mkdir(PL_DESIGNS_PATH)

zernike_names = [
    "tip", "tilt", "focus", "astig", "astig45", "coma90", "coma", "tricoma90", "tricoma", "spherical", "astig5th45", "astig5th"
] + [f"Z{i}" for i in range(12, 82)]

def is_list_or_dim1_array(x):
    return isinstance(x, list) or (isinstance(x, np.ndarray) and len(x.shape) == 1)

def rms(x):
    return np.sqrt(np.mean((x - np.mean(x)) ** 2))

def angles_relative_to_center(x, y):
    xc, yc = np.mean(x), np.mean(y)
    xd, yd = x - xc, y - yc
    return (np.arctan2(yd, xd) + 3 * np.pi / 2) % (2 * np.pi)

def nanify(phase_screen, aperture=None):
    if aperture is None:
        aperture = phase_screen
    x = copy(phase_screen)
    x = x - np.mean(x)
    x[np.where(aperture == 0)] = np.nan
    return Field(x, phase_screen.grid)

def imshow_psf(f: Field, **kwargs):
    imshow_field(np.log10(f / np.max(f)), **kwargs)

def peak_to_valley(x):
    return np.max(x) - np.min(x)

def norm(a, b):
    return np.sum(a * np.conj(b))

def corr(a, b):
    return np.abs(norm(a, b)) / np.sqrt(norm(a, a) * norm(b, b))

def input_to_2d(input_efield, input_footprint, extent):
	"""
	Takes in an input electric field as a 1D array (coordinates represented by input_footprint) and fills it in to a 2D grid for plotting.
	"""
	xl, yl = extent[0][1] - extent[0][0] + 1, extent[1][1] - extent[1][0] + 1
	xm, ym = extent[0][0], extent[1][0]
	input_efield_2d = np.zeros((xl, yl), dtype=np.complex64)
	input_efield_2d[input_footprint[0] - xm, input_footprint[1] - ym] = input_efield
	return input_efield_2d

def normalize(x):
    return x / np.sum(x)

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