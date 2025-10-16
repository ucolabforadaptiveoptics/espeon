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
    x = np.maximum(x, 0.0)
    if np.all(x == 0.0):
        return x
    return x / np.sum(x)
