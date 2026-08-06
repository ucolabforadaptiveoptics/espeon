"""Compare Kim six-port outputs for the same HCIPy-generated input PSFs."""

import hcipy as hc
import numpy as np

from espeon import PhotonicLanternOptics


pupil_grid = hc.make_pupil_grid(256, 1)
pupil = hc.make_circular_aperture(1)(pupil_grid)
zernikes = hc.make_zernike_basis(4, 1, pupil_grid)
labels = ("flat", "tip", "tilt", "focus")
phases = (0, *(0.3 * mode for mode in zernikes[1:4]))

results = {}
for backend, tag in (
    ("Lightbeam", "kim-6port-260729"),
    ("CBeam", "kim-6port-cbeam"),
):
    lantern = PhotonicLanternOptics(tag)
    propagate = hc.FraunhoferPropagator(
        pupil_grid, lantern.input_grid, focal_length=7.5
    )
    powers = []
    for phase in phases:
        pupil_wavefront = hc.Wavefront(
            pupil * np.exp(1j * phase), wavelength=1.55e-6
        )
        power = lantern.readout(propagate(pupil_wavefront))
        powers.append(power / power.sum())
    results[backend] = np.asarray(powers)

for index, label in enumerate(labels):
    print(f"{label:5}: Lightbeam {results['Lightbeam'][index].round(3)}")
    print(f"       CBeam     {results['CBeam'][index].round(3)}")

maximum_difference = np.max(np.abs(results["Lightbeam"] - results["CBeam"]))
print(f"maximum normalized port-power difference: {maximum_difference:.3f}")
assert maximum_difference < 0.05
