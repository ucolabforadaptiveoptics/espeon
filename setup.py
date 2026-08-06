from setuptools import setup, find_packages

setup(
    name="espeon",
    author='Aditya R. Sengupta',
    version="3.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "jax",
        "matplotlib",
        "hcipy",
        "h5py",
    ],
    extras_require={
        "lightbeam": [
            "numpy<2.4",
            "numba",
            "lightbeam @ git+https://github.com/aditya-sengupta/lightbeam@2eae558",
        ],
        "cbeam": [
            "cbeam @ git+https://github.com/jw-lin/cbeam.git",
        ],
    },
)
