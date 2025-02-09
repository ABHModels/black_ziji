# setup.py
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extension = Extension(
    name="rebin_spectrum_wrapper",
    sources=["conv_core/core/rebin_spectrum_wrapper.pyx"],
    include_dirs=[numpy.get_include(), 'conv_core/core/src/']  # This line is important
)

setup(
    name='Rebin Spectrum',
    ext_modules=cythonize(extension),
    zip_safe=False,
)
