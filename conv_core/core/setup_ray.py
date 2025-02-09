# setup.py
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy

extension = Extension(
    name="ray_line",
    sources=["conv_core/core/ray_line_wrap.pyx"],
    include_dirs=[numpy.get_include(), 'conv_core/core/src/'],  # This line is important
    extra_compile_args=["-fopenmp"],
    extra_link_args=["-fopenmp"],
)

setup(
    name="Ray Line",
    ext_modules=cythonize(extension),
    zip_safe=False,
)
