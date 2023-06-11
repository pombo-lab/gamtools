from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

extensions = [
    Extension('gamtools.cosegregation_internal',
        ["lib/gamtools/cosegregation_internal.pyx"],
              include_dirs=[numpy.get_include()]),
    Extension('gamtools.mirnylib_numutils_internal',
        ["lib/gamtools/mirnylib_numutils_internal.pyx"],
              include_dirs=[numpy.get_include()]),
	   ]
 
setup(
    ext_modules = cythonize(extensions, language_level = "3"),
)
