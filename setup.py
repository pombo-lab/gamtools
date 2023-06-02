from setuptools import Extension, setup
from Cython.Build import cythonize
import numpy

extensions = [
    Extension('gamtools.cosegregation_internal',
        ["lib/gamtools/cosegregation_internal.pyx"]),
    Extension('gamtools.mirnylib_numutils_internal',
        ["lib/gamtools/mirnylib_numutils_internal.pyx"])
	   ]
 
setup(
    ext_modules = cythonize(extensions),
)
