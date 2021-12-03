import os
import sys
from setuptools import setup, Extension
from setuptools.command.sdist import sdist
from setuptools.command.build_ext import build_ext

try:
    from wheel.bdist_wheel import bdist_wheel
except ImportError:
    bdist_wheel = None


class CustomBuildExtCommand(build_ext):
    """Customized setuptools build_ext command - checks numpy is installed."""
    def run(self):

        # Check numpy is installed before trying to find the location
        # of numpy headers

        try:
            import numpy
        except ImportError:
            raise ImportError('numpy need to be installed before GAMtools can be '
                              'compiled. Try installing with "pip install numpy" '
                              'before installing GAMtools.')

        self.include_dirs.append(numpy.get_include())
        self.include_dirs.append('lib/gamtools/data/include')

        build_ext.run(self)


class custom_cythonize_sdist(sdist):
    def run(self):
        # Make sure the compiled Cython files in the distribution are up-to-date
        from Cython.Build import cythonize
        cythonize([
                   "lib/gamtools/cosegregation_internal.pyx",
                   "lib/gamtools/mirnylib_numutils_internal.pyx",
        ], language_level = "3")
        cythonize(
            Extension(
                'gamtools.slice_wrapper',
                 sources=["lib/gamtools/slice_wrapper.pyx", "lib/gamtools/slice_internals.cpp"],
                 libraries=["gsl", "gslcblas"],
                 language='c++'), language='c++', language_level=3),
        sdist.run(self)


custom_classes = {
      'sdist': custom_cythonize_sdist,
      'build_ext': CustomBuildExtCommand,
    }


if not bdist_wheel is None:
    class CustomBdistWheelCommand(bdist_wheel):
        """Customized bdist_wheel command - checks numpy is installed and cythonises"""
        def run(self):

            from Cython.Build import cythonize
            cythonize([
                       "lib/gamtools/cosegregation_internal.pyx",
                       "lib/gamtools/mirnylib_numutils_internal.pyx",
            ], language_level = "3")
            cythonize(
                Extension(
                    'gamtools.slice_wrapper',
                     sources=["lib/gamtools/slice_wrapper.pyx", "lib/gamtools/slice_internals.cpp"],
                     libraries=["gsl", "gslcblas"],
                     language='c++'), language='c++', language_level=3),

            bdist_wheel.run(self)

    custom_classes['bdist_wheel'] = CustomBdistWheelCommand

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "gamtools",
    version = "2.0.0alpha7",
    author = "Rob Beagrie",
    author_email = "rob@beagrie.co.uk",
    url = "https://gam.tools",
    description = ("A package containing some utilities for analyzing GAM data."),
    license = "Apache2.0",
    package_dir = {'': 'lib'},
    packages=['gamtools', 'gamtools.qc'],
    ext_modules = [Extension('gamtools.cosegregation_internal',
                   ["lib/gamtools/cosegregation_internal.c"]),
                   Extension('gamtools.mirnylib_numutils_internal',
                      ["lib/gamtools/mirnylib_numutils_internal.c"],),
                   Extension('gamtools.slice_wrapper',
                             sources=["lib/gamtools/slice_wrapper.cpp", "lib/gamtools/slice_internals.cpp"],
                             libraries=["gsl", "gslcblas"],
                             language='c++'),
                   ],
    cmdclass = custom_classes,
    install_requires=[
      'numpy',
      'scipy',
      'pandas',
      'wrapit>=0.3.0',
      'pytest'],
    extras_require={
      ':python_version<"3.0"': ['doit==0.29.0'],
      ':python_version>="3.0"': ['doit==0.30.0'],
      ':python_version<"3.0"': ['mock'],
    },
    entry_points = {
                    'console_scripts': [
                        'gamtools = gamtools.main:main',
                        'create_empty_bedgraph = gamtools.utils:empty_bedgraph_from_cmdline',
                    ]
                   },
    long_description=read('README.md'),
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "License :: OSI Approved :: Apache Software License",
    ],
)
