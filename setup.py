import os
from setuptools import setup, Extension

try:
    from Cython.Distutils import build_ext
except:
    from setuptools.command.build_ext import build_ext
    ext_modules = [Extension('gamtools.cosegregation_internal',
                   ["lib/gamtools/cosegregation_internal.c"])]
else:
    ext_modules = [Extension('gamtools.cosegregation_internal',
                   ["lib/gamtools/cosegregation_internal.pyx"])]


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

        build_ext.run(self)

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "gamtools",
    version = "1.1.0",
    author = "Rob Beagrie",
    author_email = "rob@beagrie.com",
    url = "http://gam.tools",
    description = ("A package containing some utilities for analyzing GAM data."),
    license = "Apache2.0",
    package_dir = {'': 'lib'},
    packages=['gamtools', 'gamtools.qc'],
    ext_modules = ext_modules,
    install_requires=[
      'cython',
      'numpy',
      'scipy',
      'pandas',
      'wrapit',
      'pytest'],
    extras_require={
      ':python_version<"3.0"': ['doit==0.29.0'],
      ':python_version>="3.0"': ['doit==0.30.0'],
      ':python_version<"3.0"': ['mock'],
    },
    # Set include_dirs in a custom build_ext class so that numpy is only
    # required if we are compiling C files
    cmdclass={
          'build_ext': CustomBuildExtCommand,
      },
    entry_points = {
                    # TODO: make new EIYBrowse filetypes using IO functions in gamtools.matrix
                    #'EIYBrowse.filetypes': [
                    #    'gam_segmentation_file = gamtools.segmentation:GamSegmentationFile',
                    #],
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
