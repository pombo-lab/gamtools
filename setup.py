import os
import sys
from setuptools import setup, Extension

def parse_setuppy_commands():
    """Check the commands and respond appropriately (avoid
    parsing Cython and template files if False).
    """
    args = sys.argv[1:]

    if not args:
        # User forgot to give an argument probably, let setuptools handle that.
        return True

    info_commands = ['--help-commands', '--name', '--version', '-V',
                     '--fullname', '--author', '--author-email',
                     '--maintainer', '--maintainer-email', '--contact',
                     '--contact-email', '--url', '--license', '--description',
                     '--long-description', '--platforms', '--classifiers',
                     '--keywords', '--provides', '--requires', '--obsoletes',
                     'egg_info', 'install_egg_info', 'rotate']

    for command in info_commands:
        if command in args:
            return False
    else:
        return True


if parse_setuppy_commands():

    from Cython.Distutils import build_ext

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

else:

    from setuptools.command.build_ext import build_ext

    CustomBuildExtCommand = build_ext

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "gamtools",
    version = "2.0.0-alpha1",
    author = "Rob Beagrie",
    author_email = "rob@beagrie.com",
    url = "http://gam.tools",
    description = ("A package containing some utilities for analyzing GAM data."),
    license = "Apache2.0",
    package_dir = {'': 'lib'},
    packages=['gamtools', 'gamtools.qc'],
    ext_modules = [Extension('gamtools.cosegregation_internal',
                   ["lib/gamtools/cosegregation_internal.pyx"]),
                   Extension('gamtools.mirnylib_numutils_internal',
                      ["lib/gamtools/mirnylib_numutils_internal.pyx"],),
                   ],
    install_requires=[
      'cython',
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
