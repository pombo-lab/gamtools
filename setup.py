import os
from distutils.extension import Extension
from setuptools import setup

failed_imports = []

try:
    from Cython.Build import cythonize
except ImportError:
    failed_imports.append('cython')

try:
    import numpy
except ImportError:
    failed_imports.append('numpy')

if len(failed_imports) > 0:
    raise ImportError('{} need to be installed before GAMtools can be '
                      'compiled. Try installing them with "pip '
                      'install {}" before installing GAMtools.'.format(
                          ' and '.join(failed_imports),
                          ' '.join(failed_imports)))

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "gamtools",
    version = "0.0.1",
    author = "Rob Beagrie",
    author_email = "rob@beagrie.com",
    description = ("A package containing some utilities for analyzing GAM data."),
    license = "BSD",
    packages=['gamtools'],
    install_requires=['numpy', 'scipy', 'doit', 'pandas', 'wrapit', 'pytest'],
    include_dirs=[numpy.get_include()],
    ext_modules = [Extension('gamtools.cosegregation_internal',
                            ["gamtools/cosegregation_internal.c"])],
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
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
