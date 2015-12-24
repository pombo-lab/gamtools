import os
from distutils.core import setup
from Cython.Build import cythonize

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "GamTools",
    version = "0.0.1",
    author = "Rob Beagrie",
    author_email = "rob@beagrie.com",
    description = ("A package containing some utilities for analyzing GAM data."),
    license = "BSD",
    packages=['GamTools', 'wrapit', 'doit'],
    ext_modules = cythonize("GamTools/cosegregation_internal.pyx"),
    entry_points = {'EIYBrowse.filetypes': [
                        'gam_segmentation_file = GamTools.segmentation:GamSegmentationFile',
                    ]
                   },
    long_description=read('README'),
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
    ],
)
