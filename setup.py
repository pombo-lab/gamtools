import os
from distutils.extension import Extension
from setuptools import setup

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
    install_requires=['numpy', 'scipy', 'doit', 'pandas'],
    ext_modules = Extension('gamtools.cosegregation_internal', ["gamtools/cosegregation_internal.c"])
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
