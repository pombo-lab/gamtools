"""
==================
The qc sub-package
==================

The qc sub-package contains modules for doing quality control
of GAM datasets:

- gamtools.qc.fastqc: module for parsing fastqc output files
- gamtools.qc.merge: module for merging qc metrics from different sources
- gamtools.qc.pass_qc: module for identifying NPs which pass/fail QC
- gamtools.qc.screen: module for parsing fastqc_screen output files
- gamtools.qc.segregation: module for extracting qc stats from segregation tables.

"""

from . import screen, fastqc, segregation, merge, pass_qc
