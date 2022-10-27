############
resolution_qc
############

The ``gamtools resolution_qc`` tool is used to perform some basic QC checks
on a segregation table. This allows the user to generate segregation tables
at a variety of resolutions, then check each resolution to determine the 
minimum window size that can be used for this specific dataset.

This tool requires a segregation table, plus the thickness of each NP and
the average nuclear radius (in any units as long as both parameters are
specified in the same way).

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools resolution_qc [-x CHROM [CHROM ...]] [-g GENOME_SIZE]
                         [-p NPS_PER_SAMPLE] [-i] SEGREGATION_TABLE
                         SLICE_THICKNESS NUCLEAR_RADIUS
                   

