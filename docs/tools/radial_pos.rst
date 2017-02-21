###########
radial_pos
###########

The ``gamtools radial_pos`` tool is used to calculate chromatin radial
positioning from GAM segregation tables. Radial positioning is estimated from
the average size of NPs that contain a given chromatin region, since chromatin
that occupies a more peripheral position will be intersected by smaller, more
apical NPs (i.e. those which slice the nucleus close to the top/bottom/sides),
whereas more central chromatin can only be intersected by larger, more
equatorial NPs (i.e. those which slice the nucleus through the middle).

The size of each NP is estimated from it's genomic coverage, i.e. the number
of positive windows. NPs which contain a larger number of positive windows
are assumed to also be larger in volume.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools radial_pos [OPTIONS] -s <SEGREGATION_FILE> -o <OUTPUT_FILE> 

**Optional parameters:**

+----------------------+----------------------------------------------------------------------------------+
| Option               | Description                                                                      |
+======================+==================================================================================+
| -n, --no-blanks      | Exclude regions that were never detected from the output (for making bedgraphs)  |
+----------------------+----------------------------------------------------------------------------------+
