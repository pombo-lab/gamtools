##########
compaction
##########

The ``gamtools compaction`` tool is used to calculate chromatin compaction from
GAM segregation tables. Chromatin compaction is estimated from the number of
NPs that contain a given chromatin region, since chromatin that occupies a larger
volume will be intersected by a greater number of NPs.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools compaction [OPTIONS] -s <SEGREGATION_FILE> -o <OUTPUT_FILE> 

**Optional parameters:**

+----------------------+----------------------------------------------------------------------------------+
| Option               | Description                                                                      |
+======================+==================================================================================+
| -n, --no-blanks      | Exclude regions that were never detected from the output (for making bedgraphs)  |
+----------------------+----------------------------------------------------------------------------------+
