############
call_windows
############

The ``gamtools call_windows`` tool determines which genomic regions were present in each NP.
The input file is a tab delimited table in which the first three columns indicate the
genomic region in bed format (chrom, start, stop) and the remaining columns give the number
of reads mapping to each genomic region for each NP. For example:

::

  chrom   start    stop     NP_1   NP_2   NP_3   NP_4   NP_5
  chr19   0        50000    0      0      0      0      0
  chr19   50000    100000   1      26     1      54     0
  chr19   100000   150000   0      34     0      0      1
  chr19   150000   200000   2      16     0      0      0
  chr19   200000   250000   0      1      32     7      0
  chr19   250000   300000   3      0      0      0      1
  chr19   300000   350000   1      4      50     12     0
  chr19   350000   400000   2      3      32     1      3
  chr19   400000   450000   0      0      115    0      0

This type of coverage table can be generated quite easily using 
`bedtools multicov`_ and is generated automatically by the GAMtools
:doc:`process_nps` command.

The output file is in the same format, but each entry is either a 1
(indicating the region was present) or a 0 (indicating that it was
absent).

::

  chrom   start    stop     NP_1   NP_2   NP_3   NP_4   NP_5
  chr19   0        50000    0      0      0      0      0
  chr19   50000    100000   0      1      0      1      0
  chr19   100000   150000   0      1      0      0      0
  chr19   150000   200000   0      1      0      0      0
  chr19   200000   250000   0      0      1      0      0
  chr19   250000   300000   0      0      0      0      0
  chr19   300000   350000   0      0      1      1      0
  chr19   350000   400000   0      0      1      0      0
  chr19   400000   450000   0      0      1      0      0

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools call_windows [OPTIONS] <COVERAGE_TABLE>

**Optional parameters:**

+----------------------+----------------------------------------------------------------------------------+
| Option               | Description                                                                      |
+======================+==================================================================================+
| -d, --details-file   | Write a table of fitting parameters to this path                                 |
+----------------------+----------------------------------------------------------------------------------+
| -o, --output-file    | Output segregation file to create (or "-" to write to stdout), default is stdout |
+----------------------+----------------------------------------------------------------------------------+
| -f, --fitting-folder | Save plots for each individual curve fitting to this folder                      |
+----------------------+----------------------------------------------------------------------------------+


.. _bedtools multicov: https://bedtools.readthedocs.io/en/latest/content/tools/multicov.html
