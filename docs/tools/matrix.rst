############
matrix
############

The ``gamtools matrix`` tool is used to calculate proximity matrices from
segregation tables.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools matrix [OPTIONS] -s <SEGREGATION_FILE> -r <REGION> [<REGION> ...] 


**Optional parameters:**

+---------------------+-------------------------------------------------+
| Option              | Description                                     |
+=====================+=================================================+
| -f, --output-format | Output matrix file format (choose from: csv.gz, |
|                     | txt.gz, npz, txt, csv, png, default is txt.gz)  |
+---------------------+-------------------------------------------------+
| -t, --matrix-type   | Method used to calculate the interaction matrix |
|                     | (choose from: cosegregation, linkage, dprime,   |
|                     | default is dprime)                              |
+---------------------+-------------------------------------------------+
| -o, --output-file   | Output matrix file. If not specified, new file  |
|                     | will have the same name as the segregation file |
|                     | and an extension indicating the genomic         |
|                     | region(s) and the matrix method                 |
+---------------------+-------------------------------------------------+

**Specifying regions**

The -r/--regions parameter allows the user to specify the
specific genomic regions to calculate matrices for. If
one region is specified, a matrix is calculated for
that region against itself. If more than one region is
specified, a matrix is calculated for each region
against the other. Regions are specified using UCSC
browser syntax, i.e. "chr4" for the whole of
chromosome 4 or "chr4:100000-200000" for a sub-region
of the chromosome.

