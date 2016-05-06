############
convert
############

The ``gamtools convert`` tool is used to convert proximity matrices output by GAMtools to/from
various different formats.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools convert [OPTIONS] <INPUT_MATRIX> <OUTPUT_MATRIX> 

**Optional parameters:**

+-----------------------+----------------------------------------------------------------------------------------------------------------------------+
| Option                | Description                                                                                                                |
+=======================+============================================================================================================================+
| -i, --input-format    | Input matrix file format                                                                                                   |
+-----------------------+----------------------------------------------------------------------------------------------------------------------------+
| -o, --output-format   | Output matrix file format                                                                                                  |
+-----------------------+----------------------------------------------------------------------------------------------------------------------------+
| -t, --thresholds-file | Thresholds file. If specified, any values lower than the specified thresholds will be masked/excluded from the output file |
+-----------------------+----------------------------------------------------------------------------------------------------------------------------+
| -w, --windows-file    | File containing the genomic locations of matrix bins (only required if not specified in input matrix file).                |
+-----------------------+----------------------------------------------------------------------------------------------------------------------------+
| -r, --region          | Region covered by the input matrix (required if -w /--windows-file is specified)                                           |
+-----------------------+----------------------------------------------------------------------------------------------------------------------------+
