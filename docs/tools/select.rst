############
select
############

The ``gamtools select`` tool is used to select or exclude samples
from a `segregation table`_.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools select [OPTIONS] -s <SEGREGATION_FILE> -o <OUTPUT_FILE>
                            -n [<SAMPLE_NAME> [<SAMPLE_NAME> ...]]  

**Required parameters:**

+------------------------+--------------------------------------------------+
| Option                 | Description                                      |
+========================+==================================================+
| -s, --segregation-file | A file containing the segregation of all samples |
+------------------------+--------------------------------------------------+
| -n, --sample-names     | Names of the samples to remove                   |
+------------------------+--------------------------------------------------+
| -o, --output-file      | Output file path (or - to write to stdout)       |
+------------------------+--------------------------------------------------+


**Optional parameters:**

+--------------------+------------------------------------------------------+
| Option             | Description                                          |
+====================+======================================================+
| -d, --drop-samples | Discard the listed samples (default: discard samples |
|                    | not in the list)                                     |
+--------------------+------------------------------------------------------+


