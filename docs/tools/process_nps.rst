############
process_nps
############

The ``gamtools process_nps`` tool is used to map raw sequencing data
from a collection of NPs and call positive windows from those NPs
to generate a `segregation table`_. It can optionally also calculate
various QC metrics for each NP, generate bigwig/bed files for
visualising the raw data and calculate proximity matrices.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools process_nps [OPTIONS] -g <GENOME_FILE> <FASTQ_FILE> [<FASTQ_FILE> ...] 


**Optional parameters:**

+-----------------------+--------------------------------------------------------+
| Option                | Description                                            |
+=======================+========================================================+
| -o, --output_dir      | Write segregation, matrix etc. to this directory       |
+-----------------------+--------------------------------------------------------+
| -q, --minimum-mapq    | Filter out any mapped read with a mapping quality less |
|                       | than x (default is 20, use -q 0 for no filtering)      |
+-----------------------+--------------------------------------------------------+
| -c, --do-qc           | Perform sample quality control.                        |
+-----------------------+--------------------------------------------------------+
| -i, --bigwigs         | Make bigWig files.                                     |
+-----------------------+--------------------------------------------------------+
| -b, --bigbeds         | Make bed files of positive windows                     |
+-----------------------+--------------------------------------------------------+
| -w, --window-sizes    | One or more window sizes for calling positive windows  |
+-----------------------+--------------------------------------------------------+
| -s, --matrix-sizes    | Resolutions for which proximity matrices should be     |
|                       | produced.                                              |
+-----------------------+--------------------------------------------------------+
| --qc-window-size      | Use this window size for qc (default is median window  |
|                       | size).                                                 |
+-----------------------+--------------------------------------------------------+
| -f, --fittings_dir    | Write segregation curve fitting plots to this          |
|                       | directory                                              |
+-----------------------+--------------------------------------------------------+
| -d, --details-file    | If specified, write a table of fitting parameters to   |
|                       | this path                                              |
+-----------------------+--------------------------------------------------------+
| --additional-qc-files | Any additional qc files to filter on                   |
+-----------------------+--------------------------------------------------------+

**Parameters inherited from doit:**

``gamtools process_nps`` uses doit_ as a task dependency engine, to
determine what actions need to be performed and in which order. A number
of additional command line parameters are available that control doit's behaviour.

+----------------------+----------------------------------------------------------------+
| Option               | Description                                                    |
+======================+================================================================+
| --doit-db-file       | Doit saves information about each run in a                     |
|                      | database file. This parameter specifies the                    |
|                      | location of that database file.                                |
+----------------------+----------------------------------------------------------------+
| --doit-backend       | Doit database format. (one of                                  |
|                      | sqlite3, json, dbm. default: dbm)                              |
+----------------------+----------------------------------------------------------------+
| --doit-verbosity     | ``0`` capture (do not print) stdout/stderr from task.          |
|                      | ``1`` capture stdout only.                                     |
|                      | ``2`` do not capture anything (print everything                |
|                      | immediately). Default: 1                                       |
+----------------------+----------------------------------------------------------------+
| --doit-reporter      | Where should doit report the output from each task. One        |
|                      | of (json, console, zero, executed-only). Default: console      |
+----------------------+----------------------------------------------------------------+
| --doit-process       | Number of subprocesses (default is 0, i.e.  serial processing) |
+----------------------+----------------------------------------------------------------+
| --doit-parallel-type | Tasks can be executed in parallel in different ways:           |
|                      | ``process``: uses python multiprocessing module                |
|                      | ``thread``: uses threads. Default is process.                  |
+----------------------+----------------------------------------------------------------+

.. _doit: http://pydoit.org
