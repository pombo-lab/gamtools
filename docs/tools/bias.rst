############
bias
############

The ``gamtools bias`` tool is used to calculate the bias that could be
attributed to a particular genomic feature. Given a bedgraph file with a value
for each genomic window, the tool generates ten bins, each containing an equal
number of genomic windows, based on the feature. For example, if the bedgraph
file contains restriction site density, bin one would contain the genomic
windows with the 10% lowest site density and bin ten would contain the 10
highest. It then calculates the mean interaction frequency between windows in
each possible combination of bins, normalised for linear genomic distance.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools bias [OPTIONS] -f <FEATURE_FILE> -o <OUTPUT_PATH> -m <MATRIX> [<MATRIX> ...] 


**Required parameters:**

+---------------------+-------------------------------------------------+
| Option              | Description                                     |
+=====================+=================================================+
| -f, --feature-path  | Path to input bedgraph containing one value per |
|                     | genomic window.                                 |
+---------------------+-------------------------------------------------+
| -m, --matrix-paths  | Path to one or more interaction matrices to use |
|                     | for calculating biases.                         |
+---------------------+-------------------------------------------------+
| -o, --output-path   | Output bias matrix file. Path to use for saving |
|                     | the result of the bias calculation.             |
+---------------------+-------------------------------------------------+


