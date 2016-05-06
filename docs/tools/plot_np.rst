############
plot_np
############

The ``gamtools plot_np`` tool is used to plot the sequencing coverage and the
window calling output for a single NP over the whole genome. It can be very
useful as a quick quality control, as NPs that have been sequenced correctly
display a characteristic pattern of coverage, with many large regions (up to
whole chromosomes) displaying very low coverage.

Plotting the coverage for an NP requires three files, a bigwig file made from
the raw sequencing data, a bed file containing the positive windows and a
`genome file`_.

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools plot_np -w BIGWIG_FILE -b BED_FILE
                   -g GENOME_FILE -o OUTPUT_FILE 

