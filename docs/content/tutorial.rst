############
Tutorial
############

===========
First steps
===========

Installing **GAMtools**
-----------------------

The first step in the **GAMtools** tutorial is to make sure that **GAMtools**
is properly installed. Try to run ``gamtools --help`` and make sure that you
get the following ouput:

.. code-block:: bash

  $ gamtools --help
  usage: gamtools [-h]
                {call_windows,convert,enrichment,matrix,permute_segregation,plot_np,process_nps,select}
                ...

If this command gives you an error message, it is likely that **GAMtools** has
not been installed correctly. Please ensure you have followed the steps
outlined in the installation guide.

Downloading the tutorial data
-----------------------------

Once **GAMtools** is working correctly, you need to download some example data
to work with during the tutorial. The tutorial data is located on the
`GAMtools website <tutorialdata>`_. Download the tutorial data (e.g. by
using wget), extract it and cd into the newly created directory. The directory
should contain a folder called ``fastqs`` and a file called `clean.sh`.

.. code-block:: bash

  $ wget http:/gam.tools/tutorial_data.tar.gz
  $ tar zxvf tutorial_data.tar.gz
  $ cd gamtools_tutorial
  $ ls
  clean.sh  fastqs/

The ``fastqs`` folder contains sequencing data from 100 separate nuclear
profiles (NPs):

.. code-block:: bash

  $ ls fastqs/
  NP_001.fq.gz  NP_026.fq.gz  NP_051.fq.gz  NP_076.fq.gz
  NP_002.fq.gz  NP_027.fq.gz  NP_052.fq.gz  NP_077.fq.gz
  NP_003.fq.gz  NP_028.fq.gz  NP_053.fq.gz  NP_078.fq.gz
  NP_004.fq.gz  NP_029.fq.gz  NP_054.fq.gz  NP_079.fq.gz
  NP_005.fq.gz  NP_030.fq.gz  NP_055.fq.gz  NP_080.fq.gz
  ...
  NP_025.fq.gz  NP_050.fq.gz  NP_075.fq.gz  NP_100.fq.gz


These files are the primary raw output of a GAM experiment.
The first thing we need to do with the sequencing data is to "map" it to a
genome. The example data comes from mouse embryonic stem cells, so we need to
map it to the mouse genome, which we will do using bowtie2_. If you already
have bowtie2 and a mouse genome assembly installed and configured on your local
machine, you can skip the next step (mouse assembly mm9 is preferred, but any
other assembly should work with this tutorial).

Configuring bowtie2
-------------------

If you have not yet installed bowtie2_, please follow the
`installation instructions on the bowtie2 homepage <bowtie2install>`_. Once
you have bowtie installed, verify that everything is working correctly:

.. code-block:: bash

  $ bowtie2 --version
  /home/rob_000/bowtie2-2.2.9/bowtie2-align-s version 2.2.9
  64-bit
  Built on Windows8
  30 Apr 2016 18:13:39

We next need to provide the sequence of the mouse genome for bowtie to map
against. If you wish, you can download and configure the
`full mouse mm9 "index" <bowtie2index>`_ from the bowtie2 homepage. However,
the 100 sequencing datasets provided as part of the tutorial only contain
sequencing data from a small region of chromosome 19, so you can also
use a `special truncated index <tutorialindex>`_ containing only the
sequence of mouse chromosome 19. This will allow bowtie to run much faster
whilst using less RAM, and is perfectly sufficient for completing this tutorial.
If you wish to use the tutorial index, download it from the **GAMtools** website,
extract it to the same folder as fastqs and configure bowtie to use the new
truncated index:

.. code-block:: bash

  $ wget http://gam.tools/tutorial_index.tar.gz
  $ zxvf tutorial_index.tar.gz
  $ ls
  clean.sh fastqs/ genome/
  $ export BOWTIE2_INDEXES=$(pwd)/genome/
  $ ls $BOWTIE2_INDEXES
  genome.1.bt2  genome.3.bt2  genome.rev.1.bt2  mm9.chrom.sizes
  genome.2.bt2  genome.4.bt2  genome.rev.2.bt2

========================================================
Mapping the sequencing data and calling positive windows
========================================================

The **GAMtools** command used for mapping NP sequencing data is ``gamtools
process_nps``. The ``process_nps`` command has a lot of different parameters
and options, you can use the ``--help`` flag to get a full description of all
the available parameters. Further information about the ``process_nps`` 
command can also be found on the `process_NPs <tools/process_nps>`_ page.

.. code-block:: bash

  $ gamtools process_nps --help
  usage: gamtools process_nps [-h] -g GENOME_FILE [-o OUPUT_DIRECTORY]
                              [-f FITTINGS_DIRECTORY] [-d DETAILS_FILE] [-i]
                              [-b] [-c] [-w WINDOW_SIZE [WINDOW_SIZE ...]] [-m]
                              [-s MATRIX_SIZE [MATRIX_SIZE ...]]
                              [--qc-window-size QC_WINDOW_SIZE]
                              [--additional-qc-files [ADDITIONAL_QC_FILES [ADDITIONAL_QC_FILES ...]]]
                              [-q MINIMUM_MAPQ] [--doit-db-file DEP_FILE]
                              [--doit-backend {sqlite3,json,dbm}]
                              [--doit-verbosity {0,1,2}]
                              [--doit-reporter {json,console,zero,executed-only}]
                              [--doit-process NUM_PROCESS]
                              [--doit-parallel-type {process,thread}]
                              INPUT_FASTQ [INPUT_FASTQ ...]

For now, we can just use the default options. That means that all we need to specifiy
is a genome file (using ``-g/--genome-file``) and a list of input fastq files:

.. code-block:: bash

  $ gamtools process_nps -g genome/chr19.size fastqs/*.fq.gz

This tells **GAMtools** to use the genome file ``genome/chr19.size`` .You
will have this file if you downloaded the special truncated index. If you
are using your own mouse genome index, you will have to specify your own
genome file (which is usually named something like ``mm9.chrom.sizes``).
The next argument tells **GAMtools** to process all of the files with the
extension ".fq.gz" in the folder called "fastqs". When you run the
command, **GAMtools** will start mapping the sequencing data, and you
should see an output like this:

.. code-block:: bash

  $ gamtools process_nps -g genome/mm9.chrom.sizes fastqs/*.fq.gz
  -- Creating output directory
  .  Mapping fastq:fastqs/NP_025.fq.gz
  .  Mapping fastq:fastqs/NP_017.fq.gz
  .  Mapping fastq:fastqs/NP_065.fq.gz
  .  Mapping fastq:fastqs/NP_014.fq.gz
  .  Mapping fastq:fastqs/NP_090.fq.gz
  .  Mapping fastq:fastqs/NP_078.fq.gz
 
**GAMtools** will then proceed to map all 100 individual sequencing files to
the mouse genome. This will take around 5 minutes if you are using the
truncated index and a moderately fast computer. If you are using your own full
mouse genome index, it may take a little longer. Once it has mapped the files,
**GAMtools** will sort the mapped files, remove PCR duplicates and create an
index for fast data retrieval. The final steps are to compute the number of
reads from each NP that overlap each 50kb window in the supplied genome file,
and then to use this read coverage count to determine which of the windows was
present in the original NP.

  

.. _tutorialdata: http://gam.tools/tutorial_data.tar.gz
.. _tutorialindex: http://gam.tools/tutorial_index.tar.gz
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2
.. _bowtie2install: http://bowtie-bio.sourceforge.net/bowtie2

