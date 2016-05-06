############
enrichment
############

The ``gamtools enrichment`` tool is used to calculate the enrichment of pairwise
interactions between windows of different classes. For example, it can be used to
answer the question of whether a particular set of interactions connects windows
that contain genes with windows that contain enhancers more or less frequently
than would be expected by chance.

``gamtools enrichment`` requires two input files. The first is a tab-delimited table giving
the pairwise interactions. This table must contain the following columns:

+-----------------+-----------------------------------------------------+
| Column          | Description                                         |
+=================+=====================================================+
| **chrom**       | Name of the chromosome                              |
+-----------------+-----------------------------------------------------+
| **Pos_A**       | index of the window on the left of the interaction  |
+-----------------+-----------------------------------------------------+
| **Pos_B**       | index of the window on the right of the interaction |
+-----------------+-----------------------------------------------------+
| **interaction** | strength of the interaction                         |
+-----------------+-----------------------------------------------------+

Window indices used for **Pos_A** and **Pos_B** are 0-based, such that ``0``
would be the first window on the chromosome and ``20`` would be the 19th. An
example file might look like:

::

    chrom    Pos_A  Pos_B    interaction
    chr1     10     20       0.75
    chr2     10     20       0.50
    chr1     10     30       0.40

The second input file is a comma-delimited (csv) table giving the classes of
each window. This table can have the following columns:

+-----------+----------------------------------------------+
| Column    | Description                                  |
+===========+==============================================+
| **chrom** | Name of the chromosome                       |
+-----------+----------------------------------------------+
| **i**     | index of the window                          |
+-----------+----------------------------------------------+
| **start** | Start co-ordinate of the window *(optional)* |
+-----------+----------------------------------------------+
| **stop**  | Stop co-ordinate of the window *(optional)*  |
+-----------+----------------------------------------------+

Any additional columns will be interpreted as different classes. Each
class column should indicate whether the given window is a member of
the class with the values ``False`` or ``True``. An example classification
table might look like this:

::

  chrom,i,Enhancer,Gene
  chr1,10,True,False
  chr1,20,False,True
  chr2,20,True,False
  chr2,30,False,True

The output is a csv file in the following format:

+--------------+-------------------------------------------------------------------------------------------+
| Column       | Description                                                                               |
+==============+===========================================================================================+
| **class1**   | First class involved in interaction                                                       |
+--------------+-------------------------------------------------------------------------------------------+
| **class2**   | Second class involved in interaction                                                      |
+--------------+-------------------------------------------------------------------------------------------+
| **count**    | Number of interactions where a window in **class1** interacts with a window in **class2** |
+--------------+-------------------------------------------------------------------------------------------+
| **permuted** | Whether or not the interactions table was randomly permuted before counting               |
+--------------+-------------------------------------------------------------------------------------------+

For example:

::

  class1,class2,count,permuted
  Gene,Gene,234148,yes
  Gene,Enhancer,268228,yes
  Enhancer,Enhancer,10598,yes

===============================
Usage and option summary
===============================
**Usage**:
::

  gamtools enrichment [OPTIONS] -i <INTERACTIONS_FILE> -c <CLASSES_FILE> 

**Optional parameters:**

+---------------------+----------------------------------------------------------------------+
| Option              | Description                                                          |
+=====================+======================================================================+
| -o, --output-prefix | First part of the output file name (default is "enrichment_results") |
+---------------------+----------------------------------------------------------------------+
| -p, --permutations* | Number of times to randomly permute the input file                   |
+---------------------+----------------------------------------------------------------------+
| -n, --no-permute*   | Do not permute the input file, instead calculate observed counts     |
+---------------------+----------------------------------------------------------------------+

\* Options -p/-n are mutually exclusive, exactly one of these two options **must** be given
