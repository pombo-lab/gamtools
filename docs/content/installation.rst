############
Installation
############

**GAMtools** is intended to run in a command line environment on UNIX, LINUX
and Apple OS X operating systems. GAMtools can also be installed on Windows
using cygwin_.

The recommended way to install **GAMtools** is to use python's package manager,
``pip``. This method should ensure that all of GAMtools required dependencies
are installed automatically.

Alternatively, **GAMtools** can be installed by downloading the source code
and compiling it manually.

====================================
Installing stable releases using pip
====================================

To install the latest stable release using pip, you need to run the
following command:

.. code-block:: bash

  $ pip install gamtools

pip should automatically find and install any `mandatory dependencies`_
that are not currently installed. Additional `optional dependencies`_ are
required for full **GAMtools** functionality, but these must be
installed manually.

===============================
Installing from source (GitHub)
===============================

If you want to install the latest development version of **GAMtools**,
you will need to install from source code. First, clone the
`GAMtools repository`_ from GitHub:

.. code-block:: bash

  $ git clone https://github.com/pombo-lab/gamtools.git

Then install the downloaded package using pip:

.. code-block:: bash

  $ pip install gamtools/

Or if pip is not installed:

.. code-block:: bash

  $ cd gamtools
  $ python setup.py install

Installation using pip is the preferred method, as this will handle installing
the `mandatory dependencies`_ automatically. If **GAMtools** is installed using
``python setup.py install`` you may need to manually install
`mandatory dependencies`_ yourself.

===============
Troubleshooting
===============

GAMtools requires numpy_ and cython_ to be installed before it can compile
properly. If you are installing using ``pip``, numpy and cython should be installed
automatically, but there is a chance this might not work. If you are having issues
installing **GAMtools**, the first step is to ensure both numpy and cython are
properly installed:

.. code-block:: bash

  $ pip install cython numpy

If you are still having problems, please post a ticket on our `GitHub issues`_ page.

===============================
Mandatory dependencies
===============================

**GAMtools** depends on a number of additional python libraries, which must
be installed for it to function correctly. These libraries are normally
installed automatically during the **GAMtools** installation process.

Mandatory python dependencies
-----------------------------

  * doit_
  * numpy_
  * scipy_
  * cython_
  * pandas_
  * wrapit_

These python libraries can all be installed using pip:

.. code-block:: bash

  $ pip install doit numpy scipy cython pandas wrapit

===============================
Optional dependencies
===============================

Some features in **GAMtools** depend on additional libraries and/or
programs which are not installed automatically.

Making plots
------------

The ``gamtools matrix`` command requires some python plotting libraries to be
installed. These may also be required for the ``gamtools call_windows`` command
if the ``--fitting-folder`` flag is specified.

Optional python dependencies
----------------------------

  * matplotlib_

Working with raw sequencing data
--------------------------------

The ``gamtools process_nps`` command is used to map and process raw sequencing
data from NPs. This can require a number of additional command line programs
to be installed and configured:

Mapping and processing programs
-------------------------------

+-------------------+-------------------------------------------------------+
| Program           | Required for                                          |
+===================+=======================================================+
| Bowtie2_          | Mapping raw sequencing data.                          |
+-------------------+-------------------------------------------------------+
| samtools_         | Mapping raw sequencing data.                          |
+-------------------+-------------------------------------------------------+
| bedtools_         | Calling positive windows for an NP.                   |
+-------------------+-------------------------------------------------------+
| bedGraphToBigWig_ | Creating bigwigs (``--bigwigs`` flag)                 |
+-------------------+-------------------------------------------------------+
| bedToBigBed_      | Creating bigbeds (``--bigbeds`` flag)                 |
+-------------------+-------------------------------------------------------+
| fastqc_           | Performing dataset quality control (``--do-qc`` flag) |
+-------------------+-------------------------------------------------------+
| fastq_screen_     | Performing dataset quality control (``--do-qc`` flag) |
+-------------------+-------------------------------------------------------+

=========================
Testing your installation
=========================

To test that you have installed gamtools and all its dependencies correctly you
can run the command ``gamtools test``. If you have skipped installing any
optional dependencies, you may get a warning message saying something like "x
could not be found, and is required for y". You can safely ignore these
messages unless you need the particular gamtools functionality in the message.


.. _cygwin: https://cygwin.com
.. _GAMtools repository: https://github.com/pombo-lab/GAMtools
.. _doit: http://pydoit.org
.. _numpy: http://www.numpy.org
.. _scipy: http://www.scipy.org
.. _cython: http://cython.org
.. _pandas: http://pandas.pydata.org
.. _wrapit: https://github.com/rbeagrie/wrapit/
.. _matplotlib: http://matplotlib.org/
.. _pybedtools: https://pythonhosted.org/pybedtools/
.. _metaseq: https://pythonhosted.org/metaseq/
.. _bowtie2: http://bowtie-bio.sourceforge.net/bowtie2
.. _bedtools: https://bedtools.readthedocs.io/en/latest/index.html
.. _samtools: http://www.htslib.org/
.. _bedGraphToBigWig: http://hgdownload.cse.ucsc.edu/admin/exe/
.. _bedToBigBed: http://hgdownload.cse.ucsc.edu/admin/exe/
.. _fastqc: http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
.. _fastq_screen: http://www.bioinformatics.bbsrc.ac.uk/projects/fastq_screen/
.. _GitHub issues: https://github.com/pombo-lab/GAMtools/issues

