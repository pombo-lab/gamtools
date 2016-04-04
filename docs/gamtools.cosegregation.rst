gamtools.cosegregation module
==============================

.. _regions:

Regions
-------

Regions are :class:`pandas.DataFrame` objects, where columns represent samples
and rows represent windows. The values should be either 1 (signifying that the
given window was detected as present in the given sample) or 0 (signifying that
the window was not detected.)

.. _samples:

Samples
-------

Each column in a :ref:`region <regions>` is an integer
:class:`numpy arrays <numpy.ndarray>` of length x, representing the segregation
of a particular genomic window across x samples.



.. automodule:: gamtools.cosegregation
    :members:
    :undoc-members:
    :show-inheritance:
