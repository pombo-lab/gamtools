import numpy as np
# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).
cimport numpy as np
# We now need to fix a datatype for our arrays. I've used the variable
# DTYPE for this, which is assigned to the usual NumPy runtime
# type info object.
DTYPE = np.int
# "ctypedef" assigns a corresponding compile-time type to DTYPE_t. For
# every type in the numpy module there's a corresponding compile-time
# type with a _t-suffix.
ctypedef np.int_t DTYPE_t

def cosegregation_frequency(np.ndarray[DTYPE_t, ndim=2] samples):
    """Take a table of n columns and return the co-segregation frequencies"""


    cdef np.ndarray[DTYPE_t, ndim=2] counts = np.zeros([2, 2], dtype=DTYPE)
    cdef np.ndarray s
    cdef int x
    cdef int y

    for s in samples:

        x = s[0]
        y = s[1]
        counts[x,y] += 1

    return counts

