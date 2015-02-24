import numpy as np

# "cimport" is used to import special compile-time information
# about the numpy module (this is stored in a file numpy.pxd which is
# currently part of the Cython distribution).

cimport numpy as np
cimport cython

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

def cosegregation_frequency_3(np.ndarray[DTYPE_t, ndim=2] samples):

    """Take a table of n columns and return the co-segregation frequencies"""

    cdef np.ndarray[DTYPE_t, ndim=3] counts = np.zeros([2, 2, 2], dtype=DTYPE)

    cdef np.ndarray s

    cdef int x

    cdef int y

    cdef int z


    for s in samples:


        x = s[0]

        y = s[1]

        z = s[2]

        counts[x,y,z] += 1


    return counts


def coseg_explicit_2_cy(np.ndarray[DTYPE_t, ndim=1] seg1,
                        np.ndarray[DTYPE_t, ndim=1] seg2):

    cdef np.ndarray[DTYPE_t, ndim=2] counts = np.zeros([2, 2], dtype=DTYPE)

    cdef int i
    cdef int x
    cdef int y

    for i in range(len(seg1)):

        x = seg1[i]
        y = seg2[i]

        counts[x, y] += 1

    return counts


def coseg_explicit_3_cy(np.ndarray[DTYPE_t, ndim=1] seg1,
                        np.ndarray[DTYPE_t, ndim=1] seg2,
                        np.ndarray[DTYPE_t, ndim=1] seg3):

    cdef np.ndarray[DTYPE_t, ndim=3] counts = np.zeros([2, 2, 2], dtype=DTYPE)

    cdef int i

    cdef int x

    cdef int y

    cdef int z

    for i in range(len(seg1)):

        x = seg1[i]
        y = seg2[i]
        z = seg3[i]

        counts[x, y, z] += 1

    return counts

@cython.boundscheck(False)
@cython.wraparound(False)
def coseg_all_3(np.ndarray[DTYPE_t, ndim=2] a,
                np.ndarray[DTYPE_t, ndim=2] b,
                np.ndarray[DTYPE_t, ndim=2] c):
    
    cdef int aN = a.shape[0]
    cdef int bN = b.shape[0]
    cdef int cN = c.shape[0]
    
    cdef int M = a.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=6] temp = np.zeros([aN, bN, cN, 2, 2, 2], dtype=DTYPE)
    cdef np.ndarray[np.float64_t, ndim=3] result = np.zeros([aN, bN, cN], dtype=np.float64)

    cdef int ai
    cdef int bi
    cdef int ci
    cdef int i

    cdef int ax
    cdef int bx
    cdef int cx

    cdef int A
    cdef int B
    cdef int C
    cdef int ABC

    cdef float ABCp
    cdef float Ap
    cdef float Bp
    cdef float Cp
    cdef float Ep
    cdef float d
    
    for ai in range(aN):
        for bi in range(bN):
            for ci in range(cN):

                for i in range(M):
                    
                    ax = a[ai,i]
                    bx = b[bi,i]
                    cx = c[ci,i]
                    
                    temp[ai, bi, ci, ax, bx, cx] += 1

                ABC = temp[ai, bi, ci, 1, 1, 1]
                A = (temp[ai, bi, ci, 1, 1, 1] +
                     temp[ai, bi, ci, 1, 0, 1] +
                     temp[ai, bi, ci, 1, 1, 0] +
                     temp[ai, bi, ci, 1, 0, 0])
                B = (temp[ai, bi, ci, 1, 1, 1] +
                     temp[ai, bi, ci, 0, 1, 1] +
                     temp[ai, bi, ci, 1, 1, 0] +
                     temp[ai, bi, ci, 0, 1, 0])
                C = (temp[ai, bi, ci, 1, 1, 1] +
                     temp[ai, bi, ci, 1, 0, 1] +
                     temp[ai, bi, ci, 0, 1, 1] +
                     temp[ai, bi, ci, 0, 0, 1])

                if (ABC == 0) or (A == 0) or (B == 0) or (C==0):
                    result[ai, bi, ci] = np.NaN

                else:
                    ABCp = ABC / <float>M
                    Ap = A / <float>M
                    Bp = B / <float>M
                    Cp = C / <float>M
                    Ep = Ap * Bp * Cp

                    d = ABCp - Ep

                    if d > 0.:

                        if (Ap <= Bp) and (Ap <= Cp):

                            result[ai, bi, ci] = d / (Ap - Ep)

                        elif (Bp <= Cp) and (Bp <= Ap):

                            result[ai, bi, ci] = d / (Bp - Ep)

                        elif (Cp <= Ap) and (Cp <= Bp):

                            result[ai, bi, ci] = d / (Cp - Ep)

                    elif d < 0.:

                        result[ai, bi, ci] = d / Ep

                    else:

                        result[ai, bi, ci] = 0.



    return result
