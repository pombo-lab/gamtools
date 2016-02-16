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

@cython.boundscheck(False)
@cython.wraparound(False)
def cosegregation_2d(np.ndarray[DTYPE_t, ndim=2] a,
                     np.ndarray[DTYPE_t, ndim=2] b):

    cdef int aN = a.shape[0]
    cdef int bN = b.shape[0]

    cdef int M = a.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=4] temp = np.zeros([aN, bN, 2, 2], dtype=DTYPE)
    cdef np.ndarray[np.float64_t, ndim=2] result = np.zeros([aN, bN], dtype=np.float64)

    cdef int ai
    cdef int bi
    cdef int i

    cdef int ax
    cdef int bx

    cdef int AB

    for ai in range(aN):
        for bi in range(bN):

            for i in range(M):

                ax = a[ai,i]
                bx = b[bi,i]

                temp[ai, bi, ax, bx] += 1

            AB = temp[ai, bi, 1, 1]

            result[ai, bi] = AB

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
def cosegregation_3d(np.ndarray[DTYPE_t, ndim=2] a,
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

    cdef int ABC

    
    for ai in range(aN):
        for bi in range(bN):
            for ci in range(cN):

                for i in range(M):
                    
                    ax = a[ai,i]
                    bx = b[bi,i]
                    cx = c[ci,i]
                    
                    temp[ai, bi, ci, ax, bx, cx] += 1

                ABC = temp[ai, bi, ci, 1, 1, 1]

                result[ai, bi, ci] = ABC

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
def linkage_2d(np.ndarray[DTYPE_t, ndim=2] a,
               np.ndarray[DTYPE_t, ndim=2] b):

    cdef int aN = a.shape[0]
    cdef int bN = b.shape[0]

    cdef int M = a.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=4] temp = np.zeros([aN, bN, 2, 2], dtype=DTYPE)
    cdef np.ndarray[np.float64_t, ndim=2] result = np.zeros([aN, bN], dtype=np.float64)

    cdef int ai
    cdef int bi
    cdef int i

    cdef int ax
    cdef int bx

    cdef int A
    cdef int B
    cdef int AB

    cdef float ABp
    cdef float Ap
    cdef float Bp
    cdef float Ep
    cdef float d

    for ai in range(aN):
        for bi in range(bN):

            for i in range(M):

                ax = a[ai,i]
                bx = b[bi,i]

                temp[ai, bi, ax, bx] += 1

            AB = temp[ai, bi, 1, 1]
            A = (temp[ai, bi, 1, 1] +
                 temp[ai, bi, 1, 0])
            B = (temp[ai, bi, 1, 1] +
                 temp[ai, bi, 0, 1])

            if (A == 0) or (B == 0):
                result[ai, bi] = np.NaN

            else:
                ABp = AB / <float>M
                Ap = A / <float>M
                Bp = B / <float>M
                Ep = Ap * Bp

                d = ABp - Ep

                result[ai, bi] = d

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
def linkage_3d(np.ndarray[DTYPE_t, ndim=2] a,
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
    
    cdef int AB
    cdef int AC
    cdef int BC

    cdef int ABC

    cdef float ABCp
    cdef float Ap
    cdef float Bp
    cdef float Cp
    cdef float ABp
    cdef float BCp
    cdef float ACp
    cdef float dAB
    cdef float dAC
    cdef float dBC
    cdef float d
    
    for ai in range(aN):
        for bi in range(bN):
            for ci in range(cN):

                for i in range(M):
                    
                    ax = a[ai,i]
                    bx = b[bi,i]
                    cx = c[ci,i]
                    
                    temp[ai, bi, ci, ax, bx, cx] += 1

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

                AB = (temp[ai, bi, ci, 1, 1, 1] +
                      temp[ai, bi, ci, 1, 1, 0])
                AC = (temp[ai, bi, ci, 1, 1, 1] +
                      temp[ai, bi, ci, 1, 0, 1])
                BC = (temp[ai, bi, ci, 1, 1, 1] +
                      temp[ai, bi, ci, 0, 1, 1])

                ABC = temp[ai, bi, ci, 1, 1, 1]

                if (A == 0) or (B == 0) or (C==0):
                    result[ai, bi, ci] = np.NaN

                else:
                    ABCp = ABC / <float>M
                    Ap = A / <float>M
                    Bp = B / <float>M
                    Cp = C / <float>M
                    ABp = AB / <float>M
                    ACp = AC / <float>M
                    BCp = BC / <float>M
                    dAB = ABp - (Ap * Bp)
                    dAC = ACp - (Ap * Cp)
                    dBC = BCp - (Bp * Cp)

                    d = ABCp - (Ap * dBC) - (Bp * dAC) - (Cp * dAB) - (Ap * Bp * Cp)

                    result[ai, bi, ci] = d

    return result

@cython.boundscheck(False)
@cython.wraparound(False)
def dprime_2d(np.ndarray[DTYPE_t, ndim=2] a,
              np.ndarray[DTYPE_t, ndim=2] b):
    
    cdef int aN = a.shape[0]
    cdef int bN = b.shape[0]
    
    cdef int M = a.shape[1]

    cdef np.ndarray[DTYPE_t, ndim=4] temp = np.zeros([aN, bN, 2, 2], dtype=DTYPE)
    cdef np.ndarray[np.float64_t, ndim=2] result = np.zeros([aN, bN], dtype=np.float64)

    cdef int ai
    cdef int bi
    cdef int i

    cdef int ax
    cdef int bx

    cdef int A
    cdef int B
    cdef int AB

    cdef float ABp
    cdef float Ap
    cdef float Bp
    cdef float Ep
    cdef float d
    
    for ai in range(aN):
        for bi in range(bN):

            for i in range(M):
                
                ax = a[ai,i]
                bx = b[bi,i]
                
                temp[ai, bi, ax, bx] += 1

            AB = temp[ai, bi, 1, 1]
            A = (temp[ai, bi, 1, 1] +
                 temp[ai, bi, 1, 0])
            B = (temp[ai, bi, 1, 1] +
                 temp[ai, bi, 0, 1])

            if (A == 0) or (B == 0):
                result[ai, bi] = np.NaN

            else:
                ABp = AB / <float>M
                Ap = A / <float>M
                Bp = B / <float>M
                Ep = Ap * Bp

                d = ABp - Ep

                if d > 0.:

                    if (Ap <= Bp):

                        result[ai, bi] = d / (Ap - Ep)

                    elif (Bp <= Ap):

                        result[ai, bi] = d / (Bp - Ep)

                elif d < 0.:

                    result[ai, bi] = d / Ep

                else:

                    result[ai, bi] = 0.

    return result

