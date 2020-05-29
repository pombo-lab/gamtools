#(c) 2012 Massachusetts Institute of Technology. All Rights Reserved
#Code written by Maksim Imakaev (imakaev@mit.edu)
"""
Cythonised functions for efficiently distance normalised proximity matrices.

Taken from the mirnylib library (https://bitbucket.org/mirnylab/mirnylib/)
with permission from Maxim Imakaev and Leonid Mirny.
"""
#@PydevCodeAnalysisIgnore
import numpy as np
cimport numpy as np
cimport cython

from libc.math cimport isnan

def logbinsnew(a, b, ratio=0, N=0):
    a = int(a)
    b = int(b)
    a10, b10 = np.log10([a, b])
    if ratio != 0:
        if N != 0:
            raise ValueError("Please specify N or ratio")
        N = int(np.log(b / a) / np.log(ratio))
    elif N == 0:
        raise ValueError("Please specify N or ratio")
    data10 = np.logspace(a10, b10, N)
    data10 = np.array(np.rint(data10), dtype=int)
    data10 = np.sort(np.unique(data10))
    assert data10[0] == a
    assert data10[-1] == b
    return data10


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def observedOverExpected(matrix):
    "Calculates observedOverExpected of any contact map. Ignores NaNs"

    cdef int i, j, k, start, end, count, offset
    cdef double x, ss, meanss

    cdef np.ndarray[np.double_t, ndim=2] data = np.array(matrix, dtype=np.double, order="C")
    cdef int N = data.shape[0]
    
    _bins = logbinsnew(1, N, 1.03)
    _bins = [(0, 1)] + [(_bins[i], _bins[i+1]) for i in range(len(_bins) - 1)]
    cdef np.ndarray[np.int64_t, ndim=2] bins = np.array(_bins, dtype=np.int64, order="C")
    cdef int M = bins.shape[0]
    
    for k in range(M):
        start, end = bins[k, 0], bins[k, 1]
        ss = 0
        count = 0
        for offset in range(start, end):
            for j in range(0, N - offset):
                x = data[offset + j, j]
                if isnan(x):
                    continue
                ss += x
                count += 1
        
        meanss = ss / count
        if meanss != 0:
            for offset in range(start,end):
                for j in range(0,N-offset):
                    data[offset + j, j] /= meanss
                    if offset > 0: data[j, offset+j] /= meanss
    return data


