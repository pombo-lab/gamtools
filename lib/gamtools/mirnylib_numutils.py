"""
Functions for efficiently distance normalised proximity matrices.

Taken from the mirnylib library (https://bitbucket.org/mirnylab/mirnylib/)
with permission from Maxim Imakaev and Leonid Mirny.
"""
# pylint: disable-all
from . import mirnylib_numutils_internal


def fillDiagonal(inArray, diag, offset=0):
    "Puts diag in the offset's diagonal of inArray"
    N = inArray.shape[0]
    assert inArray.shape[1] == N
    if offset >= 0:
        inArray.flat[offset:N * (N - offset):N + 1] = diag
    else:
        inArray.flat[(-offset * N)::N + 1] = diag


def removeDiagonals(inArray, m):
    """removes up to mth diagonal in array
    m = 0: main
    m = 1: 3 diagonals, etc.
    """
    for i in range(-m, m + 1):
        fillDiagonal(inArray, 0, i)



def observedOverExpected(matrix):
    """
    Parameters
    ----------
    matrix : a symmetric contactmap
        A Hi-C contactmap to calculate observed over expected.

    Returns
    -------
        matrix : a symmetric corrected contactmap

    .. note:: This function does not work in place; it returns a copy.

    It divides each diagonal of a Hi-C contact map by its' mean.
    It also does it in a smart way: it calculates averages
    over stripes from X to X*1.05, and divides each stripe by its mean.

    It allows to avoid divergence far from the main diagonal with a very few reads.
    """
    return mirnylib_numutils_internal.observedOverExpected(matrix)
