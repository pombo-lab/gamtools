#(c) 2012 Massachusetts Institute of Technology. All Rights Reserved
#Code written by Maksim Imakaev (imakaev@mit.edu)
"""
Cythonised utilities from numutils
"""
#@PydevCodeAnalysisIgnore
import numpy as np
cimport numpy as np
cimport cython
from math import log
import subprocess
import os



ctypedef unsigned long ulong
ctypedef unsigned int uint
ctypedef unsigned short ushort
ctypedef unsigned char uchar


from libcpp.string cimport string
from libc.math cimport isnan


cdef extern from "<sstream>" namespace "std":
    cdef cppclass stringstream:
        string str()
        string str(string val)
        stringstream& operator<< (double val)
        stringstream& operator<< (string val)
        stringstream& operator<< (char* val)



cdef extern from "stdlib.h":
    long c_libc_random "random"()
    double c_libc_drandom "drand48"()


cdef extern from "stdio.h":
    int sprintf(char *str, char *format, ...)


ctypedef fused my_type:
    cython.int
    cython.float
    cython.char
    cython.long
    cython.short
    cython.double
    cython.complex
    ulong
    uint
    ushort
    uchar


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
def fakeCisImpl(np.ndarray[np.double_t, ndim = 2] data, np.ndarray[np.int64_t,ndim = 2] mask):
    cdef int N
    N = len(data)
    cdef int i,j,r,s
    for i in range(N):
        for j in range(i,N):
            if mask[i,j] == 1:
                while True:
                    r = c_libc_random() % 2
                    if (r == 0):
                        s = c_libc_random() % N
                        if mask[i,s] == 0:
                            data[i,j] = data[i,s]
                            data[j,i] = data[i,s]
                            break
                    else:
                        s = c_libc_random() % N
                        if mask[j,s] == 0:
                            data[i,j] = data[j,s]
                            data[j,i] = data[j,s]



def logbinsnew(a, b, ratio=0, N=0):
    a = int(a)
    b = int(b)
    a10, b10 = np.log10([a, b])
    if ratio != 0:
        if N != 0:
            raise ValueError("Please specify N or ratio")
        N = np.log(b / a) / np.log(ratio)
    elif N == 0:
        raise ValueError("Please specify N or ratio")
    data10 = np.logspace(a10, b10, N)
    data10 = np.array(np.rint(data10), dtype=int)
    data10 = np.sort(np.unique(data10))
    assert data10[0] == a
    assert data10[-1] == b
    return data10



@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef my_type[:] _boolIndex (my_type[:] array, uchar[:] indexes, my_type[:] output, bounds):
    cdef int N,M,j,i
    j = 0
    N = len(array)
    M = len(output)
    if bounds == True:
        for i in range(N):
            if indexes[i] == True:
                if j == M:
                    raise ValueError("Array out of bounds")
                output[j] = array[i]

                j += 1
        if j != M:
            raise ValueError("Out array is too big")
    else:
        for i in range(N):
            if indexes[i] == True:
                output[j] = array[i]
                j += 1
    return output


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def createRWCython(N, Ndim=2):
    cdef np.ndarray[np.int64_t,cast = True,ndim = 2] offsets = np.zeros((N, Ndim), dtype=int)
    cdef np.ndarray[np.int64_t,cast=True,ndim=1] rand1 = np.random.randint(0,Ndim,N)
    cdef np.ndarray[np.int64_t,cast=True,ndim=1] rand2 = np.random.randint(0,2,N)
    cdef int i
    for i in xrange(N):
        offsets[i,rand1[i]] = rand2[i] * 2 - 1
    return np.cumsum(offsets, axis=0)


def fasterBooleanIndexing(np.ndarray array, np.ndarray indexes,output = None,outLen = None, bounds = True):
    """
    A faster way of writing "output = a[my_boolean_array]"
    Is approximately 30-50 % faster, but requires pre-allocation of the output array.

    .. warning :: this function relies on the fact that you know the length of the output array,
    and either supply it as an output array, or provide outLen.
    If not supplied, it will be estimated, and function will be almost as slow as numpy indexing

    Parameters
    ----------
    array : numpy.array compatible
        array to be indexed
    indexes : numpy.array compatible
        array if indexes of the same length as array
    output : numpy.array compatible or None, optional
        output array. If not provided, please provide outLen parameter
        Note that if you create it, use numpy.empty, not numpy.zeros!
    outLen : int, optional
        length of output array that must be equal to indexes.sum()
    bounds : bool, optional
        check "on the fly" if length of output is indeed equal to indexes.sum()

    """
    if (len(array) < 10000) or (array.dtype == np.bool): return array[indexes]
    _indexes = np.asarray(indexes,dtype = bool)
    indexes = _indexes.view(dtype = np.uint8)

    cdef int N,M,i,j
    N = array.shape[0]

    if output == None:
        if outLen == None:
            outLen = indexes.sum()
        output = np.empty(outLen, dtype = array.dtype)
    else:
        output = np.asarray(output)
    M = output.shape[0]
    if indexes.shape[0] != array.shape[0]:
        raise ValueError("Length mismatch: indexes - {0}, array - {1}".format(len(indexes),len(array)))

    if array.dtype == np.float64:
        result = _boolIndex[cython.double](array,indexes,output,bounds)
    elif array.dtype == np.int64:
        result = _boolIndex[cython.long](array,indexes,output,bounds)
    elif array.dtype == np.int32:
        result = _boolIndex[cython.int](array,indexes,output,bounds)
    elif array.dtype == np.int16:
        result = _boolIndex[cython.short](array,indexes,output,bounds)
    elif array.dtype == np.float32:
        result = _boolIndex[cython.float](array,indexes,output,bounds)
    elif array.dtype == np.int8:
        result = _boolIndex[cython.char](array,indexes,output,bounds)
    elif array.dtype == np.uint64:
        result = _boolIndex[ulong](array,indexes,output,bounds)
    elif array.dtype == np.uint32:
        result = _boolIndex[uint](array,indexes,output,bounds)
    elif array.dtype == np.uint16:
        result = _boolIndex[ushort](array,indexes,output,bounds)
    elif array.dtype == np.uint8:
        result = _boolIndex[uchar](array,indexes,output,bounds)
    else:
        print array.dtype
        raise ValueError('data type not implemented')
    return output






@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
def _arrayInArray(np.ndarray array,np.ndarray filterarray):
    "Actual implementation of arrayInArray"
    cdef np.ndarray[np.uint8_t,cast = True,ndim = 1] mask = np.zeros(len(array),'bool')
    cdef np.ndarray[np.int64_t,ndim = 1] args = np.argsort(array)
    arsort = array.take(args,axis = 0)
    cdef np.ndarray[np.int64_t,ndim = 1] diffs = np.r_[0,np.nonzero(np.diff(arsort) )[0]+1,len(arsort)]  #places where sorted values of an array are changing
    values = arsort.take(diffs[:len(diffs)-1])  #values at that places
    cdef np.ndarray[np.int64_t,ndim = 1] allinds = np.searchsorted(values[:len(values)-1],filterarray)
    cdef np.ndarray[np.uint8_t,cast = True,ndim = 1] exist = values.take(allinds) == filterarray                 #check that value in filterarray exists in array
    cdef int N = len(allinds)
    cdef int i,j
    N,args,exist  #Eclipse warning remover
    for i in range(N):
        if exist[i] == 0:
            continue
        for j  in range(diffs[allinds[i]],diffs[allinds[i]+1]):
            mask[args[j]] = 1
    return mask

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.nonecheck(False)
@cython.cdivision(True)
def _arraySumByArray(array,filterarray,weightarray):
    "faster [sum(weightarray [array == i]) for i in filterarray]"
    array = np.asarray(array,dtype = float)
    cdef np.ndarray[np.float64_t,ndim = 1] _weightarray = np.asarray(weightarray,dtype = float)
    if len(array) != len(_weightarray): raise ValueError
    cdef np.ndarray[np.int64_t,ndim = 1] args = np.argsort(array)
    arsort = array.take(args)
    cdef np.ndarray[np.int64_t,ndim = 1] diffs = np.r_[0,np.nonzero(np.diff(arsort) > 0.)[0]+1,len(arsort)]
    values = arsort.take(diffs[:len(diffs)-1])
    cdef np.ndarray[np.int64_t,ndim = 1] allinds = np.searchsorted(values[:len(values)-1],filterarray)
    cdef np.ndarray[np.uint8_t,cast = True,ndim = 1] exist = values.take(allinds) == filterarray
    cdef int i,j,N
    N = len(allinds)
    cdef np.ndarray[np.float64_t,ndim = 1] ret = np.zeros(len(allinds),float)
    for i in range(N):
        if (exist[i] == 0): continue
        for j in range(diffs[allinds[i]],diffs[allinds[i]+1]):
            ret[i] += _weightarray[args[j]];


    return ret

@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
@cython.nonecheck(False)

def ultracorrectSymmetricWithVector(x,v = None,M=None,diag = -1,
                                    tolerance=1e-5):
    """Main method for correcting DS and SS read data. Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction"""
    print("Hi")
    if M == None:
        M = 599
    totalBias = np.ones(len(x),float)
    if v == None: v = np.zeros(len(x),float)  #single-sided reads
    x = np.array(x,np.double,order = 'C')
    cdef np.ndarray[np.double_t, ndim = 2] _x = x
    cdef np.ndarray[np.double_t, ndim = 1] s
    v = np.array(v,float,order = "C")
    cdef int i , j, N
    N = len(x)

    for iternum in range(M):
        s0 = np.sum(_x,axis = 1)

        mask = s0 == 0

        v[mask] = 0   #no SS reads if there are no DS reads here
        nv = v / (totalBias * (totalBias[mask==False]).mean())
        s = s0 + nv
        for dd in range(diag + 1):   #excluding the diagonal
            if dd == 0:
                s -= np.diagonal(_x)
            else:
                dia = np.array(np.diagonal(_x,dd))
                s[dd:] = s[dd:] -  dia
                s[:len(s)-dd] = s[:len(s)-dd] - dia
        s = s / np.mean(s[s0!=0])
        s[s0==0] = 1
        s -= 1
        s *= 0.8
        s += 1
        totalBias *= s

        for i in range(N):
            for j in range(N):
                _x[i,j] = _x[i,j] / (s[i] * s[j])

        if M == 599:
            if np.abs(s-1).max() < tolerance:
                #print "IC used {0} iterations".format(iternum+1)
                break


    corr = totalBias[s0!=0].mean()  #mean correction factor
    x  = x * corr * corr #renormalizing everything
    totalBias /= corr
    return x, v/totalBias, totalBias

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def ultracorrectAny(x,v = None,M=None,diag = -1, attemptInPlace=True,
                                    tolerance=1e-5):
    """Main method for correcting DS and SS read data. Possibly excludes diagonal.
    By default does iterative correction, but can perform an M-time correction"""
    if M == None:
        M = 599

    totalBias1 = np.ones(len(x),float)
    if v == None: v = np.zeros(len(x),float)  #single-sided reads
    if attemptInPlace:
        x = np.asarray(x,np.double,order = 'C')
    else:
        x = np.array(x,np.double,order = 'C')
    N1, N2 = x.shape
    cdef np.ndarray[np.double_t, ndim = 2] _x = x
    cdef np.ndarray[np.double_t, ndim = 1] s
    v = np.array(v,float,order = "C")
    cdef int i , j, N
    N = len(x)
    for iternum in xrange(M):
        s0 = np.sum(_x,axis = 1)
        mask = [s0 == 0]
        v[mask] = 0   #no SS reads if there are no DS reads here
        nv = v / (totalBias * (totalBias[mask==False]).mean())
        s = s0 + nv
        for dd in xrange(diag + 1):   #excluding the diagonal
            if dd == 0:
                s -= np.diagonal(_x)
            else:
                dia = np.diagonal(_x,dd)
                print dia
                s[dd:] = s[dd:] -  dia
                s[:-dd] = s[:-dd] - dia
        s = s / np.mean(s[s0!=0])
        s[s0==0] = 1
        s -= 1
        s *= 0.8
        s += 1
        totalBias *= s

        for i in range(N):
            for j in range(N):
                _x[i,j] = _x[i,j] / (s[i] * s[j])

        if M == 599:
            if np.abs(s-1).max() < tolerance:
                print "IC used {0} iterations".format(iternum+1)
                break


    corr = totalBias[s0!=0].mean()  #mean correction factor
    x  = x * corr * corr #renormalizing everything
    totalBias /= corr
    return x, v/totalBias, totalBias


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
def ultracorrectSymmetricByMask(x, mask, M=None, tolerance=1e-5):
    """performs iterative correction excluding some regions of a heatmap from consideration.
    These regions are still corrected, but don't contribute to the sums
    """
    support = """
    #include <math.h>
    """
    assert mask.shape == x.shape
    assert mask.dtype == np.bool
    cdef np.ndarray[np.double_t , ndim = 2] _x
    cdef np.ndarray[np.uint8_t , ndim = 2] _mask
    x = np.array(x,np.double,order = 'C')
    if not  np.abs(x-x.T).mean() / (1. * np.abs(x.mean())) < 1e-10:
        raise ValueError("Please provide symmetric matrix!")
    _x = x
    cdef np.ndarray[np.double_t,ndim = 1] sums
    cdef np.ndarray[np.double_t,ndim = 1] counts
    _mask = mask.view(dtype = np.uint8)
    cdef int i, j, N
    if M == None:
        M = 9999
    N = len(x)
    allsums = np.ones(N,float)
    cdef double ss, count
    for iteration in range(M):
        sums = np.zeros(N,float)
        counts = np.zeros(N,float)

        for i in range(N):
            for j in range(N):
                if _mask[i,j] == 1:
                    sums[i] += _x[i,j]
                    counts[i] += 1
        ss = 0
        count = 0
        for i in range(N):
            if (counts[i] > 0) and (sums[i] > 0):
                sums[i] = sums[i] / counts[i]
                ss = ss + sums[i]
                count += 1
            else:
                sums[i] = -1

        for i in range(N):
            if sums[i] == -1:
                sums[i] = 1
            else:
                sums[i] = sums[i] / (ss / count)
        for i in range(N):
            for j in range(N):
                _x[i,j] = _x[i,j] / (sums[i] * sums[j])

        allsums *= sums
        if M == 9999:
            if np.abs(sums-1).max() < tolerance:
                print "IC used {0} iterations".format(iteration+1)
                break


    return x,allsums


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


@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)

def observedOverExpectedWithMask(matrix,mask):
    "Calculates observedOverExpected of any contact map, over regions where mask==1"

    cdef int i, j, bin, start, end, count, offset
    cdef double ss, meanss

    _data = np.array(matrix, dtype = np.double, order = "C")
    N = _data.shape[0]

    _datamask = np.array(mask==1, dtype = np.double, order = "C")


    cdef np.ndarray[np.double_t, ndim = 2] data = _data

    cdef np.ndarray[np.double_t, ndim = 2] datamask = _datamask

    _bins = logbinsnew(1,N,1.03)
    _bins = [(0,1)] + [(_bins[i],_bins[i+1]) for i in xrange(len(_bins)-1)]
    _bins = np.array(_bins,dtype = np.int64, order = "C")
    cdef np.ndarray[np.int64_t, ndim = 2] bins = _bins
    cdef int M = len(bins)

    for bin in range(M):
        start = bins[bin,0]
        end  = bins[bin,1]
        ss = 0
        count = 0
        for offset in range(start,end):
            for j in range(0,N-offset):
                if datamask[offset+j,j]==1:
                    ss += data[offset+j, j]
                    count += 1
        #print start, end, count
        meanss = ss / count
        if meanss != 0:
            for offset in range(start,end):
                for j in range(0,N-offset):
                    if datamask[offset+j,j]==1:
                        data[offset + j, j] /= meanss
                        if offset > 0: data[j,offset+j] /= meanss

    return _data



@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)
@cython.cdivision(True)

def calculateExpectedWithMask(matrix,mask):
    "Calculates observedOverExpected of any contact map, over regions where mask==1"

    cdef int i, j, bin, start, end, count, offset
    cdef double ss, meanss

    _data = np.array(matrix, dtype = np.double, order = "C")
    N = _data.shape[0]

    _datamask = np.array(mask==1, dtype = np.double, order = "C")


    cdef np.ndarray[np.double_t, ndim = 2] data = _data

    cdef np.ndarray[np.double_t, ndim = 2] datamask = _datamask

    _bins = logbinsnew(1,N,1.03)
    _bins = [(0,1)] + [(_bins[i],_bins[i+1]) for i in xrange(len(_bins)-1)]
    _bins = np.array(_bins,dtype = np.int64, order = "C")
    cdef np.ndarray[np.int64_t, ndim = 2] bins = _bins
    cdef int M = len(bins)

    for bin in range(M):
        start = bins[bin,0]
        end  = bins[bin,1]
        ss = 0
        count = 0
        for offset in range(start,end):
            for j in range(0,N-offset):
                if datamask[offset+j,j]==1:
                    ss += data[offset+j, j]
                    count += 1
        #print start, end, count
        meanss = ss / count
        if meanss != 0:
            for offset in range(start,end):
                for j in range(0,N-offset):
                    if datamask[offset+j,j]==1:
                        data[offset + j, j] = meanss
                        if offset > 0: data[j,offset+j] = meanss

    return _data




@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)


def commandExists(command):
    "checks if the bash command exists"
    command = command.split()[0]
    if subprocess.call(['which', command]) != 0:
        return False
    return True

def gzipWriter(filename):
    """
    creates a writing process with gzip or parallel gzip (pigz) attached to it
    """
    filename = os.path.abspath(filename)
    with open(filename, 'wb') as outFile:
        if commandExists("pigz"):
            writer = ["pigz", "-c", "-9"]
        else:
            writer = ["gzip", "-c", "-2"]

        pwrite = subprocess.Popen(writer, stdin=subprocess.PIPE, stdout=outFile, shell=False, bufsize=-1)
    return pwrite

@cython.boundscheck(False)
@cython.nonecheck(False)
@cython.wraparound(False)

def matrixToGzippedFile(inMatrix, filename,formatString="%.3lf "):
    cdef int N = len(inMatrix)
    cdef int M = len(inMatrix[0])
    writer = gzipWriter(filename)
    cdef int i,j
    outPipe = writer.stdin
    cdef np.ndarray[np.double_t, ndim = 2] newMatrix  =  np.array(inMatrix, dtype=np.double, order="C")

    emptyString = " "*100

    cdef char* c_string  = emptyString
    cdef string s
    s.reserve(50*M)
    cdef double element
    cdef string curString
    cdef  char* curStringTemplate = formatString
    for i in xrange(N):
        s.resize(0)
        for j in xrange(M):
            element = newMatrix[i,j]
            sprintf(c_string, curStringTemplate, element)
            s.append(c_string)
        if i != N-1:
            curString = "\n"
            s.append(curString)

        outPipe.write(s)
    writer.communicate()


