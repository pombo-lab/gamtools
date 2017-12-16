cdef extern from "slice.cpp":
    int run_slice()

def slice():
    return run_slice()
