#!usr/bin/evn python
from cffi import FFI
import os
import numpy as np
import time
import gc


f = open('lu.c', 'r')
src = f.read()
f.close()

ffi = FFI()
ffi.cdef('void det_by_lu(double *y, double *B, int N);')
_C = ffi.verify(src, extra_compile_args=['-fopenmp', '-D use_openmp', '-Ofast','-march=native','-ffast-math', '-ftree-vectorizer-verbose=1'], extra_link_args=['-fopenmp'])
c_det_by_lu = _C.det_by_lu


f = open('lu_no_simd.c', 'r')
src = f.read()
f.close()

ffi = FFI()
ffi.cdef('void det_by_lu(double *y, double *B, int N);')
_C2 = ffi.verify(src)
c_det_by_lu2 = _C2.det_by_lu



x = np.random.randn(2**11, 2**11)
y = x.copy()

def run_c(A,B,y,N):
    # run c code
    np.copyto(B,A)
    
    # check that result is correct
    L = np.tril(B, -1) + np.eye(N)
    U = np.triu(B)
    #assert_almost_equal( L.dot(U),  A)
    
    pa = ffi.cast("double *", A.ctypes.data)
    pb = ffi.cast("double *", B.ctypes.data)
    gc.disable()
    st = time.time()
    
    loops = 1 + min(1000000 // (N*N), 20000)
    print 'loops', loops

    for l in range(loops):
        np.copyto(B,A)
        #c_det_by_lu(ffi.cast("double *", y.ctypes.data), ffi.cast("double *", B.ctypes.data), ffi.cast("int", N))
        c_det_by_lu(pa, pb, N)
        
    et = time.time()
    gc.enable()
   
    return  (et - st)/loops


def run_c2(A,B,y,N):
    # run c code
    np.copyto(B,A)
    
    # check that result is correct
    L = np.tril(B, -1) + np.eye(N)
    U = np.triu(B)
    #assert_almost_equal( L.dot(U),  A)
    
    pa = ffi.cast("double *", A.ctypes.data)
    pb = ffi.cast("double *", B.ctypes.data)
    gc.disable()
    st = time.time()
    
    loops = 1 + min(1000000 // (N*N), 20000)
    print 'no simd loops', loops

    for l in range(loops):
        np.copyto(B,A)
        c_det_by_lu2(pa, pb, N)
        
    et = time.time()
    gc.enable()
   
    return  (et - st)/loops


print run_c(x,y,x,2**11)
print run_c2(x,y,x,2**11)

