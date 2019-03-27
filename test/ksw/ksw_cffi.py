#!usr/bin/env python
from cffi import FFI
import os


ffi = FFI()


ffi.cdef("""
#define KSW_XBYTE  0x10000
#define KSW_XSTOP  0x20000
#define KSW_XSUBO  0x40000
#define KSW_XSTART 0x80000

struct _kswq_t;
typedef struct _kswq_t kswq_t;

typedef struct {
        int score; // best score
        int te, qe; // target end and query end
        int score2, te2; // second best score and ending position on the target
        int tb, qb; // target start and query start
} kswr_t;


        kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry);
        kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra);

        int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int h0, int *_qle, int *_tle);
        int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *_n_cigar, uint32_t **_cigar);
        kswq_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat);
""")

f =  open('./klib/ksw.c')
src = f.read()
f.close()
_C = ffi.verify(src, extra_compile_args=['-I' + './klib/'])


x = ffi.new('kswq_t *q[2]')

"""
A = ffi.new('int8_t mat[256]')
A[0:256] = [4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 0, 4, 1,  4, 4, 4, 2,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  3, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4, 
        4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4,  4, 4, 4, 4]
"""


A = ffi.new('int8_t mat[25]')

A[0:25] = [1, -1, -1, -1, 0,
           -1, 1, -1, -1, 0,
           -1, -1, 1, -1, 0,
           -1, -1, -1, 1, 0,
           0, 0, 0, 0, 0]

#s = 'A'*512
s = chr(0) * 512
n = len(s)
z = _C.ksw_align2(n, s, n, s, 5, A, 2, 1, 0x10000)
