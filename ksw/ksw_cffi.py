#!usr/bin/env python
from cffi import FFI
import os


ffi = FFI()


ffi.cdef("""
struct _kswq_t;
typedef struct _kswq_t kswq_t;

typedef struct {
        int score; // best score
        int te, qe; // target end and query end
        int score2, te2; // second best score and ending position on the target
        int tb, qb; // target start and query start
} kswr_t;


        kswr_t ksw_align(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra, kswq_t **qry);
        int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int h0, int *_qle, int *_tle);
        int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *_n_cigar, uint32_t **_cigar);
        kswq_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat);

""")

f =  open('./klib/ksw.c')
src = f.read()
f.close()
_C = ffi.verify(src, extra_compile_args=['-I' + './klib/'])
