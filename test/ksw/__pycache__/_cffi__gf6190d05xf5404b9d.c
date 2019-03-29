
#include <stdio.h>
#include <stddef.h>
#include <stdarg.h>
#include <errno.h>
#include <sys/types.h>   /* XXX for ssize_t on some platforms */

/* this block of #ifs should be kept exactly identical between
   c/_cffi_backend.c, cffi/vengine_cpy.py, cffi/vengine_gen.py
   and cffi/_cffi_include.h */
#if defined(_MSC_VER)
# include <malloc.h>   /* for alloca() */
# if _MSC_VER < 1600   /* MSVC < 2010 */
   typedef __int8 int8_t;
   typedef __int16 int16_t;
   typedef __int32 int32_t;
   typedef __int64 int64_t;
   typedef unsigned __int8 uint8_t;
   typedef unsigned __int16 uint16_t;
   typedef unsigned __int32 uint32_t;
   typedef unsigned __int64 uint64_t;
   typedef __int8 int_least8_t;
   typedef __int16 int_least16_t;
   typedef __int32 int_least32_t;
   typedef __int64 int_least64_t;
   typedef unsigned __int8 uint_least8_t;
   typedef unsigned __int16 uint_least16_t;
   typedef unsigned __int32 uint_least32_t;
   typedef unsigned __int64 uint_least64_t;
   typedef __int8 int_fast8_t;
   typedef __int16 int_fast16_t;
   typedef __int32 int_fast32_t;
   typedef __int64 int_fast64_t;
   typedef unsigned __int8 uint_fast8_t;
   typedef unsigned __int16 uint_fast16_t;
   typedef unsigned __int32 uint_fast32_t;
   typedef unsigned __int64 uint_fast64_t;
   typedef __int64 intmax_t;
   typedef unsigned __int64 uintmax_t;
# else
#  include <stdint.h>
# endif
# if _MSC_VER < 1800   /* MSVC < 2013 */
#  ifndef __cplusplus
    typedef unsigned char _Bool;
#  endif
# endif
#else
# include <stdint.h>
# if (defined (__SVR4) && defined (__sun)) || defined(_AIX) || defined(__hpux)
#  include <alloca.h>
# endif
#endif

 
#include "ksw.h"
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
        //kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, int gapo, int gape, int xtra);
        //kswr_t ksw_align3(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8 t *mat, int gapo, int gape, int xtra);

        int ksw_extend(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int h0, int *_qle, int *_tle);
        int ksw_global(int qlen, const uint8_t *query, int tlen, const uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int w, int *_n_cigar, uint32_t **_cigar);
        kswq_t *ksw_qinit(int size, int qlen, const uint8_t *query, int m, const int8_t *mat);


/*
kswr_t ksw_align2(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, int gapo, int gape, int xtra)
{
        kswq_t *q[2] = {0, 0};
        kswr_t r;
    r = ksw_align(qlen, query, tlen, target, m, AA, gapo, gape, xtra, q);
    //return r;
    //return ksw_align(qlen, query, tlen, target, m, mat, gapo, gape, xtra, q);
        free(q[0]); free(q[1]);
    return r;
};
*/

kswr_t ksw_align3(int qlen, uint8_t *query, int tlen, uint8_t *target, int m, const int8_t *mat, int gapo, int gape, int xtra)
{
        kswq_t *q[2] = {0, 0};
        kswr_t r;
    r = ksw_align(qlen, query, tlen, target, m, mat, gapo, gape, xtra, q);
    //return r;
    //return ksw_align(qlen, query, tlen, target, m, mat, gapo, gape, xtra, q);
        free(q[0]); free(q[1]);
    return r;
};

static void _cffi_check__kswr_t(kswr_t *p)
{
  /* only to generate compile-time warnings or errors */
  (void)p;
  (void)((p->score) << 1);
  (void)((p->te) << 1);
  (void)((p->qe) << 1);
  (void)((p->score2) << 1);
  (void)((p->te2) << 1);
  (void)((p->tb) << 1);
  (void)((p->qb) << 1);
}
intptr_t _cffi_layout__kswr_t(intptr_t i)
{
  struct _cffi_aligncheck { char x; kswr_t y; };
  static intptr_t nums[] = {
    sizeof(kswr_t),
    offsetof(struct _cffi_aligncheck, y),
    offsetof(kswr_t, score),
    sizeof(((kswr_t *)0)->score),
    offsetof(kswr_t, te),
    sizeof(((kswr_t *)0)->te),
    offsetof(kswr_t, qe),
    sizeof(((kswr_t *)0)->qe),
    offsetof(kswr_t, score2),
    sizeof(((kswr_t *)0)->score2),
    offsetof(kswr_t, te2),
    sizeof(((kswr_t *)0)->te2),
    offsetof(kswr_t, tb),
    sizeof(((kswr_t *)0)->tb),
    offsetof(kswr_t, qb),
    sizeof(((kswr_t *)0)->qb),
    -1
  };
  return nums[i];
  /* the next line is not executed, but compiled */
  _cffi_check__kswr_t(0);
}

void _cffi_f_ksw_align(kswr_t *r, int x0, uint8_t * x1, int x2, uint8_t * x3, int x4, int8_t const * x5, int x6, int x7, int x8, kswq_t * * x9)
{
  *r = ksw_align(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9);
}

int _cffi_f_ksw_extend(int x0, uint8_t const * x1, int x2, uint8_t const * x3, int x4, int8_t const * x5, int x6, int x7, int x8, int x9, int * x10, int * x11)
{
  return ksw_extend(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10, x11);
}

int _cffi_f_ksw_global(int x0, uint8_t const * x1, int x2, uint8_t const * x3, int x4, int8_t const * x5, int x6, int x7, int x8, int * x9, uint32_t * * x10)
{
  return ksw_global(x0, x1, x2, x3, x4, x5, x6, x7, x8, x9, x10);
}

kswq_t * _cffi_f_ksw_qinit(int x0, int x1, uint8_t const * x2, int x3, int8_t const * x4)
{
  return ksw_qinit(x0, x1, x2, x3, x4);
}

int _cffi_const_KSW_XBYTE(char *out_error)
{
  if ((KSW_XBYTE) <= 0 || (unsigned long)(KSW_XBYTE) != 65536UL) {
    char buf[64];
    if ((KSW_XBYTE) <= 0)
        sprintf(buf, "%ld", (long)(KSW_XBYTE));
    else
        sprintf(buf, "%lu", (unsigned long)(KSW_XBYTE));
    sprintf(out_error, "%s has the real value %s, not %s",
            "KSW_XBYTE", buf, "65536");
    return -1;
  }
  return 0;
}

int _cffi_const_KSW_XSTART(char *out_error)
{
  if ((KSW_XSTART) <= 0 || (unsigned long)(KSW_XSTART) != 524288UL) {
    char buf[64];
    if ((KSW_XSTART) <= 0)
        sprintf(buf, "%ld", (long)(KSW_XSTART));
    else
        sprintf(buf, "%lu", (unsigned long)(KSW_XSTART));
    sprintf(out_error, "%s has the real value %s, not %s",
            "KSW_XSTART", buf, "524288");
    return -1;
  }
  return 0;
}

int _cffi_const_KSW_XSTOP(char *out_error)
{
  if ((KSW_XSTOP) <= 0 || (unsigned long)(KSW_XSTOP) != 131072UL) {
    char buf[64];
    if ((KSW_XSTOP) <= 0)
        sprintf(buf, "%ld", (long)(KSW_XSTOP));
    else
        sprintf(buf, "%lu", (unsigned long)(KSW_XSTOP));
    sprintf(out_error, "%s has the real value %s, not %s",
            "KSW_XSTOP", buf, "131072");
    return -1;
  }
  return 0;
}

int _cffi_const_KSW_XSUBO(char *out_error)
{
  if ((KSW_XSUBO) <= 0 || (unsigned long)(KSW_XSUBO) != 262144UL) {
    char buf[64];
    if ((KSW_XSUBO) <= 0)
        sprintf(buf, "%ld", (long)(KSW_XSUBO));
    else
        sprintf(buf, "%lu", (unsigned long)(KSW_XSUBO));
    sprintf(out_error, "%s has the real value %s, not %s",
            "KSW_XSUBO", buf, "262144");
    return -1;
  }
  return 0;
}

