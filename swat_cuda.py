#!usr/bin/env python
#from __future__ import print_function
from __future__ import absolute_import
import pycuda.driver as drv
import pycuda.tools
import pycuda.autoinit
from pycuda.autoinit import context
import numpy as np
import numpy.linalg as la
from pycuda.compiler import SourceModule
from pycuda import gpuarray
import sys
from numba import njit, prange
from time import time

#prange = xrange


# a simple n bit array
# usage:
# x = np.empty(10, 'uint16')
# # set index 0 to 2
# setitem(x, 0, 2, n=5, base=0b11111)
# # set index 1 to 31
# setitem(x, 1, 31)
#
# # get index 0
# getitem(x, 0)

def setitem(y, i, c, n=5):
    p = x.strides[0] * 8 // n
    base = (1 << n) - 1
    a, b = i // p, i % p
    s = b * n
    y[a] = (y[a] & (~(base << s))) | (c << s)

def getitem(y, i, n=5):
    p = x.strides[0] * 8 // n
    base = (1 << n) - 1
    a, b = i // p, i % p
    return (y[a] >> b * n) & base

initial = time();

DIM = int(eval(sys.argv[1]))
SIM = int(eval(sys.argv[2]))


N = np.int64(DIM)
S = np.int64(SIM)
#S = np.int64(DIM)

#print(N * S / 1e9)
print N

# determin thread_per_block and block_per_grid
tpb =min(2**8, S)
#tpb =min(2**10, S//2)
#bpg = (N + S - 1) // S
#bpg = (bpg + tpb - 1) // tpb
#bpg = S // tpb
#bpg =2**15

SL = S
SL = np.int64(SL)
#SL=S
bpg = int(((N + SL - 1) // SL) **.5 +1)

print 'tpb', tpb, 'bpg', bpg

np.random.seed(42)
x = np.random.randint(1, 255//2, N, dtype='int8')
#x[:] = 0
d_x = gpuarray.to_gpu(x)


D = np.int64(S)
y = np.random.randint(1, 255//2, D, dtype='int8')
#y = x[:D*1000]
y = x[:D]

print 'y shape', y.shape
d_y = gpuarray.to_gpu(y)



#F = np.zeros(N+10**6, dtype='short')
#d_F = gpuarray.to_gpu(F)
#H = F.copy()
#d_H = gpuarray.to_gpu(H)



#y = y.copy()
#for i in xrange(0, y.shape[0], 30):
#    try:
#        y[i:i+30] = x[i+1:i+31]
#    except:
#        pass


swatc0 = """
#include <stdio.h>
#include <stdint.h>

/*
inline short max(short x, short y)
{
  return x ^ ((x ^ y) & -(x < y));
}
*/

// H is the score matrix, and its length should be 2 + ref_length
__device__ short swat(char *S0, long st0, long ed0, char *S1, long st1, long ed1, short *Hi, short *Fi, short Ge, short Go)
{
  long l0=ed0-st0, l1=ed1-st1, i=0, j=0; 
  short Hmax=0, Hij, Hi_1j_1, Hi_1j, Hij_1, Eij, Eij_1, Fij, Fi_1j, w, flag=0, tmp, x, y;
  char c0, c1;


  /*
  Hi[0] = 3;
  Fi[0] = 1;
  printf("test %d %d\\n", Hi[0], Fi[0]);
  */


  for(j=0; j<l1+1; j++)
  {
    Hi[j] = 0;
    Fi[j] = 0;
    //__syncthreads();
  }

  /*
  for(j=0; j<l1+1; j++)
  {
    printf("start_fku, %ld, %ld, Hi %d\\n", i, j, Hi[j]);
  }
  */

  //long flag=0;
  for(i=1; i<l0+1; i++)
  {

    /*
    for(j=0; j<l1+1; j++)
    {
      printf("Hi-1: %ld %ld %d \\n", i-1, j, Hi[j]);
    }
    */


    c0 = S0[i+st0-1];

    // update Fj
    for(j=1; j<l1+1; j++)
    {
      x = Hi[j]-Go; 
      y = Fi[j]-Ge;
      Fi[j] = max(x, y);
    }

    // 1st round of update Hj
    for(j=0; j<l1; j++)
    {

      c1 = S1[j+st1];
      w = (c0 == c1)?2:-2;
      Hi[j] += w;

    }

    for(j=l1; j>0; j--)
    {
      x = Hi[j-1];
      y = Fi[j];
      Hi[j] = max(x, y);
    }
    Hi[0] = 0;

    // 2nd round of update Hj
    Eij_1 = 0;
    for(j=1; j<l1+1; j++)
    {
      //Eij = max(Eij_1 - Ge, Hi[j-1]-Go);

      x = Eij_1 - Ge;
      y = Hi[j-1]-Go;
      Eij = max(x, y);

      //Hi[j] = max(max(Eij, Hi[j]), 0);
      x = Hi[j];
      y = max(Eij, x);
      Hi[j] = max(y, 0);

      Eij_1 = Eij;
      x = Hi[j];
      
      Hmax = max(Hmax, x);

    }

  }

  return Hmax;

}

__global__ void calls(char *qry, long l0, char *refs, long N, long S)
{
  __shared__ char S0[8200];
  __shared__ short Fi[8200], Hi[8200];
  const long tid = threadIdx.x + blockIdx.x * blockDim.x;
  const long start=tid*S, end=tid*S+S;

  if(start<N)
  {
    for(long i=0; i<l0; i++)
    {
      S0[i] = qry[i];
    }
    short Hmax = swat(S0, 0, l0, refs, start, end, Hi, Fi, 11, 2);
    //swat(S0, 0, l0, refs, start, end, Hi, Fi, 11, 1);
    //printf("%d|%d|%d|%d|%d|\\n", l0, start, end, Hmax);
    printf("%ld %ld %d cuda\\n", start, end, Hmax);

  }
}
"""


swatc1 = """
#include <stdio.h>
#include <stdint.h>

/*
inline short max(short x, short y)
{
  return x ^ ((x ^ y) & -(x < y));
}
*/

// H is the score matrix, and its length should be 2 + ref_length
__device__ short swat(char *S0, long st0, long ed0, char *S1, long st1, long ed1, short *Hi, short *Fi, short Ge, short Go)
{
  long l0=ed0-st0, l1=ed1-st1, i=0, j=0; 
  short Hmax=0, Hij, Hi_1j_1, Hi_1j, Hij_1, Eij, Eij_1, Fij, Fi_1j, w, flag=0, tmp, x, y;
  char c0, c1;



  for(j=0; j<l1+1; j++)
  {
    Hi[j+st1] = 0;
    Fi[j+st1] = 0;
  }

  //long flag=0;
  for(i=1; i<l0+1; i++)
  {

    c0 = S0[i+st0-1];

    // update Fj
    for(j=1; j<l1+1; j++)
    {
      x = Hi[j+st1]-Go; 
      y = Fi[j+st1]-Ge;
      Fi[j+st1] = max(x, y);
    }

    // 1st round of update Hj
    for(j=0; j<l1; j++)
    {

      c1 = S1[j+st1];
      w = (c0 == c1)?2:-2;
      Hi[j+st1] += w;

    }

    for(j=l1; j>0; j--)
    {
      x = Hi[j-1+st1];
      y = Fi[j+st1];
      Hi[j+st1] = max(x, y);
    }
    Hi[st1] = 0;

    // 2nd round of update Hj
    Eij_1 = 0;
    for(j=1; j<l1+1; j++)
    {
      //Eij = max(Eij_1 - Ge, Hi[j-1]-Go);

      x = Eij_1 - Ge;
      y = Hi[j-1+st1]-Go;
      Eij = max(x, y);

      //Hi[j] = max(max(Eij, Hi[j]), 0);
      x = Hi[j+st1];
      y = max(Eij, x);
      Hi[j+st1] = max(y, 0);

      Eij_1 = Eij;
      x = Hi[j+st1];
      
      Hmax = max(Hmax, x);

    }

  }

  return Hmax;

}

__global__ void calls(char *qry, long l0, char *refs, long N, long S, short *Fi, short *Hi)
{
  __shared__ char S0[8200];
  //__shared__ short Fi[8200], Hi[8200];
  const long tid = threadIdx.x + blockIdx.x * blockDim.x;
  const long start=tid*S, end=tid*S+S;

  if(start<N)
  {
    for(long k=0;k<l0;k+=S)
    {
      //for(long i=0; i<l0; i++)
      for(long i=0; i<S; i++)
      {
        S0[i] = qry[i+k];
      }
      printf("idx %ld %ld %ld %ld\\n", k, l0, S, tid);
      __syncthreads();

      short Hmax = swat(S0, 0, l0, refs, start, end, Hi, Fi, 11, 2);
      //swat(S0, 0, l0, refs, start, end, Hi, Fi, 11, 1);
      //printf("%d|%d|%d|%d|%d|\\n", l0, start, end, Hmax);
      //printf("%ld %ld %d cuda\\n", start, end, Hmax);
      __syncthreads();

    }
  }
}
"""


swatc="""
#include "cuda_runtime.h"
#include "device_launch_parameters.h" 
#include <stdio.h>


__global__ void swat(char *qry, long N, char *S1, long M, long Go, long Ge)
{
    __shared__ char S0[6000];
    __shared__ short F[6000], H0[6000], H1[6000];
    short x, y, H1tid, H1i, H0i_old;
    long I, j;
    char c0, c1;
    const long tid = threadIdx.x + blockIdx.x * blockDim.x;
    if(tid < N)
    {
        S0[tid] = qry[tid];
        F[tid] = 0;
        H0[tid] = 0;
        H1[tid] = 0;
    }
    __syncthreads();


    for(long J=0; J<M; J++)
    {
        //break;
        if(tid>0)
        {
            //printf("H' finish %d %d\\n", tid, J);
            //printf("H1 finish %ld %ld N:%ld M:%ld\\n", tid, J, N, M);
            // update F
            x = H0[tid]-Go;
            y = F[tid]-Ge;
            F[tid] = max(x, y);
            //__syncthreads();
            //printf("F finish\\n");

            // update H'
            x = H0[tid-1] + S0[tid-1]==S1[J]?2:-2;
            y = F[tid];
            x = max(x, 0);
            H1[tid] = max(x, y);
            //__syncthreads();

            // update H
            H0[tid] = (H1[tid]-Go>H1[tid+1])?tid:-1;
            //__syncthreads();

            if(H0[tid]==tid)
            {
                H1tid = H1[tid];
                for(long i=tid+1; i<N; i++)
                {
                    H1i = H1tid - Go - (i-tid-1) * Ge;
                    H0i_old = H0[i];
                    if(H1i>H1[i])
                    {
                        H0[i] = tid;

                    }else{
                        break;
                    }
                    if(H0i_old != -1)
                    {
                        break;
                    }
                }
            }
            //__syncthreads();

            /*
            I = tid;
            j = H0[tid];
            while(j!=-1)
            {
                printf("j is %d\\n", j);
                if(I == j)
                {
                    H1[tid] = I < tid and  H1[I] - Go - (tid-I-1) * Ge or H1[I];
                    break;
                }else{
                    I = j;
                    j = H0[I];
                }
            }
            __syncthreads();
            H0[tid] =  H1[tid];
            */
        }
    }
}
"""



#mod = SourceModule(swatc0)
swatc = open('swat_cuda.cu', 'r').read()
#swatc = swatc % {'CONST': S}
#print(swatc)
#mod = SourceModule(swatc, options=['-rdc=true','-lcudart','-lcudadevrt','--machine=64'], arch='sm_52')
mod = SourceModule(swatc)
#swat_cuda = mod.get_function("calls")
#swat_cuda = mod.get_function("swat_print")
#swat_cuda = mod.get_function("swat_strip")
swat_cuda = mod.get_function("swat_strip_warp")


#swat_cuda(d_y, D, d_x, N, S, block=(tpb, 1, 1), grid=(bpg, 1))


#raise SystemExit()

GO = np.int16(11)
GE = np.int16(1)

#swat_cuda(d_y, D, d_x, N, GO, GE, block=(tpb, 1, 1), grid=(bpg, 1))
#raise SystemExit()
print 'cpu time', time() - initial
initial = time()


start = drv.Event()
end = drv.Event()
start.record()
sq = 512
#for itr in xrange((sq+D-1)//D):
if 1:
    #swat_cuda(d_y, D, d_x, N, S, block=(tpb, 1, 1), grid=(bpg, 1))
    #swat_cuda(d_y, D*1000, d_x, N, S, d_H, d_F, block=(tpb, 1, 1), grid=(bpg, 1))
    #tpb=np.int64(min(blk*2, 2**10))
    #bpg = np.int64(N / tpb)
    print 'D', D, 'N', N,  'SL', SL, 'tpb', tpb, 'bpg', bpg, bpg*bpg
    #//tpb = np.int64(tpb // 4)
    #swat_cuda(d_y, D, d_x, N, GO, GE, block=(tpb, 1, 1), grid=(bpg, bpg))
    #swat_cuda(d_y, tpb, d_x, N, GO, GE, block=(tpb, 1, 1), grid=(bpg, bpg))
    swat_cuda(d_y, D, d_x, np.int64(N), GO, GE, SL, block=(tpb, 1, 1), grid=(bpg, bpg))
    #swat_cuda(d_y, D, d_x, np.int64(N), GO, GE, SL, block=(tpb, 1, 1), grid=(bpg, 1))

    #context.synchronize()

end.record()
end.synchronize()
secs = start.time_till(end)*1e-3
print("SourceModule time and first three results:", secs)

print 'cpu', time() - initial

@njit
def swat(x, y, S, Go=11, Ge=2):
    n = x.size
    m = y.size
    H = np.zeros((n+1, S+1), dtype=np.int32)
    E = H.copy()
    F = H.copy()

    output = np.empty((m//S+1, 2), dtype=np.int32)
    for idx in xrange(0, m, S):
        H[:, 0] = 0
        H[0, :] = 0
        E[:, 0] = 0
        E[0, :] = 0
        F[:, 0] = 0
        F[0, :] = 0
        Hmax = 0
        for J in xrange(idx, idx+S):
            for I in xrange(n):
                #i, j = I + 1, J + 1
                ci, cj = x[I], y[J]
                i = I+1
                j = J + 1 - idx
                #print(i, j, E.shape, H.shape)
                E[i, j] = max(E[i, j-1]-Ge, H[i, j-1]-Go)
                F[i, j] = max(F[i-1, j]-Ge, H[i-1, j]-Go)
                H[i, j] = max(E[i, j], F[i, j], H[i-1, j-1] + (ci==cj and 2 or -2), 0)
                Hmax = max(H[i, j], Hmax) 
        output[idx//S, 0] = Hmax
        output[idx//S, 1] = idx

    return output




@njit
def swat_m(x, y, S, Go=11, Ge=2):
    n = x.size
    m = y.size
    H = np.zeros(S+1, dtype=np.int32)
    F = H.copy()

    output = np.empty((m//S+1, 2), dtype=np.int32)
    for idx in xrange(0, m, S):
        H[:] = 0
        F[:] = 0
        Hmax = 0
        for J in xrange(idx, idx+S):
            for I in xrange(n):
                F[I+1] = max(H[I+1]-Go, F[I+1]-Ge)

            for I in xrange(n):
                #i, j = I + 1, J + 1
                ci, cj = x[I], y[J]
                H[I] += (ci==cj and 2 or -2)

            for I in xrange(n, 0, -1):
                H[I] = max(H[I-1], F[I])
            H[0] = 0

            Eij_1 = 0
            for I in xrange(n):
                Eij = max(Eij_1-Ge, H[I]-Go)
                H[I+1] = max(Eij, H[I+1], 0)
                Hmax = max(H[I+1], Hmax) 

        output[idx//S, 0] = Hmax
        output[idx//S, 1] = idx

    return output


@njit(parallel=True)
def swat_vec(x, y, S, Go=11, Ge=2):
    qst = 0
    qed = x.size
    m = y.size
    H0 = np.zeros(S+1, dtype=np.int32)
    H1 = H0.copy()
    F = H0.copy()

    output = np.empty((m//S+1, 2), dtype=np.int32)
    for sst in xrange(0, m, S):
        sed = sst + S
        H0[:] = 0
        F[:] = 0
        Hmax = 0
        Hmax_I = 0
        for J in xrange(sst, sed):
            #print('#############################################################################')
            #print('H0', H0)
            #print('0F', F)
            #print('aFO', F - Ge)
            #print('bFO', H0 - Go)

            # update F, can be parallelized
            #for i in prange(qst+1, qed+1):
            #    F[i] = max(H0[i]-Go, F[i]-Ge)

            for I in prange(qst+1, qed+1):
                F[I] = max(H0[I]-Go, F[I]-Ge)

            #print('1F', F)

            # 1st round of updating H1, can be parallelized
            H1[0] = 0
            #print('0H0', H0)
            for I in prange(qst, qed):
                ci, cj = x[I], y[J]

                #if I==J:
                #    print ('iH0', H0, 'I, J', I, J)

                H1[I+1] = max(H0[I] + (ci==cj and 2 or -2), F[I+1], 0)
                #if Hmax < H1[I+1]:
                #    Hmax = H1[I+1]
                #    Hmax_I = I+1
                Hmax = max(Hmax, H1[I+1])

            #print('1H1', H1)

            #print('sst', J//S, 'J', J, 'H1max', Hmax, 'H1max_idx', Hmax_I, 'H0max', H0max, 'H0max_idx', H0max_I)
            # determine gopen position, can be parallelized
            #print('1H0', H0)
            for I in prange(qst, qed):
                if H1[I] - Go > H1[I+1]:
                    H0[I] = I
                else:
                    H0[I] = -1

            H0[qed] = -1
            # need sync
            #print('2H0', H0)


            # find step to update H1
            for I in prange(qst, qed):
                if H0[I] == I:
                    H1I = H1[I]
                    for i in xrange(I+1, qed+1):
                        H1i = H1I - Go - (i-I-1) * Ge
                        H0i_old = H0[i]
                        if H1i > H1[i]:
                            H0[i] = I
                        else:
                            break

                        if H0i_old != -1:
                            break

            #print('2H1', H1)
            #print('3H0', H0)
            # can be parallelized
            for i in prange(qst+1, qed+1):
                I = i
                j = H0[i]
                while j != -1:
                    if i == j:
                        #H1[I] = H1[i]
                        H1[I] = i < I and  H1[i] - Go - (I-i-1) * Ge or H1[i]
                        break
                    else:
                        i = j
                        j = H0[i]

            # swatch H0, H1
            #H0[:] = H1
            H0, H1 = H1, H0
            #print('4H0', H0)

        #Hmax = H0.max()
        output[sst//S, 0] = Hmax
        output[sst//S, 1] = sst

    return output



# fast version of swat_vec
@njit(parallel=True)
def swat_vec_fast(x, y, S, Go=11, Ge=2):
    qst = 0
    qed = x.size
    m = y.size
    H0 = np.zeros(S+1, dtype=np.int32)
    H1 = H0.copy()
    F = H0.copy()

    output = np.empty((m//S+1, 2), dtype=np.int32)
    for sst in xrange(0, m, S):
        sed = sst + S
        H0[:] = 0
        F[:] = 0
        Hmax = 0
        Hmax_I = 0
        for J in xrange(sst, sed):
            for I in prange(qst+1, qed+1):
                F[I] = max(H0[I]-Go, F[I]-Ge)

            # 1st round of updating H1, can be parallelized
            H1[0] = 0
            for I in prange(qst, qed):
                ci, cj = x[I], y[J]

                H1[I+1] = max(H0[I] + (ci==cj and 2 or -2), F[I+1], 0)
                Hmax = max(Hmax, H1[I+1])

            # determine gopen position, can be parallelized
            for I in prange(qst, qed):
                if H1[I] - Go > H1[I+1]:
                    H0[I] = I
                else:
                    H0[I] = -1

            H0[qed] = -1
            # need sync


            # find step to update H1
            for I in prange(qst, qed):
                if H0[I] == I:
                    H1I = H1[I]
                    for i in xrange(I+1, qed+1):
                        H1i = H1I - Go - (i-I-1) * Ge
                        if H1i <= H1[i]:
                            break
                        elif H0[i] == i:
                            H1[i] = H1i
                            H0[i] = I
                            break
                        else:
                            H0[i] = I

            # can be parallelized
            for i in prange(qst+1, qed+1):
                j = H0[i]
                k = H0[j]
                go = j == k and Go or Ge 
                H1i = j < i and  H1[j] - go - (i-j-1) * Ge or H1[i]


            '''
            # can be parallelized
            for i in prange(qst+1, qed+1):
                H1[i] = H1[H0[i]]
            '''

            H0, H1 = H1, H0

        output[sst//S, 0] = Hmax
        output[sst//S, 1] = sst

    return output




raise SystemExit()
st = time()
#output = swat(y, x, S)
output1 = swat_m(y, x, S)
print 'single thread', time() - st

st = time()
output2 = swat_vec_fast(y, x, S)
print 'multiple thread', time() - st

for i in xrange(output2.shape[0]):
    #print(i, output[i, :], x.size, y.size)
    s, idx = output1[i, :]
    s2, idx2 = output2[i, :]

    print(idx, idx+S, s, 'cpu', s2, 'equal', s==s2)



