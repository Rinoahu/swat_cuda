#!usr/bin/env python
from __future__ import print_function
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

# the best performance
src1 = """
#include <stdio.h>
__global__ void multiply_them(float *dest, float *a, float *b, long N, long D, long S, long M)
{
  __shared__ float sdata[%(S0)s];
  //const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int start = tid * S, end = start + S;
  int i, j, k;
  float ai, bi;

  for(k=0;k<M;k++)
  {
    //printf("|%(S1)s|", k);
  
    for(j=0;j<S;j++)
    {
      sdata[j] = b[j];
      __syncthreads();
    }

    //if(i < N)
    for(i=start; i < end; i++)
    {
      ai = a[i];
      for(j=0; j<S; j++)
      {
        //bi = ai < dest[j]?ai:sdata[j];
        dest[i] += ai * sdata[j];
      }
    }
  }
  __syncthreads();

  /*if(tid %(S2)s 4096 == 0)
  {
    printf("|%(S1)s|", (tid));
  }
  */
}
"""



src2 = """
#include <stdio.h>
__global__ void multiply_them(float *dest, float *a, float *b, long N, long D, long S, long M)
{
  __shared__ float sdata[%(S0)s];
  //float sdata[%(S0)s];
  //const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int start = tid * S, end = start + S;
  int i, j, k;
  float ai, bi, tmp;

  for(k=0;k<M;k++)
  {
    //printf("|%(S1)s|", k);
  
    for(j=0;j<S;j++)
    {
      sdata[j] = b[j];
      __syncthreads();
    }

    //if(i < N)
    for(i=start; i < end; i++)
    {
      ai = a[i];
      for(j=0; j<S; j++)
      {
        bi = ai < sdata[j]?ai:sdata[j];
        tmp = dest[i] + ai * sdata[j];
        dest[i] = tmp;
      }
    }
  }
  __syncthreads();

  /*if(tid %(S2)s 4096 == 0)
  {
    printf("|%(S1)s|", (tid));
  }
  */
}
"""


src3 = """
#include <stdio.h>
__global__ void multiply_them(float *dest, float *a, float *b, long N, long D, long S, long M)
{
  //__shared__ float sdata[%(S0)s];
  float sdata[%(S0)s];
  //const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int start = tid * S, end = start + S;
  int i, j, k;
  float ai, bi, tmp;

  for(k=0;k<M;k++)
  {
    //printf("|%(S1)s|", k);
  
    for(j=0;j<S;j++)
    {
      sdata[j] = b[j];
      __syncthreads();
    }

    //if(i < N)
    for(i=0; i < N; i++)
    {
      ai = a[i];
      sdata[tid] = ai < sdata[tid]?ai:sdata[tid];
      __syncthreads();
      if(tid > 0)
      {
        sdata[tid] = sdata[tid-1] < sdata[tid]?sdata[tid-1]:sdata[tid];
      }
      __syncthreads();

    }
  }
  __syncthreads();

  /*if(tid %(S2)s 4096 == 0)
  {
    printf("|%(S1)s|", (tid));
  }
  */
}
"""

src4 = """
#include <stdio.h>
__global__ void multiply_them(float *dest, float *a, float *b, long N, long D, long S, long M)
{
  __shared__ float sdata[%(S0)s];
  //float sdata[%(S0)s];
  //const int i = threadIdx.x + blockDim.x * blockIdx.x;
  const int tid = threadIdx.x + blockIdx.x * blockDim.x;
  const int start = tid * S, end = start + S;
  int i, j, k;
  float ai, bi, tmp, ei, fi, bj;

  for(k=0;k<M;k++)
  {
    //printf("|%(S1)s|", k);
  
    for(j=0;j<S;j++)
    {
      sdata[j] = b[j];
      __syncthreads();
    }

    //if(i < N)
    for(j=0; j<S; j++)
    {
      bj = b[j];
      ei = 0;
      fi = 0;
      for(i=start; i < end; i++)
      {
        ai = a[i];
        tmp = (i!=0)?(dest[i-1]+((ai==bj)?2:-2)):0;
        ai = dest[i];
        tmp = max(tmp, ai);
        ei -=2;
        tmp = max(tmp, ei);
        dest[i] = tmp;
      }
    }
  }
  __syncthreads();

  /*if(tid %(S2)s 4096 == 0)
  {
    printf("|%(S1)s|", (tid));
  }
  */
}
"""



import sys
N = int(eval(sys.argv[1]))
N = np.int64(N)
S = np.int64(2**13)

print(N, S)


# determin thread_per_block and block_per_grid
tpb = 32
bpg = (N + S - 1) // S
bpg = (bpg + tpb - 1) // tpb
#bpg = 8

x = np.random.randn(N).astype(np.float32)
d_x = gpuarray.to_gpu(x)

z = np.zeros_like(x)
d_z = gpuarray.to_gpu(z)


D = np.int64(S)
y = np.random.randn(D).astype(np.float32)
d_y = gpuarray.to_gpu(y)


mod = SourceModule(src4%{'S0': D, 'S1': "%d", 'S2': "%"})
#mod = SourceModule(src3%{'S0': D, 'S1': "%d", 'S2': "%"})
multiply_them = mod.get_function("multiply_them")


M = np.int64(1)

d_z.fill(0)
start = drv.Event()
end = drv.Event()
start.record()
for itr in xrange(1):
    multiply_them(d_z, d_x, d_y, N, D, S, M, block=(tpb, 1, 1), grid=(bpg, 1))
    context.synchronize()
end.record()
end.synchronize()
secs = start.time_till(end)*1e-3

print("SourceModule time and first three results:", secs)




