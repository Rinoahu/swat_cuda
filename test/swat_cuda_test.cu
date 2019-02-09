#include "cuda_runtime.h"
#include "device_launch_parameters.h" 
#include <stdio.h>

// set the "thread per block" to 128
__global__ void swat(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    __shared__ char S0[4097];
    __shared__ int32_t F[4097], H0[4097], H1[4097];
    int16_t x, y;
    //char c0, c1;
    const long tid=threadIdx.x, bid=blockIdx.x, bdm=blockDim.x;
    //long tid=threadIdx.x, bid=blockIdx.x, bdm=blockDim.x;
    //const long idx=tid + bid * bdm;
    long th=2048, step=4096;

    //const long tid = threadIdx.x + blockIdx.x * blockDim.x;

    //printf("hello tid %ld bid %ld M %ld \n", tid, bid, M);
    //__syncthreads();

    // load all data into shr array;
    for(long i=tid; i<N;i+=th)
    {
        S0[tid] = qry[tid];
        F[tid] = 0;
        H0[tid] = 0;
        H1[tid] = 0;
    }
    //__syncthreads();

}

