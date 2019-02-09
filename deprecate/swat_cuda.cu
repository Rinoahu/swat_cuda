#include "cuda_runtime.h"
#include "device_launch_parameters.h" 
#include <stdio.h>



// set the "thread per block" to 128
__global__ void swat_print(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    __shared__ char S0[4097], S1[4097];
    __shared__ short F[4097], H0[4097], H1[4097];
    //__shared__ int H[4097];
    short x, y;
    //char c0, c1;
    //const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, dmx=blockDim.x, dmy=blockDim.y, bid=by+bx*dmy;
    const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, bdmx=blockDim.x, bid=bx+by*gridDim.x;
    const long th=1024, step=1024*4, bth=256, B=32, W=(N+31)/B;

    /*
    if(tid==0)
    {
        printf("bx %ld bid %ld M %ld dmx %ld\n", bx, bid, M/step, dmx);
    }
    __syncthreads();
    */

    if(bid<M/step)
    {
        // load qry and ref into S0, S1;
        #pragma unroll(4)
        for(long i=tid; i<N; i+=th)
        {
            S0[i] = qry[i];
            S1[i] = refs[bid+i];
            F[i] = 0;
            H0[i] = 0;
            H1[i] = 0;
        }
        __syncthreads();

        //return;
        for(long j=0; j<step; j++)
        {
            char S1j = S1[j];
            // update F
            H0[N] = 0;
            #pragma unroll(4)
            for(long i=tid, i1; i<N; i+=th)
            //for(long i0=0, i1, i, end=(N+th-1)/th; i0<end; i0++)
            {
                //i = tid + i0*th;
                // update F
                i1 = i+1;
                x = H0[i1]-Go;
                y = F[i1]-Ge;
                F[i1] = max(x, y);
                //x = x ^ ((x ^ y) & -(x < y));
                F[i1] = x;

                // update H'
                x = H0[i] + 4*(S0[i]==S1j) - 2;
                y = F[i1];
                x = max(x, y);
                x = max(x, 0);
                H1[i1] = x;
                __syncthreads();

                //continue;
                // use prefix scan to update H1
                // find update start site
                x = (H1[i]-Go>H1[i1]);
                H0[i] = x*i-(!x);

                if(H0[i]==i)
                {
                    short H1i = H1[i];
                    #pragma unroll(2)
                    for(i1=i+1; i1<N+1; i1++)
                    {
                        short H1i1 = H1i - Go - (i1-i-1) * Ge;
                        short H0i1_old = H0[i1];

                        if(H1i1>H1[i1])
                        {
                            H0[i1] = i;
                        }else{
                            break;
                        }

                        if(H0i1_old != -1)
                        {
                            break;
                        }
                    }
                }

                // update H 
                long i0 = i;
                i1 = H0[i0];
                while(i1!=-1)
                {
                    //printf("start %ld %ld %ld %ld %ld %ld %ld %ld I J\n", start, N, idx, tid, bid, bdm, I, j);
                    if(i0==i1)
                    {
                        H1[i] = (i1==i)?H1[i]: H1[i]-Go-(i1-i-1)*Ge;
                        break;
                    }else{
                        i0 = i1;
                        i1 = H0[i0];
                    }
                }
                //continue;
            }
        }
    }
}


