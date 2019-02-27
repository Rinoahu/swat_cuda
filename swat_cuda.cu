#include "cuda_runtime.h"
#include "device_launch_parameters.h" 
#include <stdio.h>

// one sequenc per thread
__global__ void swat6(char *S0, long N, char *S1, long M, short Go, short Ge)
{
    __shared__ char s0[8097], s1[4100];
    __shared__ short F[8097], H[8097], E, x, y, Fi, Hi, i, I, i1;
    char c1;
    const long tid=threadIdx.x, B=32, W=(N+31)/B;
    for(long i0=0; i0<N; i0++)
    {
        i = i0%B + i0 / B * 33;
        i1 = (i0+1)%B + (i0+1) / B * 33;
        //i = i0;
        //i1 = i0+1;

        s0[i] = S0[i0];
        F[i1] = 0;
        H[i1] = 0;
    }
    F[0] = 0;
    H[0] = 0;
    for(long j=0; j<M; j++)
    {
        c1=S1[j];
        E = 0;
        for(long i0=0; i0<N; i0++)
        {

            //i = i0;
            //i1 = i0+1;
            // update F
            x = F[i1] - Ge;
            y = H[i1] - Go;
            Fi = max(x, y);
            F[i1] = Fi;

            // update E
            x = E - Ge;
            y = H[i] - Go;
            E = max(x, y);

            // update H 
            Hi = H[i] + (c1 == s0[i]) * 4 - 2;
            Hi = max(Hi, E);
            Hi = max(Hi, Fi);
            H[i1] = max(Hi, 0); 
        }

        /*
        E=0;
        for(long i=0; i<N; i++)
        {
            // update H
            x = E - Ge;
            y = H[i] - Go;
            E = max(x, y);
            x = H[i+1];
            H[i+1] = max(x, E);
        }
        */
    }
}



// set the "thread per block" to 128
__global__ void swat_print(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    //int N1=N+1;
    const int D=4097;

    char S0[D], S1[D], S3[500000];
    //__shared__ char S0[D], S1[D];
    //short F[D], H0[D], H1[D];
    __shared__ short F[D], H0[D]; 
    __shared__ int H1[D];

    //extern __shared__ char S0[], S1[];
    //extern __shared__ short F[], H0[], H1[];


    //__shared__ int H[4097];
    short x, y, z;
    //const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, dmx=blockDim.x, dmy=blockDim.y, bid=by+bx*dmy;
    const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, bdmx=blockDim.x, bid=bx+by*gridDim.x;
    const long th=bdmx, step=N, bth=256;

    if(tid>N)
    {
        return;
    }

    if(bid<M/step)
    {
        // load qry and ref into S0, S1;
        #pragma unroll 4
        for(long i=tid; i<N; i+=th)
        {
            S0[i] = qry[i];
            //S1[i] = refs[bid+i];
            F[i] = 0;
            H0[i] = 0;
            H1[i] = 0;
        }

        #pragma unroll 4 
        for(long i=tid, start=bid*step; i<step; i+=th)
        {
            S1[i]=refs[i+start];
        }

        __syncthreads();

        #pragma unroll 4 
        for(long j=0; j<step; j++)
        {

            #pragma unroll 4
            for(long tmp=20; tmp>=0; tmp--)
            {
                S1[j] = refs[j+tmp];
            }

            char S1j = S1[j];

            #pragma unroll 
            for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            //long i=tid, i1=tid+1;
            {
                short z = H0[i];
                __syncthreads();

                // update F
                x = H0[i1]-Go;
                y = F[i1]-Ge;
                F[i1] = max(x, y);
                F[i1] = x;

                // update H'
                //x = S0[i]!=S1j;
                //x = H0[i] + 4-x*2;
                //x = H0[i] + 2;
                //x = H0[i] + 4*(S0[i]!=S1j) - 2;

                //x = H0[i] + (S0[i]!=S1j)?-2:2;
                x = z + (S0[i]!=S1j)?-2:2;
                y = F[i1];
                x = max(x, y);
                x = max(x, 0);
                H1[i1] = x;
                H0[i] = -1;
            }

            /*
            #pragma unroll
            for(short tmp=1; tmp<N; tmp*=2)
            if(tid>=tmp)
            {
                H0[tid] = max(H0[tid], H1[tid-tmp] - tmp);
                __syncthreads();

            }
            continue;
            __syncthreads();
            */

            /*
            short H=H1[tid];
            #pragma unroll
            for(short tmp=1; tmp<32; tmp*=2)
            {
                short lane = tid % 32;
                H = lane>tmp?__shfl(H, lane-tmp, 32):H;
            }
            continue;
            */

            for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            {
                // find affect range of each elem
                x = H1[i] - Go;
                //#pragma unroll 4
                for(long k=i1; k<N; k++)
                {
                    if(x<=H1[k])
                    {
                        break;
                    }
                    else{
                        atomicMax(&(H1[k]), x);
                        //H1[k] = max(H1[k], x);
                        x -= Ge;
                    }
                }
            }
            continue;



            for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            {
                // find affect range of each elem
                x = H1[i] - Go;
                //#pragma unroll 4
                for(long k=i1; k<N; k++)
                {
                    if(x<=H1[k])
                    {
                        break;
                    }
                    else{
                        H0[k] = 1;
                        x -= Ge;
                    }
                }
            }
            __syncthreads();
            //continue;

            for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            {
                if(H0[i]==-1&&H0[i1]==1)
                {
                    x = H1[i1] - Go;
                    #pragma unroll 4
                    for(long k=i1; k<N; k++)
                    {
                        if(H0[k]!=1)
                        {
                            break;
                        }else{
                            H1[k] = x;
                            x -= Ge;
                        }

                    }
                }

            //continue;
            H0[i1] = H1[i1];
            if(i1==N-2)
            {
                S1[j] = H0[i1];
            }
            }
        }
    }
}






// set the "thread per block" to 128
__global__ void swat_print0(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    //int N1=N+1;
    const int D=257, D1=256;

    char S0[D], S1[D];
    //__shared__ char S0[D], S1[D];
    //short F[D], H0[D], H1[D];
    __shared__ short F[D], H0[D], H1[D];

    //extern __shared__ char S0[], S1[];
    //extern __shared__ short F[], H0[], H1[];


    //__shared__ int H[4097];
    short x, y;
    //const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, dmx=blockDim.x, dmy=blockDim.y, bid=by+bx*dmy;
    const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, bdmx=blockDim.x, bid=bx+by*gridDim.x;
    const long th=bdmx, step=N, bth=256;

    if(tid>N)
    {
        return;
    }

    if(bid<M/step)
    {
        // load qry and ref into S0, S1;
        #pragma unroll 4
        for(long i=tid; i<N; i+=th)
        {
            S0[i] = qry[i];
            //S1[i] = refs[bid+i];
            F[i] = 0;
            H0[i] = 0;
            H1[i] = 0;
        }

        #pragma unroll 4 
        for(long i=tid, start=bid*step; i<step; i+=th)
        {
            S1[i]=refs[i+start];
        }

        __syncthreads();

        #pragma unroll 4 
        for(long j=0; j<step; j++)
        {

            #pragma unroll 4
            for(long tmp=20; tmp>=0; tmp--)
            {
                S1[j] = refs[j+tmp];
            }

            char S1j = S1[j];

            #pragma unroll 
            for(long I=tid, I1=tid+1; I<N; I+=th, I1+=th)
            //for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            //long i=tid, i1=tid+1;
            {
                long i=I%D, i1=I1%D;
                // update F
                x = H0[i1]-Go;
                y = F[i1]-Ge;
                F[i1] = max(x, y);
                F[i1] = x;

                // update H'
                //x = S0[i]!=S1j;
                //x = H0[i] + 4-x*2;
                //x = H0[i] + 2;
                //x = H0[i] + 4*(S0[i]!=S1j) - 2;
                x = H0[i] + (S0[i]!=S1j)?-2:2;
                y = F[i1];
                x = max(x, y);
                x = max(x, 0);
                H1[i1] = x;
                H0[i] = -1;
            }
            __syncthreads();

            for(long I=tid, I1=tid+1; I<N; I+=th, I1+=th)
            //for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            {
                long i=I%D1, i1=i+1;
                // find affect range of each elem
                x = H1[i] - Go;
                //#pragma unroll 4
                //for(long k=i1; k<N; k++)
                for(long K=I1; K<N; K++)
                {
                    long k=K%D1;
                    if(x<=H1[k])
                    {
                        break;
                    }
                    else{
                        H0[k] = 1;
                        x -= Ge;
                    }
                }
            }
            __syncthreads();
            //continue;

            for(long I=tid, I1=tid+1; I<N; I+=th, I1+=th)
            //for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            {
                long i=I%D1, i1=i+1;
                if(H0[i]==-1&&H0[i1]==1)
                {
                    x = H1[i1] - Go;
                    #pragma unroll 4
                    //for(long k=i1; k<N; k++)
                    for(long K=I1; K<N; K++)
                    {
                        long k=K%D1;
                        if(H0[k]!=1)
                        {
                            break;
                        }else{
                            H1[k] = x;
                            x -= Ge;
                        }

                    }
                }

            //continue;
            H0[i1] = H1[i1];
            if(i1==(N-2)%D1)
            {
                S1[j] = H0[i1];
            }
            }
        }
    }
}



