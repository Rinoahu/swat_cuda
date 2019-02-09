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
__global__ void swat1(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    __shared__ char S0[4097], S1[4097];
    __shared__ short F[4097], H0[4097], H1[4097];
    //__shared__ int H[4097];
    short x, y;
    //char c0, c1;
    const long tid=threadIdx.x, bid=blockIdx.x, bdm=blockDim.x;
    //long tid=threadIdx.x, bid=blockIdx.x, bdm=blockDim.x;
    //const long idx=tid + bid * bdm;
    //const long th=1024, step=1024*16, bth=256;
    const long th=1024, step=1024*8, bth=256, B=32, W=(N+31)/B;
    //const long tid = threadIdx.x + blockIdx.x * blockDim.x;
    //printf("hello tid %ld bid %ld M %ld \n", tid, bid, M);
    //__syncthreads();

    // load all data into shr array;
    for(long i=tid; i<N; i+=th)
    {
        S0[i] = qry[i];
        F[i] = 0;
        H0[i] = 0;
        H1[i] = 0;
    }

    __syncthreads();
    //return;

    //for(long start=bid*step; start<N; start+=step)
    //for(long start=bid*step; start<M; start+=step * 32)
    for(long bi=bid, bend=(M+step-1)/step; bi<bend; bi+=bth)
    {
        /*
        if(tid==0)
        {
            //printf("tid %ld bid %ld %ld\n",tid, bid, bi);
            printf("bi %ld %ld\n", bi, N);
        }
        */

        //printf("hello tid %ld bid %ld bdm %ld \n", tid, bid, bdm);
        //printf("start %ld end %ld N %ld tid %ld bid %ld bdm %ld \n", start, start+step, N, tid, bid, bdm);
        // get ref seq
        //for(long j=start, end=start+step; j<end; j++)
        //for(long tmp=0;tmp<8;tmp++)
        //{
        //S1[tid] = refs[bi*step+tid];
        //__syncthreads();
        //continue;
        //long flag = 0;

        /*
        long offset=bi*step;
        for(long i=tid; i<N; i+=th)
        {
            S1[i] = refs[i+offset];
        }
        __syncthreads();
        */

        //for(long j=0; j<step; j++)
        for(long j=bi*step, end=(bi+1)*step; j<end; j++)
        {
            char S1j = refs[j];
            //char S1j = S1[j];
            /*
            for(long i=tid; i<N; i+=th)
            {
                H0[i]+=1;
                //H[i];
                __syncthreads();
                /*
                if(tid==0)
                {
                    printf("tid %ld j %ld bi %ld bid %ld th %ld N %ld bend %ld start %ld end %ld \n", i, bi*step, bi, bid, th, N, bend, bi*step, end);
                }
                flag += 1;
            }
            */
            //__syncthreads();
            //continue;
            // update F
            H0[N] = 0;
            for(long i=tid; i<N; i+=th)
            {
                // update F
                x = H0[i+1]-Go;
                y = F[i+1]-Ge;
                //F[i+1] = max(x, y);
                x = x ^ ((x ^ y) & -(x < y));
                F[i+1] = x;


                // update H'
                //x = H0[i] + S0[i]==refs[j]?2:-2;
                //x = H0[i] + S0[i]==S1j?2:-2;
                x = H0[i] + 4*(S0[i]==S1j) - 2;
                y = F[i+1];
                //x = max(x, 0);
                //H1[i+1] = max(x, y);
                x = x ^ ((x ^ y) & -(x < y));
                x = x ^ x  & -x;
                H1[i+1] = x;
                __syncthreads();

                // find update start site
                x = (H1[i]-Go>H1[i+1]);
                H0[i] = x*i-(!x);

                if(H0[i]==i)
                {
                    short H1i = H1[i];
                    #pragma unroll
                    for(long i1=i+1; i1<N+1; i1++)
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
                long i1 = H0[i0];
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



            }

            continue;
            //printf("hello\n");
            H0[N] = 0;
            for(long i=tid; i<N+1; i+=th)
            {
                // find the potential point that influence its next and set the id to id
                //x = (H1[i]-Go>H1[i+1]);
                //H0[i] = x*i-(!x);

                /*
                if(i<N)
                {
                    H0[i] = (H1[i]-Go>H1[i+1])?i:-1;
                }else{
                    H0[i] = -1;
                }
                */

            }

            __syncthreads();
            continue;

            // find the range of each potential point
            for(long i=tid; i<N; i+=th)
            {
                if(H0[i]==i)
                {
                    short H1i = H1[i];
                    for(long i1=i+1; i1<N+1; i1++)
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
            }

            // update H 
            for(long i=tid; i<N; i+=th)
            {
                long i0 = i;
                long i1 = H0[i0];
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
        }
        /*
        if(tid==0)
        {
            printf("tid %ld bid %ld flag %ld\n", tid, bid, flag);
        }
        __syncthreads();
        */

        //}
    }
}
}


__global__ void swat_print0(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, gdx=gridDim.x, bid=bx+by*gdx;

    if(tid==0)
    {
        printf("x %ld y %ld bid %ld\n", bx, by, bid);
        __syncthreads();
    }

    /*
    if(threadIdx.x==1023)
    {
        printf("Hello from block %d, %d thread %d\n", blockIdx.x, blockIdx.y, threadIdx.x);
    }
    */
}
 
// set the "thread per block" to 128
__global__ void swat_print1(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    const int D=4097;
    //__shared__ char S0[D];
    __shared__ char S0[D], S1[D];
    __shared__ short F[D], H0[D], H1[D];
    //__shared__ int H[4097];
    short x, y;
    //char c0, c1;
    //const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, dmx=blockDim.x, dmy=blockDim.y, bid=by+bx*dmy;
    const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, bdmx=blockDim.x, bid=bx+by*gridDim.x;
    const long th=bdmx, step=1024, bth=256;

    /* 
    if(tid==0)
    {
        //printf("bx %ld bid %ld M %ld dmx %ld\n", bx, bid, M/step, dmx);
        printf("th %ld \n", th);

    }
    */

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
        for(long i=tid; i<step; i+=th)
        {
            S1[i] = refs[bid*step+i];
        }

        __syncthreads();

        //return;
        #pragma unroll 4 
        for(long j=0; j<step; j++)
        //for(long j=bid*step, end=j+step; j<end; j++)
        {
            char S1j = S1[j];
            //char S1j = refs[j];

            __syncthreads();
            //continue;

            // update F
            #pragma unroll 4
            for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            //for(long i0=0, i1, i, end=N/th; i0<end; i0++)
            //long i0 = 0, i1, i;
            //if(tid<N)
            {
                //i = tid + i0*th;
                // update F
                //i1 = i+1;
                x = H0[i1]-Go;
                y = F[i1]-Ge;
                F[i1] = max(x, y);
                //x = x ^ ((x ^ y) & -(x < y));
                F[i1] = x;
                //continue;

                // update H'
                x = H0[i] + 4*(S0[i]==S1j) - 2;
                y = F[i1];
                x = max(x, y);
                x = max(x, 0);
                H1[i1] = x;

                // use prefix scan to update H1
                // find update start site
                x = (H1[i]-Go>H1[i1]);
                H0[i] = x*i-(!x);
                //H0[i] =  H1[i]-Go>H1[i1]?i:-1;


                if(i==0)
                {
                    H0[N] = -1;
                }

                __syncthreads();

                if(H0[i1-1] = -1)
                {
                    H0[i1] = -1;
                }

                __syncthreads();

                //continue;
  
                /* 
                long I1 = i1;
                short H1i = H1[i];
                //pragma unroll
                for(i1=i+1; i1<N+1; i1++)
                {
                    short H1i1 = H1i - Go - (i1-i-1) * Ge;
                    short H0i1_old = H0[i1];

                    if(H1i1<=H1[i1])
                    {
                        break;
                    }else{
                        H0[i1] = i;
                    }

                    if(H0i1_old != -1)
                    {
                        break;
                    }
                    
                    if(H1i1<=H1[i1])
                    {
                        break;
                    }else if(H0[i1]==i1)
                    {
                        H0[i1] = i;
                        break;
                    }else
                    {
                        H0[i1] = i;
                    }
                }
                */

                if(H0[i]==i)
                {
                    short H1i = H1[i];
                    #pragma unroll 4 
                    for(i1=i+1; i1<N+1; i1++)
                    {
                        short H1i1 = H1i - Go - (i1-i-1) * Ge;
                        short H0i1_old = H0[i1];
                        //H0i1_old = H0[i1];

                        if(H1i1<=H1[i1])
                        {
                            break;
                        }else{
                            H0[i1] = i;
                        }

                        if(H0i1_old != -1)
                        {
                            break;
                        }

                        //H1i1 -= Ge;

                        /*
                        if(H1i1<=H1[i1])
                        {
                            break;
                        }else if(H0[i1]==i1)
                        {
                            H0[i1] = i;
                            break;
                        }else
                        {
                            H0[i1] = i;
                        }
                        */
                    }
                }

                __syncthreads();
                continue;
                // update H 
                long a = i1, b = H0[a];
                while(b!=-1)
                {
                    //printf("start %ld %ld %ld %ld %ld %ld %ld %ld I J\n", start, N, idx, tid, bid, bdm, I, j);
                    if(a==b)
                    {
                        H1[i1] = (i1==i)?H1[i1]: H1[b]-Go-(i1-i-1)*Ge;
                        break;
                    }else{
                        a = b;
                        b = H0[a];
                    }
                }
                //H0[I1] = H1[I1];
                H0[i] = H1[i];
                //continue;
            }
        }
    }
}







// set the "thread per block" to 128
__global__ void swat_print(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    const int D=4097;
    //__shared__ char S0[D];
    __shared__ char S0[D], S1[D];
    __shared__ short F[D], H0[D], H1[D];
    //__shared__ int H[4097];
    short x, y;
    //const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, dmx=blockDim.x, dmy=blockDim.y, bid=by+bx*dmy;
    const long tid=threadIdx.x, bx=blockIdx.x, by=blockIdx.y, bdmx=blockDim.x, bid=bx+by*gridDim.x;
    const long th=bdmx, step=1024, bth=256;

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
        for(long i=tid; i<step; i+=th)
        {
            S1[i] = refs[bid*step+i];
        }

        __syncthreads();

        #pragma unroll 4 
        for(long j=0; j<step; j++)
        {
            char S1j = S1[j];

            #pragma unroll 
            //for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
            long i=tid, i1=tid+1;
            {
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

            //for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
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

            //for(long i=tid, i1=tid+1; i<N; i+=th, i1+=th)
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
            }
            if(i1==N-2)
            {
                S1[j] = H0[i1];
            }
        }
    }
}


__global__ void swat2(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    __shared__ char S0[4097];
    __shared__ int16_t F[4097], H[4097];
    int16_t x, y;
    //char c0, c1;
    const long tid=threadIdx.x, bid=blockIdx.x, bdm=blockDim.x;
    //long tid=threadIdx.x, bid=blockIdx.x, bdm=blockDim.x;
    //const long idx=tid + bid * bdm;
    const long th=32, step=4096, bth=1024, blk=bth*th;

    //const long tid = threadIdx.x + blockIdx.x * blockDim.x;

    //printf("hello tid %ld bid %ld M %ld \n", tid, bid, M);
    //__syncthreads();

    // load all data into shr array;
   __syncthreads();

    for(long i=tid; i<N;i+=th)
        {
            S0[i] = qry[i];
        }

    __syncthreads();

    //for(long start=bid*step; start<N; start+=step)
    //for(long start=bid*step; start<M; start+=step * 32)
    for(long bi=tid, bend=(M+step-1)/step; bi<bend; bi+=blk)
    {
        ///*
        if(tid==0)
        {
            printf("bid %ld bdm %ld\n", bid, bdm);
        }
        //*/
        for(long i=0; i<N;i++)
        {
            F[i] = 0;
            H[i] = 0;
        }

        for(long j=bi*step, end=(bi+1)*step; j<end; j++)
        {
            char S1j = refs[j];
            // update F
            short E=0;
            for(long i=0; i<N; i++)
            {
                // update F
                x = H[i+1]-Go;
                y = F[i+1]-Ge;
                //F[i+1] = max(x, y);
                x = x ^ ((x ^ y) & -(x < y));
                F[i+1] = x;

                y = H[i] + 4*(S0[i]==S1j) - 2;
                x = x ^ ((x ^ y) & -(x < y));
                x = x ^ x  & -x;
                y = max(E-Ge, x-Go);
                y = max(0, y);
                H[i+1] = max(x, y);

            }
        }
    }
}


// set the "thread per block" to 128
__global__ void swat3(char *qry, long N, char *refs, long M, short Go, short Ge)
{
    __shared__ char S0[4097];
    __shared__ int16_t F[4097], H0[4097], H1[4097];
    int16_t x, y;
    //char c0, c1;
    const long tid=threadIdx.x, bidx=blockIdx.x, bidy=blockIdx.y, bdm=blockDim.x;
    //long tid=threadIdx.x, bid=blockIdx.x, bdm=blockDim.x;
    //const long idx=tid + bid * bdm;
    const long th=1024, step=4096, bth=1024;
    const long bid=bidx+bidy*bdm;
    //const long tid = threadIdx.x + blockIdx.x * blockDim.x;

    //__syncthreads();

    // load all data into shr array;
    for(long i=tid; i<N;i+=th)
    {
        S0[tid] = qry[tid];
        F[tid] = 0;
        H0[tid] = 0;
        H1[tid] = 0;
    }
    __syncthreads();

    //for(long start=bid*step; start<N; start+=step)
    //for(long start=bid*step; start<M; start+=step * 32)
    //for(long bi=bid, bend=(M+step-1)/step; bi<bend; bi+=bth)
    //{
        ///*
        if(tid==0)
        {
            //printf("tid %ld bid %ld %ld\n",tid, bid, bi);
            printf("bid %ld bdm %ld bid %ld M %ld\n", bidx, bdm, bid, M/step);

        }
        //*/

        //printf("hello tid %ld bid %ld bdm %ld \n", tid, bid, bdm);
        //printf("start %ld end %ld N %ld tid %ld bid %ld bdm %ld \n", start, start+step, N, tid, bid, bdm);
        // get ref seq
        //for(long j=start, end=start+step; j<end; j++)
        //for(long tmp=0;tmp<8;tmp++)
        //{
        if(bid>=M/step)
        {
            return;
        }
        for(long j=bid*step, end=(bid+1)*step; j<end; j++)
        {

            char S1j = refs[j];
            // update F
            for(long i=tid; i<N; i+=th)
            {
                // update F
                x = H0[i+1]-Go;
                y = F[i+1]-Ge;
                //F[i+1] = max(x, y);
                x = x ^ ((x ^ y) & -(x < y));
                F[i+1] = x;


                // update H'
                //x = H0[i] + S0[i]==refs[j]?2:-2;
                //x = H0[i] + S0[i]==S1j?2:-2;
                x = H0[i] + 4*(S0[i]==S1j) - 2;
                y = F[i+1];
                //x = max(x, 0);
                //H1[i+1] = max(x, y);
                x = x ^ ((x ^ y) & -(x < y));
                x = x ^ x  & -x;
                H1[i+1] = x;

            }

            //continue;
            //printf("hello\n");
            H0[N] = 0;
            for(long i=tid; i<N+1; i+=th)
            {
                // find the potential point that influence its next and set the id to id
                x = (H1[i]-Go>H1[i+1]);
                H0[i] = x*i-(!x);

                /*
                if(i<N)
                {
                    H0[i] = (H1[i]-Go>H1[i+1])?i:-1;
                }else{
                    H0[i] = -1;
                }
                */

            }

            __syncthreads();
            //continue;

            // find the range of each potential point
            for(long i=tid; i<N; i+=th)
            {
                if(H0[i]==i)
                {
                    short H1i = H1[i];
                    for(long i1=i+1; i1<N+1; i1++)
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
            }

            // update H 
            for(long i=tid; i<N; i+=th)
            {
                long i0 = i;
                long i1 = H0[i0];
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
            }
            __syncthreads();
        }
        //}
    //}
}


