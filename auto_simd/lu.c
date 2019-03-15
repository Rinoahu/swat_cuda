#include "omp.h"

// lu.c
inline int _(int row, int col, int rows){
    return row*rows + col;
}

void det_by_lu(double *y, double *x, int N){
    int i,j,k;
    *y = 1.;

    for(k = 0; k < N; ++k){
        *y *= x[_(k,k,N)];
        for(i = k+1; i < N; ++i){
            x[_(i,k,N)] /= x[_(k,k,N)];
        }
        for(i = k+1; i < N; ++i){
            
            #pragma omp simd
	    for(j = k+1; j < N; ++j){
                x[_(i,j,N)] -= x[_(i,k,N)] * x[_(k,j,N)];
            }
        }
    }
}

