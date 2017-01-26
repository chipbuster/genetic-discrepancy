/*
 * --------------------------------------------------------------------------
 * BLISLAB 
 * --------------------------------------------------------------------------
 * Copyright (C) 2016, The University of Texas at Austin
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met:
 *  - Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  - Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 *  - Neither the name of The University of Texas nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *
 * bl_dgemm_ref.c
 *
 *
 * Purpose:
 * implement reference mkl using GEMM (optional) in C.
 *
 * Todo:
 *
 *
 * Modification:
 *
 * 
 * */


#include <omp.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>

#include <bl_dgemm.h>
#include <bl_dgemm_ref.h>
#include <mkl.h>

#ifdef USE_BLAS
/* 
 * dgemm prototype
 *
 */ 
//void dgemm(char*, char*, int*, int*, int*, double*, double*, 
//        int*, double*, int*, double*, double*, int*);
extern void dgemm_(char*, char*, int*, int*, int*, double*, double*, 
        int*, double*, int*, double*, double*, int*);
#endif

unsigned int calcFitness(const float* P, const float* B,
                         unsigned int n, unsigned int M, int d, unsigned int m,
                         float* fitness, float *max_D, unsigned int *max_pos) {

  double total_g = 0;
  int print_count = 0;
#pragma omp parallel for reduction(+:total_g)
    for (unsigned int i = 0; i < M; ++i) {
      int SSum = 0;
      for (unsigned int j = 0; j < n; ++j) {
        bool Csub = true;
        for (int k = 0; k < d; ++k) {
           Csub &= P[j*d + k] <= B[i*d + k];
        }
        if (Csub) SSum++;
      }
    }
    return SSum;   

}

unsigned int calcFitnessFullMat(const float* P, const float* B,
                         unsigned int n, unsigned int M, int d, unsigned int m,
                         float* fitness, float *max_D, unsigned int *max_pos) {

// P = n x d
// B = M x d
// Equivalent problem = P * B.T
// m = 0 --> M

  double total_g = 0;
  int print_count = 0;
#pragma omp parallel for reduction(+:total_g)
    for (unsigned int i = 0; i < M; ++i) {
      int SSum = 0;
      for (unsigned int j = 0; j < n; ++j) {
        bool Csub = true;
        for (int k = 0; k < d; ++k) {
           Csub += P[j*d + k] <= B[i*d + k];
        }
        if (Csub) SSum++;
      }
    }
    return SSum;

}

void bl_dgemm_ref(
        int    m,
        int    n,
        int    k,
        double *XA,
        int    lda,
        double *XB,
        int    ldb,
        double *XC,
        int    ldc
        )
{

    // XA = m x k
    // XB = k x n
    // XC = m x n

    // Local variables.
    int    i, j, p;
    double beg, time_collect, time_dgemm, time_square;
    double alpha = 1.0, beta = 1.0;

    // Sanity check for early return.
    if ( m == 0 || n == 0 || k == 0 ) return;

    // Reference GEMM implementation.
    beg = omp_get_wtime();

    #pragma omp parallel for private( i, p )
    for ( j = 0; j < n; j ++ ) {
        for ( i = 0; i < m; i ++ ) {
            for ( p = 0; p < k; p ++ ) {
                XC[ j * ldc + i ] += XA[ p * lda + i ] * XB[ j * ldb + p ];
            }
        }
    }
    time_dgemm = omp_get_wtime() - beg;

}

