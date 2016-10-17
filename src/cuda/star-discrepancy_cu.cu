#include <cmath>
#include <cstdio>
#include <cfloat>
#include <iostream>
#include <omp.h>

#include <cublas_v2.h>
#include <thrust/device_vector.h>

#include "timer.h"

#define CPU_CHECK

INIT_TIMER(cu);

/*
   Star discrepancy is defined as follows:


   For a given pointset, P, of size N and dimension d:
   WLOG, assume that the convex hull of P has a volume of 1

                ( |                                            |                         )
                | | S | \forall s \in S, s_i <= p_i for i=1..d |                         |
                | |                                            |                         |
   max_{p\in P} | ----------------------------------------------  -  \prod p_i (%volume) |
                |                                                                        |
                |                      N                                                 |
                (                                                                        )

   This can be computed as follows
   (let MMX be the max value at any dimension)

   max_disc = 0
   for i = 1..N:
     this_count = 0
     for j = 1..M:
       continue if i == j

       inside = true
       for k = 1..d:
         if S[i][k] > S[j][k]
           inside = false
           break

       if inside:
         this_count++

     this_vol = 1
     for k = 1..d:
       this_vol *= s[i][j]/MMX

     max_disc = max(max_disc, this_count/N - this_vol)

   Taken from the CUDA sample NVIDIA_CUDA-7.0_Samples/0_Simple/matrixMul/matrixMul.cu

   Performs a matrix multiplication of SSums = S*S
   SSums needs to be a float for the cublas library

   wA = d, wB = n
 */
template<int BLOCK_SIZE> __global__ void 
matrixMul(
    float *S, float *SSums,
    unsigned int n, unsigned int d) {
  // Note: Need to have SSums initialized to 0 before this function is called

  __shared__ float shared_counts[BLOCK_SIZE][BLOCK_SIZE];

  // Block index
  int bx = blockIdx.x;
  int by = blockIdx.y;

  // Thread index
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  // global_x is the column number for the C matrix, global_y is the row number
  // global_y is also the final sample to add results to SSums
  int global_x = BLOCK_SIZE * bx + tx,
      global_y = BLOCK_SIZE * by + ty;

  // Index of the first sub-matrix of S processed by the block
  int aBegin = d * BLOCK_SIZE * by;
  int bBegin = d * BLOCK_SIZE * bx;

  // Index of the last sub-matrix of S processed by the block
  int aEnd   = aBegin + d - 1;

  // Step size used to iterate through the sub-matrices
  int aStep  = BLOCK_SIZE;
  int bStep  = BLOCK_SIZE;

  // Csub is used to store the element of the block sub-matrix
  // that is computed by the thread
  bool Csub = true;

  // Loop over all the sub-matrices of A and B
  // required to compute the block sub-matrix
  for (int a = aBegin, b = bBegin; 
       a <= aEnd; 
       a += aStep, b += bStep) {

    // Declaration of the shared memory array As used to
    // store the sub-matrix of A
    // Should be an nxd matrix (or smaller)
    __shared__ float S1[BLOCK_SIZE][BLOCK_SIZE];

    // Declaration of the shared memory array Bs used to
    // store the sub-matrix of B
    __shared__ float S2[BLOCK_SIZE][BLOCK_SIZE];

    // Load the matrices from device memory
    // to shared memory; each thread loads
    // one element of each matrix
    S1[ty][tx] = S[a + d * ty + tx];
    S2[ty][tx] = S[b + d * ty + tx];

    // Synchronize to make sure the matrices are loaded
    __syncthreads();

    /*
       if (bx * BLOCK_SIZE + tx == 16 && by * BLOCK_SIZE + ty == 9984) {
       for (int i = 0; i < BLOCK_SIZE; ++i) {
       printf("%d: ", bx * BLOCK_SIZE + i);
       for (int j = 0; j < d; ++j) {
    //printf("% 4d,% 4d ", S1[ty][j], S2[tx][j]);
    printf("% 4d,% 4d ", S1[i][j], S2[i][j]);
    }
    printf("\n");
    }
    }
     */

    // Multiply the two matrices together;
    // each thread computes one element
    // of the block sub-matrix
#pragma unroll

    for (int k = 0; k < BLOCK_SIZE && k < d; ++k)
    {
      // ty corresponds to the sample# in S1 and tx corresponds to sample# in S2
      // bitwise and of 
      Csub &= S2[tx][k] <= S1[ty][k];
    }

    // Synchronize to make sure that the preceding
    // computation is done before loading two new
    // sub-matrices of A and B in the next iteration
    __syncthreads();
  }

  // Write the block sub-matrix to device memory;
  // each thread writes one element
  //int c = n * BLOCK_SIZE * by + BLOCK_SIZE * bx;


  // If it's not entirely within, and it's not the same x and y, add one to the
  // value of SSums.
  if (Csub && global_x != global_y && global_y < n && global_x < n) {
    shared_counts[ty][tx] = 1;
  } else {
    shared_counts[ty][tx] = 0;
  }

  /*
     if (global_x != global_y)
     C[c + n * ty + tx] = Csub ? 1.0 : 0.0;
     else
     C[c + n * ty + tx] = 0;
   */
  //C[c + n * ty + tx] = -1;


  // Do a reduction on the shared_counts matrix, sum across the rows, so keep ty (S1) constant
#pragma unroll
  for (int i = 2; i <= BLOCK_SIZE; i *= 2) {
    if (tx % i == 0) {
      shared_counts[ty][tx] += shared_counts[ty][tx + i/2];
    }
  }
  // Only have first proc perform the atomicAdd
  if (tx == 0) {
    if (global_y < n) {
      atomicAdd(&SSums[global_y], shared_counts[ty][0]);
    }
  }

  //if (Csub) 
  //  printf("for i=%d and j=%d (%d=%d), Csub =%c\n", ty, tx, c + n * ty + tx, 8561, Csub?'Y':'N');
  //if (global_x == 16 && global_y == 9984) 
  //  printf("**** for i=%d and j=%d (%d=%d), Csub =%c errd on %d\n", global_y, global_x, c + n * ty + tx, 8561, Csub?'Y':'N', Cerrd);
}

/**
 * Cuda function that performs SSums = A * B for two different matrices
 */
template <int BLOCK_SIZE> __global__ void
matrixMul_Diff(float *SSums, float *A, float *B,
               int nA, int nB, int d, int m)
{
  // Needs to be float so we can add it to SSums
  __shared__ float shared_counts[BLOCK_SIZE][BLOCK_SIZE];

  // Block index
  int bx = blockIdx.x;
  int by = blockIdx.y;

  // Thread index
  int tx = threadIdx.x;
  int ty = threadIdx.y;

  // global_A is the column number for the C matrix, global_y is the row number
  // global_B is also the final sample to add results to SSums
  //int global_A = BLOCK_SIZE * by + ty;
  int global_B = BLOCK_SIZE * bx + tx;

  // Index of the first sub-matrix of A processed by the block
  int aBegin = d * BLOCK_SIZE * by;
  // Index of the first sub-matrix of B processed by the block
  int bBegin = d * BLOCK_SIZE * bx;

  // Index of the last sub-matrix of A processed by the block
  int aEnd   = aBegin + d - 1;

  // Step size used to iterate through the sub-matrices of A
  int aStep  = BLOCK_SIZE;
  // Step size used to iterate through the sub-matrices of B
  int bStep  = BLOCK_SIZE;

  // Csub is used to store the element of the block sub-matrix
  // that is computed by the thread
  bool Csub = true;
  // Initialize this value to zero.
  shared_counts[tx][ty] = 0;

  // Loop over all the sub-matrices of A and B
  // required to compute the block sub-matrix
  for (int a = aBegin, b = bBegin;
       a <= aEnd;
       a += aStep, b += bStep)
  {

    // Declaration of the shared memory array As used to
    // store the sub-matrix of A
    __shared__ float S1[BLOCK_SIZE][BLOCK_SIZE];

    // Declaration of the shared memory array Bs used to
    // store the sub-matrix of B
    __shared__ float S2[BLOCK_SIZE][BLOCK_SIZE];

    // Load the matrices from device memory
    // to shared memory; each thread loads
    // one element of each matrix
    //As[ty][tx] = A[a + wA * ty + tx];
    //Bs[ty][tx] = B[b + wB * ty + tx];
    //S1[ty][tx] = A[a + d * ty + tx];
    //S2[ty][tx] = B[b + d * ty + tx];
    float aval, bval;
    int ib = b + d * ty + tx;
    if (ib < nB * d) {
      //printf("B at %u,%u is %u %u %u\n",
      //     ib/d, ib%d, B[ib-1], B[ib], B[ib+1]);
      bval = B[ib];
    } else {
      bval = -1;
    }
    int ia = a + d * ty + tx;
    if (ia < nA * d) {
      //printf("A at %u,%u is %u %u %u\n",
      //       ia/d, ia%d, A[ia-1], A[ia], A[ia+1]);
      aval = A[ia];
    } else {
      aval = m + 1;
    }
    S1[ty][tx] = aval;
    S2[ty][tx] = bval;

    // Synchronize to make sure the matrices are loaded
    __syncthreads();

    // Multiply the two matrices together;
    // each thread computes one element
    // of the block sub-matrix
#pragma unroll

    for (int k = 0; k < BLOCK_SIZE && k < d; ++k)
    {
      /*
      //Csub += As[ty][k] * Bs[k][tx];
      if (S1[tx][k] == -1 || S2[ty][k] == -1) {
        int ix = a + d * tx + k;
        int iy = b + d * ty + k;
        //printf("[%d,%d/%d,%d] Invalid at %u,%u:%u,%u for %d,%d\n", 
        //       tx,ty,bx,by, ix/d,ix%d, iy/d,iy%d, S1[tx][k],S2[ty][k]);
        Csub = false;
      }
      else
      */
        Csub &= S1[ty][k] <= S2[tx][k];
    }

    // Synchronize to make sure that the preceding
    // computation is done before loading two new
    // sub-matrices of A and B in the next iteration
    __syncthreads();
  }
  
  //// Write the block sub-matrix to device memory;
  //// each thread writes one element
  //int c = wB * BLOCK_SIZE * by + BLOCK_SIZE * bx;
  //C[c + wB * ty + tx] = Csub;

  //int c = nA * BLOCK_SIZE * by + BLOCK_SIZE * bx;

  // If it's entirely within the bounds of both matrices, add one to the value of SSums.
  //shared_counts[ty][tx] += (Csub && global_A < nA && global_B < nB);
  //shared_counts[ty][tx] += (Csub && global_B < nA && global_A < nB);
  shared_counts[tx][ty] = Csub ? 1 : 0;
  

  // Do a reduction on the shared_counts matrix, sum across the rows,
  // so keep ty (S1) constant
#pragma unroll
  for (int i = 2; i <= BLOCK_SIZE; i *= 2) {
    __syncthreads();
    if (ty % i == 0) {
      shared_counts[tx][ty] += shared_counts[tx][ty + i/2];
    }
  }
  /*
  if (shared_counts[tx][0] > 0) {
    printf("found %.0f at %d (%d) bs=%d,bx=%d,by=%d,tx=%d,ty=%d\n",
           shared_counts[tx][0], global_B, global_A,
           BLOCK_SIZE,bx,by,tx,ty);
  }
  */
  // Only have first proc perform the atomicAdd
  if (ty == 0) {
    if (global_B < nB) {
      atomicAdd(&SSums[global_B], shared_counts[tx][0]);
    }
  }
}

void matMulCPU_Diff(float *SSums, float *A, float *B,
                    unsigned int nA, unsigned int nB, unsigned int d) {
  for (int i = 0; i < nB; ++i) {
    int SSum = 0;
    for (int j = 0; j < nA; ++j) {
      bool Csub = true;
      for (int k = 0; k < d; ++k) {
        Csub &= A[j*d + k] <= B[i*d + k];
      }

      //C[i*n + j] = Csub ? 1.f : 0.f;
      //if (Csub) printf("A:%d,B:%d is 1\n", j, i);
      if (Csub) SSum++;
    }

    SSums[i] = SSum;
  }
}

void matMulCPU(float *SSums, float *S,
               unsigned int n, unsigned int d) {
  for (int i = 0; i < n; ++i) {
    int SSum = 0;
    for (int j = 0; j < n; ++j) {
      if (i == j) {
        // C[i*n + j] = 0.f;
        continue;
      }
      bool Csub = true;
      for (int k = 0; k < d; ++k) {
        Csub &= S[j*d + k] <= S[i*d + k];
      }

      /*
         if (i == 9984 && j == 16) {
         printf("\n\ni (% 6d):", i);
         for (int pp = 0; pp < d; pp++) {
         printf(" % 3d", S[i*d + pp]);
         }
         printf("\nj (% 6d):", j);
         for (int pp = 0; pp < d; pp++) {
         printf(" % 3d", S[j*d + pp]);
         }
         printf("\n");
         }
       */

      //C[i*n + j] = Csub ? 1.f : 0.f;
      if (Csub) SSum++;
    }

    SSums[i] = SSum;
  }
}

#define HANDLE_ERROR( err ) ( HandleError( err, __FILE__, __LINE__ ) )
static void HandleError( cudaError_t err, const char *file, int line ) {
  if (err != cudaSuccess) {
    fprintf(stderr, "%s in %s at line %d\n",
            cudaGetErrorString( err ), file, line );
    exit( EXIT_FAILURE );
  }
}

void checkCudaError(const char *msg) {
  // make the host block until the device is finished with foo
  cudaThreadSynchronize();

  cudaError_t err = cudaGetLastError();
  if( cudaSuccess != err) {
    fprintf(stderr, "Cuda error: %s: %s (%d).\n", 
            msg, cudaGetErrorString(err), err);
    exit(EXIT_FAILURE);
  }
}

/*
 * Calculates the 'fitness score' of each box in B, according to step (ii) and 
 * (iii) of Shah paper.
 *
 * Output will be in fitness vector, which should be of size M
 * fitness_cu needs to be a float for cublas
 *
 * Will return the number of fitness that are not 0.
 */
unsigned int DStarFitness(float *P, float *P_cu, 
                   float *B, float *B_cu,
                   unsigned int n, unsigned int M,
                   unsigned int d, unsigned int m,
                   float *fitness, float* fitness_cu,
                   float &max_D) {
  // Define the blocksize and the gridsize
  int block_size = 32;
  dim3 threads(block_size, block_size);
  int grid_size_y = n < block_size ? 1 : ceil((float)n/block_size);
  int grid_size_x = M < block_size ? 1 : ceil((float)M/block_size);
  dim3 grid(grid_size_x, grid_size_y);
  printf("Block size is %dx%d and grid is %dx%d\n", block_size,block_size, grid_size_x,grid_size_y);

  START_TIMER(cu);
  HANDLE_ERROR( cudaMemset(fitness_cu, 0, M*sizeof(float)) );
  // Also need to copy this vector over.
  HANDLE_ERROR( cudaMemcpy(B_cu, B, M*d*sizeof(float), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(P_cu, P, n*d*sizeof(float), cudaMemcpyHostToDevice) );
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"cudacp");

  START_TIMER(cu);
  matrixMul_Diff<32><<< grid, threads >>>(fitness_cu, P_cu, B_cu, n, M, d, m);
  checkCudaError("matrixMul fitness_Diff");
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"GPU execute");

#ifdef CPU_CHECK
  float* fitness_dev = new float[M];
  HANDLE_ERROR( cudaMemcpy(fitness_dev, fitness_cu, M*sizeof(float), cudaMemcpyDeviceToHost) );
  float* fitness_test = new float[M];
  START_TIMER(cu);
  matMulCPU_Diff(fitness_test, P, B, n, M, d);
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"CPU execute");

  for(int i = 0; i < M; ++i) {
    //for (int j = 0; j < n; ++j) {
    //    if (C[i*n + j] != C_dev[i*n + j]) {
    //      printf("Error, difference at %d,%d: %.0f\n", i, j, C[i*n+j]);
    //    }
    //}
    if (fitness_dev[i] != fitness_test[i]) {
      printf("Error, difference at %d: (CPU)%.0f vs (GPU)%.0f\n", i, fitness_test[i], fitness_dev[i]);
    }
  }
  fflush(stdout);
#endif

  START_TIMER(cu);
  HANDLE_ERROR ( cudaMemcpy(fitness, fitness_cu, M*sizeof(float), cudaMemcpyDeviceToHost) );

  int print_count = 0;

  // Find the g(B_i; P)
  double total_g = 0;
  for (int i = 0; i < M; ++i) {
    // compute the volume
    double volume = 1;
    for (int k = 0; k < d; ++k) {
      // Need to normalize by dividing by m
      volume *= (double)B[i*d + k]/m;
    }

    // Number of points in A divided by the volume
    double disc = std::abs(fitness[i]/n - volume);
    if ( (fitness[i] > 0 && print_count++ < 5) || i < 3) {
      fprintf(stderr,"%d: %e, %.0f: %f\n", i, volume, fitness[i], disc);
    }
    /*
       if (fitness[i] > 0) {
       printf("%d: %f, %f: %f\n", i, volume, fitness[i], disc);
       }
     */
    fitness[i] = disc;
    total_g += disc;
  }
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"calculating fitness");

  //Now, calculate the fitness, according to (iii) of Shah
  unsigned int nr = 0;
  max_D = 0;
  for (unsigned int i = 0; i < M; ++i) {
    if (fitness[i] > max_D) {
      max_D = fitness[i];
    }
    if (fitness[i] > 1.0/M*total_g) {
      fitness[i] = 1.0 / pow(1 - fitness[i], 2.0);
      nr++;
    } else {
      fitness[i] = 0.0;
    }
  }

  fprintf(stderr, "Found nr of %u\n", nr);
  return nr;
}

// Star discrepancy, defined by the number of points in A that are within points
// in B. For the Alon paper, B is the smaller matrix and A is the original
// sample points.
float DStarBoth(float *A, float *A_cu, 
                float *B, float *B_cu,
                unsigned int nA, unsigned int nB,
                unsigned int d, unsigned int m) {
  // Define the blocksize and the gridsize
  int block_size = 32;
  dim3 threads(block_size, block_size);
  int grid_size = nA < block_size ? 1 : ceil((float)nA/block_size);
  dim3 grid(grid_size, grid_size);
  printf("Block size is %dx%d and grid is %dx%d\n", block_size,block_size, grid_size,grid_size);

  // Needs to be a float for the cublas library
  float *SSums_cu;
  START_TIMER(cu);
  cudaMalloc(&SSums_cu, nB*sizeof(float));
  HANDLE_ERROR( cudaMemset(SSums_cu, 0, nB*sizeof(float)) );
  // Also need to copy this vector over.
  HANDLE_ERROR( cudaMemcpy(B_cu, B, nB*d*sizeof(float), cudaMemcpyHostToDevice) );
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"cudaMemset");

  START_TIMER(cu);
  matrixMul_Diff<32><<< grid, threads >>>(SSums_cu, A_cu, B_cu, nA, nB, d, m);
  checkCudaError("matrixMul SSums_Diff");
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"GPU execute");

#ifdef CPU_CHECK
  float* SSums_dev = new float[nB];
  HANDLE_ERROR( cudaMemcpy(SSums_dev, SSums_cu, nB*sizeof(float), cudaMemcpyDeviceToHost) );
  float* SSums_test = new float[nB];
  matMulCPU_Diff(SSums_test, A, B, nA, nB, d);

  for(int i = 0; i < nB; ++i) {
    //for (int j = 0; j < n; ++j) {
    //    if (C[i*n + j] != C_dev[i*n + j]) {
    //      printf("Error, difference at %d,%d: %.0f\n", i, j, C[i*n+j]);
    //    }
    //}
    if (SSums_dev[i] != SSums_test[i]) {
      printf("Error, difference at %d: (CPU)%.0f vs (GPU)%.0f\n", i, SSums_test[i], SSums_dev[i]);
    }
  }
#endif

  START_TIMER(cu);
  float* SSums = new float[nB];
  HANDLE_ERROR ( cudaMemcpy(SSums, SSums_cu, nB*sizeof(float), cudaMemcpyDeviceToHost) );

  int print_count = 0;

  // Find the max
  double d_max = 0;
  for (int i = 0; i < nB; ++i) {
    // compute the volume
    double volume = 1;
    for (int k = 0; k < d; ++k) {
      // Need to normalize by dividing by m
      volume *= (double)B[i*d + k]/m;
    }

    // Number of points in A divided by the volume
    double disc = std::abs(SSums[i]/nA - volume);
    if (SSums[i] > 1 && print_count++ < 31) {
      printf("%d: %e, %f: %f\n", i, volume, SSums[i], disc);
    }
    /*
       if (SSums[i] > 0) {
       printf("%d: %f, %f: %f\n", i, volume, SSums[i], disc);
       }
     */
    d_max = std::max(d_max, disc);
  }
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"Calculating dmax");

  return d_max;
}

float DStarSingle(float *S, float *S_cu, unsigned int n, unsigned int d, unsigned int m) {
  // Define the blocksize and the gridsize
  int block_size = 32;
  dim3 threads(block_size, block_size);
  int grid_size = n < block_size ? 1 : ceil((float)n/block_size);
  dim3 grid(grid_size, grid_size);
  printf("Block size is %dx%d and grid is %dx%d\n", block_size,block_size, grid_size,grid_size);

  START_TIMER(cu);
  //unsigned int *S_cu;
  //cudaMalloc(&S_cu, n*d*sizeof(unsigned int));
  //HANDLE_ERROR( cudaMemcpy(S_cu, S, n*d*sizeof(unsigned int), cudaMemcpyHostToDevice) );
  float *SSums_cu;
  cudaMalloc(&SSums_cu, n*sizeof(float));
  HANDLE_ERROR( cudaMemset(SSums_cu, 0, n*sizeof(float)) );
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"cudaMalloc");

  START_TIMER(cu);
  matrixMul<32><<< grid, threads >>>(S_cu, SSums_cu, n, d);
  checkCudaError("matrixMul SSums");
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"GPU execute");

#ifdef CPU_CHECK
  float* SSums_dev = new float[n];
  HANDLE_ERROR( cudaMemcpy(SSums_dev, SSums_cu, n*sizeof(float), cudaMemcpyDeviceToHost) );
  float* SSums_test = new float[n];
  matMulCPU(SSums_test, S, n, d);

  for(int i = 0; i < n; ++i) {
    //for (int j = 0; j < n; ++j) {
    //    if (C[i*n + j] != C_dev[i*n + j]) {
    //      printf("Error, difference at %d,%d: %.0f\n", i, j, C[i*n+j]);
    //    }
    //}
    if (SSums_dev[i] != SSums_test[i]) {
      printf("Error, difference at %d: (CPU)%.0f vs (GPU)%.0f\n", i, SSums_test[i], SSums_dev[i]);
    }
  }
#endif

  START_TIMER(cu);
  float* SSums = new float[n];
  HANDLE_ERROR ( cudaMemcpy(SSums, SSums_cu, n*sizeof(float), cudaMemcpyDeviceToHost) );

  // Find the max
  double d_max = 0;
  for (int i = 0; i < n; ++i) {
    // compute the volume
    double volume = 1;
    for (int k = 0; k < d; ++k) {
      // Need to normalize by dividing by m
      volume *= (double)S[i*d + k]/m;
    }

    double disc = std::abs(SSums[i]/n - volume);
    if (SSums[i] > 0 && i < 31) {
      printf("%d: %f, %f: %f\n", i, volume, SSums[i], disc);
    }
    /*
       if (SSums[i] > 0) {
       printf("%d: %f, %f: %f\n", i, volume, SSums[i], disc);
       }
     */
    d_max = std::max(d_max, disc);
  }
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"Calculating dmax");

  return d_max;
}

void deviceSetup(float **P_cu, float **B_cu, float **fitness_cu,
                 float *P, float *B,
                 unsigned int n, unsigned int M, unsigned int d) {
  fprintf(stderr, "Device setup with n=%u, M=%u, d=%u\n", n, M, d);
  START_TIMER(cu);
  cudaMalloc(P_cu, n*d*sizeof(float));
  cudaMalloc(B_cu, M*d*sizeof(float));
  cudaMalloc(fitness_cu, M*sizeof(float));
  HANDLE_ERROR( cudaMemcpy(*P_cu, P, n*d*sizeof(float), cudaMemcpyHostToDevice) );
  HANDLE_ERROR( cudaMemcpy(*B_cu, B, M*d*sizeof(float), cudaMemcpyHostToDevice) );
  STOP_TIMER(cu);
  PRINT_TIMER(cu,"cudaMemcpy");
}

void deviceDestruct(float* P_cu, float *B_cu, float* fitness_cu) {
  cudaFree(P_cu);
  cudaFree(B_cu);
  cudaFree(fitness_cu);
}
