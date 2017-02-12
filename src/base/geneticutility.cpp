#include <iostream>
#include <cstdio>
#include <climits>
#include "geneticutility.h"


/* Calculates the number of points that lie inside or outside a given box.
 * All matrices are assumed to be 1D arrays in row-major order.
 *
 *
 * pts: n x d -- Matrix of n points, of d dimensions each. Each point lives in
 *      a row of the given matrix.
 * bxs: M x d -- Matrix of x boxes. Each box has one corner at the origin, the
 *      opposite corner lies at the coordinate given by the d-dim point.
 * inout: M x d -- Matrix of booleans. Tells the program whether to use strict
 *        inequality or not.
 * res: M x n -- Matrix of outcomes. The value in res[i,j] will be nonzero if
 *       pts[j] lies inside bxs[i]
 * n, M, d: unsigned -- The dimensions mentioned above.
 */



// Caller is responsible for initializing memory
void calcPointInsideBoxInfo(const float * const pts, const float * const bxs,
                            const bool * const inout,
                            unsigned * res,
                            unsigned n, unsigned M, unsigned d){

// Apparently the user cannot be trusted to initialize memory correctly.
// This might be removed for raw speed testing.
for(int i = 0; i < M * n; i++){
   res[i] = UINT_MAX;
}

#pragma omp parallel for
for (unsigned int i = 0; i < M; ++i) {
  for (unsigned int j = 0; j < n; ++j) {
    for (int k = 0; k < d; ++k) {
      if (inout[i*d + k]) {
        res[i * n + j] &= pts[j*d + k] < bxs[i*d + k];
      } else {
        res[i * n + j] &= pts[j*d + k] <= bxs[i*d + k];
      }
    }
  }
 }
}


/* Reads given input file and returns a new array containing the points
 * contained in that file. Also modifies the input pointers to give a correct
 * set of dimensions for the user.
 *
 * User is responsible for deleting the returned array with delete[].
 */
float* readInputFile(const char* filename, unsigned& m, unsigned& d, unsigned& n){
  FILE* inf = fopen(filename, "r");

  if (!inf) {
    perror("Could not open input file");
    throw -2; //Not optimal, but since it's rarely called, probably okay?
  }

  fscanf(inf, "%u %u %u", &m, &d, &n);
  float* P = new float[n*d];
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < d; ++j) {
      unsigned int temp;
      fscanf(inf, "%u", &temp);
      P[i*d + j] = temp;
    }
  }
  fclose(inf);

  return P;
}
