#include <iostream>
#include <cstdio>
#include "geneticutility.h"


/* Calculates the number of points that lie inside or outside a given box.
 * All matrices are assumed to be 1D arrays in row-major order.
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
void accel_calcPointInsideBoxInfo(const float * const pts, const float * const bxs,
                            const bool * const inout,
                            unsigned * res,
                            unsigned n, unsigned M, unsigned d){

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
