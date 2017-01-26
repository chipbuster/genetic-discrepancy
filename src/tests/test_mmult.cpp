#include <iostream>
#include <random>
#include <cmath>
#include <cstdlib>

#include "geneticutility.h"
#include "accelmult.h"

using namespace std;

// pts is n * d, bxs is M x d, inout is M * d
// Fill these matrices with randomly generated values
void initMatrices(float* pts, float* bxs, bool* inout,
                 unsigned n, unsigned M, unsigned d){

  std::mt19937 gen;
  // gen.seed(0);

  std::uniform_real_distribution<float> realDist(0.0, 1.0);
  std::uniform_int_distribution<unsigned> boolDist(0,1);

  // Fill pts
  for(unsigned i = 0; i < n * d; i++){
    pts[i] = realDist(gen);
  }
  // Fill pts
  for(unsigned i = 0; i < M * d; i++){
    bxs[i] = realDist(gen);
  }
    // Fill pts
  for(unsigned i = 0; i < M * d; i++){
    inout[i] = boolDist(gen);
  }
}

/* Calculate whether the matrices meaningfully differ, down to a tolerance
 * Both matrices should be m x n in row-major order */
bool matricesDiffer(const unsigned* const matA, const unsigned* const matB,
                    unsigned m, unsigned n){

  for(unsigned i = 0; i < m; i++){
    for(unsigned j = 0; j < n; j++){
      bool diff = matA[i * m + j] != matB[i * m + j];
      if (diff){
        cerr << "[FAIL]: Matrices differ by more than specified amount" << endl;
        cerr << "        at position (" << i << " , " << j << ")" << endl;

        return true;
      }
    }
  }

  return false;
}

int main(int argc, char** argv){
  if (argc != 2){
    cerr << "Usage: " << argv[0] << " numtests" << endl;
    return -1;
  }

  int numTrials = atoi(argv[1]);

  // Matrix dimensions
  unsigned n = 500;
  unsigned m = 500;
  unsigned d = 100;

  // Inputs to the algorithm
  float pts[n * d];
  float bxs[m * d];
  bool inout[m * d];

  // Outputs, one for each version
  unsigned resUnaccel[n * m];
  unsigned resAccel[n * m];

  for(int i = 0; i < numTrials; i++){
    initMatrices(pts, bxs, inout, n, m, d);

    accel_calcPointInsideBoxInfo(pts, bxs, inout, resAccel, n, m, d);
    calcPointInsideBoxInfo(pts, bxs, inout, resUnaccel, n, m, d);

    bool failed = matricesDiffer(resAccel, resUnaccel, m, n);
    if(failed) return 1;

  }

  cout << "[PASS]: Matrix multiplication function" << endl;

  return 0;
}
