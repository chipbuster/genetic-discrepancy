#include <algorithm>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cfloat>
#include <iostream>
#include <omp.h>
#include <random>
#include <sstream>

#include "timer.h"

#define CLONE_PERC      0.5     // percent of r boxes to clone (nc)
#define CROSSOVER_PERC  0.6667  // percent of remaining (M - nc) boxes to clone
#define CROSSOVER_PROB  0.5     // probability to crossover a given dimension
#define MUTATION_PROB   0.80    // probability to mutate a given dimension
#define STOPPING_SAME   0.0001
#define MAX_RUNS        50
#define MUT_PROB_ONE    0.8     // probability of mutating to a 1

#define NUM_M           1000
#define MIN_RUNS        100

// Used for mic pragmas
#define ALLOC alloc_if(1) free_if(0)
#define FREE alloc_if(0) free_if(1)
#define REUSE alloc_if(0) free_if(0)

// Comment this line to copy every time.
#define KEEP_POINTERS

//#define DEBUG
//#define GEN_DEBUG
INIT_TIMER(cpp);

// global pointers for mic
__declspec(target(mic))
  float* g_P; // Pointer to original points
__declspec(target(mic))
  float* g_B; // Pointer to points used it fitness function
__declspec(target(mic))
  float* g_fitness; // Pointer to fitness scoring function

/**
 * Randomly generates points according to Lemma 1.1 of M Shah paper.
 */
void genPoints(float *inp, float *outp,
               unsigned int inN, unsigned int outM, unsigned int d) {
  // Set up random number generator
  std::random_device rd;
  std::uniform_int_distribution<uint> dis(0, inN - 1);
  // Not sure what this generator is, but it seems to be the default
  std::mt19937 gen(rd());

  for (unsigned int i = 0; i < outM; i++) {
    for (unsigned int j = 0; j < d; j++) {
      // Get random x cord
      unsigned int ii = dis(rd);
      outp[i * d + j] = inp[ii * d + j];
    }
  }
}

unsigned int calcFitness(unsigned int n, unsigned int M, int d, unsigned int m,
                         float *max_D, unsigned int *max_pos) {

  double total_g = 0;
  int print_count = 0;
#ifdef KEEP_POINTERS
#pragma offload target(mic) \
  nocopy(g_P:length(n*d) REUSE ) \
  in(g_B:length(M*d) REUSE ) \
  in(d) in(m) \
  out(g_fitness:length(M) REUSE )
#else
#pragma offload target(mic) \
  in(g_P:length(n*d)) \
  in(g_B:length(M*d)) \
  in(d) in(m) \
  out(g_fitness:length(M))
#endif
#pragma omp parallel for reduction(+:total_g)
    for (unsigned int i = 0; i < M; ++i) {
      int SSum = 0;
      for (unsigned int j = 0; j < n; ++j) {
        bool Csub = true;
        for (int k = 0; k < d; ++k) {
          Csub &= g_P[j*d + k] <= g_B[i*d + k];
        }
        if (Csub) SSum++;
      }
      float volume = 1;
      for (int k = 0; k < d; ++k) {
        volume *= g_B[i*d + k]/m;
      }
      float disc = fabs(SSum/n - volume);
#ifdef DEBUG
      if ( (SSum > 0 && print_count++ < 5) || i < 3) {
        fprintf(stderr,"%d: %e, %d: %f\n", i, volume, SSum, disc);
      }
#endif
      g_fitness[i] = disc;
    }

    //Now, calculate the fitness, according to (iii) of Shah
    *max_pos = 0;
    unsigned int nr = 0;
    for (unsigned int i = 0; i < M; ++i) {
      if (g_fitness[i] > *max_D) {
        *max_D = g_fitness[i];
      }
      if (g_fitness[i] > 1.0/M*total_g) {
        g_fitness[i] = 1.0 / pow(1 - g_fitness[i], 2.0);
        nr++;
      } else {
        g_fitness[i] = 0.0;
      }
      if (g_fitness[i] > g_fitness[*max_pos]) {
        *max_pos = i;
      }
    }

    return nr;
}

void moveToNew(float *from, float *to,
               unsigned int fr_idx, unsigned int to_idx, unsigned int d) {
  for (unsigned int i = 0; i < d; ++i) {
    to[to_idx * d + i] = from[fr_idx * d + i];
  }
}

/**
 * Need to swap both the values in S and the values in fitness
 */
void swap(float *S, float *f, //unsigned int *indices,
          unsigned int f_idx, unsigned int t_idx, unsigned int d) {
  // Swap the indices
  //indices[t_idx] = f_idx;
  //indices[f_idx] = t_idx;
  // Swap the fitness values
  float temp_f = f[t_idx];
  f[t_idx] = f[f_idx];
  f[f_idx] = temp_f;
  for (unsigned int i = 0; i < d; i++) {
    float temp = S[t_idx * d + i];
    S[t_idx * d + i] = S[f_idx * d + i];
    S[f_idx * d + i] = temp;
  }
}

void setupMic(unsigned int n, unsigned int M, int d) {
#pragma offload_transfer target(mic) \
  in(g_P:length(n*d) ALLOC ) \
  nocopy(g_B:length(M*d) ALLOC ) \
  nocopy(g_fitness:length(M) ALLOC )
}

void cleanMic(unsigned int n, unsigned int M, int d) {
#pragma offload_transfer target(mic) \
  nocopy(g_P:length(n*d) FREE ) \
  nocopy(g_B:length(M*d) FREE ) \
  nocopy(g_fitness:length(M) FREE )
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s <file> <nthreads> <outfile> [M]\n", argv[0]);
    return -1;
  }
 
  int nthreads = atoi(argv[2]);
	omp_set_dynamic(0);
	omp_set_num_threads(nthreads);
	printf("Running with %d threads\n", omp_get_max_threads());

  START_TIMER(cpp);
  FILE* inf = fopen(argv[1], "r");
  //FILE* inf = fopen("C:/Users/Nathan Clement/Documents/gpu/gpu_code/samples.m1k.d128.n10k.uniform", "r");
  //FILE* inf = fopen("C:/Users/Nathan Clement/Documents/gpu/gpu_code/samples.m1k.d32.n1m.uniform", "r");
  //FILE* inf = fopen("C:/Users/Nathan Clement/Documents/gpu/gpu_code/samples.m1000.d10.n10k.prgn", "r");

  if (!inf) {
    perror("Could not open input file");
    return -2;
  }
  // All values read from the input file.
  unsigned int m, d, n;
  fscanf(inf, "%d %d %d", &m, &d, &n);
  // Use the global P pointer
  g_P = new float[n*d];

  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < d; ++j) {
      unsigned int temp;
      fscanf(inf, "%d", &temp);
      g_P[i*d + j] = temp;
    }
  }
  fclose(inf);
  STOP_TIMER(cpp);
  fprintf(stderr,"Number of samples: %u, with dimension: %u, and max: %u\n", n, d, m);
  PRINT_TIMER(cpp,"Reading samples in");

  START_TIMER(cpp);
  float* BStar = new float[(2*n+1)*d];
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < d; ++j) {
      BStar[2*i*d + j] = g_P[i*d + j];
      // Also need to add the non-inclusive edge.
      BStar[(2*i+1)*d + j] = g_P[i*d + j] - FLT_EPSILON;
    }
  }
  for (unsigned int j = 0; j < d; ++j) {
    BStar[2*n*d + j] = m;
  }
  STOP_TIMER(cpp);
  PRINT_TIMER(cpp,"generating BStar");

  // Create a uniform generator from 0-1
  std::uniform_real_distribution<double> unif_01(0.0,1.0);
  std::default_random_engine gen;
  // Create uniform sampler for number of dimensions.
  std::uniform_int_distribution<uint> unif_d(0, d - 1);
  // Generate a sample from the base gene pool, which is 2*n+1
  std::uniform_int_distribution<uint> unif_N(0, 2*n);

  // Value of M from Shah, the initial population of the gene pool.
  unsigned int M = NUM_M;
  if (argc == 5) {
    M = atoi(argv[4]);
  }

  //Create B from the input P
  // use the global B pointer
  g_B = new float[M*d];
  float* B_next = new float[M*d];
  genPoints(BStar, g_B, 2*n+1, M, d);
  // Initial guess for D* is 0
  float best_DStar = 0, prev_DStar = 0;
  bool keep_going = true;
  unsigned int num_runs = 0;

  // Use the global fitness pointer
  g_fitness = new float[M];

  unsigned int loop_no = 0;

#ifdef KEEP_POINTERS
  // Need to setup the mic
  setupMic(n, M, d);
#endif

  while (keep_going || 1) {
    //unsigned int* indices = new unsigned int[M];
    //float ds = DStarSingle(P, P_cu, n, d, m);
    //float ds = DStarBoth(P, P_cu, P, P_cu, n, n, d, m);
#ifdef GEN_DEBUG
    fprintf(stderr, "Running DStarFitness with n=%u, M=%u, d=%u\n", n, M, d);
#endif
    //DStarFitness(P, P_cu, B, B_cu, n, M, d, m, fitness, fitness_cu, this_DStar);

    // Do the calculations!!
    float this_DStar = 0;
    unsigned int nr = 0;
    unsigned int max_f = 0;
    START_TIMER(cpp);
    nr = calcFitness(n, M, d, m, &this_DStar, &max_f);
    STOP_TIMER(cpp);
    PRINT_TIMER(cpp, "calculating fitness");

    fprintf(stderr, "Found max_f at position %u and nr of %u\n", max_f, nr);
//#ifdef GEN_DEBUG
    fprintf(stderr, "Successful! Found nr of %u and discrepancy of %e\n", nr, this_DStar);
//#endif

    if (nr == 0) {
      // Bad sample. Generate more points.
      genPoints(BStar, g_B, 2*n+1, M, d);
#ifdef GEN_DEBUG
      fprintf(stderr,"After %u consecutive (%u total) runs, the predicted star discrepancy is %f (%f%% change)\n",
             num_runs, loop_no++, best_DStar, 0);
#endif
      continue;
    }

    /*
    // Let's sort this array
    for (int i = 0; i < M; i++) {
    indices[i] = i;
    }
    std::sort(indices, indices+M, sort_fitness(fitness));
    */

    START_TIMER(cpp);
    unsigned int curr_count = 0;
    // Step (vi) - clone nc boxes from r
    unsigned int nc = CLONE_PERC*nr;
    // Clone max with Pr 1
    moveToNew(g_B, B_next, max_f, 0, d);
    //swap(B, fitness, max_f, 0, d);
    ++curr_count;
    // Set the best amount so far.
    best_DStar = std::max(best_DStar, this_DStar);
    // Clone the additional nc - 1 boxes from R (boxes in B that don't have f=0)
    // Clone with uniform probability of 1/nr
    float minprob = 1.0/nr;
    for (unsigned int i = 0; i < M && curr_count < nc; ++i) {
      if (g_fitness[i] == 0.f) continue;
      float r = unif_01(gen);
      if (r < minprob) {
        moveToNew(g_B, B_next, i, curr_count, d);
        //swap(B, fitness, i, curr_count, d);
        ++curr_count;
      }
    }
    // (vii) Crossover nx = 67% of (M - nc)
    unsigned int nx = CROSSOVER_PERC * (M - curr_count);
    std::uniform_int_distribution<uint> unif_nr(0, nr - 1);
    // Crossover happens two at a time.
    for (unsigned int i = 0; i < nx; i += 2) {
      // Get the random positions for crossover
      uint r1 = unif_nr(gen);
      uint r2 = unif_nr(gen);
      //fprintf(stderr, "Looking for positions %u and %u (out of %u, M=%u, nr=%u, max=%u)\n",
      //        r1, r2, nx, M, nr, max_pos);
      // Need indices for each of these, but skip over the zero-valued ones.
      uint idx1 = 0, idx2 = 0;
      uint valid_i = 0;
      // Keep going until both are found.
      uint found = 0;
      for (unsigned int j = 0; j < M && found < 2; j++) {
        if (g_fitness[j] > 0.f) {
          if (valid_i == r1) {
            idx1 = j;
            found++;
          }
          if (valid_i == r2) {
            idx2 = j;
            found++;
          }
          valid_i++;
        }
      }
      //fprintf(stderr, "[%u=%u/%u] Found positions %u=%u and %u=%u (found=%u, valid_i=%u)\n",
      //        curr_count, i,nx, r1,idx1, r2,idx2, found, valid_i);
      assert(found == 2 && "Could not find both values!!");

      // Add these two to the new list, but do crossover with them first.
      //swap(B, fitness, curr_count, idx1, d);
      //swap(B, fitness, curr_count+1, idx2, d);
      moveToNew(g_B, B_next, idx1, curr_count, d);
      moveToNew(g_B, B_next, idx2, curr_count+1, d);

      // Crossover randomly with probability CROSSOVER_PROB
      for (unsigned int j = 0; j < d; ++j) {
        float tocross = unif_01(gen);
        if (tocross < CROSSOVER_PROB) {
          unsigned int temp = B_next[curr_count*d + j];
          B_next[curr_count*d + j] = B_next[(curr_count+1)*d + j];
          B_next[(curr_count+1)*d + j] = temp;
        }
      }
      curr_count += 2;
    }

    // Now, need to just do mutations on the remaining M - nc - nx boxes.
    for (; curr_count < M; ++curr_count) {
      for (unsigned int i = 0; i < d; ++i) {
        float tomut = unif_01(gen);
        if (tomut < MUTATION_PROB) {
          tomut = unif_01(gen);
          if (tomut > MUT_PROB_ONE) {
            B_next[curr_count*d + i] = m;
          } else {
            unsigned int base_idx = unif_N(gen);
            B_next[curr_count*d + i] = BStar[base_idx*d + i];
          }
        }
      }
    }

    float perc_change = (best_DStar - prev_DStar) / best_DStar;
    if ( perc_change < STOPPING_SAME) {
      if (num_runs++ > MAX_RUNS && loop_no > MIN_RUNS) {
        keep_going = false;
      } else {
        keep_going = true;
      }
    } else {
      num_runs = 0;
    }

    fprintf(stderr,"After %u consecutive (%u total) runs, the predicted star discrepancy is %f (%e) (%f%% change)\n",
           num_runs, loop_no++, best_DStar, best_DStar, perc_change*100);
    prev_DStar = best_DStar;

#ifdef DEBUG
    std::stringstream ss;
    ss << argv[3] << "_" << loop_no;
    // Save to a file, if provided
    FILE* ofs = fopen(ss.str().c_str(), "w");
    fprintf(ofs, "%u %u %u\n", m, d, M);
    for (int i = 0; i < M; i++) {
      for (int j = 0; j < d; j++) {
        fprintf(ofs, "%u ", B_next[i * d + j]);
      }
      fprintf(ofs, "\n");
    }
    fclose(ofs);
    fprintf(stderr, "-- Just wrote to %s\n", ss.str().c_str());
#endif
    
    memcpy(g_B, B_next, M*d*sizeof(unsigned int));
    STOP_TIMER(cpp);
    PRINT_TIMER(cpp, "total time for mutations");
  }

  fprintf(stderr,"Star discrepancy after %u generations is %f\n", loop_no, best_DStar);

  // Save to a file, if provided
  if (argc > 2) {
    FILE* ofs = fopen(argv[3], "w");
    fprintf(ofs, "%u %u %u\n", m, d, M);
    for (unsigned int i = 0; i < M; i++) {
      for (unsigned int j = 0; j < d; j++) {
        fprintf(ofs, "%f ", B_next[i * d + j]);
      }
      fprintf(ofs, "\n");
    }
    fclose(ofs);
  }
  
#ifdef KEEP_POINTERS
  // Need to destruct the mic
  cleanMic(n, M, d);
#endif

  return 0;
}
