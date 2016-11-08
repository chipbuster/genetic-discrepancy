#include <algorithm>
#include <cassert>
#include <cmath>
#include <cfloat>
#include <cstdio>
#include <cstring>
#include <ctime> // time(NULL)
#include <cfloat> // max float
#include <iostream>
#include <omp.h> // parallelization
#include <random> // random numbers
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

//#define DEBUG
//#define GEN_DEBUG
INIT_TIMER(cpp);

// Global random generator.
std::mt19937 gen;

// Print out some stuff about the random numbers we generate.
#define RANDOM_CHECK 0

/**
 * Randomly generates points according to Lemma 1.1 of M Shah paper.
 */
void genPoints(float *inp, bool *inout, float *outp, bool *out_inout,
               unsigned int inN, unsigned int outM, unsigned int d) {
  // Set up random number generator
  std::uniform_int_distribution<uint> dis(0, inN - 1);

  if (RANDOM_CHECK) printf("rand genPoints:");
  for (unsigned int i = 0; i < outM; i++) {
    for (unsigned int j = 0; j < d; j++) {
      // Get random x cord
      unsigned int ii = dis(gen);
      if (RANDOM_CHECK) printf(" %u", ii);
      outp[i * d + j] = inp[ii * d + j];
      out_inout[i * d + j] = inout[ii * d + j];
    }
  }
  if (RANDOM_CHECK) printf("\n");
}

/**
 * Calculate fitness, according to Shah genetic algorithm:
 *
 * @param P input points, n x d
 * @param B input potential boxes, m x d
 * @param inout boolean values for whether B[i][j] should be inclusive or
 *          exclusive, m x d
 * @param n dim of P
 * @param M dim of B (and inout)
 * @param d matching dim of P and B
 * @param m the maximum entry of P or B
 * @param fitness output vector of size m that has fitness values for each m
 *          input boxes
 * @param max_D output maxiumum discrepancy value
 * @param max_pos output maximum fitness location
 */
unsigned int calcFitness(const float* P, const float* B, const bool* inout,
                         unsigned int n, unsigned int M, int d, unsigned int m,
                         float* fitness, float *max_D, unsigned int *max_pos) {

  double total_g = 0;
  unsigned int nr = 0;
#pragma omp parallel
{
#pragma omp for reduction(+:total_g)
    for (unsigned int i = 0; i < M; ++i) {
      int SSum = 0;
      for (unsigned int j = 0; j < n; ++j) {
        bool Csub = true;
        for (int k = 0; k < d; ++k) {
          if (inout[i*d + k]) {
            Csub &= P[j*d + k] < B[i*d + k];
          } else {
            Csub &= P[j*d + k] <= B[i*d + k];
          }
        }
        if (Csub) SSum++;
      }
      float volume = 1;
      for (int k = 0; k < d; ++k) {
        volume *= B[i*d + k]/m;
      }
      float disc = std::abs(SSum/n - volume);
      fitness[i] = disc;
      total_g += disc;
    }

#pragma omp barrier
    int local_max = -1;
    float local_max_D = -1;

    //Now, calculate the fitness, according to (iii) of Shah
    nr = 0;
#pragma omp for reduction(+:nr)
    for (unsigned int i = 0; i < M; ++i) {
      // Save the discrepancy before we change it to fitness.
      if (local_max_D < 0 || fitness[i] > local_max_D) {
        local_max_D = fitness[i];
      }
      if (fitness[i] > 1.0/M*total_g) {
        fitness[i] = 1.0 / pow(1.f - fitness[i], 2.0);
        nr++;
      } else {
        fitness[i] = 0.0;
      }
      if (local_max == -1 || fitness[i] > fitness[local_max]) {
        local_max = i;
      }
    }
#pragma omp single
{
    *max_pos = 0;
    *max_D = 0;
} // end omp single
#pragma omp barrier

#pragma omp critical
{
      if (fitness[local_max] > fitness[*max_pos]) {
        *max_pos = local_max;
      }
      if (local_max_D > *max_D) {
        *max_D = local_max_D;
      }
} // end critical

} // end omp parallel region
    return nr;
}

void moveToNew(float *from, bool *from_inout, float *to, bool *to_inout,
               unsigned int fr_idx, unsigned int to_idx, unsigned int d) {
  for (unsigned int i = 0; i < d; ++i) {
    to[to_idx * d + i] = from[fr_idx * d + i];
    to_inout[to_idx * d + i] = from_inout[fr_idx * d + i];
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

int main(int argc, char* argv[]) {
  if (argc < 2) {
    fprintf(stderr, "usage: %s <file> <outfile> [M]\n", argv[0]);
    return -1;
  }

 
  printf("Running with %d threads\n", omp_get_max_threads());

  // Set up random generator
  // Can set the random seed by setting RANDOM_SEED as an environment variable.
  if (const char* env_p = std::getenv("RANDOM_SEED")) {
    unsigned seed = atoi(env_p);
    printf("Using random seed %u\n", seed);
    gen = std::mt19937(seed);
  } else {
    std::random_device rd;
    gen = std::mt19937(rd());
  }

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
  float* P = new float[n*d];
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < d; ++j) {
      unsigned int temp;
      fscanf(inf, "%d", &temp);
      P[i*d + j] = temp;
    }
  }
  fclose(inf);
  STOP_TIMER(cpp);
  fprintf(stderr,"Number of samples: %u, with dimension: %u, and max: %u\n", n, d, m);
  PRINT_TIMER(cpp,"Reading samples in");

  START_TIMER(cpp);
  float* BStar = new float[(2*n+1)*d];
  bool* BStar_inout = new bool[(2*n+1)*d];
  for (unsigned int i = 0; i < n; ++i) {
    for (unsigned int j = 0; j < d; ++j) {
      BStar[2*i*d + j] = P[i*d + j];
      BStar_inout[2*i*d + j] = true;
      // Also need to add the non-inclusive edge.
      BStar[(2*i+1)*d + j] = P[i*d + j];
      BStar_inout[(2*i+1)*d + j] = false;
    }
  }
  for (unsigned int j = 0; j < d; ++j) {
    BStar[2*n*d + j] = m;
  }
  STOP_TIMER(cpp);
  PRINT_TIMER(cpp,"generating BStar");

  // Create a uniform generator from 0-1
  std::uniform_real_distribution<double> unif_01(0.0,1.0);
  // Create uniform sampler for number of dimensions.
  std::uniform_int_distribution<uint> unif_d(0, d - 1);
  // Generate a sample from the base gene pool, which is 2*n+1
  std::uniform_int_distribution<uint> unif_N(0, 2*n);

  // Value of M from Shah, the initial population of the gene pool.
  unsigned int M = NUM_M;
  if (argc == 4) {
    M = atoi(argv[3]);
  }

  //Create B from the input P
  float* B = new float[M*d];
  bool* B_inout = new bool[M*d];
  float* B_next = new float[M*d];
  bool* B_next_inout = new bool[M*d];
  genPoints(BStar, BStar_inout, B, B_inout, 2*n+1, M, d);
  // Initial guess for D* is 0
  float best_DStar = 0, prev_DStar = 0;
  bool keep_going = true;
  unsigned int num_runs = 0;

  float* fitness = new float[M];

  unsigned int loop_no = 0;

  while (keep_going) {
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
    nr = calcFitness(P, B, B_inout, n, M, d, m, fitness, &this_DStar, &max_f);
    STOP_TIMER(cpp);
    PRINT_TIMER(cpp, "calculating fitness");

    fprintf(stderr, "Found max_f at position %u and nr of %u\n", max_f, nr);
//#ifdef GEN_DEBUG
    fprintf(stderr, "Successful! Found nr of %u and discrepancy of %e\n", nr, this_DStar);
//#endif

    if (nr == 0) {
      // Bad sample. Generate more points.
      genPoints(BStar, BStar_inout, B, B_inout, 2*n+1, M, d);
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
    moveToNew(B, B_inout, B_next, B_next_inout, max_f, 0, d);
    //swap(B, fitness, max_f, 0, d);
    ++curr_count;
    // Set the best amount so far.
    best_DStar = std::max(best_DStar, this_DStar);
    // Clone the additional nc - 1 boxes from R (boxes in B that don't have f=0)
    // Clone with uniform probability of 1/nr
    float minprob = 1.0/nr;
    if (RANDOM_CHECK) printf("clone prob:");
    for (unsigned int i = 0; i < M && curr_count < nc; ++i) {
      if (fitness[i] == 0.f) continue;
      float r = unif_01(gen);
      if (RANDOM_CHECK) printf(" %f", r);
      if (r < minprob) {
        moveToNew(B, B_inout, B_next, B_next_inout, i, curr_count, d);
        //swap(B, fitness, i, curr_count, d);
        ++curr_count;
      }
    }
    if (RANDOM_CHECK) printf("\n");

    // (vii) Crossover nx = 67% of (M - nc)
    unsigned int nx = CROSSOVER_PERC * (M - curr_count);
    std::uniform_int_distribution<uint> unif_nr(0, nr - 1);
    if (RANDOM_CHECK) printf("crossover prob:\n");
    // Crossover happens two at a time.
    for (unsigned int i = 0; i < nx; i += 2) {
      // Get the random positions for crossover
      uint r1 = unif_nr(gen);
      uint r2 = unif_nr(gen);
      if (RANDOM_CHECK) printf("  %d,%d: ", r1, r2);
      //fprintf(stderr, "Looking for positions %u and %u (out of %u, M=%u, nr=%u, max=%u)\n",
      //        r1, r2, nx, M, nr, max_pos);
      // Need indices for each of these, but skip over the zero-valued ones.
      uint idx1 = 0, idx2 = 0;
      uint valid_i = 0;
      // Keep going until both are found.
      uint found = 0;
      for (unsigned int j = 0; j < M && found < 2; j++) {
        if (fitness[j] > 0.f) {
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
      moveToNew(B, B_inout, B_next, B_next_inout, idx1, curr_count, d);
      moveToNew(B, B_inout, B_next, B_next_inout, idx2, curr_count+1, d);

      if (RANDOM_CHECK) printf("[");
      // Crossover randomly with probability CROSSOVER_PROB
      for (unsigned int j = 0; j < d; ++j) {
        float tocross = unif_01(gen);
        if (RANDOM_CHECK) printf("%f, ", tocross);
        if (tocross < CROSSOVER_PROB) {
          unsigned int temp = B_next[curr_count*d + j];
          B_next[curr_count*d + j] = B_next[(curr_count+1)*d + j];
          B_next[(curr_count+1)*d + j] = temp;
          bool temp_b = B_next_inout[curr_count*d + j];
          B_next_inout[curr_count*d + j] = B_next_inout[(curr_count+1)*d + j];
          B_next_inout[(curr_count+1)*d + j] = temp_b;
        }
      }
      if (RANDOM_CHECK) printf("]\n");
      curr_count += 2;
    }

    // Now, need to just do mutations on the remaining M - nc - nx boxes.
    if (RANDOM_CHECK) printf("mutations:");
    for (; curr_count < M; ++curr_count) {
      for (unsigned int i = 0; i < d; ++i) {
        float tomut = unif_01(gen);
        if (RANDOM_CHECK) printf(" %f", tomut);
        if (tomut < MUTATION_PROB) {
          tomut = unif_01(gen);
          if (RANDOM_CHECK) printf(",%f", tomut);
          if (tomut > MUT_PROB_ONE) {
            if (RANDOM_CHECK) printf(",NA");
            B_next[curr_count*d + i] = m;
            B_next_inout[curr_count*d + i] = true;
          } else {
            unsigned int base_idx = unif_N(gen);
            if (RANDOM_CHECK) printf(",%u", base_idx);
            printf("%u ", base_idx);
            B_next[curr_count*d + i] = BStar[base_idx*d + i];
            B_next_inout[curr_count*d + i] = BStar_inout[base_idx*d + i];
          }
        } else {
          if (RANDOM_CHECK) printf(",NA,NA");
        }
      }
    }
    if (RANDOM_CHECK) printf("\n");
    printf("\n");

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
    ss << argv[2] << "_" << loop_no;
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

    memcpy(B, B_next, M*d*sizeof(unsigned int));
    memcpy(B_inout, B_next_inout, M*d*sizeof(bool));
    STOP_TIMER(cpp);
    PRINT_TIMER(cpp, "total time for mutations");
  }

  fprintf(stderr,"Star discrepancy after %u generations is %f\n", loop_no, best_DStar);

  // Save to a file, if provided
  if (argc > 2) {
    FILE* ofs = fopen(argv[2], "w");
    fprintf(ofs, "%u %u %u\n", m, d, M);
    for (unsigned int i = 0; i < M; i++) {
      for (unsigned int j = 0; j < d; j++) {
        fprintf(ofs, "%f ", B_next[i * d + j]);
      }
      fprintf(ofs, "\n");
    }
    fclose(ofs);
  }

  return 0;
}
