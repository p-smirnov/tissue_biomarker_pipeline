/*
// Fast permutations for rCI using a naive matrix based approach.
*/

#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <R.h>
#include <Rinternals.h>

// #include "xoroshiro128+.h"


/* This is xoroshiro128+ 1.0, our best and fastest small-state generator
   for floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than
   xoroshiro128**. It passes all tests we are aware of except for the four
   lower bits, which might fail linearity tests (and just those), so if
   low linear complexity is not considered an issue (as it is usually the
   case) it can be used to generate 64-bit outputs, too; moreover, this
   generator has a very mild Hamming-weight dependency making our test
   (http://prng.di.unimi.it/hwd.php) fail after 5 TB of output; we believe
   this slight bias cannot affect any application. If you are concerned,
   use xoroshiro128** or xoshiro256+.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. If you have
   a 64-bit seed, we suggest to seed a splitmix64 generator and use its
   output to fill s. 

   NOTE: the parameters (a=24, b=16, b=37) of this version give slightly
   better results in our test than the 2016 version (a=55, b=14, c=36).
*/

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}


// static uint64_t s[2];

uint64_t next(uint64_t *s) {
  const uint64_t s0 = s[0];
  uint64_t s1 = s[1];
  const uint64_t result = s0 + s1;

  s1 ^= s0;
  s[0] = rotl(s0, 24) ^ s1 ^ (s1 << 16); // a, b
  s[1] = rotl(s1, 37); // c

  return result;
}

// Code to trim random numbers from :https://stackoverflow.com/questions/822323/how-to-generate-a-random-int-in-c
// maxindex here is exclusive of the right edge
uint64_t generate_random_index(uint64_t *state, uint64_t maxindex){
   // uint64_t maxindex = llround(maxindexd);
   
   if ((maxindex-1) == UINT64_MAX){
      return next(state);
   } else {
    // Supporting larger values for n would requires an even more
    // elaborate implementation that combines multiple calls to rand()
    // assert (maxindex <= UINT64_MAX)

    // Chop off all of the values that would cause skew...
    uint64_t end = UINT64_MAX / (maxindex); // truncate skew
    // assert (end > 0);
    end *= maxindex;

    // ... and ignore results from rand() that fall above that limit.
    // (Worst case the loop condition should succeed 50% of the time,
    // so we can expect to bail out of this loop pretty quickly.)
    uint64_t r;
    while ((r = next(state)) >= end);

    return r % maxindex;
  }
}



void printVec(uint64_t *list, uint64_t N){
  for(uint64_t i = 0; i < N; i ++){
    printf("%lld, ", list[i]);
  }
}

// Using "inside out" Fisher Yates https://en.wikipedia.org/wiki/Fisher–Yates_shuffle
void sampleIdx(uint64_t N, uint64_t *permPointer, uint64_t *state){
  uint64_t j;
  permPointer[0] = 0;
  // printf("%lld", N);

  for(uint64_t i = 0; i <= N-1; i++){
      j = generate_random_index(state, i+1);
      if(j != i){
        permPointer[i] = permPointer[j];
      }
      permPointer[j] = i;
  }

}
  


double runPerm(double *out, double *xvec, double *yvec, double obsCor, int *dsOneHot, uint64_t R, uint64_t N, int *DSsize, int numDS, uint64_t *state){

  double pval;
  double currCor;
  uint64_t totalSeenLarger = 0;
  uint64_t i = 0;
  // uint64_t j = 0;
  uint64_t *permIdx = malloc(N * sizeof(uint64_t));

  double *dsDenom = malloc(numDS*sizeof(double));
  double *dsSum = malloc(numDS*sizeof(double));
  double *yShuffled = malloc(N * sizeof(double));




  while(i < R){

    sampleIdx(N, permIdx, state);
    // printVec(permIdx, N);
    // printf("\n");
    
    for(int j = 0; j < numDS; j++){
      dsDenom[j] = 0;
      dsSum[j] = 0;
    }

    for(int j = 0; j < N; j++){
      yShuffled[j] = yvec[permIdx[j]];
        // yShuffled[j] = yvec[j];
    }

    for(int j = 0; j < N; j++){

      dsSum[dsOneHot[j]] += yShuffled[j];
    }

    // Think about optimizing the left hand side here for consecutive access?
    // Maybe counting up until we switch datasets is more efficent? (maybe not, since that creates branching)
    for(int j = 0; j < N; j++){
      dsDenom[dsOneHot[j]] += pow((DSsize[dsOneHot[j]] * yShuffled[j] - dsSum[dsOneHot[j]]), 2);
    }

    currCor = 0;

    for(int j = 0; j < N; j++){
      currCor += xvec[j] * sqrt(DSsize[dsOneHot[j]] - 1) * (DSsize[dsOneHot[j]] * yShuffled[j] - dsSum[dsOneHot[j]])/sqrt(dsDenom[dsOneHot[j]]);
    }
    // printf("%f\n", currCor);

    if(fabs(currCor) >= fabs(obsCor)){
      totalSeenLarger++;
    }    

    out[i] = currCor;
    // out[i] = permIdx[0];
    i++;
  }

  // TODO:: Figure out how the 0 based indexing affects calculation of the pvalues
  // printf("%ld\n", totalSeenLarger);

  free(permIdx);
  free(dsDenom);
  free(dsSum);
  free(yShuffled);

  if(totalSeenLarger==0){
    pval = 1/(((double)R )+ 1);
  } else {
    pval = (((double)totalSeenLarger+1))/(((double)R+1));
  }
  
  return(pval);
  // return(currCor);

}


SEXP metaPermC(SEXP pin_x,
               SEXP pin_y,
               SEXP pobsCor,
               SEXP pdsOneHot,
               SEXP pDSsize,
               SEXP pnumDS,
               SEXP pR,
               SEXP pn,
               SEXP pseed){
  
  double Ndouble = *REAL(pn);

  double Rdouble = *REAL(pR);  
  double obsCor = *REAL(pobsCor);

  double temp;

  uint64_t N = (uint64_t) Ndouble;
  uint64_t R = (uint64_t) Rdouble;

  SEXP pout = PROTECT(allocVector(REALSXP,R));
  
  double *out = REAL(pout);
  
  double *seed = REAL(pseed);
  uint64_t *state = (uint64_t*) seed;

  int numDS = *INTEGER(pnumDS);

  temp = runPerm(out, REAL(pin_x), REAL(pin_y), obsCor, INTEGER(pdsOneHot), R, N, INTEGER(pDSsize), numDS, state);

  UNPROTECT(1);
  
  return pout;
  
}


