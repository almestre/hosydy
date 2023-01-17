// Rng class definition

#ifndef RNG_H
#define RNG_H

#include <random>
#include <cstdint> // uint32_t and uint64_t types
#include <string>

class Rng {
public:
   // Default constructor
   // Rng();


   static void Rng_init(); // use thread_local here in a multithreading context
   static void Rng_save(); // use thread_local here in a multithreading context
   static uint32_t random_uint32();
   static uint64_t random_uint64();
   static double unif_01();

// We don't use thread_local in a distribution if we plan to use different parameter values among calls
   double uniform_real(double,double); // Parameters: min and max
   int uniform_int(int,int); // Parameters: min and max
   bool bernoulli(double); // Parameters: probability of success
   int binomial(int,double); // Parameters: number of trials and probability of success
   int poisson(double); // Parameters: mean
   double normal(double,double); // Parameters: mean and standard deviation
   int negative_binomial(double,double); // Parameters: mean and dispersion parameter*
      // *Alternative parameterization
      // see https://en.wikipedia.org/wiki/Negative_binomial_distribution

// To get a random 64-bit integer with k 1-bits:
   uint64_t random_k_bits(int);

// To get a vector with randomized indexes from 0 to n-1:
   std::vector<int> randIndexVect( int );

private:
   static std::mt19937 Rng_engine; // Random Number Algorithm; use thread_local here in a multithreading context

   // Utility functions
   uint64_t k_bit_helper(int, int, uint64_t, uint64_t);

   };

   #endif // RNG_H
