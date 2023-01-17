// Rng class definition

#include <numeric> // function std::iota
#include<random>
#include<algorithm>
#include <fstream>
#include <array>
#include "Rng.h"
#include <boost/math/distributions/negative_binomial.hpp> // Boost version of negative binomial dist*
#include "Simul.h"
//  *The negative binomial of Boost has a extended definition where the parameter r can take on a positive real value (i.e. Pólya distribution)

using namespace std;

// ---Static data members---
// Static data members should be defined and initialized in the cpp. file
mt19937 Rng::Rng_engine;  // use thread_local here in a multithreading context

// ---Static member functions---
void Rng::Rng_init() {
  std::string replid = std::to_string( Simul::getReplID() );
  ifstream frand("RNG_internal_state/RNGis" + Simul::getScenID() + "_" + replid + ".dat");
   if (frand.is_open())
    {
      frand >> Rng_engine;
    }
   else
   {
//      array<int,624> seed_data; // the state size of mt19937 is 624
      array<int,mt19937::state_size> seed_data; // the state size of mt19937 is 624
      random_device r;
      generate_n(seed_data.data(), seed_data.size(), ref(r));
      seed_seq seq(begin(seed_data), end(seed_data));
      Rng_engine.seed(seq);
   }
   frand.close();
}

void Rng::Rng_save() {
  std::string replid = std::to_string( Simul::getReplID() );
   ofstream frand("RNG_internal_state/RNGis" + Simul::getScenID() + "_" + replid + ".dat");
   if (frand.is_open())
      frand << Rng_engine;
   frand.close();
}

uint32_t Rng::random_uint32() {
   uniform_int_distribution<uint32_t> dist{0, 0xFFFFFFFF}; // 0xFFFFFFFF = UINT32_MAX
   return dist(Rng_engine);
}

uint64_t Rng::random_uint64() {
   uniform_int_distribution<uint64_t> dist{0, 0xFFFFFFFFFFFFFFFF}; // 0xFFFFFFFFFFFFFFFF = UINT64_MAX
   return dist(Rng_engine);
}

double Rng::unif_01() {
uniform_real_distribution<double> dist{0,1};
   return dist(Rng_engine);
}

// ---Non-static member functions---
double Rng::uniform_real(double min, double max) {
   uniform_real_distribution<double> dist{min, max};
   return dist(Rng_engine);
}

int Rng::uniform_int(int  min, int max) {
   uniform_int_distribution<int> dist{min, max};
   return dist(Rng_engine);
}

bool Rng::bernoulli(double prob) {
   bernoulli_distribution dist{prob};
   return dist(Rng_engine);
}

int Rng::binomial(int trials, double prob) {
   binomial_distribution<> dist{trials, prob};
   return dist(Rng_engine);
}

int Rng::poisson(double mean) {
   poisson_distribution<> dist{mean};
   return dist(Rng_engine);
}

double Rng::normal(double mean, double stdev) {
   normal_distribution<> dist{mean, stdev};
   return dist(Rng_engine);
}

int Rng::negative_binomial(double mean, double theta) {
   // Here we use an alternative parameterization based on mean and theta
   // E.g. see https://stat.ethz.ch/R-manual/R-devel/library/stats/html/NegBinomial.html
   double r = theta;
   double p = theta / ( theta + mean );
   // r and p are the standard parameters of a negative binomial
   // The r parameter is real (i.e. Pólya distribution)
   boost::math::negative_binomial_distribution<> dist(r,p); // Uses default `RealType = double`
   // Applying inverse random sampling
   // See https://en.wikipedia.org/wiki/Inverse_transform_sampling
   double randnumber = unif_01();
   return ( boost::math::quantile(dist, randnumber) );
}

// To get a random 64-bit integer with k 1-bits:
uint64_t Rng::random_k_bits(int k) {
  return k_bit_helper(64, k, 1, 0);
}

// To get a vector with randomized indexes from 0 to n-1:
std::vector<int> Rng::randIndexVect( int n ) {
   std::vector<int> vecRandPop(n);
   std::iota (std::begin(vecRandPop), std::end(vecRandPop), 0);
   std::shuffle( begin(vecRandPop), end(vecRandPop), Rng_engine );
   return( vecRandPop );
}

   // Utility functions

uint64_t Rng::k_bit_helper(int n, int k, uint64_t bit, uint64_t accum) {
  if (!(n && static_cast<uint64_t>(k)))
    return accum;
  if (k > uniform_int(0,n-1))
    return k_bit_helper(n - 1, k - 1, bit + bit, accum + bit);
  else
    return k_bit_helper(n - 1, k, bit + bit, accum);
}

