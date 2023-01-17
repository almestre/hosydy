// Definition of Population<Symbiont> class member functions
// Specialized class' member functions are defnied outside the class definition's body

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric> // funstion std::iota
#include <cstdint> // uint32_t and uint64_t types
#include <cstdlib> // exit() function
#include<vector> // C++ standard vector class template
#include<stack> // C++ standard stack class template
#include <string> // C++ standard string class
#include "Population.h" // Organism class definition

using namespace std;

// Population<Symbiont> member function definitions
// Full specialization is not a template: no "template <>" prefix in member function definitions

   // Constructor
Population<Symbiont>::Population(const size_t PopVecSize, int n, uint64_t hid, uint64_t id): Pop(std::vector<Symbiont>(PopVecSize)), N(n), PatchID(hid), ID(id) {} // Here, the patch is a host

   // Non-static member functions

void Population<Symbiont>::newIndFromSource( Population<Host>& hpop, const SourcePatch& continent, Rng& rng ) { // local immigrant (not adapted locally)
   // Create new individual and get its index
   int indexInd = createInd(hpop);
   // Sample a phenotype value from the continent
   double newPhen = rng.normal( continent.getOptPhen(), continent.getStdPhenS() );
   // Calculate number of 1-bits (k) that would produce the newPhen
   double l = static_cast<double>(Symbiont::getL());
   double alpha2 = Symbiont::getAlpha() + Symbiont::getAlpha();
   double kd = l + newPhen/alpha2;
   int k = static_cast<int>(round( kd ));
   // Get a random 64-bit integer with k 1-bits
   uint64_t newGen = rng.random_k_bits(k);
   // Initialize sex, genotype and phenotype of the new individual:
   Pop[indexInd].setSex( ( Rng::unif_01() > 0.5 )? 'f' : 'm' );
   Pop[indexInd].setGen(newGen);
   Pop[indexInd].Phen_init();
}

void Population<Symbiont>::newLocAdIndFromSource( Population<Host>& hpop, const SourcePatch& continent, Rng& rng ) { // local resident (locally adapted)
   // Create new individual and get its index
   int indexInd = createInd(hpop);
   // Get the phenotype value of the host (i.e. the local optimum)
   Host* ptrHost = hpop.getInd ( getPatchID() ); // hpop is the name of the host population
   double LocOptPhen = ptrHost->getPhen();
   ptrHost = nullptr;
   // Sample a locally adapted phenotype value from the continent (based on the local optimum)
   double newPhen = rng.normal( LocOptPhen, continent.getLocStdPhenS() );
   // Calculate number of 1-bits (k) that would produce the newPhen
   double l = static_cast<double>(Symbiont::getL());
   double alpha2 = Symbiont::getAlpha() + Symbiont::getAlpha();
   double kd = l + newPhen/alpha2;
   int k = static_cast<int>(round( kd ));
   // Get a random 64-bit integer with k 1-bits
   uint64_t newGen = rng.random_k_bits(k);
   // Initialize sex, genotype and phenotype of the new individual:
   Pop[indexInd].setSex( ( Rng::unif_01() > 0.5 )? 'f' : 'm' );
   Pop[indexInd].setGen(newGen);
   Pop[indexInd].Phen_init();
}

void Population<Symbiont>::newBorn( Population<Host>& hpop, Rng& rng, Gamete& gamete1, Gamete& gamete2) {
   // Create new individual and get its index
   int indexInd = createInd(hpop);
   // Create a new genotype from the parental gametes
   uint64_t newGen = GenfromGametes( gamete1, gamete2 );
   // Apply mutation
   if ( Symbiont::getmutRate() > 0 ) {
      newGen = genMutation( newGen, rng );
   }
   // Initialize sex, genotype and phenotype of the new individual:
   Pop[indexInd].setSex( ( Rng::unif_01() > 0.5 )? 'f' : 'm' );
   Pop[indexInd].setGen(newGen);
   Pop[indexInd].Phen_init();
}

int Population<Symbiont>::newImmigrant( Population<Host>& hpop) {
   int indexInd = createInd(hpop);
   return(indexInd);
}

void Population<Symbiont>::removeInd( Population<Host>& hpop , int index) {
   if ( index < N-1 ) {
   std::swap(Pop[index],Pop[N-1]);
   --N;
   // Update Nsymbiont of the host harbouring this symbiont population
   Host* ptrHost = hpop.getInd ( getPatchID() ); // hpop is the name of the host population
   --(*ptrHost);
   ptrHost = nullptr;
   }
   else {
   --N;
   // Update Nsymbiont of the host harbouring this symbiont population
   Host* ptrHost = hpop.getInd ( getPatchID() ); // hpop is the name of the host population
   --(*ptrHost);
   }
}

void Population<Symbiont>::setPop(const std::vector<Symbiont>& pop) { Pop = pop; }
const std::vector<Symbiont>& Population<Symbiont>::getPop() const {return Pop;}

Symbiont* Population<Symbiont>::getPtrToInd( int indexInd ) {return &Pop[indexInd];}

void Population<Symbiont>::setN(int n) {N = n;}
int Population<Symbiont>::getN() const {return N;}

void Population<Symbiont>::setPatchID(uint64_t hid) {PatchID = hid;}
uint64_t Population<Symbiont>::getPatchID() const {return PatchID;}

void Population<Symbiont>::setID(uint64_t id) {ID = id;}
uint64_t Population<Symbiont>::getID() const {return ID;}

void Population<Symbiont>::immigrFromSource( const SourcePatch& continent, Rng& rng, Population<Host>& hpop) {
   if ( (continent.getSPrev() > Rng::unif_01() ) || (continent.getSPrev() == 1) ) {  // Determine whether the population is not empty
      // Sample the number of individuals in the population (for a non-empty population)
      int n_ind = rng.negative_binomial( continent.getSAb(), continent.getSTheta() );
      while ( n_ind == 0 ) { // To sample an n_ind > 0
         n_ind = rng.negative_binomial( continent.getSAb(), continent.getSTheta() );
      }
      if ( n_ind > static_cast<int>( getPop().size() ) ) { // To prevent from Pop vector overflow
         n_ind = static_cast<int>( getPop().size() );
      }
      // Estimate the number of locally adapted symbionts
      int n_locad = rng.binomial( n_ind, continent.getPSLocAd() );
      // Create the new source individuals (both globally and locally adapted individuals)
         // Locally adapted individuals
      if ( n_locad > 0 ) {
         for ( int counter = 0; counter < n_locad; ++counter ) {
            newLocAdIndFromSource( hpop, continent, rng );
         }
      }
      if ( (n_ind - n_locad) > 0 ) {
         // Globally adapted individuals
         for ( int counter = 0; counter < (n_ind - n_locad); ++counter ) {
            newIndFromSource( hpop, continent, rng );
         }
      }
   }
}

void Population<Symbiont>::popReproduction ( Population<Host>& hpop, Rng& rng) {
   // Produce the gamete pool:
   gpool gamPool = produceGametePool( hpop, rng );
   // Reset population
   N=0;
   Host* ptrHost = hpop.getInd( getPatchID() );
   ptrHost->setNsymbiont(0);
   ptrHost = nullptr;
   // Produce newborn individuals from the gamete pool:
   vector<int> randgF = rng.randIndexVect( gamPool.first.size() );
   vector<int> randgM = rng.randIndexVect( gamPool.second.size() );
//cout << "\n\nrandgF vector: \n\t";
//for (int i:randgF) {cout<<i<<" ";}
//cout << "\n\nrandgM vector: \n\t";
//for (int i:randgM) {cout<<i<<" ";}
   int NnewBorn = min( gamPool.first.size(), gamPool.second.size() );
//cout << "\n\nNnewBorn = " << NnewBorn;
//cout << "\n\nSize of female gamPool = " << gamPool.first.size();
//cout << "\n\nSize of male gamPool = " << gamPool.second.size();
   if ( NnewBorn > 0 ) {
      for ( int i = 0; i < NnewBorn; ++i ) {
         newBorn( hpop, rng, gamPool.first[randgF[i]], gamPool.second[randgM[i]] );
      }
   }
}

void Population<Symbiont>::printIndividuals () const {
   for (int counter = 0; counter < N; ++counter) {
   Pop[counter].printIndividual();
   }
}

void Population<Symbiont>::printPopulation () const {
   std::cout << PopulationtoString() << std::endl;
}

double Population<Symbiont>::getAvSPhen() const {
   if (N>0) {
      double SumPhen = 0;
      for (int counter = 0; counter < N; ++counter) {
         SumPhen += Pop[counter].getPhen();
      }
      return ( SumPhen / static_cast<double>(N) );
   }
   else { return ( 0 ); }
}

double Population<Symbiont>::getSSPhen() const {
   if (N>0) {
      double mean = getAvSPhen();
      double var = 0;
      for (int counter = 0; counter < N; ++counter) {
         var += pow( ( Pop[counter].getPhen() - mean ), 2) ;
      }
      return ( var );
   }
   else { return ( 0 ); }
}

double Population<Symbiont>::getSumPhen() const {
   if (N>0) {
      double SumPhen = 0;
      for (int counter = 0; counter < N; ++counter) {
         SumPhen += Pop[counter].getPhen();
      }
      return ( SumPhen );
   }
   else { return ( 0 ); }
}

double Population<Symbiont>::getDevSPhen( Population<Host>& hpop ) const {
   if (N>0) {
      // Obtain the average symbiont phenotype within the host
      double avphensymb = getAvSPhen();
      // Obtain the host phenotype
      Host* ptrHost = hpop.getInd( getPatchID() );
      double hostphen = ptrHost->getPhen();
      ptrHost = nullptr;
      // Notice that here we do not have to divide by hostphen because our optimal phenotype for the island is equal to zero
      return ( avphensymb - hostphen );
   }
   else { return ( 0 ); }
}

void Population<Symbiont>::getAlFreq( vector<int>& alfreq ) const {
   if ( N > 0 ) {
      for ( int counter = 0; counter < N; ++counter ) {
         pair<uint32_t, uint32_t> pairgen = pair32Int( Pop[counter].getGen() );
        // Lower "haploid" part of the int representing the genotipe
         sumAlToAlFreq( pairgen.first, alfreq );
         sumAlToAlFreq( pairgen.second, alfreq );
      }
   }
}

int Population<Symbiont>::getSumHetPop () const {
   int sumhetpop = 0;
   for ( int i=0; i < N ; ++i ) {
      sumhetpop += Pop[i].sumHetLocInd();
   }
   return(sumhetpop);
}

   // Utility functions

std::string Population<Symbiont>::PopulationtoString() const {
   std::ostringstream output;
   output << "N = " << getN()
   << "\nHost ID: " << getPatchID()
   << "\t" << "Index: " << static_cast<uint32_t>((getPatchID() << 32) >> 32)
   << "\t" << "Version: " << static_cast<uint32_t>(getPatchID() >> 32)
   << "\nID: " << getID()
   << "\t" << "Index: " << static_cast<uint32_t>((getID() << 32) >> 32)
   << "\t" << "Version: " << static_cast<uint32_t>(getID() >> 32) << "\n\n";
   return output.str();
}

int Population<Symbiont>::createInd( Population<Host>& hpop) {
   if (static_cast<uint32_t>(N) < Pop.size()) {
   int n = N;
   // Update N
   ++N;
   // Update Nsymbiont of the host harbouring this symbiont population
   Host* ptrHost = hpop.getInd ( getPatchID() ); // hpop is the name of the host population
   ++(*ptrHost);
   ptrHost = nullptr;
   // Return the index of the new symbiont
   return(n);
   }
   else {
   std:: cout << "SPop vector is full";
   exit(1);
   }
}

uint64_t Population<Symbiont>::GenfromGametes(const Gamete& gamete1, const Gamete& gamete2) const {
   // (uint64_t) mostSignificantWord << 32 | leastSignificantWord
   uint64_t newGen = (uint64_t) gamete1.getHaplGen() << 32 | gamete2.getHaplGen();
   return(newGen);
}

Population<Symbiont>::gpool Population<Symbiont>::produceGametePool( const Population<Host>& hpop, Rng& rng) const {
   vector<Gamete> FemaleGam;
   vector<Gamete> MaleGam;
   FemaleGam.reserve(Host::getKsymbiont());
   MaleGam.reserve(Host::getKsymbiont());
   for ( int i = 0; i < N; ++i ) {
      vector<Gamete> GamInd = Pop[i].produceGametes( PatchID, hpop, rng );
      if ( Pop[i].getSex() == 'f' ) {
         FemaleGam.insert( FemaleGam.end(), GamInd.begin(),GamInd.end());
//cout << "Added to female pool: " << GamInd.size() << " gametes\n";
//for (const auto& i:GamInd) {i.printHaplGen();}
      }
      else {
         MaleGam.insert( MaleGam.end(), GamInd.begin(),GamInd.end());
//cout << "Added to male pool: " << GamInd.size() << " gametes\n";
//for (const auto& i:GamInd) {i.printHaplGen();}
      }
   }
   return (pair<vector<Gamete>,vector<Gamete>>( FemaleGam, MaleGam ));
}

void Population<Symbiont>::sumAlToAlFreq ( uint32_t haplgen , vector<int>& alfrq ) const {
   const uint32_t SHIFT{8 * sizeof(uint32_t) - 1};
   const uint32_t MASK{static_cast<const uint32_t>(1 << SHIFT)};
   for (uint32_t i{1}; i <= SHIFT + 1; ++i) {
      if (haplgen & MASK) {++alfrq[i-1];}
      haplgen <<= 1; // shift haplgen left by 1
   }
}

pair<uint32_t, uint32_t> Population<Symbiont>::pair32Int (uint64_t value) const {
    return pair<uint32_t, uint32_t>((value << 32) >> 32, value >> 32);
}

uint64_t Population<Symbiont>::genMutation (uint64_t gen, Rng& rng) const {
   // Determine the number of mutations that occurred in the genotype (nmut)
   int nal = Symbiont::getL() + Symbiont::getL(); // nal is the number of alleles per genotype
   int nmut = rng.binomial( nal , Symbiont::getmutRate() );
   if ( nmut > 0) { // if any mutation occurred
      // Get a vector with randomized indexes from 0 to nal-1:
      std::vector<int> rand_al = rng.randIndexVect ( nal );
      // Toggle nmut random bits in the genotype using rand_al
      for ( int i=0; i < nmut ; ++i ) {
         gen ^= 1ULL << rand_al[ i ];
      }
   }
   return (gen);
}
