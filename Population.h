// Population template class definition
// Class template member functions are defined within the class definition's body (no interface-implementation separation)

#ifndef POPULATION_H
#define POPULATION_H

#include <sstream>
#include <iostream>
#include <numeric> // function std::iota
#include <cmath>
#include <cstdint> // uint32_t and uint64_t types
#include <cstdlib> // exit() function
#include<vector> // C++ standard vector class template
#include<stack> // C++ standard stack class template
#include <string> // C++ standard string class
#include "Host.h" // Organism class definition
#include "Symbiont.h" // Organism class definition
#include "SourcePatch.h"
#include "Gamete.h"
#include "Param.h"

// Forward declarations:
template<typename T>
class Metapopulation;

/**** Primary template (member definitions in the header) ****/
template<typename T>
class Population {

public:
   // Within-class type definitions:
   typedef std::pair<std::vector<HostGamete>,std::vector<HostGamete>> hgpool; // gamete pool

   explicit Population (const size_t = Param::getHPopVecSize(), const std::stack<uint32_t>& = std::stack<uint32_t>(), int = 0, uint64_t = 0, uint64_t = 0);

   uint64_t newIndFromSource( const SourcePatch&, Rng& );
   uint64_t newBorn( Rng&, Gamete&, Gamete& );
   T* getInd(uint64_t);
   const T* getConstInd(uint64_t) const; // variant of getInd for access to individual info without modifying it
   void removeInd(uint64_t);

   void setPop(const std::vector<T>&);
   const std::vector<T>& getPop() const;

   void setindList(const std::vector<int>&);
   std::vector<int> getindList() const;
   void initindList();

   void setfreeList(const std::stack<uint32_t>&);
   std::stack<uint32_t> getfreeList() const;
   void initfreeList();

   void setN(int);
   int getN() const;

   void setPatchID(uint64_t);
   uint64_t getPatchID() const;

   void setID(uint64_t);
   uint64_t getID() const;

   void initPop( Patch& patch );

   void initImmigrFromSource( const SourcePatch&, Rng&, Metapopulation<Population<Symbiont>>& );
   void immigrFromSource( const SourcePatch&, Rng&, Metapopulation<Population<Symbiont>>& );
   void pulseImmigrFromSource( const SourcePatch&, Rng&, Metapopulation<Population<Symbiont>>& );

   void popReproduction ( const Patch&, Rng&, Metapopulation<Population<Symbiont>>& );
   void popMortality ( Rng&, Metapopulation<Population<Symbiont>>& );

   void printIndividuals() const;
   void printiList() const;
   void printPopulation() const;

   double getSPrev() const; // Get symbiont prevalence
   double getAvSAbOcH() const; // Get average symbiont abundance per occupied host
   double getAvSAb() const; // Get average symbiont abundance (including empty hosts)
   double getVarSAb() const; // Get variance of symbiont abundances (including empty hosts)
   double getStdSAb() const; // Get standard deviation of symbiont abundances (including empty hosts)
   double getVmrSAb() const; // Get the variance-to-mean ratio of symbiont abundances
   double getCvSAb() const; // Get the coefficient of variation of symbiont abundances
   double getAvHPhen() const; // Get the average host phenotype
   double getStdHPhen() const; // Get the average symbiont phenotype
   std::vector<int> getAlFreq() const; // Get a vector with the frequencies of allele "1" (vector size equals the number of loci). Here, frequency is the total number of alleles "1" present in the population for a given locus (considering the two alleles per locus belonging to each individual in the population).
   double getHexp () const; // Get the expected heterozigosity of the host population based on Hardy-Weinberg equilibrium
   double getHobs () const; // Get the observed heterozigosity of the host population


private:
   std::vector<T> Pop; // Vector of individuals representing the population
   std::vector<int> indList; // Indirection list to access individuals by id
   std::stack<uint32_t> freeList; // List of free id to use when new individuals are created
   int N; // Number of individuals in the population
   uint64_t PatchID; // id of the patch inhabited by the population (in case of a symbiont population the patch is a host individual)
   uint64_t ID; // Population id

   // Utility functions
   std::string PopulationtoString() const;
   uint64_t createInd();
   uint64_t GenfromGametes(const Gamete&, const Gamete&) const;
   hgpool produceGametePool( const Patch& , Rng& ) const;
   void sumAlToAlFreq ( uint32_t , std::vector<int>& ) const;
   std::pair<uint32_t, uint32_t> pair32Int (uint64_t) const;
   int sumHetPop () const;
   uint64_t genMutation (uint64_t, Rng&) const;
};

/**** Population class template specialization for Symbiont (member definitions in cpp file) ****/

template<>  // empty angle brackets
class Population<Symbiont> { // template name with arguments in angle brackets

// Within-class type definitions:
typedef std::pair<std::vector<Gamete>,std::vector<Gamete>> gpool; // gamete pool

public:
   explicit Population (const size_t = Param::getSPopVecSize(), int = 0, uint64_t = 0, uint64_t = 0);

   void newIndFromSource( Population<Host>&, const SourcePatch&, Rng& );
   void newLocAdIndFromSource( Population<Host>&, const SourcePatch&, Rng& );
   void newBorn( Population<Host>&, Rng&, Gamete& gamete1, Gamete& gamete2);
   int newImmigrant( Population<Host>& );
   void removeInd(Population<Host>&, int);

   void setPop(const std::vector<Symbiont>&);
   const std::vector<Symbiont>& getPop() const;

   Symbiont* getPtrToInd( int );

   void setN(int);
   int getN() const;

   void setPatchID(uint64_t); // Here, the patch is a host
   uint64_t getPatchID() const; // Here, the patch is a host

   void setID(uint64_t);
   uint64_t getID() const;

   void immigrFromSource( const SourcePatch&, Rng&, Population<Host>& );

   void popReproduction ( Population<Host>&, Rng& );

   void printIndividuals() const;
   void printPopulation() const;

   double getAvSPhen() const; // Get the average symbiont phenotype
   double getSSPhen() const; // Get the average symbiont phenotype
   double getSumPhen() const; // Get the sum of all symbiont phenotypes
   double getDevSPhen( Population<Host>& ) const; // Get the average of the deviation of the mean symbiont phenotype from the host phenotype, relative to the host phenotype (to calculate EvC for symbionts)
   void getAlFreq( std::vector<int>& ) const;
   int getSumHetPop () const; // get the sum of heterozigotic loci in the population

private:
   std::vector<Symbiont> Pop; // Vector of individuals representing the population
   int N; // Number of individuals in the population
   uint64_t PatchID; // id of the host inhabited by the symbiont population
   uint64_t ID; // id of the symbiont population

   // Utility functions
   std::string PopulationtoString() const;
   int createInd( Population<Host>& );
   uint64_t GenfromGametes(const Gamete&, const Gamete&) const;
   gpool produceGametePool( const Population<Host>&, Rng&) const;
   void sumAlToAlFreq ( uint32_t , std::vector<int>& ) const;
   std::pair<uint32_t, uint32_t> pair32Int (uint64_t) const;
   uint64_t genMutation (uint64_t, Rng&) const;
};


// Including Metapopulation class template full definition after the Population class definitions, and before their implementation
#include "Metapopulation.h"

// Primary template member function definitions (same header as the class template definition)

   // Constructor
template<typename T>
   Population<T>::Population(const size_t PopVecSize, const std::stack<uint32_t>& frl, int n, uint64_t pid, uint64_t id): Pop(std::vector<T>(PopVecSize)), indList(std::vector<int>(PopVecSize)), freeList(frl), N(n), PatchID(pid), ID(id) {}

   // Non-static member functions

template<typename T>
   uint64_t Population<T>::newIndFromSource( const SourcePatch& continent, Rng& rng ) {
   // Create new individual and get its id
   uint64_t IDInd = createInd();
   T* ind = getInd(IDInd);
   // Sample a phenotype value from the continent
   double newPhen = rng.normal( continent.getOptPhen(), continent.getStdPhenH() );
   // Calculate number of 1-bits (k) that would produce the newPhen
   double l = static_cast<double>(T::getL());
   double alpha2 = T::getAlpha() + T::getAlpha();
   double kd = l + newPhen/alpha2;
   int k = static_cast<int>(round( kd ));
   // Get a random 64-bit integer with k 1-bits
   uint64_t newGen = rng.random_k_bits(k);
   // Initialize sex, genotype and phenotype of the new individual:
   ind->setSex( ( Rng::unif_01() > 0.5 )? 'f' : 'm' );
   ind->setGen(newGen);
   ind->Phen_init();
   ind = nullptr;
   return IDInd;
}

template<typename T>
   uint64_t Population<T>::newBorn( Rng& rng, Gamete& gamete1, Gamete& gamete2) {
   // Create new individual and get its ID
   uint64_t IDInd = createInd();
   T* ind = getInd(IDInd);
   // Create a new genotype from the parental gametes
   uint64_t newGen = GenfromGametes( gamete1, gamete2 );
   // Apply mutation
   if ( T::getmutRate() > 0 ) {
      newGen = genMutation( newGen, rng );
   }
   // Initialize sex, genotype and phenotype of the new individual:
   ind->setSex( ( Rng::unif_01() > 0.5 )? 'f' : 'm' );
   ind->setGen(newGen);
   ind->Phen_init();
   // Return the new born's ID:
   return IDInd;
}

template<typename T>
   T* Population<T>::getInd(uint64_t id) {
   uint32_t index = indList[ (id & 0xFFFFFFFF) ]; // 'id & 0xFFFFFFFF' is the lower 32-bit part of the ID (i.e. the index)
   if (Pop[index].getID() == id) {
      return &Pop[index];
      }
   else {return nullptr;}
}

template<typename T>
   const T* Population<T>::getConstInd(uint64_t id) const {
   uint32_t index = indList[ (id & 0xFFFFFFFF) ]; // 'id & 0xFFFFFFFF' is the lower 32-bit part of the ID (i.e. the index)
   if (Pop[index].getID() == id) {
      return &Pop[index];
      }
   else {return nullptr;}
}

template<typename T>
   void Population<T>::removeInd(uint64_t id ) {
//std::cout << "\n\n";
//printiList();
   uint32_t index = indList[ (id & 0xFFFFFFFF) ]; // 'id & 0xFFFFFFFF' is the lower 32-bit part of the ID (i.e. the index)
//std::cout << "\nindex part of ID= " << (id & 0xFFFFFFFF);
//std::cout << "\nIndex obtained from the indlist = " << indList[ (id & 0xFFFFFFFF) ];
//std::cout << "\nindex = " << index;
   if (Pop[index].getID() == id) {
      // Update ID
      Pop[index].setID ( ((Pop[index].getID()) & 0xFFFFFFFF) | ((( Pop[index].getID() >> 32 ) + 1 ) << 32) );
      // Reset NSymbiont and SPopID
      Pop[index].setNsymbiont ( 0 );
//      Pop[index].setSPopID ( 0 );
      // Get the ID of the last living individual in the population vector (to update indlist)
      uint32_t LastAliveID = Pop[N-1].getID();
      // Swap the removed individual with the last living individual in the populaiton vector
      std::swap(Pop[index],Pop[N-1]);
      // Update the indirection list
      std::swap(indList[(id & 0xFFFFFFFF)], indList[(LastAliveID & 0xFFFFFFFF)]);
      // Return index part of the removed id to the freelist
      freeList.push( id & 0xFFFFFFFF );
//std::cout << "\n";
//printiList();
//std::cout << "\n\n\n";
      --N;
      }
}

template<typename T>
   void Population<T>::setPop(const std::vector<T>& pop) { Pop = pop; }
template<typename T>
   const std::vector<T>& Population<T>::getPop() const {return Pop;}

template<typename T>
   void Population<T>::setindList(const std::vector<int>& inl) {indList = inl;}
template<typename T>
   std::vector<int> Population<T>::getindList() const {return indList;}
template<typename T>
   void Population<T>::initindList() {
   std::iota (std::begin(indList), std::end(indList), 0);
   }

template<typename T>
   void Population<T>::setfreeList(const std::stack<uint32_t>& frl) {freeList = frl;}
template<typename T>
   std::stack<uint32_t> Population<T>::getfreeList() const {return freeList;}
template<typename T>
   void Population<T>::initfreeList() {
   for (uint32_t n = indList.size() ; n > 0 ; --n) {freeList.push(n - 1);}
   }

template<typename T>
   void Population<T>::setN(int n) {N = n;}
template<typename T>
   int Population<T>::getN() const {return N;}

template<typename T>
   void Population<T>::setPatchID(uint64_t pid) {PatchID = pid;}
template<typename T>
   uint64_t Population<T>::getPatchID() const {return PatchID;}

template<typename T>
   void Population<T>::setID(uint64_t id) {ID = id;}
template<typename T>
   uint64_t Population<T>::getID() const {return ID;}

template<typename T>
void Population<T>::initPop( Patch& patch ) {
   // Initialize indirection list and free list:
   initindList();
   initfreeList();
   // Initialize Patch_ID and patch.HPopID
   setPatchID( patch.getID() );
   patch.setHPopID( getID() );
}

template<typename T> // Function for initial migration event
void Population<T>::initImmigrFromSource( const SourcePatch& continent, Rng& rng, Metapopulation<Population<Symbiont>>& smpop ) {
   // Calculate the number of new immigrants from the source population
   int n_ind = rng.poisson( continent.getHNmigrants() );
   while ( n_ind == 0 ) { // To sample an n_ind > 0
      n_ind = rng.poisson( continent.getHNmigrants() );
   }
   if ( n_ind > static_cast<int>(getPop().size()) ) { // To prevent from Pop vector overflow
      n_ind = static_cast<int>(getPop().size());
   }
   // Create the new source individuals
   for ( int counter = 0; counter < n_ind; ++counter ) {
      int64_t newHostID = newIndFromSource( continent, rng );
      smpop.newPopFromSource( continent, rng, *this, newHostID );
   }
}

template<typename T> // Function for periodic migration
void Population<T>::immigrFromSource( const SourcePatch& continent, Rng& rng, Metapopulation<Population<Symbiont>>& smpop ) {
   // Calculate the number of new immigrants from the source population
   int n_ind = rng.poisson( continent.getHNmigrants() );
   // Notice that here we do not use: while ( n_ind == 0 ) { n_ind = rng.poisson( continent.getHNmigrants() );}
   // This is because we allow n_ind = 0 for an event of periodic migration
   if ( n_ind > static_cast<int>(getPop().size()) ) { // To prevent from Pop vector overflow
      n_ind = static_cast<int>(getPop().size());
   }
   // Create the new source individuals
   if ( n_ind > 0 ) { // because in this function n_ind can be 0
      for ( int counter = 0; counter < n_ind; ++counter ) {
         int64_t newHostID = newIndFromSource( continent, rng );
         smpop.newPopFromSource( continent, rng, *this, newHostID );
      }
   }
}

template<typename T> // Function for pulse migration
void Population<T>::pulseImmigrFromSource( const SourcePatch& continent, Rng& rng, Metapopulation<Population<Symbiont>>& smpop ) {

   if ( rng.bernoulli( continent.getPrPulseMigr() ) ) { // If the immigration pulse occurs
      // Calculate the number of new immigrants from the source population
      int n_ind = rng.poisson( continent.getHNmigrants() );
      // Force n_ind > 0 (a successful event of pulse migration always has at least 1 host immigrant)
      while ( n_ind == 0 ) { // To sample an n_ind > 0
         n_ind = rng.poisson( continent.getHNmigrants() );
      }
      if ( n_ind > static_cast<int>(getPop().size()) ) { // To prevent from Pop vector overflow
         n_ind = static_cast<int>(getPop().size());
      }
      // Create the new source individuals
      for ( int counter = 0; counter < n_ind; ++counter ) {
         int64_t newHostID = newIndFromSource( continent, rng );
         smpop.newPopFromSource( continent, rng, *this, newHostID );
      }
   }  // End if the immigration pulse occurs
}

template<typename T>
void Population<T>::popReproduction ( const Patch& island, Rng& rng, Metapopulation<Population<Symbiont>>& smpop) {
   if (N>0) { // If the population is not empty
      // Produce the gamete pool:
      hgpool gamPool = produceGametePool( island, rng );
      // Produce newborn individuals from the gamete pool:
      std::vector<int> randgF = rng.randIndexVect( gamPool.first.size() );
      std::vector<int> randgM = rng.randIndexVect( gamPool.second.size() );
//std::cout << "\n\nrandgF vector: \n\t";
//for (int i:randgF) {std::cout<<i<<" ";}
//std::cout << "\n\nrandgM vector: \n\t";
//for (int i:randgM) {std::cout<<i<<" ";}
      int NnewBorn = std::min( gamPool.first.size(), gamPool.second.size() );
//std::cout << "\n\nNnewBorn = " << NnewBorn;
//std::cout << "\n\nSize of female gamPool = " << gamPool.first.size();
//std::cout << "\n\nSize of male gamPool = " << gamPool.second.size();
      if ( NnewBorn > 0 ) { // If we have newborns
         for ( int i = 0; i < NnewBorn; ++i ) { // For each newborn
            // Create a newborn host and get its ID
            uint64_t newbornID = newBorn( rng, gamPool.first[randgF[i]], gamPool.second[randgM[i]] );
            // Create an empty symbiont population in the newborn host and gets its ID
            uint64_t nbSpopID = smpop.newBornPop(*this, newbornID);
            // Get IDs of parental hosts
            uint64_t fparentID = gamPool.first[randgF[i]].getHostID();
            uint64_t mparentID = gamPool.second[randgM[i]].getHostID();
            // Vertical transmission from female parent
            smpop.verTrans( rng, *this, nbSpopID, fparentID );
            // Vertical transmission from male parent
            smpop.verTrans( rng, *this, nbSpopID, mparentID );
         } // End for each newborn
      } // End if we have newborns
   } // End if the population is not empty
}

template<typename T>
void Population<T>::popMortality ( Rng& rng, Metapopulation<Population<Symbiont>>& smpop ) {
   // Get host mortality rate
   double d = Host::getd();
   // Determine the number of dead hosts (Ndead)
   int Ndead = rng.binomial( N, d );
//std::cout << "\nNdead = " << Ndead;
   // Generate Ndead random deaths of a host and its associated symbiont population
   for ( int counter = 0; counter < Ndead; ++counter ) { // for each host death
//   for ( int counter = 0; counter < Ndead; ++counter ) { // for each host death
      // Determine which host dies
      int indexDeadHost = rng.uniform_int( 0, N-1 );
//std::cout << "\nN = " << N;
//std::cout << "\nindexDeadHost = " << indexDeadHost;
      // Get the ID of the host and its associated symbiont population
      int64_t deadHostID = Pop[ indexDeadHost ].getID();
      int64_t deadSPopID = Pop[ indexDeadHost ].getSPopID();
//std::cout << "\nindex part of deadHostID= " << (deadHostID & 0xFFFFFFFF);
      // Remove the host and its associated symbiont population
      removeInd(deadHostID);
      smpop.removePop(deadSPopID);
   } // end for each host death
}

template<typename T>
void Population<T>::printIndividuals () const {
   for (int counter = 0; counter < N; ++counter) {
   Pop[counter].printIndividual();
   }
}

template<typename T>
void Population<T>::printiList () const {
   for (auto const& item : indList) {
   std::cout << item << ' ';
   }
   std::cout << std::endl;
}

template<typename T>
void Population<T>::printPopulation () const {
   std::cout << PopulationtoString() << std::endl;
}

template<typename T>
double Population<T>::getSPrev () const {
   if (N>0) {
      int OccupHosts = 0;
      for (int counter = 0; counter < N; ++counter) {
         if (Pop[counter].getNsymbiont() > 0) { ++OccupHosts; }
      }
      return ( static_cast<double>(OccupHosts)/static_cast<double>(N) );
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getAvSAbOcH() const {
   if (N>0) {
      int OccupHosts = 0;
      int TotalSAb = 0;
      for (int counter = 0; counter < N; ++counter) {
         if (Pop[counter].getNsymbiont() > 0) {
            ++OccupHosts;
            TotalSAb += Pop[counter].getNsymbiont();
         }
      }
      if ( OccupHosts>0 ) { return ( static_cast<double>(TotalSAb) / static_cast<double>(OccupHosts) ); }
      else { return ( 0 ); }
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getAvSAb() const {
   if (N>0) {
      int TotalSAb = 0;
      for (int counter = 0; counter < N; ++counter) {
         TotalSAb += Pop[counter].getNsymbiont();
      }
      return ( static_cast<double>(TotalSAb) / static_cast<double>(N) );
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getVarSAb() const {
   if (N>0) {
      double mean = getAvSAb();
      if ( mean > 0 ) {
         double var = 0;
         for (int counter = 0; counter < N; ++counter) {
            var += pow( ( static_cast<double>(Pop[counter].getNsymbiont()) - mean ), 2) ;
         }
         return ( var );
      }
      else { return ( 0 ); }
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getStdSAb() const {
   if (N>0) {
      double mean = getAvSAb();
      if ( mean > 0 ) {
         double var = 0;
         for (int counter = 0; counter < N; ++counter) {
            var += pow( ( static_cast<double>(Pop[counter].getNsymbiont()) - mean ), 2) ;
         }
         return ( sqrt(var) );
      }
      else { return ( 0 ); }
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getVmrSAb() const {
   if (N>0) {
      double mean = getAvSAb();
      if (mean > 0) {
         double var = getVarSAb();
         return ( var/mean );
      }
      else { return ( 0 ); }
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getCvSAb() const {
   if (N>0) {
      double mean = getAvSAb();
      if (mean > 0) {
         double std = getStdSAb();
         return ( std/mean );
      }
      else { return ( 0 ); }
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getAvHPhen() const {
   if (N>0) {
      double SumPhen = 0;
      for (int counter = 0; counter < N; ++counter) {
         SumPhen += Pop[counter].getPhen();
      }
      return ( SumPhen / static_cast<double>(N) );
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getStdHPhen() const {
   if (N>0) {
      double mean = getAvHPhen();
      double var = 0;
      for (int counter = 0; counter < N; ++counter) {
         var += pow( ( Pop[counter].getPhen() - mean ), 2) ;
      }
      return ( sqrt(var) );
   }
   else { return ( 0 ); }
}

template<typename T>
double Population<T>::getHexp () const {
   std::vector<int> nAl1 = getAlFreq();
   double numerator_Hexp = 0;
   double nal = static_cast<double>(N) * 2; // nal is the total number of alleles per locus in the population
   double nloci = static_cast<double>( T::getL() );
   for ( int i = 0; i < static_cast<int>( nAl1.size() ) ; ++i ) {
      double nal1 = static_cast<double>(nAl1[i]);
      numerator_Hexp += 2 * nal1 * ( nal - nal1 );
   }
   return ( numerator_Hexp / ( nloci * nal * nal ) );
}

template<typename T>
std::vector<int> Population<T>::getAlFreq() const {
   std::vector<int> alfreq( Host::getL() );
   if ( N > 0 ) {
      for ( int counter = 0; counter < N; ++counter ) {
         std::pair<uint32_t, uint32_t> pairgen = pair32Int( Pop[counter].getGen() );
         sumAlToAlFreq( pairgen.first, alfreq );
         sumAlToAlFreq( pairgen.second, alfreq );
      }
   }
   return alfreq;
}

template<typename T>
double Population<T>::getHobs () const {
   double sumhetpop = static_cast<double>( sumHetPop() );
   double n = static_cast<double>( N );
   double l = static_cast<double>( T::getL() );
   return ( sumhetpop / ( n * l ) );
}

   // Utility functions

template<typename T>
std::string Population<T>::PopulationtoString() const {
   std::ostringstream output;
   output << "N = " << getN()
   << "\nPatch ID: " << getPatchID()
   << "\nID: " << getID()
   << "\t" << "Index: " << static_cast<uint32_t>((getID() << 32) >> 32)
   << "\t" << "Version: " << static_cast<uint32_t>(getID() >> 32) << "\n\n";
   return output.str();
}

template<typename T>
   uint64_t Population<T>::createInd() {
   if (!freeList.empty()) {
      uint32_t newID = freeList.top(); // freeList elements are uint32_t
      freeList.pop();
      Pop[N].setID (static_cast<uint64_t>(newID) | ((( Pop[N].getID() >> 32 ) + 1 ) << 32)); // set the ID's lower32 part (ie. the index) equal to newID, and increment by 1 the upper32 part (i.e. the version)
      uint64_t newID64 = Pop[N].getID();
      indList[newID] = N;
      ++N;
      return newID64;
      }
   else {
   std::cout << "The free list of Host Population is empty";
   exit(1);
   }
}

template<typename T>
   uint64_t Population<T>::GenfromGametes(const Gamete& gamete1, const Gamete& gamete2) const {
   // (uint64_t) mostSignificantWord << 32 | leastSignificantWord
   uint64_t newGen = (uint64_t) gamete1.getHaplGen() << 32 | gamete2.getHaplGen();
   return(newGen);
}

template<typename T> // Function for host populations
typename Population<T>::hgpool Population<T>::produceGametePool( const Patch& island, Rng& rng) const {
   std::vector<HostGamete> FemaleGam;
   std::vector<HostGamete> MaleGam;
   FemaleGam.reserve(island.getKhost());
   MaleGam.reserve(island.getKhost());
   for ( int i = 0; i < N; ++i ) {
      std::vector<HostGamete> GamInd = Pop[i].produceGametes( island, N, rng );
      if ( Pop[i].getSex() == 'f' ) {
         FemaleGam.insert( FemaleGam.end(), GamInd.begin(),GamInd.end());
//std::cout << "Added to female pool: " << GamInd.size() << " gametes\n";
//for (const auto& i:GamInd) {i.printHaplGen();}
      }
      else {
         MaleGam.insert( MaleGam.end(), GamInd.begin(),GamInd.end());
//std::cout << "Added to male pool: " << GamInd.size() << " gametes\n";
//for (const auto& i:GamInd) {i.printHaplGen();}
      }
   }
   return (std::pair<std::vector<HostGamete>,std::vector<HostGamete>>( FemaleGam, MaleGam ));
}

template<typename T>
void Population<T>::sumAlToAlFreq ( uint32_t haplgen , std::vector<int>& alfrq ) const {
   const uint32_t SHIFT{8 * sizeof(uint32_t) - 1};
   const uint32_t MASK{static_cast<const uint32_t>(1 << SHIFT)};
   for (uint32_t i{1}; i <= SHIFT + 1; ++i) {
      if (haplgen & MASK) {++alfrq[i-1];}
      haplgen <<= 1; // shift haplgen left by 1
   }
}

template<typename T>
std::pair<uint32_t, uint32_t> Population<T>::pair32Int (uint64_t value) const {
    return std::pair<uint32_t, uint32_t>((value << 32) >> 32, value >> 32);
}

template<typename T>
int Population<T>::sumHetPop () const {
   int sumhetpop = 0;
   for ( int i=0; i < N ; ++i ) {
      sumhetpop += Pop[i].sumHetLocInd();
   }
   return(sumhetpop);
}

template<typename T>
uint64_t Population<T>::genMutation (uint64_t gen, Rng& rng) const {
   // Determine the number of mutations that occurred in the genotype (nmut)
   int nal = T::getL() + T::getL(); // nal is the number of alleles per genotype
   int nmut = rng.binomial( nal , T::getmutRate() );
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

#endif // POPULATION_H
