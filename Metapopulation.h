// Metapopulation template class definition
// Class template member functions are defined within the class definition's body (no interface-implementation separation)

#ifndef METAPOPULATION_H
#define METAPOPULATION_H

#include <sstream>
#include <iostream>
#include <numeric> // function std::iota
#include <cstdint> // uint32_t and uint64_t types
#include <cstdlib> // exit() function
#include<vector> // C++ standard vector class template
#include<stack> // C++ standard stack class template
#include <string> // C++ standard string class
#include "Host.h" // Organism class definition
#include "Symbiont.h" // Organism class definition
#include "SourcePatch.h"
#include "Gamete.h"
//#include "Population.h"
#include "Param.h"

// Forward declarations:
template<typename T>
class Population;

/**** Primary template (member definitions in the header) ****/
template<typename T>
class Metapopulation {
public:
   Metapopulation (const size_t = Param::getHPopVecSize(), const std::stack<uint32_t>& = std::stack<uint32_t>(), int = 0, uint64_t = 0);

// Definitions of operators for Pears. coef. calculation (definitions of their implementation are at the end of the header)
friend std::vector<double> operator-(const std::vector<double>&, double); // Used for Pearson cor. coef.
friend std::vector<double> operator*(const std::vector<double>&, const std::vector<double>&); // Used for Pearson cor. coef.

   void initMetapop( ); // version for symbiont metapop

      // Overloaded functions
   void newPopFromSource( const SourcePatch&, Rng&, Patch& ); // version for host metapop
   void newPopFromSource( const SourcePatch&, Rng&, Population<Host>&, uint64_t ); // version for symbiont metapop

   uint64_t newBornPop(Population<Host>&, uint64_t); // Creates an empty symbiont population from a new born host

   void metapopReproduction( Population<Host>&, Rng& ); // version for symbiont metapop

   T* getPop(uint64_t);
   void removePop(uint64_t);

   void horizTrans( Rng&, Population<Host>& );
   void verTrans( Rng&, Population<Host>&, int64_t, int64_t );

   void setMetapop(const std::vector<T>&);
   const std::vector<T>& getMetapop() const;

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

   void printPopulations() const;
   void printiList() const;
   void printMetapopulation() const;

   double getGAvPhen() const; // Get grand mean, i.e. overall mean of all the individuals in the metapopulation
   double getMsbPhen() const; // Get symbiont phenotype variation between hosts (mean square between groups)
   double getMswPhen() const; // Get symbiont phenotype variation within hosts (mean square within groups)
   double getHSPhenCor( Population<Host>& ) const; // Get correlation between the phenotype of a host individual and the average phenotype of the symbionts it harbours
   double getEvC( Population<Host>& ) const; // Get the evolutionary change (EvC) for symbionts
   std::vector<int> getAlFreq() const;
   double getHexp() const; // Get the expected heterozigosity of the symbiont metapopulation based on Hardy-Weinberg equilibrium
   double getHobs() const; // Get the observed heterozigosity of the symbiont metapopulation


private:
   std::vector<T> Metapop; // Vector of populations representing the metapopulation
   std::vector<int> indList; // Indirection list to access populations by id
   std::stack<uint32_t> freeList; // List of free id to use when new populations are created
   int N; // Number of populations in the metapopulation
   uint64_t PatchID; // id of the patch inhabited by the metapopulation; in the case of a symbiont metapopulation, the patch is a host population. This PatchID will be useful if we create a host metapopulation.

   // Utility functions
   std::string MetapopulationtoString() const;
      // Overloaded function
   uint64_t createPop(); // version for host metapop
   uint64_t createPop( Population<Host>& , uint64_t ); // version for symbiont metapop
//   double getCor( const std::vector<double>& , const std::vector<double>&, int ) const; // function that returns correlation coefficient
   double sum(const std::vector<double>&) const; // Used for Pearson cor. coef.
   double mean(const std::vector<double>&) const; // Used for Pearson cor. coef.
   double sqsum(const std::vector<double>&) const; // Used for Pearson cor. coef.
   double stdev(const std::vector<double>&) const; // Used for Pearson cor. coef.
   double getPearsCoef(const std::vector<double>&, const std::vector<double>&) const; // Used for Pearson cor. coef.

};

// Including Population class template full definition after the Metapopulation class definition, and before its implementation
#include "Population.h"

   // Constructor
template<typename T>
Metapopulation<T>::Metapopulation(const size_t MaxMetapSize, const std::stack<uint32_t>& frl, int n, uint64_t pid): Metapop(std::vector<T>(MaxMetapSize)), indList(std::vector<int>(MaxMetapSize)), freeList(frl), N(n), PatchID(pid) {}

   // Non-static member functions

template<typename T> // version for symbiont metapop
void Metapopulation<T>::initMetapop( ) {
   // Initialize indirection list and free list:
   initindList();
   initfreeList();
   // If we create a host metapopulation, PatchID should be initialised here
}

template<typename T> // Overloaded function for host metapopulations
void Metapopulation<T>::newPopFromSource( const SourcePatch& continent, Rng& rng, Patch& patch ) {
   // Create new population and get its id
   uint64_t popID = createPop();
   T* pop = getPop(popID);
   pop->initFromSource();
   // Set the pointer to null
   pop = nullptr;
}

template<typename T> // Overloaded function for symbiont metapopulations
void Metapopulation<T>::newPopFromSource( const SourcePatch& continent, Rng& rng, Population<Host>& hpop, uint64_t hostid ) {
   // Create new symbiont population and get its ID
   uint64_t popID = createPop(hpop,hostid);
   // Get the new symbiont population
   Population<Symbiont>* pop = getPop(popID);
   // Generate a pulse of immigration from the source population
   pop->immigrFromSource( continent, rng, hpop );
   // Set the pointer to null
   pop = nullptr;
}

template<typename T> // Function to create an empty symbiont population from a new born host
   uint64_t Metapopulation<T>::newBornPop(Population<Host>& hpop, uint64_t newborn_id) {
   // Create a new (empty) symbiont population associated with the new born host and get its ID
   uint64_t sPopID = createPop(hpop,newborn_id);
   // Return the ID of the new symbiont population:
   return sPopID;
}

template<typename T> // Overloaded function for symbiont metapopulations
   void Metapopulation<T>::metapopReproduction( Population<Host>& hpop, Rng& rng) {
   for ( int counter = 0; counter < N; ++counter ) { // For each population
   Metapop[counter].popReproduction ( hpop, rng );
   } // End for each population
}

template<typename T>
T* Metapopulation<T>::getPop(uint64_t id) {
   uint32_t index = indList[ (id & 0xFFFFFFFF) ]; // 'id & 0xFFFFFFFF' is the lower 32-bit part of the ID (i.e. the index)
   if (Metapop[index].getID() == id) {
      return &Metapop[index];
      }
   else {return nullptr;}
}

template<typename T>
void Metapopulation<T>::removePop(uint64_t id) {
//std::cout << "\n\nSMetapoplist before one death\n";
//printiList();
   uint32_t index = indList[ (id & 0xFFFFFFFF) ]; // 'id & 0xFFFFFFFF' is the lower 32-bit part of the ID (i.e. the index)
   if (Metapop[index].getID() == id) {
      // Update ID
      Metapop[index].setID ( ((Metapop[index].getID()) & 0xFFFFFFFF) | ((( Metapop[index].getID() >> 32 ) + 1 ) << 32) );
      // Reset N and PatchID
      Metapop[index].setN ( 0 );
//      Metapop[index].setPatchID ( 0 );
      // Get the ID of the last living individual in the population vector (to update indlist)
      uint32_t LastAliveID = Metapop[N-1].getID();
      // Swap the removed population with the last living population in the Metapop vector
      std::swap(Metapop[index],Metapop[N-1]);
      // Update the indirection list
      std::swap(indList[(id & 0xFFFFFFFF)],indList[(LastAliveID & 0xFFFFFFFF)]);
//std::cout << "\n\nSMetapoplist after one death\n";
//printiList();
//std::cout << "\nIndex pushed" << index;
//std::cout << "\n\n\n";
      --N;
      freeList.push(id & 0xFFFFFFFF);
      }
}

template<typename T> // Function to implement horizontal transmission in a symbiont metapopulation
void Metapopulation<T>::horizTrans( Rng& rng, Population<Host>& hpop) {
   // Get emigration rate
   if ( N > 1 ) { // if we have more than one symbiont population
      double e = Symbiont::getEht();
      // Create a vector with randomly ordered population indexes to get populations in random order
      std::vector<int> vecRandPop = rng.randIndexVect( N );
      for ( int counter = 0; counter < N; ++counter ) { // for each population
         // Get index and population size of the source population
         int indexSrcPop = vecRandPop[counter];
         int n = Metapop[indexSrcPop].getN();
         if (n>0) { // If the population is not empty
            // Determine number of emigrants from the source population
            int Nemigrants = rng.binomial( n, e );
//std::cout << "\nN_emigrants = " << Nemigrants << " in loop " << counter;
//std::cout << "\n\tSource population index: " << indexSrcPop;
            // Distribute n random emigrants (where n=Nemigrants) among other randomly selected populations
            for ( int counter2 = 0; counter2 < Nemigrants; ++counter2 ) { // for each emigrant in a population
               // Determine which individual emigrates
               int indexEmigrant = rng.uniform_int( 0, n-1 );
//std::cout << "\n\tEmigrant index: " << indexEmigrant;
               // Determine the population of destination
               int indexDestPop = rng.uniform_int( 0, N-1 );
               // Make sure we do not select the source population (warning: potential infinite loop if N=1!)
               while ( indexDestPop == indexSrcPop ) { indexDestPop = rng.uniform_int( 0, N-1 ); }
               // Create a new immigrant in the population of destination
               int indexNewIm = Metapop[indexDestPop].newImmigrant( hpop );
//std::cout << "\n\tDestination population index: " << indexDestPop;
//std::cout << "\n\tNew immigrant index: " << indexNewIm;
               // Swap data between target locations of the source and destination vectors
               Symbiont* ptrSrcLoc = Metapop[indexSrcPop].getPtrToInd( indexEmigrant );
               Symbiont* ptrDestLoc = Metapop[indexDestPop].getPtrToInd( indexNewIm );
               std::swap(
                  *(ptrSrcLoc), // Source
                  *(ptrDestLoc) // Destination
               );
               // Remove emigrant from the source population:
               Metapop[indexSrcPop].removeInd( hpop, indexEmigrant );
               // Update size of the source population (n)
               n = Metapop[indexSrcPop].getN();
               // Set pointers to null:
               ptrSrcLoc = nullptr;
               ptrDestLoc = nullptr;
            } // End for each emigrant in a population
         }  // End if the population is not empty
      } // End for each population
   } // End if we have more than one symbiont population
}

template<typename T> // Function to implement vertical transmission from a parent to a newborn
void Metapopulation<T>::verTrans( Rng& rng, Population<Host>& hpop, int64_t nbSpop_id, int64_t parent_id ) {
   // Get emigration rate
   double e = Symbiont::getEvt();
   // Get host parent and symbiont populations involved
   Host* ParentPtr = hpop.getInd(parent_id);
   Population<Symbiont>* nbSpopPtr = getPop( nbSpop_id );
   Population<Symbiont>* pSpopPtr = getPop( ParentPtr->getSPopID() );
//std::cout << "\nParent\n";
//ParentPtr->printIndividual();
//printiList();
//printMetapopulation();
//printPopulations();
//std::cout << "\nSPopNB\n";
//nbSpopPtr->printPopulation();
//std::cout << "\nSPopP\n";
//pSpopPtr->printPopulation();
   // Set pointers to null:
   ParentPtr = nullptr;
   if (pSpopPtr->getN() > 0) { // if the parent harbours symbionts
      // Determine number of emigrants from the parent to the newborn
      int Nemigrants = rng.binomial( pSpopPtr->getN(), e);
//std::cout << "\n\ne =" << e;
//std::cout << "\nNemigrants =" << Nemigrants;
//std::cout << "\nN =" << pSpopPtr->getN();
      if ( Nemigrants > 0) { // If we have immigrants
         // Transfer n random emigrants (where n=Nemigrants) from fparent to newborn
         for ( int counter = 0; counter < Nemigrants; ++counter ) { // for each emigrant in the parent's population
            // Determine which individual emigrates
            int indexEmigrant = rng.uniform_int( 0, pSpopPtr->getN()-1 );
//std::cout << "\n\tEmigrant index: " << indexEmigrant;
            // Create a new immigrant in the population of destination (i.e. newborn)
            int indexNewIm = nbSpopPtr->newImmigrant( hpop );
            // Swap data between target locations of the source and destination vectors
            Symbiont* ptrSrcLoc = pSpopPtr->getPtrToInd( indexEmigrant );
            Symbiont* ptrDestLoc = nbSpopPtr->getPtrToInd( indexNewIm );
            std::swap(
               *(ptrSrcLoc), // Source
               *(ptrDestLoc) // Destination
            );
            // Remove emigrant from the source population:
            pSpopPtr->removeInd( hpop, indexEmigrant );
            // Set pointers to null:
            ptrSrcLoc = nullptr;
            ptrDestLoc = nullptr;
         } // End for each emigrant in the parent's population
      } // End if we have immigrants
   }  // End if the parent harbours symbionts
   // Set pointers to null:
   nbSpopPtr = nullptr;
   pSpopPtr = nullptr;
}

template<typename T>
void Metapopulation<T>::setMetapop(const std::vector<T>& metapop) { Metapop = metapop; }
template<typename T>
const std::vector<T>& Metapopulation<T>::getMetapop() const {return Metapop;}

template<typename T>
void Metapopulation<T>::setindList(const std::vector<int>& inl) {indList = inl;}
template<typename T>
std::vector<int> Metapopulation<T>::getindList() const {return indList;}
template<typename T>
void Metapopulation<T>::initindList() {
   std::iota (std::begin(indList), std::end(indList), 0);
   }

template<typename T>
void Metapopulation<T>::setfreeList(const std::stack<uint32_t>& frl) {freeList = frl;}
template<typename T>
std::stack<uint32_t> Metapopulation<T>::getfreeList() const {return freeList;}
template<typename T>
void Metapopulation<T>::initfreeList() {
   for (uint32_t n = indList.size() ; n > 0 ; --n) {freeList.push(n - 1);}
}

template<typename T>
void Metapopulation<T>::setN(int n) {N = n;}
template<typename T>
int Metapopulation<T>::getN() const {return N;}

template<typename T>
void Metapopulation<T>::setPatchID(uint64_t pid) {PatchID = pid;}
template<typename T>
uint64_t Metapopulation<T>::getPatchID() const {return PatchID;}

template<typename T>
void Metapopulation<T>::printPopulations () const {
   for (int counter = 0; counter < N; ++counter ) {
   Metapop[counter].printPopulation();
   }
}

template<typename T>
void Metapopulation<T>::printiList () const {
   for (auto const& item : indList) {
   std::cout << item << ' ';
   }
   std::cout << std::endl;
}

template<typename T>
void Metapopulation<T>::printMetapopulation () const {
   std::cout << MetapopulationtoString() << std::endl;
}

template<typename T>
double Metapopulation<T>::getGAvPhen() const {
   if (N>0) {
      double TotalSumPhen = 0;
      int TotalN = 0;
      for ( int counter=0; counter < N; ++counter ) {
         TotalSumPhen += Metapop[counter].getSumPhen();
         TotalN += Metapop[counter].getN();
      }
      if ( TotalN > 0 ) {
         return( TotalSumPhen / static_cast<double>(TotalN) );
      }
      else { return(0); }
   }
   else { return(0); }
}

template<typename T>
double Metapopulation<T>::getMsbPhen() const {
   if (N>0) {
      double SSB = 0;
      double gmean = getGAvPhen();
      for ( int counter=0; counter < N; ++counter ) {
         SSB += Metapop[counter].getN() * pow ( ( Metapop[counter].getAvSPhen() - gmean ), 2) ;
      }
      if ( (N - 1) > 0 ) {
         return( SSB / static_cast<double>( N - 1 ) );
      }
      else { return(0); }
   }
   else { return(0); }
}

template<typename T>
double Metapopulation<T>::getMswPhen() const {
   if (N>0) {
      double SSW = 0;
      int TotalN = 0;
      for ( int counter=0; counter < N; ++counter ) {
         SSW += Metapop[counter].getSSPhen();
         TotalN += Metapop[counter].getN();
//std::cout << "\nSSW = " << SSW;
      }
      if ( (TotalN > 0) && ( TotalN > N  ) ) {
         return( SSW / static_cast<double>( TotalN - N ) );
      }
      else { return(0); }
   }
   else { return(0); }
}

template<typename T>
double Metapopulation<T>::getHSPhenCor( Population<Host>& hpop ) const {
   int NHocc = 0;
   std::vector<double> hphen;
   hphen.reserve(N);
   std::vector<double> sphen;
   hphen.reserve(N);
   for ( int counter = 0; counter < N; ++counter ) {
      if ( Metapop[counter].getN() > 0 ) {
         sphen.push_back(Metapop[counter].getAvSPhen());
         Host* ptrHost = hpop.getInd ( Metapop[counter].getPatchID() );
         hphen.push_back(ptrHost->getPhen());
         ptrHost = nullptr;
         ++NHocc;
      }
   }
   return ( getPearsCoef ( hphen, sphen ) );
}

template<typename T>
double Metapopulation<T>::getEvC( Population<Host>& hpop ) const {
   int NHocc = 0;
   double sumDevSPhen = 0;
   for ( int counter = 0; counter < N; ++counter ) {
      if ( Metapop[counter].getN() > 0 ) {
         sumDevSPhen += Metapop[counter].getDevSPhen( hpop );
         ++NHocc;
      }
   }
   return ( sumDevSPhen / NHocc );
}

template<typename T>
std::vector<int> Metapopulation<T>::getAlFreq() const {
   std::vector<int> alfreq( Symbiont::getL() );
   if (N>0) {
      for ( int counter = 0; counter < N ; ++counter ) {
         Metapop[counter].getAlFreq( alfreq );
      }
   }
   return( alfreq );
}

template<typename T>
double Metapopulation<T>::getHexp() const {
   std::vector<int> nal1( Symbiont::getL() );
   int ntotal = 0;
   if (N>0) {
      for ( int counter = 0; counter < N ; ++counter ) {
         Metapop[counter].getAlFreq( nal1 );
         ntotal += Metapop[counter].getN();
      }
   }
   double numerator_Hexp = 0;
   double nal = static_cast<double>(ntotal) * 2; // nal is the total number of alleles per locus in the
   double nloci = static_cast<double>( Organism::getL() );
   for ( int i = 0; i < static_cast<int>( nal1.size() ) ; ++i ) {
      int nal1_int = nal1[i];
      double nal1 = static_cast<double>(nal1_int);
      numerator_Hexp += 2 * nal1 * ( nal - nal1 );
   }
   return ( numerator_Hexp / ( nloci * nal * nal ) );
}

template<typename T>
double Metapopulation<T>::getHobs () const {
   int smhtpp = 0;
   int ntotal = 0;

   for ( int i=0; i < N ; ++i ) {
      smhtpp += Metapop[i].getSumHetPop();
      ntotal += Metapop[i].getN();
   }
   double sumhetpop = static_cast<double>( smhtpp );
   double n = static_cast<double>( ntotal );
   double l = static_cast<double>( Organism::getL() );

   return ( sumhetpop / ( n * l ) );
}

   // Utility functions

template<typename T>
std::string Metapopulation<T>::MetapopulationtoString() const {
   std::ostringstream output;
   output << "N = " << getN()
   << "\nPatch ID: " << getPatchID() << "\n\n";
   return output.str();
}

template<typename T> // Overloaded function for host metapopulations
uint64_t Metapopulation<T>::createPop() {
   if (!freeList.empty()) {
      uint32_t newID = freeList.top(); // freeList elements are uint32_t
      freeList.pop();
      Metapop[N].setID (static_cast<uint64_t>(newID) | ((( Metapop[N].getID() >> 32 ) + 1 ) << 32)); // set the ID's lower32 part (ie. the index) equal to newID, and increment by 1 the upper32 part (i.e. the version)
      uint64_t newID64 = Metapop[N].getID();
      indList[newID] = N;
      ++N;
      return newID64;
      }
   else {
   exit(1);
   }
}

template<typename T> // Overloaded function for symbiont metapopulations
uint64_t Metapopulation<T>::createPop( Population<Host>& hpop, uint64_t hostid ) {
//std::cout << "\n\nSMetapoplist before new Spop born\n";
//printiList();
   if (!freeList.empty()) {
      uint32_t newID = freeList.top(); // freeList elements are uint32_t
      freeList.pop();
      // Initialize PatchID with the host ID
      Metapop[N].setPatchID( hostid );
      // Set ID of the new symbiont population
      Metapop[N].setID (static_cast<uint64_t>(newID) | ((( Metapop[N].getID() >> 32 ) + 1 ) << 32)); // set the ID's lower32 part (ie. the index) equal to newID, and increment by 1 the upper32 part (i.e. the version)
      // Update Host.SPopID
      Host* ptrHost = hpop.getInd ( hostid ); // hpop is the name of the host population
      ptrHost->setSPopID(Metapop[N].getID());
      ptrHost = nullptr;
      uint64_t newID64 = Metapop[N].getID();
      indList[newID] = N;
//std::cout << "\n\nSMetapoplist after new Spop born\n";
//printiList();
//std::cout << "\n\n\n";
      ++N;
      return newID64;
      }
   else {
   exit(1);
   }
}

template<typename T> // Used for Pearson cor. coef.
double Metapopulation<T>::sum(const std::vector<double>& x) const {
   double s = 0;
   for (int i = 0; i < static_cast<int>(x.size()); i++)
      {
         s += x[i];
      }
   return s;
}

template<typename T> // Used for Pearson cor. coef.
double Metapopulation<T>::mean(const std::vector<double>& x) const {
   return sum(x) / static_cast<double>(x.size());
}

template<typename T> // Used for Pearson cor. coef.
double Metapopulation<T>::sqsum(const std::vector<double>& x) const {
   double s = 0;
   for (int i = 0; i < static_cast<int>(x.size()); i++)
      {
         s += pow(x[i], 2);
      }
   return s;
}

template<typename T> // Used for Pearson cor. coef.
double Metapopulation<T>::stdev(const std::vector<double>& x) const {
   double N = static_cast<double>(x.size());
   return pow(sqsum(x) / N - pow(sum(x) / N, 2), 0.5);
}

template<typename T>
double Metapopulation<T>::getPearsCoef(const std::vector<double>& x, const std::vector<double>& y) const {
   return sum((x - mean(x))*(y - mean(y))) / (static_cast<double>(x.size())*stdev(x)* stdev(y));
}


// Defining implementation of operators for vectors (outside the class)
// Used for Pearson cor. coef.
inline std::vector<double> operator-(const std::vector<double>& x, double y) { // inline used to avoid multiple definitions in the same header
   std::vector<double> retvect;
   retvect.reserve(x.size());
   for (int i = 0; i < static_cast<int>(x.size()); i++)
      {
         retvect.push_back(x[i] - y);
      }
   return retvect;
}

inline std::vector<double> operator*(const std::vector<double>& x, const std::vector<double>& y) { // inline used to avoid multiple definitions in the same header
   std::vector<double> retvect;
   retvect.reserve(x.size());
   for (int i = 0; i < static_cast<int>(x.size()) ; i++)
      {
         retvect.push_back(x[i] * y[i]);
      }
   return retvect;
}

#endif // METAPOPULATION_H
