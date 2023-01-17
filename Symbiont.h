// Symbiont derived class definition

#ifndef SYMBIONT_H
#define SYMBIONT_H

#include <cstdint> // uint32_t and uint64_t types
#include <string> // C++ standard string class
#include "Organism.h" // Organism class definition
#include "Gamete.h" // Organism class definition
#include "Rng.h" // Population class definition


// Forward declarations (#include in the cpp file)
template<typename T>
class Population;
class Host;

class Symbiont : public Organism {

public:

   explicit Symbiont (char = '\0', uint64_t = 0, double = 0.0); // constructor

   void printIndividual() const; // Print data members of the symbiont individual (function with extended functionality)

   std::vector<Gamete> produceGametes(int64_t, const Population<Host>&, Rng& ) const; // Generate a vector of gametes produced by the symbiont in the current reproductive event

   // Static member functions (they cannot be declared as constant!)
   static double getRmax(); // Get Rmax
   static void setRmax( double ); // Set Rmax
   static double getVs(); // Get Vs
   static void setVs( double ); // Set Vs
   static double getEht(); // Get Eht
   static void setEht( double ); // Set Eht
   static double getEvt(); // Get Evt
   static void setEvt( double ); // Set Evt
   static double getmutRate(); // Get mutRate
   static void setmutRate( double ); // Set mutRate

private:
   static double Rmax; // Intrinsic population growth rate
   static double Vs; // Width of stabilizing selection
   static double Eht; // Per capita per generation emigr. rate for horiz. transm.
   static double Evt; // Per capita per generation emigr. rate for vert. transm.
   static double mutRate; // Per allele, per generation mutation rate

   // Utility functions

   std::string toString() const; // Get string representation of the object data members (function with extended functionality)
   Gamete createOneGamete() const; // Create one gamete
   int calculateNGametes(int64_t, const Population<Host>&, Rng&) const; // Calculate the number of gametes produced by the host in the current reprodutive event
   };

   #endif // SYMBIONT_H
