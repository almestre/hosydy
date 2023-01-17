// Host derived class definition

#ifndef HOST_H
#define HOST_H

#include <cstdint> // uint32_t and uint64_t types
#include <string> // C++ standard string class
#include "Organism.h" // Organism class definition
#include "HostGamete.h" // HostGamete class definition
#include "Rng.h" // Rng class definition
#include "Patch.h" // Rng class definition

class Host : public Organism {

public:

   explicit Host (char = '\0', uint64_t = 0, double = 0, int = 0, uint64_t = 0, uint64_t = 0); // constructor

   void setNsymbiont(int); // Set number of symbionts inhabiting the host
   int getNsymbiont() const; // Get sex ("m" for male, : "f" for male)

   void setID(uint64_t); // Set unique identifier for the individual host
   uint64_t getID() const; // Get unique identifier for the individual host
   
   void setSPopID(uint64_t); // Set ID for the symbiont population associated to this host
   uint64_t getSPopID() const; // Get ID for the symbiont population associated to this host 
   
   void printID() const; // Print ID (function with extended functionality)
   void printIndividual() const; // Print host individual data (function with extended functionality)

   std::vector<HostGamete> produceGametes(const Patch&, int, Rng& ) const; // Generate a vector of gametes produced by the host in the current reproductive event

   // Overloaded operators

   Host& operator++(); // overlodaded prefix increment operator (add 1 symbiont)
   Host& operator--(); // overlodaded prefix decrement operator (remove 1 symbiont)
   Host& operator+=(int); // overlodaded += operator (add symbionts)
   Host& operator-=(int); // overlodaded -= operator (remove symbionts)

   // Static member functions (they cannot be declared as constant!)

   static int getKsymbiont(); // Get Ksymbiont
   static void setKsymbiont( int ); // Set Ksymbiont
   static double getRmax(); // Get Rmax
   static void setRmax( double ); // Set Rmax
   static double getVs(); // Get Vs
   static void setVs( double ); // Set Vs
   static double getlambda(); // Get lambda
   static void setlambda( double ); // Set lambda
   static double getb(); // Get b
   static void setb( double ); // Set b
   static double getd(); // Get d
   static void setd( double ); // Set d
   static double getmutRate(); // Get mutRate
   static void setmutRate( double ); // Set mutRate

private:

   int Nsymbiont; // Number of symbionts inhabiting the host
   uint64_t ID; // Host ID
   uint64_t SPopID; // ID of the symbiont population inhabiting the host

   // Static data

   static int Ksymbiont; // Carrying capacity of a symbiont population inhabiting the host
   static double Rmax; // Intrinsic population growth rate
   static double Vs; // Width of stabilizing selection
   static double lambda; // strength of effects of symbiont load
   static double b; // growth rate of the interaction effect size
   static double d; // host death rate per capita per reproductive event
   static double mutRate; // Per allele, per generation mutation rate

   // Utility functions

   void displayBits(uint32_t) const; // Display bits of a 32-bit uinteger
   std::pair<uint32_t, uint32_t> pair32Int (uint64_t) const;  // Create a pair of 32-bit uintegers from a 64-bit uinteger
   std::string toString() const; // Get string representation of the object data members (function with extended functionality)
   HostGamete createOneGamete() const; // Create one gamete
   int calculateNGametes(const Patch&, int, Rng&) const; // Calculate the number of gametes produced by the host in the current reprodutive event

   };

   #endif // HOST_H
