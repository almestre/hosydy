// Organism base class definition

#ifndef ORGANISM_H
#define ORGANISM_H

#include <cstdint> // uint32_t and uint64_t types
#include <string>

class Organism {

public:

   explicit Organism(char = '\0', uint64_t = 0, double = 0); // constructor

   void setSex(char); // Set sex (i.e., gender: "m" for male, : "f" for male)
   char getSex() const; // Get sex ("m" for male, : "f" for male)

   void setGen(uint64_t); // Set genotype
   uint64_t getGen() const; // Get genotype

   void setPhen(double); // Set phenotype value
   double getPhen() const; // Get phenotype value
   void Phen_init(); // Initialize phenotype

   void printGen() const; // Print genotype
   std::string toString() const; // Get string representation of the object data members

   uint32_t createOneHaplGen() const; // Creates a new haploid genotip by free recombination
   
   int sumHetLocInd () const; // calculates the sum of the heterozigotic loci of the individual

   // Static member functions (they cannot be declared as constant!)

   static int getL(); // Get number of bi-allelic loci per genotype
   static void setL( int ); // Set number of bi-allelic loci per genotype
   static double getAlpha(); // Get effect size for each allele
   static void setAlpha( double ); // Set effect size for each allele

private:

   char Sex; // gender: "m" for male, : "f" for male
   uint64_t Gen; // genotype
   double Phen; // phenotype value

   // Static data
   static int L; // number of bi-allelic loci per genotype
   static double Alpha; // effect size for each allele

   // Utility functions
   int popcount64b(uint64_t) const; // Count 1-bits in a 64-bit integer
   double sumAlpha(int nbites) const; // Used as part of the function that calculate Phen
   void displayBits(uint32_t) const; // Display bits of a 32-bit uinteger
   std::pair<uint32_t, uint32_t> pair32Int (uint64_t) const; // Create a pair of 32-bit uintegers from a 64-bit uinteger 
   uint32_t freeRecombination(uint32_t,uint32_t) const; // Free-recombination algorithm

   };

   #endif // ORGANISM_H
