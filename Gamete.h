// Gamete base class definition

#ifndef GAMETE_H
#define GAMETE_H

#include <cstdint> // uint32_t and uint64_t types
#include <string>

class Gamete {
public:
   explicit Gamete(uint32_t = 0);

   void setHaplGen(uint32_t);
   uint32_t getHaplGen() const;

   void printHaplGen() const;

private:
   uint32_t HaplGen; // genotype

   // Utility functions
   void displayBits(uint32_t) const;

   };

   #endif // GAMETE_H
