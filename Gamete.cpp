// Organism base class member function definitions

#include <cstdint> // uint32_t and uint64_t types
#include <string>
#include <bitset>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "Gamete.h"
using namespace std;

// constructor
Gamete::Gamete(uint32_t hgen): HaplGen(hgen) {}

// Set HaplGen
void Gamete::setHaplGen(uint32_t hgen) {HaplGen = hgen;}

// Get HaplGen
uint32_t Gamete::getHaplGen() const {return HaplGen;}

// Print HaplGen
void Gamete::printHaplGen() const {
   cout << "Haploid genotype: " << endl;
   displayBits(HaplGen);
}

// ---Utility functions---

void Gamete::displayBits(uint32_t value) const {
   const uint32_t SHIFT{8 * sizeof(uint32_t) - 1};
   const uint32_t MASK{static_cast<const uint32_t>(1 << SHIFT)};

   cout << setw(10) << value << " = ";

   // display bits
   for (uint32_t i{1}; i <= SHIFT + 1; ++i) {
   cout << (value & MASK ? '1' : '0');
   value <<= 1; // shift value left by 1

   if (i % 8 == 0) { // output a space after 8 bits
      cout << ' ';
      }
   }

   cout << endl;
}
