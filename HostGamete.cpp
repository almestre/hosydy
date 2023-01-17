// Organism base class member function definitions

#include <cstdint> // uint32_t and uint64_t types
#include <string>
#include <bitset>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "HostGamete.h"
using namespace std;

// constructor
HostGamete::HostGamete(uint32_t hgen, uint64_t hid): Gamete (hgen) {
   setHostID(hid);
}

// Set HostID
void HostGamete::setHostID(uint64_t hid) {HostID = hid;}

// Get HostID
uint64_t HostGamete::getHostID() const {return HostID;}

// Print HostID
void HostGamete::printHostID() const {
   auto pairHostID{pair32Int(HostID)};
   cout << "Host ID: " << endl;
   cout <<  "Index\t";
   displayBits(pairHostID.first);
   cout <<  "Version\t";
   displayBits(pairHostID.second);
}

// ---Utility functions---

void HostGamete::displayBits(uint32_t value) const {
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

pair<uint32_t, uint32_t> HostGamete::pair32Int (uint64_t value) const {
    return pair<uint32_t, uint32_t>((value << 32) >> 32, value >> 32);
}

