// Patch class member function definitions

#include <iostream>
#include <string>
#include <iomanip>
#include <stdexcept>
#include <sstream>
#include "Patch.h"
using namespace std;

// constructor
Patch::Patch( int kh, double op, uint64_t id, uint64_t hpopid): Khost(kh), OptPhen(op), ID(id), HPopID(hpopid) {}

// Non-static member functions


void Patch::setKhost(int Kh) {Khost = Kh;}
int Patch::getKhost() const {return Khost;}

void Patch::setOptPhen(double op) {OptPhen = op;}
double Patch::getOptPhen() const {return OptPhen;}

void Patch::setID(uint64_t id) {ID = id;}
uint64_t Patch::getID() const {return ID;}


void Patch::setHPopID(uint64_t hpopid) {HPopID = hpopid;}
uint64_t Patch::getHPopID() const {return HPopID;}

void Patch::printPatch() {
   cout << toString() << endl;
}

   // Utility functions

string Patch::toString() {
   ostringstream output;
   output << fixed << setprecision(2);
   output << "K host: " << getKhost()
      << "\nOptimal phenotype = " << getOptPhen()
      << "\nID = " << getID()
      << "\t" << "Index: " << static_cast<uint32_t>((getID() << 32) >> 32)
      << "\t" << "Version: " << static_cast<uint32_t>(getID() >> 32)
      << "\nID = " << getHPopID()
      << "\t" << "Index: " << static_cast<uint32_t>((getHPopID() << 32) >> 32)
      << "\t" << "Version: " << static_cast<uint32_t>(getHPopID() >> 32);
   return output.str();
}
