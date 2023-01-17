// Organism base class member function definitions

#include <cstdint> // uint32_t and uint64_t types
#include <string>
#include <bitset>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "Organism.h"
#include "Rng.h"
#include "Param.h"

using namespace std;

// ---Static member functions---

int Organism::getL() {return L;}
void Organism::setL( int l ) { L = l; }

double Organism::getAlpha() {return Alpha;}
void Organism::setAlpha( double alpha ) { Alpha = alpha; }

// ---Static data members---

int Organism::L = Param::getL();
double Organism::Alpha = Param::getAlpha();

// constructor
Organism::Organism(char sx, uint64_t gn, double ph): Sex(sx), Gen(gn), Phen(ph) {}

// ---Non-static member functions---

void Organism::setSex(char sx) {Sex = sx;}
char Organism::getSex() const {return Sex;}


void Organism::setGen(uint64_t gn) {Gen = gn;}
uint64_t Organism::getGen() const {return Gen;}


void Organism::setPhen(double ph) {Phen = ph;}
double Organism::getPhen() const {return Phen;}

void Organism::Phen_init() {
   int Nbites1 = popcount64b(Gen);
   int Nbites0 = 64 - Nbites1;
   double phen = sumAlpha(Nbites1) - sumAlpha(Nbites0);
   setPhen(phen);
}

void Organism::printGen() const {
   auto pairGen{pair32Int(Gen)};
   cout << "Individual genotype: " << endl;
   displayBits(pairGen.first);
   displayBits(pairGen.second);
}

string Organism::toString() const {
   ostringstream output;
   output << fixed << setprecision(2);
   output << "Gender: " << getSex()
      << "\nGenotype: " << std::bitset<64>(getGen())
      << "\nPhenotype = " << getPhen()
      << "\nL = " << getL()
      << "\nAlpha = " << getAlpha();
   return output.str();
}

uint32_t Organism::createOneHaplGen() const {
   pair<uint32_t, uint32_t> Hplgens = pair32Int(Gen);
   return (freeRecombination(Hplgens.first,Hplgens.second));
}

int Organism::sumHetLocInd () const {
   pair<uint32_t, uint32_t> Hplgens = pair32Int(Gen);
   const uint32_t SHIFT{8 * sizeof(uint32_t) - 1};
   const uint32_t MASK{static_cast<const uint32_t>(1 << SHIFT)};

   int NLocHet = 0; // NLocHet is the number of heterozigotic loci

   for (uint32_t i{1}; i <= SHIFT + 1; ++i) { // for each locus
      // sum the two allele values of the locus
      int sumAlLoc = (Hplgens.first & MASK ? 1 : 0) + (Hplgens.second & MASK ? 1 : 0);

      if (sumAlLoc == 1) { // If the locus is heterozigotic
         ++NLocHet;
      }
     Hplgens.first <<= 1; // shift value left by 1
     Hplgens.second <<= 1; // shift value left by 1
   }
   return(NLocHet);
}

// ---Utility functions---

int Organism::popcount64b(uint64_t x) const {
   x -= (x >> 1) & 0x5555555555555555;             //put count of each 2 bits into those 2 bits
   x = (x & 0x3333333333333333) + ((x >> 2) & 0x3333333333333333); //put count of each 4 bits into those 4 bits
   x = (x + (x >> 4)) & 0x0f0f0f0f0f0f0f0f;        //put count of each 8 bits into those 8 bits
   x += x >>  8;  //put count of each 16 bits into their lowest 8 bits
   x += x >> 16;  //put count of each 32 bits into their lowest 8 bits
   x += x >> 32;  //put count of each 64 bits into their lowest 8 bits
   return x & 0x7f;
}

double Organism::sumAlpha(int nbites) const {
   double alphaSum = 0;
   for (int i = 1; i <= nbites; ++i) { alphaSum += Alpha; }
   return(alphaSum);
}

void Organism::displayBits(uint32_t value) const {
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

pair<uint32_t, uint32_t> Organism::pair32Int (uint64_t value) const {
    return pair<uint32_t, uint32_t>((value << 32) >> 32, value >> 32);
}

uint32_t Organism::freeRecombination(uint32_t Hplgen1,uint32_t Hplgen2) const {
   uint32_t Randuint32 = Rng::random_uint32();
   // Step 1
   uint32_t A = Hplgen1 ^ Hplgen2;
   // Step2
   uint32_t B = Hplgen1 | A;
   // Step3
   uint32_t C = ( A & ( Randuint32 ^ 0xFFFFFFFF ) );
   // Step4
   uint32_t result = B ^ C;
   return(result);
}

