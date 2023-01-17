// Host derived class member function definitions

#include <cstdint> // uint32_t and uint64_t types
#include <string>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include "Host.h" // Host class definition
#include "Param.h"
using namespace std;

// ---Static member functions---
int Host::getKsymbiont() {return Ksymbiont;}
void Host::setKsymbiont( int ks ) { Ksymbiont = ks; }
double Host::getRmax() {return Rmax;}
void Host::setRmax( double rm ) { Rmax = rm; }
double Host::getVs() {return Vs;}
void Host::setVs( double vs ) { Vs = vs; }
double Host::getlambda() {return lambda;}
void Host::setlambda( double lmbd ) { lambda = lmbd; }
double Host::getb() {return b;}
void Host::setb( double binput ) { b = binput; }
double Host::getd() {return d;}
void Host::setd( double dinput ) { d = dinput; }
double Host::getmutRate() {return mutRate;}
void Host::setmutRate( double mtr ) { mutRate = mtr; }

// ---Static data members---
int Host::Ksymbiont = Param::getKsymbiont();
double Host::Rmax = Param::getRmaxH();
double Host::Vs = Param::getVsH();
double Host::lambda = Param::getlambda();
double Host::b = Param::getb();
double Host::d = Param::getd();
double Host::mutRate = Param::getmutRateH();

// constructor
Host::Host(char sx, uint64_t gn, double ph, int Ns, uint64_t id, uint64_t spopid)
   : Organism (sx, gn, ph) {
   setNsymbiont(Ns);
   setID(id);
   setSPopID(spopid);
}

void Host::setNsymbiont(int Ns) {Nsymbiont = Ns;}
int Host::getNsymbiont() const {return Nsymbiont;}

void Host::setID(uint64_t id) {ID = id;}
uint64_t Host::getID() const {return ID;}

void Host::setSPopID(uint64_t spopid) {SPopID = spopid;}
uint64_t Host::getSPopID() const {return SPopID;}

void Host::printID() const {
   auto pairID{pair32Int(ID)};
   cout << "Host ID: " << endl;
   cout <<  "Index\t";
   displayBits(pairID.first);
   cout <<  "Version\t";
   displayBits(pairID.second);
}

void Host::printIndividual() const {
    cout << toString() << "\n" << endl;
 }

vector<HostGamete> Host::produceGametes(const Patch& island, int n, Rng& rng) const {
   int Ngametes = calculateNGametes(island, n, rng);
   vector<HostGamete> vgametes;
   vgametes.reserve(Ngametes);
   for (int i = 0; i < Ngametes; ++i) {
      vgametes.push_back( createOneGamete() );
   }
   return(vgametes);
}

// ---Overloaded Operators---

Host& Host::operator++() {
   ++Nsymbiont;
   return *this;
   }

Host& Host::operator--() {
   --Nsymbiont;
   return *this;
   }

Host& Host::operator+=(int addedSymbionts) {
   Nsymbiont += addedSymbionts;
   return *this;
   }

Host& Host::operator-=(int removedSymbionts) {
   Nsymbiont -= removedSymbionts;
   return *this;
   }

// ---Utility functions---

void Host::displayBits(uint32_t value) const {
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

pair<uint32_t, uint32_t> Host::pair32Int (uint64_t value) const {
    return pair<uint32_t, uint32_t>((value << 32) >> 32, value >> 32);
}

// Get string representation of the object data members
string Host::toString() const {
   ostringstream output;
   output << Organism::toString()
      << "\nNsymbiont: " << getNsymbiont()
      << "\nK symbiont: " << getKsymbiont()
      << "\nRmax: " << getRmax()
      << "\nVs: " << getVs()
      << "\nlambda: " << getlambda()
      << "\nb: " << getb()
      << "\nID: " << getID()
      << "\t" << "Index: " << static_cast<uint32_t>((getID() << 32) >> 32)
      << "\t" << "Version: " << static_cast<uint32_t>(getID() >> 32)
      << "\nSPopID: " << getSPopID()
      << "\t" << "Index: " << static_cast<uint32_t>((getSPopID() << 32) >> 32)
      << "\t" << "Version: " << static_cast<uint32_t>(getSPopID() >> 32) << "\n\n";
   return output.str();
}

HostGamete Host::createOneGamete() const {
   HostGamete gamete;
   int32_t hplgen = createOneHaplGen();
   gamete.setHaplGen( hplgen );
   gamete.setHostID( getID() );
   return(gamete);
}

int Host::calculateNGametes(const Patch& island, int n, Rng& rng) const {
   // Getting data from hpop and Patch
   double N = static_cast<double>(n);
   double K = static_cast<double>( island.getKhost() );
   double OptPhen = island.getOptPhen();
   // Calculation of the expected mean number of gamete produced
      double term1 = Rmax * ( 1 - N / K );
      double term2 = ( getPhen() - OptPhen ) * ( getPhen() - OptPhen );
      double term3 = Vs + Vs;
      double term4 = lambda * ( Ksymbiont / ( 1 + Ksymbiont * exp ( -b * Nsymbiont ) ) );
      double R = term1 -  term2 / term3 + term4 ;
      // Sampling number of gametes
      double mean = (exp(R) + exp(R)) * d;
//cout << "\nmean = " << mean;
//cout << "\nterm1 = " << term1;
return ( rng.poisson(mean) );
}

