// Symbiont derived class member function definitions

#include <cstdint> // uint32_t and uint64_t types
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <sstream>
#include <cmath>
#include "Symbiont.h" // Symbiont class definition
#include "Host.h" // Forward declaration in the header
#include "Population.h" // Forward declaration in the header
#include "Param.h"

using namespace std;

// ---Static member functions---

double Symbiont::getRmax() {return Rmax;}
void Symbiont::setRmax( double rm ) { Rmax = rm; }

double Symbiont::getVs() {return Vs;}
void Symbiont::setVs( double vs ) { Vs = vs; }

double Symbiont::getEht() {return Eht;}
void Symbiont::setEht( double eht ) { Eht = eht; }

double Symbiont::getEvt() {return Evt;}
void Symbiont::setEvt( double evt ) { Evt = evt; }

double Symbiont::getmutRate() {return mutRate;}
void Symbiont::setmutRate( double mtr ) { mutRate = mtr; }

// ---Static data members---
double Symbiont::Rmax = Param::getRmaxS();
double Symbiont::Vs = Param::getVsS();
double Symbiont::Eht = Param::getEht();
double Symbiont::Evt = Param::getEvt();
double Symbiont::mutRate = Param::getmutRateS();

// constructor
Symbiont::Symbiont(char sx, uint64_t gn, double ph)
   : Organism (sx, gn, ph) {}

void Symbiont::printIndividual() const {
   cout << Organism::toString() << "\n" << endl;
}

vector<Gamete> Symbiont::produceGametes(int64_t hostid, const Population<Host>& hpop, Rng& rng)  const {
   int Ngametes = calculateNGametes(hostid, hpop, rng);
//cout << "\nNgametes = " << Ngametes << "\n";
   vector<Gamete> vgametes;
   vgametes.reserve(Ngametes);
   for (int i = 0; i < Ngametes; ++i) {
      vgametes.push_back( createOneGamete() );
   }
   return(vgametes);
}

// ---Utility functions---

Gamete Symbiont::createOneGamete() const {
   Gamete gamete;
   int32_t hplgen = createOneHaplGen();
   gamete.setHaplGen( hplgen );
   return(gamete);
}

int Symbiont::calculateNGametes(int64_t hostid, const Population<Host>& hpop, Rng& rng) const {
   // Getting data from Host
   const Host* ptrHost = hpop.getConstInd( hostid ); // hpop is the name of the host population
   double N = static_cast<double>(ptrHost->getNsymbiont());
   double K = static_cast<double>( Host::getKsymbiont() );
   double HostPhen = ptrHost->getPhen();
   ptrHost = nullptr;
   // Calculation of the expected mean number of gamete produced
      double term1 = Rmax * ( 1 - N / K );
      double term2 = ( getPhen() - HostPhen ) * ( getPhen() - HostPhen );
      double term3 = Vs + Vs;
      double R = term1 - term2 / term3;
      // Sampling number of gametes
      double mean = exp(R) + exp(R);
   return ( rng.poisson(mean) ); // rng.poisson(mean)
}


