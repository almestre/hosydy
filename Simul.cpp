// Implementation of Simul class

#include <string>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Param.h"
#include "Output.h"
#include "Rng.h"
#include "Patch.h"
#include "SourcePatch.h"
#include "Organism.h"
#include "Symbiont.h"
#include "Host.h"
#include "Simul.h"
#include "Population.h"
#include "Metapopulation.h"

#include <cstdint> // uint32_t and uint64_t types
#include <string>

using namespace std;

// ---Static member functions---

void Simul::runSimul() {
   // Input parameters from JSON file:
//   string inputFile("input.JSON"); // old version
//   Param::inputParamFromJsonFile( inputFile ); // old version
   Param::inputParamFromJsonFile( ); // to run with linux terminal
   Param::initParam();
 // Initialize scenario
  // Create a rng object
   Rng rng;
   rng.Rng_init();
   rng.Rng_save(); // place after Rng_init to save the initial state
  // Create the continent and the island:
   SourcePatch continent;
   Patch island;
  // Create and initialise population of hosts associated to the island:
   Population<Host> hpop;
   hpop.initPop( island );
  // Create a symbiont metapopulation:
   Metapopulation<Population<Symbiont>> smpop;
   smpop.initMetapop();
  // First immigration from source:   
   hpop.initImmigrFromSource( continent, rng, smpop );
  // Create output files and print headers into the files:
   Output::createOutputFiles();
   Output::printHeadersToFiles();
 // Run cycles
   int Ncycles = NYears*NHReprPerYear*NStepsPerHRepr;
   for (int counter = 0; counter < Ncycles; ++counter) {
      Output::setCurrSimStep( counter );
      smpop.horizTrans( rng, hpop );
      smpop.metapopReproduction( hpop, rng );
      if (counter % NStepsPerHRepr == ( NStepsPerHRepr - 1 )) { // First host reproductive cycle at step 11
         hpop.popReproduction(island, rng, smpop);
         hpop.popMortality(rng,smpop);
         if ( counter < Ncycles-1200 ) { hpop.immigrFromSource( continent, rng, smpop ); } // last 100 years without host migration
      }
//      if ( counter >= Ncycles-2400 ) { // output only for the last 200 years without host migration
         if ( counter % NStepsPerHRepr == NStepsPerHRepr - 1 ) { Output::printDataToFiles( hpop, smpop ); } // Output frequency = 1 year (starting from the end of step 0)
//      }
   }
}

void Simul::setNStepsPerHRepr( int nsphr ) { NStepsPerHRepr = nsphr; }
int Simul::getNStepsPerHRepr() {return NStepsPerHRepr;}

void Simul::setNHReprPerYear( int nhrpy ) { NHReprPerYear = nhrpy; }
int Simul::getNHReprPerYear() {return NHReprPerYear;}

void Simul::setNYears( int nyears ) { NYears = nyears ; }
int Simul::getNYears() {return NYears;}

void Simul::setScenID( string scenid ) { ScenID = scenid ; }
string Simul::getScenID() {return ScenID;}

void Simul::setReplID( int replid ) { ReplID = replid ; }
int Simul::getReplID() {return ReplID;}

// ---Static data members---
int Simul::NStepsPerHRepr = Param::getNStepsPerHRepr();
int Simul::NHReprPerYear = Param::getNHReprPerYear();
int Simul::NYears = Param::getNYears();
string Simul::ScenID = Param::getScenID();
int Simul::ReplID = Param::getReplID();

// constructor

Simul::Simul() {}
