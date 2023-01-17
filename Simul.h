// Simul class definition

#ifndef SIMUL_H
#define SIMUL_H

#include <cstdint> // uint32_t and uint64_t types
#include <string>

class Simul {
public:
   explicit Simul();

   // Static member functions

   static void runSimul(); // Run the simulation

   static void setNStepsPerHRepr( int ); // Set NStepsPerHRepr
   static int getNStepsPerHRepr(); // Get NStepsPerHRepr
   
   static void setNHReprPerYear( int ); // Set NHReprPerYear
   static int getNHReprPerYear(); // Get NHReprPerYear
   
   static void setNYears( int ); // Set NYears
   static int getNYears(); // Get NYears

   static void setScenID( std::string ); // Set ScenID
   static std::string getScenID(); // Get ScenID

   static void setReplID( int ); // Set NYears
   static int getReplID(); // Get NYears

private:
   
   static int NStepsPerHRepr; // Number of symbiont cycles (i.e. steps) per host reproductive event
   static int NHReprPerYear; // Number of host reproductive events per year 
   static int NYears; // Number of years simulated
   static std::string ScenID; // ID of the simulation scenario
   static int ReplID; // ID of the simulation replicate
   };

#endif // SIMUL_H
