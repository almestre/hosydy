// #### HOSYDY_test61
// Date 18/08/2021
// Code for testing the new Organism class' member function "int sumHetLocInd ()";

#include <typeinfo>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <vector>
#include <string>
#include <chrono>

#include "Rng.h"
#include "SourcePatch.h"
#include "Patch.h"
#include "Population.h"
#include "Metapopulation.h"
#include "Gamete.h"
#include "Param.h"
#include "Simul.h"

using namespace std;
using namespace std::chrono;

int main() {

auto start = high_resolution_clock::now();

   // Init parameters from JSON file:
//   string inputFile("input.JSON");
//   Param::inputParamFromJsonFile( inputFile );
//   Param::inputParamFromJsonFile( );
//   Param::initParam();
   // Run the simulation:
Simul::runSimul();

auto stop = high_resolution_clock::now();

auto duration = duration_cast<seconds>(stop - start);
cout << duration.count() << " seconds" << endl;
}
