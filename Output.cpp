// Implementation of Output class

#include <cstdint> // uint32_t and uint64_t types
#include <string>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Output.h"
#include "Param.h"
#include "Metapopulation.h"
#include "Population.h"
#include "Organism.h"
#include "Host.h"
#include "Symbiont.h"
#include "Simul.h"

using namespace std;

// ---Static member functions---

void Output::createOutputFiles(  ) {
  std::string replid = std::to_string( Simul::getReplID() );
  output1.open("output_full/output" + Simul::getScenID() + "_" + replid + ".csv");
  if( !output1 ) { // file couldn't be opened
      cerr << "Error: file output1 could not be opened" << endl;
      exit(1);
   }
// 25/08/21: we don't need MAFS at this moments (heterozigosity metrics will be enough)
//  output2.open("output/outputMAFH" + Simul::getScenID() + "_" + replid + ".csv");
//  if( !output2 ) { // file couldn't be opened
//      cerr << "Error: file output2 could not be opened" << endl;
//      exit(1);
//   }
//  output3.open("output/outputMAFS" + Simul::getScenID() + "_" + replid + ".csv");
//  if( !output3 ) { // file couldn't be opened
//      cerr << "Error: file output3 could not be opened" << endl;
//      exit(1);
//   }
}

void Output::printHeadersToFiles(  ) {
   printHeader1( output1 );
//   printHeader2( output2 );
//   printHeader2( output3 );
}

void Output::printDataToFiles( Population<Host>& hpp, const Metapopulation<Population<Symbiont>>& smpp ) {
   printOutput1( output1, hpp, smpp );
}

void Output::printAlFreqToFiles( Population<Host>& hpp, const Metapopulation<Population<Symbiont>>& smpp ) {
   printOutput2( output2, hpp );
   printOutput3( output3, smpp );
}

void Output::printHeader1( ofstream& outf ) {

   outf << "Scenario ID" << ","
      << "Replicate ID" << ","
      << "Year" << ","
      << "Nhost" << ","
      << "SPrev" << ","
      << "AvSAb" << ","
      << "VmrSAb" << ","
      << "CvSAb" << ","
      << "AvHPhen" << ","
      << "StdHPhen" << ","
      << "Hobs_host" << ","
      << "Hexp_host" << ","
      << "AvSPhen" << ","
      << "MsbSPhen" << ","
      << "MswSPhen" << ","
      << "HSPhenCor" << ","
      << "Hobs_symb" << ","
      << "Hexp_symb" << ","
      << "EvCS" << "\n";
}

void Output::printHeader2( ofstream& outf ) {
   outf << "Scenario ID" << ","
      << "Replicate ID" << ","
      << "Year" << ",";
   for (int i = 1; i < ( Param::getL() ); ++i ) {
      outf << "Al" << i << ",";
   }
   outf << "Al" << Param::getL() << "\n"; // last entry
}

void Output::printOutput1( ofstream& outf, Population<Host>& hpp, const Metapopulation<Population<Symbiont>>& smpp ) {
   outf << Simul::getScenID() << ","
      << Simul::getReplID() << ","    
      << (CurrSimStep/12)+1 << ","
      << hpp.getN() << ","
      << hpp.getSPrev() << ","
      << hpp.getAvSAbOcH() << ","
      << hpp.getVmrSAb() << ","
      << hpp.getCvSAb() << ","
      << hpp.getAvHPhen() << ","
      << hpp.getStdHPhen() << ","
      << hpp.getHobs() << ","
      << hpp.getHexp() << ","
      << smpp.getGAvPhen() << ","
      << smpp.getMsbPhen() << ","
      << smpp.getMswPhen() << ","
      << smpp.getHSPhenCor( hpp ) << ","
      << smpp.getHobs() << ","
      << smpp.getHexp() << ","
      << smpp.getEvC( hpp ) << "\n";
}

void Output::printOutput2( ofstream& outf, const Population<Host>& hpp ) {
   std::vector<int> alfreq = hpp.getAlFreq();
   outf << Simul::getScenID() << ","
        << Simul::getReplID() << ","    
        << (CurrSimStep/12)+1 << ",";
   for (int counter = 0; counter < (static_cast<int>(alfreq.size())-1); ++counter ) {
      outf << alfreq[counter] << ",";
   }
      outf << alfreq.back() << "\n"; // last entry
}

void Output::printOutput3( ofstream& outf, const Metapopulation<Population<Symbiont>>& smpp ) {
   std::vector<int> alfreq = smpp.getAlFreq();
   outf << Simul::getScenID() << ","
        << Simul::getReplID() << ","    
        << (CurrSimStep/12)+1 << ",";
   for (int counter = 0; counter < (static_cast<int>(alfreq.size())-1); ++counter ) {
      outf << alfreq[counter] << ",";
   }
      outf << alfreq.back() << "\n"; // last entry
}

void Output::setOutFreq( int outfreq ) { OutFreq = outfreq ; }
int Output::getOutFreq() {return OutFreq;}

void Output::setCurrSimStep( int cursimstep ) { CurrSimStep = cursimstep ; }
int Output::getCurrSimStep() {return CurrSimStep;}

// ---Static data members---
int Output::OutFreq = Param::getOutFreq();
int Output::CurrSimStep = 0;
ofstream Output::output1;
ofstream Output::output2;
ofstream Output::output3;

// ---Constructor---

Output::Output() {}
