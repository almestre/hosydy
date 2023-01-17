// Output class definition

/* Variables included in output1.csv:
   1. Nhost = Host population size
   2. SPrev = Symbiont prevalence
   3. AvSAb = Mean symbiont abundance (i.e. average number of symbionts per occupied host)
   4. VmrSAb = Variance-to-mean ratio of symbiont abundances
   5. CvSAb = Coefficient of variation of symbiont abundances
   6. AvHPhen = Mean host phenotype value
   7. StdHPhen = Standard deviation of host phenotype values
   8. AvSPhen = Global Mean symbiont phen. val.
   9.SPhenMsB = Mean square between hosts
   10.SPhenMsW = Symbiont phenotypic variance within hosts (mean square within groups)
   11. HSPhenCor = Correlation between individual host phenotype and the mean phenotype of the symbionts inhabiting the host

   Variables included in output2.csv: Host allele frequencies
   Variables included in output3.csv: Symbiont allele frequencies

*/

#ifndef OUTPUT_H
#define OUTPUT_H

#include <cstdint> // uint32_t and uint64_t types
#include <string>
#include <fstream>

// Forward declarations:
class Host;
class Symbiont;
template<typename T>
class Population;
template<typename T>
class Metapopulation;


class Output {
public:

   // Within-class type definitions:
//   typedef Population<Host> HPop; // gamete pool
//   typedef Metapopulation<Population<Symbiont>> SMpop; // gamete pool

   Output(); // constructor

   // Static member functions

   static void createOutputFiles();
   static void printHeadersToFiles();
   static void printDataToFiles( Population<Host>&, const Metapopulation<Population<Symbiont>>& );
   static void printAlFreqToFiles( Population<Host>&, const Metapopulation<Population<Symbiont>>& );

   static void printHeader1( std::ofstream& );
   static void printHeader2( std::ofstream& );

   static void printOutput1( std::ofstream&, Population<Host>&, const Metapopulation<Population<Symbiont>>& );
   static void printOutput2( std::ofstream&, const Population<Host>& );
   static void printOutput3( std::ofstream&, const Metapopulation<Population<Symbiont>>& );

   static void setOutFreq( int ); // Set OutFreq
   static int getOutFreq(); // Get OutFreq

   static void setCurrSimStep( int ); // Set OutFreq
   static int getCurrSimStep(); // Get OutFreq

private:

   static int OutFreq; // Frequency of output generation (measured in symbiont cycles)
   static int CurrSimStep; // Current simulation step (measured in symbiont cycles)
   static std::ofstream output1;
   static std::ofstream output2;
   static std::ofstream output3;
   };

   #endif // OUTPUT_H
