// Param base class definition

#ifndef PARAM_H
#define PARAM_H

class Param {
public:
   Param();

   // Static member functions

   static void inputParamFromJsonFile( ); // Input parameters from a JSON file ( get the path from an input file using linux command lines )
   static void inputParamFromJsonFile( std::string& ); // Input parameters from a JSON file

   static void initParam(); // Initialize parameters represented by static data members in Organism, Host and Symbiont classes

   static void setL( int ); // Set L
   static int getL(); // Get L

   static void setAlpha( double ); // Set Alpha
   static double getAlpha(); // Get Alpha

   static void setKsymbiont( int ); // Set Ksymbiont
   static int getKsymbiont(); // Get Ksymbiont

   static void setRmaxH( double ); // Set RmaxH
   static double getRmaxH(); // Get RmaxH

   static void setVsH( double ); // Set VsH
   static double getVsH(); // Get VsH

   static void setlambda( double ); // Set lambda
   static double getlambda(); // Get lambda

   static void setb( double ); // Set b
   static double getb(); // Get b

   static void setd( double ); // Set d
   static double getd(); // Get d

   static void setRmaxS( double ); // Set RmaxS
   static double getRmaxS(); // Get RmaxS

   static void setVsS( double ); // Set VsS
   static double getVsS(); // Get VsS

   static void setEht( double ); // Set Eht
   static double getEht(); // Get  Eht

   static void setEvt( double ); // Set Evt
   static double getEvt(); // Get Evt

   static void setKhost( int ); // Set Khost
   static int getKhost(); // Get Khost

   static void setOptPhen( double ); // Set OptPhen
   static double getOptPhen(); // Get OptPhen

   static void setOptPhenS( double ); // Set OptPhenS
   static double getOptPhenS(); // Get OptPhenS

   static void setStdPhenH( double ); // Set StdPhenH
   static double getStdPhenH(); // Get StdPhenH

   static void setStdPhenS( double ); // Set StdPhenS
   static double getStdPhenS(); // Get StdPhenS

   static void setLocStdPhenS( double ); // Set LocStdPhenS
   static double getLocStdPhenS(); // Get LocStdPhenS

   static void setPrPulseMigr( double ); // Set PrPulseMigr
   static double getPrPulseMigr(); // Get PrPulseMigr

   static void setHNmigrants( double ); // Set HNmigrants
   static double getHNmigrants(); // Get HNmigrants

   static void setPSLocAd( double ); // Set PSLocAd
   static double getPSLocAd(); // Get PSLocAd

   static void setSPrev( double ); // Set SPrev
   static double getSPrev(); // Get SPrev

   static void setSAb( double ); // Set SAb
   static double getSAb(); // Get SAb

   static void setSTheta( double ); // Set STheta
   static double getSTheta(); // Get STheta

   static void setHPopVecSize( size_t ); // Set HPopVecSize
   static size_t getHPopVecSize(); // Get HPopVecSize

   static void setSPopVecSize( size_t ); // Set SPopVecSize
   static size_t getSPopVecSize(); // Get SPopVecSize

   static void setNStepsPerHRepr( int ); // Set NStepsPerHRepr
   static int getNStepsPerHRepr(); // Get NStepsPerHRepr

   static void setNHReprPerYear( int ); // Set NHReprPerYear
   static int getNHReprPerYear(); // Get NHReprPerYear

   static void setNYears( int ); // Set NYears
   static int getNYears(); // Get NYears

   static void setOutFreq( int ); // Set OutFreq
   static int getOutFreq(); // Get OutFreq

   static void setScenID( std::string ); // Set ScenID
   static std::string getScenID(); // Get ScenID

   static void setReplID( int ); // Set ReplID
   static int getReplID(); // Get ReplID

   static double getmutRateH(); // Get mutRateH
   static void setmutRateH( double ); // Set mutRateH

   static double getmutRateS(); // Get mutRateS
   static void setmutRateS( double ); // Set mutRateS

private:

   static int L; // Number of bi-allelic loci per genotype (the same for both host and symbiont)
   static double Alpha; // Effect size for each allele (the same for both host and symbiont)
   static int Ksymbiont; // Carrying capacity of a symbiont population inhabiting a host
   static double RmaxH; // Host intrinsic population growth rate
   static double VsH; // Width of stabilizing selection in hosts
   static double lambda; // Strength of effects of symbiont load on host fitness
   static double b; // Growth rate of the interaction effect size of symbiont load on host fitness
   static double d; // Host death rate per capita per reproductive event
   static double RmaxS; // Symbiont intrinsic population growth rate
   static double VsS; // Width of stabilizing selection
   static double Eht; // Per capita per generation emigr. rate for symbiont horiz. transm.
   static double Evt; // Per capita per generation emigr. rate for symbiont vert. transm.
   static int Khost; // Patch carrying capacity for hosts
   static double OptPhen; // Optimal phenotype for hosts in the patch (island)
   static double OptPhenS; // Optimal phenotype for hosts in the source patch (continent)
   static double StdPhenH; // Standard deviation of host phenotype values in the source population
   static double StdPhenS; // Global standard deviation of symbiont phenotype values in the source metapopulation
   static double LocStdPhenS; // Local standard deviation of phenotypes for symbionts that are locally adapted to an immigrant host
   static double PrPulseMigr;  // Probability that a migration pulse occurs at the end of a host cycle
   static double HNmigrants; // Average number of migrants generated by the source per migratory pulse
   static double PSLocAd; // Probability that a symbiont inhabiting a host immigrant is locally adapted
   static double SPrev; // Average prevalence of symbionts in host migrants from the source
   static double SAb; // Average abundance of symbionts in host migrants from the source (without considering empty hosts)
   static double STheta; // Dispersion paramater of symbiont abundances
   static size_t HPopVecSize; // Size of the vector that stores a host population; also represents the size of the vector that stores a symbiont metapopulation
   static size_t SPopVecSize; // Size of vectors that store a symbiont population
   static int NStepsPerHRepr; // Number of symbiont cycles (i.e. steps) per host reproductive event
   static int NHReprPerYear; // Number of host reproductive events per year
   static int NYears; // Number of years simulated
   static int OutFreq; // Frequency of output generation (in symbiont cycles)
   static std::string ScenID; // ID of the simulation (combination of characters and numbers)
   static int ReplID; // ID of the simulation
   static double mutRateH; // Per allele, per generation mutation rate in hosts
   static double mutRateS; // Per allele, per generation mutation rate in symbionts
   };

   #endif // PARAM_H
