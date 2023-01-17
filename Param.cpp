// Implementation of Param class

#include <string>
#include <iomanip>
#include <stdexcept>
#include <iostream>
#include <fstream>
#include <sstream>
#include "Param.h" // Host class definition
#include "json.hpp"
#include "Organism.h"
#include "Symbiont.h"
#include "Host.h"
#include "Simul.h"
#include "Output.h"

using json = nlohmann::json;
using namespace std;

// ---Static member functions---

void Param::inputParamFromJsonFile( ) {
   // Read input parameter values from a JSON file
   json inputData;
   string scenid;
   string replid;
   cin >> scenid; // get the ID of the simulation scenario (run with linux command line)
   cin >> replid; // get the ID of the simulation replicate (run with linux command line)
   setScenID( scenid );
   setReplID( stoi(replid) ); // the function stoi() converts a string into an int
   ifstream file ( "input/input" + scenid + ".JSON" );
   file >> inputData;
   // Set parameter values
//   setL( inputData[ "L" ].get<int>() );
//   setAlpha( inputData[ "Alpha" ].get<double>() );
//   setKsymbiont( inputData[ "Ksymbiont" ].get<int>() );
//   setRmaxH( inputData[ "RmaxH" ].get<double>() );
   setVsH( inputData[ "VsH" ].get<double>() ) ;
   setlambda( inputData[ "lambda" ].get<double>() );
//   setb( inputData[ "b" ].get<double>() );
//   setd( inputData[ "d" ].get<double>() );
//   setRmaxS( inputData[ "RmaxS" ].get<double>() );
   setVsS( inputData[ "VsS" ].get<double>() );
   setEht( inputData[ "Eht" ].get<double>() );
//   setEvt( inputData[ "Evt" ].get<double>() );
//   setKhost( inputData[ "Khost" ].get<int>() );
   setOptPhen( inputData[ "OptPhen" ].get<double>() );
//   setOptPhenS( inputData[ "OptPhenS" ].get<double>() );
//   setStdPhenH( inputData[ "StdPhenH" ].get<double>() );
//   setStdPhenS( inputData[ "StdPhenS" ].get<double>() );
//   setLocStdPhenS( inputData[ "LocStdPhenS" ].get<double>() );
//   setPrPulseMigr( inputData[ "PrPulseMigr" ].get<double>() );
   setHNmigrants( inputData[ "HNmigrants" ].get<double>() );
//   setPSLocAd( inputData[ "PSLocAd" ].get<double>() );
   setSPrev( inputData[ "SPrev" ].get<double>() );
//   setSAb( inputData[ "SAb" ].get<double>() );
//   setSTheta( inputData[ "STheta" ].get<double>() );
//   setHPopVecSize( inputData[ "HPopVecSize" ].get<size_t>() );
//   setSPopVecSize( inputData[ "SPopVecSize" ].get<size_t>() );
//   setNStepsPerHRepr( inputData[ "NStepsPerHRepr" ].get<int>() );
//   setNHReprPerYear( inputData[ "NHReprPerYear" ].get<int>() );
//   setNYears( inputData[ "NYears" ].get<int>() );
//   setOutFreq( inputData[ "OutFreq" ].get<int>() );
//   setmutRateH( inputData[ "mutRateH" ].get<double>() );
//   setmutRateS( inputData[ "mutRateS" ].get<double>() );
}

void Param::inputParamFromJsonFile( string& path ) {
   // Read input parameter values from a JSON file
   json inputData;
   ifstream file ( path );
   file >> inputData;
   // Set parameter values
//   setL( inputData[ "L" ].get<int>() ); // default
//   setAlpha( inputData[ "Alpha" ].get<double>() ); // default
//   setKsymbiont( inputData[ "Ksymbiont" ].get<int>() ); // default
//   setRmaxH( inputData[ "RmaxH" ].get<double>() ); // default
   setVsH( inputData[ "VsH" ].get<double>() ) ;
   setlambda( inputData[ "lambda" ].get<double>() );
//   setb( inputData[ "b" ].get<double>() ); // default
//   setd( inputData[ "d" ].get<double>() ); // default
//   setRmaxS( inputData[ "RmaxS" ].get<double>() ); // default
   setVsS( inputData[ "VsS" ].get<double>() ); // default
   setEht( inputData[ "Eht" ].get<double>() );
//   setEvt( inputData[ "Evt" ].get<double>() ); // default
//   setKhost( inputData[ "Khost" ].get<int>() ); // default
   setOptPhen( inputData[ "OptPhen" ].get<double>() );
//   setOptPhenS( inputData[ "OptPhenS" ].get<double>() ); // default
//   setStdPhenH( inputData[ "StdPhenH" ].get<double>() ); // default
//   setStdPhenS( inputData[ "StdPhenS" ].get<double>() ); // default
//   setLocStdPhenS( inputData[ "LocStdPhenS" ].get<double>() );
//   setPrPulseMigr( inputData[ "PrPulseMigr" ].get<double>() );
   setHNmigrants( inputData[ "HNmigrants" ].get<double>() );
//   setPSLocAd( inputData[ "PSLocAd" ].get<double>() ); // default
   setSPrev( inputData[ "SPrev" ].get<double>() );
//   setSAb( inputData[ "SAb" ].get<double>() ); // default
//   setSTheta( inputData[ "STheta" ].get<double>() ); // default
//   setHPopVecSize( inputData[ "HPopVecSize" ].get<size_t>() ); // default
//   setSPopVecSize( inputData[ "SPopVecSize" ].get<size_t>() ); // default
//   setNStepsPerHRepr( inputData[ "NStepsPerHRepr" ].get<int>() ); // default
//   setNHReprPerYear( inputData[ "NHReprPerYear" ].get<int>() ); // default
//   setNYears( inputData[ "NYears" ].get<int>() ); // default
//   setOutFreq( inputData[ "OutFreq" ].get<int>() ); // default
//   setmutRateH( inputData[ "mutRateH" ].get<double>() ); // default
//   setmutRateS( inputData[ "mutRateS" ].get<double>() ); // default
}

void Param::initParam() {
   Organism::setL( getL() );
   Organism::setAlpha( getAlpha() );
   Host::setKsymbiont( getKsymbiont() );
   Host::setRmax( getRmaxH() );
   Host::setVs( getVsH() );
   Host::setlambda( getlambda() );
   Host::setb( getb() );
   Host::setd( getd() );
   Symbiont::setRmax( getRmaxS() );
   Symbiont::setVs( getVsS() );
   Symbiont::setEht( getEht() );
   Symbiont::setEvt( getEvt() );
   Simul::setNStepsPerHRepr( getNStepsPerHRepr() );
   Simul::setNHReprPerYear( getNHReprPerYear() );
   Simul::setNYears( getNYears() );
   Simul::setScenID( getScenID() );
   Simul::setReplID( getReplID() );
   Output::setOutFreq( getOutFreq() );
   Host::setmutRate( getmutRateH() );
   Symbiont::setmutRate( getmutRateS() );
}

void Param::setL( int l ) { L = l; }
int Param::getL() {return L;};

void Param::setAlpha( double alpha ) { Alpha = alpha ; }
double Param::getAlpha() {return Alpha;}

void Param::setKsymbiont( int ks ) { Ksymbiont = ks ; }
int Param::getKsymbiont() {return Ksymbiont;}

void Param::setRmaxH( double rmh ) { RmaxH = rmh ; }
double Param::getRmaxH() {return RmaxH;}

void Param::setVsH( double vsh ) { VsH = vsh ; }
double Param::getVsH() {return VsH;}

void Param::setlambda( double lmbd ) { lambda = lmbd ; }
double Param::getlambda() {return lambda;}

void Param::setb( double binput ) { b = binput ; }
double Param::getb() {return b;}

void Param::setd( double dinput ) { d = dinput ; }
double Param::getd() {return d;}

void Param::setRmaxS( double rms ) { RmaxS = rms ; }
double Param::getRmaxS() {return RmaxS;}

void Param::setVsS( double vss ) { VsS = vss ; }
double Param::getVsS() {return VsS;}

void Param::setEht( double eht ) { Eht = eht ; }
double Param::getEht() {return Eht;}

void Param::setEvt( double evt ) { Evt = evt ; }
double Param::getEvt() {return Evt;}

void Param::setKhost( int kh ) { Khost = kh ; }
int Param::getKhost() {return Khost;}

void Param::setOptPhen( double op ) { OptPhen = op ; }
double Param::getOptPhen() {return OptPhen;}

void Param::setOptPhenS( double ops ) { OptPhenS = ops ; }
double Param::getOptPhenS() {return OptPhenS;}

void Param::setStdPhenH( double stdph ) { StdPhenH = stdph ; }
double Param::getStdPhenH() {return StdPhenH;}

void Param::setStdPhenS( double stdps ) { StdPhenS = stdps ; }
double Param::getStdPhenS() {return StdPhenS;}

void Param::setLocStdPhenS( double locstdps ) { LocStdPhenS = locstdps ; }
double Param::getLocStdPhenS() {return LocStdPhenS;}

void Param::setPrPulseMigr( double hnm ) { PrPulseMigr = hnm ; }
double Param::getPrPulseMigr() {return PrPulseMigr;}

void Param::setHNmigrants( double hnm ) { HNmigrants = hnm ; }
double Param::getHNmigrants() {return HNmigrants;}

void Param::setPSLocAd( double pslocad ) { PSLocAd = pslocad ; }
double Param::getPSLocAd() {return PSLocAd;}

void Param::setSPrev( double sprev ) { SPrev = sprev ; }
double Param::getSPrev() {return SPrev;}

void Param::setSAb( double sab ) { SAb = sab ; }
double Param::getSAb() {return SAb;}

void Param::setSTheta( double stheta ) { STheta = stheta ; }
double Param::getSTheta() {return STheta;}

void Param::setHPopVecSize( size_t hpops ) { HPopVecSize = hpops ; }
size_t Param::getHPopVecSize() {return HPopVecSize;}

void Param::setSPopVecSize( size_t spops ) { SPopVecSize = spops ; }
size_t Param::getSPopVecSize() {return SPopVecSize;}

void Param::setNStepsPerHRepr( int nsphr ) { NStepsPerHRepr = nsphr ; }
int Param::getNStepsPerHRepr() {return NStepsPerHRepr;}

void Param::setNHReprPerYear( int nhrpy ) { NHReprPerYear = nhrpy ; }
int Param::getNHReprPerYear() {return NHReprPerYear;}

void Param::setNYears( int nyears ) { NYears = nyears ; }
int Param::getNYears() {return NYears;}

void Param::setScenID( string scenid ) { ScenID = scenid ; }
string Param::getScenID() {return ScenID;}

void Param::setReplID( int replid ) { ReplID = replid ; }
int Param::getReplID() {return ReplID;}

void Param::setOutFreq( int outfreq ) { OutFreq = outfreq ; }
int Param::getOutFreq() {return OutFreq;}

void Param::setmutRateH( double mtrh ) { mutRateH = mtrh; }
double Param::getmutRateH() {return mutRateH;}

void Param::setmutRateS( double mtrs ) { mutRateS = mtrs; }
double Param::getmutRateS() {return mutRateS;}

   // Static data members
int Param::L = 32;
double Param::Alpha = 0.15625;
int Param::Ksymbiont = 200;
double Param::RmaxH = 1;
double Param::VsH = 20;
double Param::lambda = 0;
double Param::b = 0.05;
double Param::d = 0.25;
double Param::RmaxS = 1;
double Param::VsS = 20;
double Param::Eht = 0;
double Param::Evt = 0.2;
int Param::Khost = 1000;
double Param::OptPhen =  -5;
double Param::OptPhenS =  -5;
double Param::StdPhenH = 1;
double Param::StdPhenS = 2;
double Param::LocStdPhenS = 0.2;
double Param::PrPulseMigr = 1;
double Param::HNmigrants = 100;
double Param::PSLocAd = 1;
double Param::SPrev = 1;
double Param::SAb = 180;
double Param::STheta = 10;
size_t Param::HPopVecSize = 4000; // 1200 for parasite and commensal; 2000 for mutualist
size_t Param::SPopVecSize = 1000;
int Param::NStepsPerHRepr = 12;
int Param::NHReprPerYear = 1;
int Param::NYears = 1100;
string Param::ScenID = "";
int Param::ReplID = 0;
int Param::OutFreq = 120;
double Param::mutRateH = 0;
double Param::mutRateS = 0;

// constructor

Param::Param() {}


