// Patch class definition

#ifndef PATCH_H
#define PATCH_H

// #include <string> // C++ standard string class

#include "Param.h"

class Patch {
public:
   explicit Patch( int=Param::getKhost(), double=Param::getOptPhen(), uint64_t = 0, uint64_t = 0); // constructor

   // Non-static member functions

   void setKhost(int); // set patch carrying capacity for hosts
   int getKhost() const; // get patch carrying capacity for hosts

   void setOptPhen(double); // set patch optimal phenotype for hosts
   double getOptPhen() const; // get patch optimal phenotype for hosts

   void setID(uint64_t); // set patch id
   uint64_t getID() const; // get patch id

   void setHPopID(uint64_t); // set id of the host population inhabiting the patch
   uint64_t getHPopID() const; // get id of the host population inhabiting the patch

   void printPatch(); // Print patch data members

private:

   int Khost; // patch carrying capacity for hosts
   double OptPhen; // patch optimal phenotype for hosts
   uint64_t ID; // patch id
   uint64_t HPopID; // id of the host population inhabiting the patch 

   // Utility functions
   std::string toString() ; // Print string representation of the object data members

   };



   #endif
