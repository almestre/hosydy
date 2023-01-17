// HostGamete derived class definition

#ifndef HOSTGAMETE_H
#define HOSTGAMETE_H

#include <cstdint> // uint32_t and uint64_t types
#include <string>
#include "Gamete.h"

class HostGamete : public Gamete  {
public:
   explicit HostGamete(uint32_t = 0, uint64_t = 0);

   void setHostID(uint64_t);
   uint64_t getHostID() const;

   void printHostID() const;

private:
   uint64_t HostID; // genotype

   // Utility functions
   void displayBits(uint32_t) const;
   std::pair<uint32_t, uint32_t> pair32Int (uint64_t) const;

   };

   #endif // HOSTGAMETE_H
