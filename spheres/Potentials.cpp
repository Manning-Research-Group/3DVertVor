#include <iostream>
#include "Potentials.h"

#if USE_NETCDF
const std::string HarmonicPotential::NetCdfName("HarmonicPotential");
std::string HarmonicPotential::getNetCdfNameAndParameters() const { 
  std::string name(NetCdfName); 
  name += ":0x1";
  return name;
}

void HarmonicPotential::setParametersFromNetCdfString(const std::string &str) {
}

const std::string PowerLawPotential::NetCdfName("PowerLawPotential");
std::string PowerLawPotential::getNetCdfNameAndParameters() const { 
  std::string name(NetCdfName); 
  name += ":0x1";
  return name;
}

void PowerLawPotential::setParametersFromNetCdfString(const std::string &str) {
}

#endif