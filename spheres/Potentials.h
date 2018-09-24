#ifndef POTENTIALS_H
#define POTENTIALS_H

#include "math.h"
#include <string>

class HarmonicPotential {
public:
  HarmonicPotential() {};
#if USE_NETCDF
  static const std::string NetCdfName;
  std::string getNetCdfNameAndParameters() const;
  void setParametersFromNetCdfString(const std::string &str);
#endif
  double energy(double relativeOverlap) { return 0.5*relativeOverlap*relativeOverlap; }
  double firstDerivative(double relativeOverlap) { return relativeOverlap; }
  double secondDerivative(double relativeOverlap) { return 1.0; }
};

class PowerLawPotential {
public:
  PowerLawPotential(double exponent) : _exponent(exponent) {};
#if USE_NETCDF
  static const std::string NetCdfName;
  std::string getNetCdfNameAndParameters() const;
  void setParametersFromNetCdfString(const std::string &str);
#endif
  double energy(double relativeOverlap) { return exp(_exponent*log(relativeOverlap))/_exponent; }
  double firstDerivative(double relativeOverlap) { return exp((_exponent-1)*log(relativeOverlap)); }
  double secondDerivative(double relativeOverlap) { return (_exponent-1)*exp((_exponent-2)*log(relativeOverlap)); }
private:
  double _exponent;
};

#endif /* POTENTIALS_H */

