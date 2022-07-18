#ifndef CGMINIMIZER_SEARCHDIRECTIONUPDATES_H
#define CGMINIMIZER_SEARCHDIRECTIONUPDATES_H

#include <sstream>
#include <Eigen/Dense>
#include "MinimizerWithDerivative.h"

class FletcherReeves {
public:
  static double firstResetInterval;
  static double subsequentResetInterval;

  FletcherReeves(const Eigen::VectorXd &gradient) : _NumberOfDofs(gradient.rows()), 
      _nextResetInterval(firstResetInterval*_NumberOfDofs), _iterationsSinceLastReset(0), 
      _p(-gradient) {}
  void update(double alpha, const Eigen::VectorXd &newGradient, const Eigen::VectorXd &oldGradient) {
    // sufficient progress?
    double beta = newGradient.squaredNorm()/oldGradient.squaredNorm();
    if(_iterationsSinceLastReset>_nextResetInterval) {
      _p = -newGradient;
      _iterationsSinceLastReset = 0;
      _nextResetInterval = subsequentResetInterval*_NumberOfDofs;
    } else {
      _p *= beta;
      _p -= newGradient;
    }
    ++_iterationsSinceLastReset;
  }
  std::string iterationAdditionalLogString() const {
    std::stringstream ss;
    ss << ",  iterationsSinceLastReset = " << _iterationsSinceLastReset;
    return ss.str();
  }
  Eigen::VectorXd searchDirectionUnitVector() const { return _p/_p.norm(); }
private:
  const long _NumberOfDofs;
  long _nextResetInterval;
  long _iterationsSinceLastReset;
  Eigen::VectorXd _p;
};


class PolakRibiere {
public:
  static double firstResetInterval;
  static double subsequentResetInterval;

  PolakRibiere(const Eigen::VectorXd &gradient) : _NumberOfDofs(gradient.rows()), 
      _nextResetInterval(firstResetInterval*_NumberOfDofs), _iterationsSinceLastReset(0), 
      _p(-gradient) {}
  void update(double alpha, const Eigen::VectorXd &newGradient, const Eigen::VectorXd &oldGradient) {
    // sufficient progress?
    double beta = (newGradient-oldGradient).dot(newGradient)/oldGradient.squaredNorm();
    if(_iterationsSinceLastReset>_nextResetInterval) {
      _p = -newGradient;
      _iterationsSinceLastReset = 0;
      _nextResetInterval = subsequentResetInterval*_NumberOfDofs;
    } else {
      _p *= beta;
      _p -= newGradient;
    }
    ++_iterationsSinceLastReset;
  }
  std::string iterationAdditionalLogString() const {
    std::stringstream ss;
    ss << ",  iterationsSinceLastReset = " << _iterationsSinceLastReset;
    return ss.str();
  }
  Eigen::VectorXd searchDirectionUnitVector() const { return _p/_p.norm(); }
private:
  const long _NumberOfDofs;
  long _nextResetInterval;
  long _iterationsSinceLastReset;
  Eigen::VectorXd _p;
};


class MemorylessBfgs {
public:
  MemorylessBfgs(const Eigen::VectorXd &gradient) : _p(-gradient), _searchDirection(_p/_p.norm()) {}
  void update(double alpha, const Eigen::VectorXd &newGradient, const Eigen::VectorXd &oldGradient) {
    // dx = alpha*searchDirection, but we spare this multiplication here, because it cancels out in most places
    Eigen::VectorXd dg = newGradient - oldGradient;
    double alphaA=0, B=0;
    double dxdgOverAlpha = _searchDirection.dot(dg);
    if(dxdgOverAlpha!=0) {
      B = _searchDirection.dot(newGradient)/dxdgOverAlpha;
      alphaA = dg.dot(newGradient)/dxdgOverAlpha - (alpha + dg.squaredNorm()/dxdgOverAlpha)*B;
    }
    _p = -newGradient + alphaA*_searchDirection + B*dg;
    _searchDirection = _p/_p.norm();
  }
  std::string iterationAdditionalLogString() const { return ""; }
  const Eigen::VectorXd &searchDirectionUnitVector() const { return _searchDirection; }
private:
  Eigen::VectorXd _p;
  Eigen::VectorXd _searchDirection;
};


class MemorylessBfgsSlow {
public:
  MemorylessBfgsSlow(const Eigen::VectorXd &gradient) : _p(-gradient), _searchDirection(_p/_p.norm()) {}
  void update(double alpha, const Eigen::VectorXd &newGradient, const Eigen::VectorXd &oldGradient) {
    Eigen::VectorXd dx = alpha*_searchDirection;
    Eigen::VectorXd dg = newGradient - oldGradient;
    double A=0, B=0;
    double dxdg = dx.dot(dg);
    if(dxdg!=0) {
      B = dx.dot(newGradient)/dxdg;
      A = dg.dot(newGradient)/dxdg - (1 + dg.squaredNorm()/dxdg)*B;
    }
    _p = -newGradient + A*dx + B*dg;
    _searchDirection = _p/_p.norm();
  }
  std::string iterationAdditionalLogString() const { return ""; }
  const Eigen::VectorXd &searchDirectionUnitVector() const { return _searchDirection; }
private:
  Eigen::VectorXd _p;
  Eigen::VectorXd _searchDirection;
};

#endif /* CGMINIMIZER_SEARCHDIRECTIONUPDATES_H */

