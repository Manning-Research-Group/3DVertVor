#ifndef CONJUGATEDGRADIENT_H
#define CONJUGATEDGRADIENT_H

#include <vector>
#include <functional>

class MinimizerWithDerivative {
public:
  MinimizerWithDerivative() : _logInterval(-1) {}
  void setLogInterval(int logInterval) { _logInterval=logInterval; }
  void clearDofs() { _dofs.clear(); _gradient.clear(); }
  void addDof(double &dof, const double &derivative) { _dofs.push_back(&dof); _gradient.push_back(&derivative); }
  void addDof(double *dof, const double *derivative) { _dofs.push_back(dof);  _gradient.push_back(derivative);  }
  virtual bool minimize(std::function<double()> functionValue, 
                        std::function<void()> updateGradient, 
                        std::function<double()> updateGradientAndFunctionValue,
                        std::function<void(int)> iterationCallback = [](int iteration){},
                        std::function<bool()> abort = [](){ return false; },
                        std::function<bool()> callbackAfterSuccessfulStep = [](){ return true; }  // true="go on", false="abort minimization"
                      ) = 0;

protected:
  int _logInterval;
  std::vector<double*> _dofs;
  std::vector<const double*> _gradient;
};

#endif /* CONJUGATEDGRADIENT_H */

