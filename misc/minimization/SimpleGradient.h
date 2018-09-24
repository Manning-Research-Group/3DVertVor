#ifndef SIMPLEGRADIENT_H
#define SIMPLEGRADIENT_H

#include "MinimizerWithDerivative.h"

class SimpleGradientMinimizer : public MinimizerWithDerivative {
public:
  SimpleGradientMinimizer() : _initialTimeStep(1e-3), _timeStepIncreaseFactor(1.1), _timeStepDecreaseFactor(0.5), _minTimeStep(1e-14), 
          _gradientTolerancePerDof(1e-4), _energyDifferenceGrace(0), _maxIterationsPerDof(100), _timeStep(_initialTimeStep), _errorCode(Undefined) {};
  void setInitialTimeStep(double initialTimeStep) { _initialTimeStep=initialTimeStep; }
  void setTimeStepIncreaseFactor(int timeStepIncreaseFactor) { _timeStepIncreaseFactor=timeStepIncreaseFactor; }
  void setTimeStepDecreaseFactor(int timeStepDecreaseFactor) { _timeStepDecreaseFactor=timeStepDecreaseFactor; }
  void setMinTimeStep(int minTimeStep) { _minTimeStep=minTimeStep; }
  void setGradientTolerancePerDof(double gradientTolerancePerDof) { _gradientTolerancePerDof=gradientTolerancePerDof; }
  void setEnergyDifferenceGrace(double energyDifferenceGrace) { _energyDifferenceGrace=energyDifferenceGrace; }
  void setMaxIterationsPerDof(int maxIterationsPerDof) { _maxIterationsPerDof=maxIterationsPerDof; }

  
  virtual bool minimize(std::function<double()> functionValue, 
                        std::function<void()> updateGradient, 
                        std::function<double()> updateGradientAndFunctionValue,
                        std::function<void(int)> iterationCallback = [](int iteration){},
                        std::function<bool()> abort = [](){ return false; },
                        std::function<bool()> callbackAfterSuccessfulStep = [](){ return true; }  // true="go on", false="abort minimization"
                      );
  double timeStep() const { return _timeStep; }
  
  enum ErrorCode {
    Success=0, ExceededMaxIterations, HitMinTimeStep, GotNan, Aborted, Undefined
  };
  ErrorCode errorCode() const { return _errorCode; }
  
private:
  double _initialTimeStep;
  double _timeStepIncreaseFactor;
  double _timeStepDecreaseFactor;
  double _minTimeStep;
  double _gradientTolerancePerDof;
  double _energyDifferenceGrace;
  int _maxIterationsPerDof;
  
  double _timeStep;
  ErrorCode _errorCode;
};

#endif /* SIMPLEGRADIENT_H */

