#ifndef FIREMINIMIZER_H
#define FIREMINIMIZER_H

#include "MinimizerWithDerivative.h"

class FireMinimizer  : public MinimizerWithDerivative {
public:
  enum Status {
    NotYetStarted,
    Success,        
    MaxIterationsReached,
    TimeStepCutoff,
    Aborted
  };
  
  FireMinimizer() :
    _status(NotYetStarted),
    _maxIterationsPerDof(1e5),
    _maxGradientComponentCutoff(1e-12),
    _numberOfDissipativeStepsUntilAdaptation(5),
    _initialTimeStep(1e-4),
    _adaptationTimeStepRelativeIncrease(1.1),
    _adaptationTimeStepRelativeDecrease(0.5),
    _adaptationMaximalTimeStep(10.*_initialTimeStep),
    _minimalTimeStepCutoff(1e-12),
    _initialAlpha(0.1),
    _adaptationAlphaRelativeIncrease(0.99),
    _warnLargeDisplacement(-1.0)
  {}
  Status status() const { return _status; }

  void setMaxIterationsPerDof(int v) { _maxIterationsPerDof=v; }
  void setMaxGradientComponentCutoff(double v) { _maxGradientComponentCutoff=v; }
  void setNumberOfDissipativeStepsUntilAdaptation(int v) { _numberOfDissipativeStepsUntilAdaptation=v; }
  void setInitialTimeStep(double v) { _initialTimeStep=v; }
  void setAdaptationTimeStepRelativeIncrease(double v) { _adaptationTimeStepRelativeIncrease=v; }
  void setAdaptationTimeStepRelativeDecrease(double v) { _adaptationTimeStepRelativeDecrease=v; }
  void setAdaptationMaximalTimeStep(double v) { _adaptationMaximalTimeStep=v; }
  void setMinimalTimeStepCutoff(double v) { _minimalTimeStepCutoff=v; }
  void setInitialAlpha(double v) { _initialAlpha=v; }
  void setAdaptationAlphaRelativeIncrease(double v) { _adaptationAlphaRelativeIncrease=v; }
  void setWarnLargeDisplacement(double v) { _warnLargeDisplacement=v; }
    
  virtual bool minimize(std::function<double()> functionValue, 
                        std::function<void()> updateGradient, 
                        std::function<double()> updateGradientAndFunctionValue,
                        std::function<void(int)> iterationCallback = [](int iteration){},
                        std::function<bool()> abort = [](){ return false; },
                        std::function<bool()> callbackAfterSuccessfulStep = [](){ return true; }  // true="go on", false="abort minimization"
                      );
  
private:
  bool converged() const;
  
  mutable Status _status;
  
  int _maxIterationsPerDof;  // ignore if negative
  double _maxGradientComponentCutoff;
  int _numberOfDissipativeStepsUntilAdaptation;
  double _initialTimeStep;
  double _adaptationTimeStepRelativeIncrease;
  double _adaptationTimeStepRelativeDecrease;
  double _adaptationMaximalTimeStep;
  double _minimalTimeStepCutoff;
  double _initialAlpha;
  double _adaptationAlphaRelativeIncrease;
  double _warnLargeDisplacement;  // ignore if negative
};

#endif /* FIREMINIMIZER_H */

