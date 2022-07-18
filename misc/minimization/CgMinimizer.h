#ifndef CGMINIMIZER_H
#define CGMINIMIZER_H

#include <Eigen/Dense>
#include "MinimizerWithDerivative.h"

template <typename LineMinimizer, typename SearchDirectionUpdate>
class CgMinimizer : public MinimizerWithDerivative {
public:
  CgMinimizer(LineMinimizer &lineMinimizer) : MinimizerWithDerivative(), 
          _lineMinimizer(lineMinimizer), _gradientTolerancePerDof(1e-12), _maxIterationsPerDof(100), _initialStepSizePerDof(1e-3), _errorCode(Undefined) {}
  void setGradientTolerancePerDof(double gradientTolerancePerDof) { _gradientTolerancePerDof=gradientTolerancePerDof; }
  void setMaxIterationsPerDof(int maxIterationsPerDof) { _maxIterationsPerDof=maxIterationsPerDof; }
  void setInitialStepSizePerDof(double initialStepSizePerDof) { _initialStepSizePerDof=initialStepSizePerDof; }
  
  virtual bool minimize(std::function<double()> functionValue, 
                      std::function<void()> updateGradient, 
                      std::function<double()> updateGradientAndFunctionValue,
                      std::function<void(int)> iterationCallback = [](int iteration){},
                      std::function<bool()> abort = [](){ return false; },
                      std::function<bool()> callbackAfterSuccessfulStep = [](){ return true; }  // true="go on", false="abort minimization"
                    );
                    
  enum ErrorCode {
    Success=0, ExceededMaxIterations, Aborted, Undefined
  };
  ErrorCode errorCode() const { return _errorCode; }
                    
private:
  LineMinimizer &_lineMinimizer;
  double _gradientTolerancePerDof;
  int _maxIterationsPerDof;
  double _initialStepSizePerDof;

  ErrorCode _errorCode;
};

#endif /* CGMINIMIZER_H */

