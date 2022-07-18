#include <iostream>
#include <math.h>
#include <Eigen/Dense>

#include "SimpleGradient.h"

bool SimpleGradientMinimizer::minimize(std::function<double()> functionValue, std::function<void()> updateDerivatives, std::function<double()> updateDerivativesAndFunctionValue,
                        std::function<void(int)> iterationCallback, std::function<bool()> abort, std::function<bool()> callbackAfterSuccessfulStep) {
  const double TotalForceSqCutoff = _dofs.size()*_gradientTolerancePerDof*_gradientTolerancePerDof;
  const int MaxNumberOfIterations = _dofs.size()*_maxIterationsPerDof;
  _timeStep = _initialTimeStep;

  double value = updateDerivativesAndFunctionValue();
  if(abort()) {
    _errorCode = Aborted;
    return false;
  }
  if(std::isnan(value)) {
    std::cerr << "SimpleGradient::minimize: got nan as value!" << std::endl;
    _errorCode = GotNan;
    return false;
  }
  Eigen::VectorXd dofs = Eigen::VectorXd::Zero(_dofs.size());
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(_dofs.size());
  for(unsigned int i=0; i<_dofs.size(); ++i) {
//    std::cout << "DEBUG: dof " << i << std::endl;
    
    dofs(i) = *_dofs[i];
    if(std::isnan(dofs(i))) {
      std::cerr << "SimpleGradient::minimize: got nan as initial state component!" << std::endl;
      _errorCode = GotNan;
      return false;
    }
    gradient(i) = *_gradient[i];
    if(std::isnan(gradient(i))) {
      std::cerr << "SimpleGradient::minimize: got nan as gradient component!" << std::endl;
      _errorCode = GotNan;
      return false;
    }
  }
  
  Eigen::VectorXd lastDofs = dofs;
  double lastValue = value;
  Eigen::VectorXd lastGradient = gradient;
  
  int iteration = 0;
  while(gradient.squaredNorm()>TotalForceSqCutoff) {
    if(iteration>MaxNumberOfIterations) {
      std::cerr << "SimpleGradient::minimize: reached max iterations!" << std::endl;
      _errorCode = ExceededMaxIterations;
      return false;
    }
    
    // step
    dofs -= _timeStep*gradient;
    for(unsigned int i=0; i<_dofs.size(); ++i) {
      (*_dofs[i]) = dofs(i);
    }
    value = updateDerivativesAndFunctionValue();
    if(abort()) {
      _errorCode = Aborted;
      return false;
    }
    if(std::isnan(value)) {
      std::cerr << "SimpleGradient::minimize: got nan as value!" << std::endl;
      _errorCode = GotNan;
      return false;
    }
    for(unsigned int i=0; i<_dofs.size(); ++i) {
      gradient(i) = *_gradient[i];
      if(std::isnan(gradient(i))) {
        std::cerr << "SimpleGradient::minimize: got nan as gradient component!" << std::endl;
        _errorCode = GotNan;
        return false;
      }
    }
    
    // log
    if((_logInterval>0) && (iteration%_logInterval==0)) {
      std::cout << "iteration " << iteration << ",  value = " << value 
              << ",  total derivative norm = " << gradient.norm()
              << ",  time step = " << _timeStep 
              << ",  scalar product = " << gradient.dot(lastGradient)/sqrt(gradient.squaredNorm()*lastGradient.squaredNorm()) 
              << ",  delta value = " << value-lastValue
              << std::endl;      
      iterationCallback(iteration);
    }
    
    // check if step was ok
    if((gradient.dot(lastGradient)<0) || (value>lastValue+_energyDifferenceGrace)) {
      // problem in last step
      _timeStep *= _timeStepDecreaseFactor;
      if(_timeStep<_minTimeStep) {
        std::cerr << "SimpleGradient::minimize: hit minimal time step!" << std::endl;
        _errorCode = HitMinTimeStep;
        return false;
      }
      // revert to last position
      dofs = lastDofs;
      value = lastValue;
      gradient = lastGradient;
    } else {
      // this was a succesful step
      if(!callbackAfterSuccessfulStep()) {
        _errorCode = Aborted;
        return false;
      }
      // last step ok
      _timeStep *= _timeStepIncreaseFactor;
      // save as last position
      lastDofs = dofs;
      lastValue = value;
      lastGradient = gradient;
    }

    ++iteration;
  }
  _errorCode = Success;
  return true;
}
