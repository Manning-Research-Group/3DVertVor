#ifndef CGMINIMIZER_INLINE_H
#define CGMINIMIZER_INLINE_H

#include <iostream>
#include <math.h>

#include "CgMinimizer.h"

template <typename LineMinimizer, typename SearchDirectionUpdate>
bool CgMinimizer<LineMinimizer, SearchDirectionUpdate>::minimize(std::function<double()> functionValue, 
                    std::function<void()> updateGradient, 
                    std::function<double()> updateGradientAndFunctionValue,
                    std::function<void(int)> iterationCallback,
                    std::function<bool()> abort,
                    std::function<bool()> callbackAfterSuccessfulStep
                  ) {
  const double TotalGradientSqCutoff = _dofs.size()*_gradientTolerancePerDof*_gradientTolerancePerDof;
  const int MaxNumberOfIterations = _dofs.size()*_maxIterationsPerDof;
  const double InitialStepSize = _initialStepSizePerDof*sqrt(_dofs.size());

  updateGradient();
  if(abort()) { _errorCode = Aborted; return false; }
  Eigen::VectorXd position = Eigen::VectorXd::Zero(_dofs.size());
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(_dofs.size());
  for(unsigned int i=0; i<_dofs.size(); ++i) {
    position(i) = *_dofs[i];
    gradient(i) = *_gradient[i];
  }
  
  SearchDirectionUpdate searchDirectionUpdate(gradient);
  double lastStepSizeTimesSlope = InitialStepSize*fabs(searchDirectionUpdate.searchDirectionUnitVector().dot(gradient));
  
  // log
  if(_logInterval>0) {
    std::cout << "initial state: " 
            << "  gradient norm = " << gradient.norm()
            << std::endl;      
  }
  
  int iteration = 0;
  while(gradient.squaredNorm()>TotalGradientSqCutoff) {
    if(iteration>MaxNumberOfIterations) {
      std::cerr << "CgMinimizer::minimize: reached max iterations!" << std::endl;
      _errorCode = ExceededMaxIterations;
      return false;
    }

    // line minimization
    Eigen::VectorXd searchDirection = searchDirectionUpdate.searchDirectionUnitVector();
    double slope = searchDirection.dot(gradient);
    if(slope>0) {
      std::cout << "CgMinimizer::minimize: positive slope before line minimization, flipping direction!" << std::endl;
      searchDirection *= -1;
      slope *= -1;
    }

    double alpha = _lineMinimizer.minimize(lastStepSizeTimesSlope, slope, 
            [&](double alpha)->double { 
              for(unsigned int i=0; i<_dofs.size(); ++i) *_dofs[i] = position(i) + alpha*searchDirection(i);
              return functionValue();
            },
            [&](double alpha)->double { 
              for(unsigned int i=0; i<_dofs.size(); ++i) *_dofs[i] = position(i) + alpha*searchDirection(i);
              updateGradient();
              double derivative = 0.0;
              for(unsigned int i=0; i<_dofs.size(); ++i) derivative += *_gradient[i]*searchDirection(i);
              return derivative;
            },
            [&](double alpha, double &derivative)->double { 
              for(unsigned int i=0; i<_dofs.size(); ++i) *_dofs[i] = position(i) + alpha*searchDirection(i);
              double value = updateGradientAndFunctionValue();
              derivative = 0.0;
              for(unsigned int i=0; i<_dofs.size(); ++i) derivative += *_gradient[i]*searchDirection(i);
              return value;
            });

    // update stuff
    lastStepSizeTimesSlope = alpha*fabs(slope);
    position += alpha*searchDirection;
    for(unsigned int i=0; i<_dofs.size(); ++i) *_dofs[i] = position(i);
    updateGradient();
    if(abort()) { _errorCode = Aborted; return false; }
    Eigen::VectorXd newGradient = Eigen::VectorXd::Zero(_dofs.size());
    for(unsigned int i=0; i<_dofs.size(); ++i) newGradient(i) = *_gradient[i];
    
    searchDirectionUpdate.update(alpha, newGradient, gradient);
    
    // log
    if((_logInterval>0) && (iteration%_logInterval==0)) {
      std::cout << "iteration " << iteration 
              << ",  gradient norm = " << newGradient.norm()
              << ",  lastStepSize = " << alpha
              << searchDirectionUpdate.iterationAdditionalLogString()
              << ",  lm its = " << _lineMinimizer.iterations()
              << std::endl;      
      iterationCallback(iteration);
    }
    
    gradient = newGradient;    
    ++iteration;
  }
  _errorCode = Success;
  return true;                      
}

#endif /* CGMINIMIZER_INLINE_H */

