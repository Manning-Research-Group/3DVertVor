#include <iostream>
#include <Eigen/Dense>

#include "FireMinimizer.h"
#include "MinimizerWithDerivative.h"

bool FireMinimizer::converged() const {
  for(const double * const c : _gradient) if(*c>_maxGradientComponentCutoff) return false;
  return true;
}

bool FireMinimizer::minimize(std::function<double()> functionValue, 
        std::function<void()> updateGradient, 
        std::function<double()> updateGradientAndFunctionValue,
        std::function<void(int)> iterationCallback, 
        std::function<bool()> abort, 
        std::function<bool()> callbackAfterSuccessfulStep) {
  const long MaxIterations = _maxIterationsPerDof * _dofs.size();

  // initialize gradient and velocity
  Eigen::VectorXd grad=Eigen::VectorXd::Zero(_dofs.size());
  Eigen::VectorXd v=Eigen::VectorXd::Zero(_dofs.size());
  
  // initial conditions; 
  double deltaT = _initialTimeStep;
  double alpha = _initialAlpha;

  // number of successive dissipative steps
  int nDissipativeSteps = 0;

  //Algorithm loop
  for(int iteration=0; (_maxIterationsPerDof<0)||(iteration<MaxIterations); ++iteration) {
    // *** Run MD integration ***
    // calculate v(t+delta_t/2) = v(t) + 0.5*a(t)*delta_t
    v -= 0.5*deltaT*grad; //a = f = -g
    
    // DEBUG
    double oldValue = functionValue();
    
    // move dofs:  x(t+delta_t) = x(t) + v(t+delta_t/2)*delta_t 
    for(unsigned int i=0; i<_dofs.size(); ++i) {
      double displacement = deltaT*v(i);
      if((_warnLargeDisplacement>0) && (fabs(displacement)>_warnLargeDisplacement)) {
        std::cout << "Warning, step[" << i << "] = " << displacement << std::endl;
      }
      *_dofs[i] += displacement;
    }
    
    // callback after step:
    if(!callbackAfterSuccessfulStep()) {
      _status = Aborted; 
      return false;
    }
    
    // calculate and copy the gradient
    updateGradient();
    for(unsigned int i=0; i<_gradient.size(); ++i) grad(i) = *_gradient[i];
    
    //calculate v(t+delta_t) = v(t+delta_t/2) + 0.5*a(t+delta_t)*delta_t
    v -= 0.5*deltaT*grad;
    
    // *** Done. *** //

    //Convergence test and report
    if(converged()) {
      _status = Success; 
      return true;
    }

    // TAS, modified by MM:
    if (deltaT<_minimalTimeStepCutoff) {
      _status = TimeStepCutoff;
      return false;
    }

    //Calculate power
//    double P = -grad.dot(v);
    double P = oldValue - functionValue();

    if((_logInterval>0) && (iteration%_logInterval==0)) {
      double maxGradient = -1;
      for(unsigned int i=0; i<grad.rows(); ++i) if(grad(i)>maxGradient) maxGradient=grad(i);
      double curValue = functionValue();
      std::cout << "iteration " << iteration << ",  value = " << functionValue() << ",  max gradient component = " << maxGradient << ",  delta t = " << deltaT << ",  power = " << P << std::endl;      
      std::cout << "diff:       " << curValue-oldValue << std::endl;
      std::cout << "from power: " << -P*deltaT << std::endl;      
      iterationCallback(iteration);
    }
    
    //set v = (1-alpha)v + alpha|v|fhat.
    double coeff1 = 1.-alpha;
    double coeff2 = alpha*sqrt(v.squaredNorm()/grad.squaredNorm());
    v = coeff1*v - coeff2*grad;

    // Was this a dissipative step?
    if(P>0) {
      // Enough dissipative steps? Then increase delta_t and decrease alpha.
      if(nDissipativeSteps > _numberOfDissipativeStepsUntilAdaptation) {
        deltaT = deltaT*_adaptationTimeStepRelativeIncrease;
        if(_adaptationMaximalTimeStep<deltaT) deltaT = _adaptationMaximalTimeStep;
        alpha *= _adaptationAlphaRelativeIncrease;
      }
      ++nDissipativeSteps;
    } else {
      // Last step wasn't dissipative. Then decrease delta_t, set v=0, and set alpha = alpha_start.
      deltaT *= _adaptationTimeStepRelativeDecrease;
      v.setZero();
      alpha = _initialAlpha;
      nDissipativeSteps = 0;
    }
  }

  _status = MaxIterationsReached;
  return false;
}

