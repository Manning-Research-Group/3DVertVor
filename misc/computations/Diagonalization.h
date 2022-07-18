#ifndef DIAGONALIZATION_H
#define DIAGONALIZATION_H

#include <Eigen/Dense>

#include "Hessian.h"

class Diagonalization : public Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> {
public:
  Diagonalization(const Hessian &h) : _h(h) { compute(_h, Eigen::ComputeEigenvectors); }
  const bool valid() const { return info() == Eigen::Success; }
  
  double modulusWithRespectToDof(int dof, const double cutoffEVal=-1, double cutoffEVecComponentSq=-1) const {
    if((cutoffEVal>=0) && (cutoffEVecComponentSq<0)) {
      cutoffEVecComponentSq = cutoffEVal;
    }    
    double sum = 0.0;
    for(int i=0; i<eigenvalues().rows(); ++i) {
      double omegaSq = eigenvalues()(i);
      double element = eigenvectors()(dof, i);
      if(omegaSq>=cutoffEVal) {
        sum += element*element/omegaSq;
      } else if(element*element>=cutoffEVecComponentSq) {
        return 0.0;
      }
    }
    return 1.0/sum;
  }

private:
  const Hessian &_h;
};

#endif /* DIAGONALIZATION_H */

