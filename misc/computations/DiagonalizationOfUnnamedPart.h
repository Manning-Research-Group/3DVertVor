#ifndef DIAGONALIZATIONOFUNNAMEDPART_H
#define DIAGONALIZATIONOFUNNAMEDPART_H

#include <Eigen/Dense>

#include "HessianWithNamedDofs.h"

template<typename NamedDofEnum, NamedDofEnum ...NamedDofs>
class DiagonalizationOfUnnamedPart : public Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd>  {
public:
  typedef HessianWithNamedDofs<NamedDofEnum, NamedDofs...>  H;
  
  DiagonalizationOfUnnamedPart(const H &h) : _h(h) { 
    compute(_h.topLeftCorner(_h.numberOfUnnamedDofs(), _h.numberOfUnnamedDofs()), Eigen::ComputeEigenvectors); 
  }
  const bool valid() const { return info() == Eigen::Success; }
  
  template <NamedDofEnum CP> double modulusWithRespectTo(const double eigenValueCutoff=-1) const {
    double modulus = _h.template nn<CP,CP>();
    for(int i=0; i<eigenvalues().rows(); ++i) {
      double omegaSq = eigenvalues()(i);
      if(omegaSq>=eigenValueCutoff) {
        double projection = eigenvectors().col(i).transpose() * _h.col(_h.template namedDofToIndex<CP>()).head(_h.numberOfUnnamedDofs());
        modulus -= projection*projection/omegaSq;
      }
    }
    return modulus;
  }
  
private:
  const H &_h;
};

#endif /* DIAGONALIZATIONOFUNNAMEDBLOCK_H */

