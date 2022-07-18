#ifndef DIAGONALIZATIONOFEXTENDEDHESSIAN_H
#define DIAGONALIZATIONOFEXTENDEDHESSIAN_H

#include "Diagonalization.h"
#include "HessianWithNamedDofs.h"

class Diagonalization;

template<typename NamedDofEnum, NamedDofEnum ...NamedDofs>
class DiagonalizationOfExtendedHessian : public Diagonalization  {
public:
  typedef HessianWithNamedDofs<NamedDofEnum, NamedDofs...>  H;
  DiagonalizationOfExtendedHessian(const H &h) : Diagonalization(h), _h(h) {}
  template <NamedDofEnum CP> double modulusWithRespectTo(const double cutoffEVal=-1, double cutoffEVecComponentSq=-1) const {
    return modulusWithRespectToDof(_h.template namedDofToIndex<CP>(), cutoffEVal, cutoffEVecComponentSq);
  }
  
private:
  const H &_h;
};

#endif /* DIAGONALIZATIONOFFULLHESSIAN_H */

