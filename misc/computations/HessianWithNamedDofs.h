#ifndef HESSIANWITHCONTROLPARAMETERS_H
#define HESSIANWITHCONTROLPARAMETERS_H

#include <climits>
#include "Hessian.h"

template<typename NamedDofEnum, NamedDofEnum ...NamedDofs> class DiagonalizationOfUnnamedPart;
template<typename NamedDofEnum, NamedDofEnum ...NamedDofs> class DiagonalizationOfExtendedHessian;

template<typename NamedDofEnum, NamedDofEnum ...NamedDofs>
class HessianWithNamedDofs : public Hessian {
public:
  HessianWithNamedDofs(unsigned int numberOfUnnamedDofs) : Hessian(numberOfUnnamedDofs + sizeof...(NamedDofs)) {}
  int totalNumberOfDofs() const { return numberOfDofs(); }
  int numberOfUnnamedDofs() const { return numberOfDofs() - sizeof...(NamedDofs); }
  int numberOfNamedDofs() const { return sizeof...(NamedDofs); }
  template <NamedDofEnum NamedDof> bool contains() const { return recPositionAndFound<NamedDof, NamedDofs...>::found; }
  template <NamedDofEnum NamedDof> int namedDofToIndex() const { return numberOfUnnamedDofs()+positionOfParam<NamedDof>::value; }
  
  double uu(int unnamedDof1, int unnamedDof2) const { return (*this)(unnamedDof1, unnamedDof2); }
  template <NamedDofEnum NamedDof> double nu(int unnamedDof) const { return (*this)(namedDofToIndex<NamedDof>(), unnamedDof); }
  template <NamedDofEnum NamedDof> double un(int unnamedDof) const { return (*this)(unnamedDof, namedDofToIndex<NamedDof>()); }
  template <NamedDofEnum NamedDof1, NamedDofEnum NamedDof2> double nn() const { return (*this)(namedDofToIndex<NamedDof1>(), namedDofToIndex<NamedDof2>()); }
  
  double &uu(int unnamedDof1, int unnamedDof2) { return (*this)(unnamedDof1, unnamedDof2); }
  template <NamedDofEnum NamedDof> double &nu(int unnamedDof) { return (*this)(namedDofToIndex<NamedDof>(), unnamedDof); }
  template <NamedDofEnum NamedDof> double &un(int unnamedDof) { return (*this)(unnamedDof, namedDofToIndex<NamedDof>()); }
  template <NamedDofEnum NamedDof1, NamedDofEnum NamedDof2> double &nn() { return (*this)(namedDofToIndex<NamedDof1>(), namedDofToIndex<NamedDof2>()); }

  DiagonalizationOfUnnamedPart<NamedDofEnum, NamedDofs...> diagonalizeUnnamedPart() const { return DiagonalizationOfUnnamedPart<NamedDofEnum, NamedDofs...>(*this); }
  DiagonalizationOfExtendedHessian<NamedDofEnum, NamedDofs...> diagonalizeExtendedHessian() const { return DiagonalizationOfExtendedHessian<NamedDofEnum, NamedDofs...>(*this); }
  
private:
  template <NamedDofEnum NDof, NamedDofEnum ...OtherNDofs> struct recPositionAndFound;
  template <NamedDofEnum NDof, NamedDofEnum FirstNDof, NamedDofEnum ...OtherNDofs> struct recPositionAndFound<NDof, FirstNDof, OtherNDofs...> {
    typedef recPositionAndFound<NDof, OtherNDofs...> previousStep;
    const static bool found = (NDof==FirstNDof) || previousStep::found;
    const static int position = (NDof==FirstNDof)?0:(previousStep::position + 1); 
  };
  template <NamedDofEnum NDof> struct recPositionAndFound<NDof> {
    const static bool found = false;
    const static int position = INT_MIN; 
  };
  template <NamedDofEnum NDof> struct positionOfParam {
    typedef recPositionAndFound<NDof, NamedDofs...> positionAndFound;
    static_assert(positionAndFound::found, "Control parameter not part of dynamical matrix!");
    const static int value = positionAndFound::position; 
  };
};

#endif /* DYNAMICALMATRIXWITHCONTROLPARAMETERS_H */

