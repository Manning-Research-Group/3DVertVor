#ifndef TESSELLATION_INLINE_H
#define TESSELLATION_INLINE_H

#include "Tessellation.h"

inline void Tessellation::addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(DynamicalMatrixOld &dm, const Cell *c1, const Cell *c2, const Matrix3x3 &derivative) const {
//  if((c1!=_cells[1]) && (c2!=_cells[1])) {
    int i1=_cellToDynamicalMatrixIndex.at(c1);
    int i2=_cellToDynamicalMatrixIndex.at(c2);
    dm.addDerivativePositionPosition(i1+0, i2+0, derivative.xx());
    dm.addDerivativePositionPosition(i1+1, i2+0, derivative.yx());
    dm.addDerivativePositionPosition(i1+2, i2+0, derivative.zx());
    dm.addDerivativePositionPosition(i1+0, i2+1, derivative.xy());
    dm.addDerivativePositionPosition(i1+1, i2+1, derivative.yy());
    dm.addDerivativePositionPosition(i1+2, i2+1, derivative.zy());
    dm.addDerivativePositionPosition(i1+0, i2+2, derivative.xz());
    dm.addDerivativePositionPosition(i1+1, i2+2, derivative.yz());
    dm.addDerivativePositionPosition(i1+2, i2+2, derivative.zz());
//  }
}

inline void Tessellation::addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(DynamicalMatrixOld &dm, int cp, const Cell *c, const Vector3D &derivative) const {
  int i=_cellToDynamicalMatrixIndex.at(c);
  dm.addDerivativeControlParameterPosition(cp, i+0, derivative.x());
  dm.addDerivativeControlParameterPosition(cp, i+1, derivative.y());
  dm.addDerivativeControlParameterPosition(cp, i+2, derivative.z());
}

template<int L> double Tessellation::bondOrientationalOrder() const {
  std::vector<std::complex<double> > sumQ;
  int numCellNeighborships=0;
  for(int m=-L;m<=L;++m) sumQ.push_back(std::complex<double>(0,0));
  for(Cell *c : cells()) {
    for(DirectedFace *f : c->faces()) {
      SphericalHarmonicsForVector<L> sph(f->cellCenterConnectionVector());
      for(int m=-L;m<=L;++m) {
        sumQ[L+m] += sph.value(L, m);
      }
      ++numCellNeighborships;
    }
  }
  double sumAbsSq = 0.0;
  for(int m=-L;m<=L;++m) {
    sumAbsSq += norm(sumQ[L+m]);
  }
  return sqrt(4*M_PI/(2*L+1)*sumAbsSq)/numCellNeighborships;
}  

// for convenience
template<int L> double bondOrientationalOrder(const Tessellation &t) { return t.bondOrientationalOrder<L>(); }

#endif /* TESSELLATION_INLINE_H */

