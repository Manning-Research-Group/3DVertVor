#ifndef MATRICES_H
#define MATRICES_H

#include <unordered_map>
#include <vector>

#include "computations/EigenDenseBackend.h"
#include "computations/DynamicalMatrix.h"
#include "computations/CompatibilityMatrix.h"

#include "core/Tessellation.h"

enum class PeriodicBoxDof {
  PureShearXy, SimpleShearYx
};

enum class SpringType {
  Surface, Volume
};

extern const std::vector<PeriodicBoxDof> DefaultBoxDofs;
extern const std::vector<SpringType> DefaultSprings;

// enums may not have std::hash defined, so I have to do it here
namespace std {
  template<> struct hash<PeriodicBoxDof> { size_t operator()(const PeriodicBoxDof& e) const { return static_cast<std::size_t>(e); } };
  template<> struct hash<SpringType> { size_t operator()(const SpringType& e) const { return static_cast<std::size_t>(e); } };
}

template<class Backend=EigenDenseBackend>
class Matrices {
public:
  Matrices(const Tessellation &tessellation, const std::vector<PeriodicBoxDof> &boxDofs=DefaultBoxDofs, const std::vector<SpringType> &springs=DefaultSprings);
  
  void computeCompatibilityMatrix();
  const CompatibilityMatrix<Backend> &compatibilityMatrix() const { return _compatibilityMatrix; }
  void computeDynamicalMatrix();
  const DynamicalMatrix<Backend> &dynamicalMatrix() const { return _dynamicalMatrix; }
  
  int dofIndexForBoxDof(PeriodicBoxDof boxDof) const { return _boxDofToMatrixIndex.at(boxDof); }
  int dofIndexForCell(const Cell *c) const { return 3*_cellToCellIndex.at(c);}
  int springIndexForCell(const Cell *c, SpringType spring) const { return _cellSpringToMatrixIndex.at(c).at(spring); }
  double modulusWrtBoxDof(PeriodicBoxDof boxDof, const double cutoffEVal=-1, double cutoffEVecComponentSq=-1) const { return _dynamicalMatrix.modulusWithRespectToDof(dofIndexForBoxDof(boxDof), cutoffEVal, cutoffEVecComponentSq); }
  
private:
  void computeCellContributionsToDynamicalMatrix(const Cell *c);
  void computeVertexContributionsToDynamicalMatrix(const VertexOfCell *v);
  void computeEdgeContributionsToDynamicalMatrix(const DirectedEdgeOfCell *e);
  void computeFaceContributionsToDynamicalMatrix(const DirectedFace *f);
  void addVertexPositionDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const DirectedEdgeOfCell *e, const Matrix3x3x3 &vertexPositionDerivative);
  void addFaceAreaVectorDerivativeToDynamicalMatrix(const DirectedFace* f1, const DirectedFace* f2, const DirectedFace* f, const Matrix3x3x3& areaVectorDerivative);
  void addCellEnergyDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const Matrix3x3 &derivative);

  void addToCompatibilityMatrixCell(int springIndex, const Cell *c, const Vector3D &derivative);
  void addToCompatibilityMatrixBoxDof(int springIndex, int boxDofIndex, double derivative);
  void addToDynamicalMatrixCellCell(const Cell *c1, const Cell *c2, const Matrix3x3 &derivative);
  void addToDynamicalMatrixBoxDofCell(int boxDofIndex, const Cell *c, const Vector3D &derivative);
  void addToDynamicalMatrixBoxDofBoxDof(int boxDofIndex1, int boxDofIndex2, double derivative);
  
  const Tessellation &_Tessellation;
  const PeriodicBox &_Box;
  const int _FirstBoxDof;
  std::unordered_map<const Cell*,int> _cellToCellIndex;
  std::unordered_map<PeriodicBoxDof,int> _boxDofToMatrixIndex;
  std::unordered_map<const Cell*, std::unordered_map<SpringType,int> > _cellSpringToMatrixIndex;
  CompatibilityMatrix<Backend> _compatibilityMatrix;
  DynamicalMatrix<Backend> _dynamicalMatrix;
};



#endif /* MATRICES_H */

