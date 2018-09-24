#ifndef MATRICES_INLINE_H
#define MATRICES_INLINE_H

#include "core/DirectedEdgeOfCell-inline.h"

#include "Matrices.h"

template<class Backend>
Matrices<Backend>::Matrices(const Tessellation& t, const std::vector<PeriodicBoxDof> &boxDofs, const std::vector<SpringType> &springs) :
  _Tessellation(t),
  _Box(t.box()),
  _FirstBoxDof(3*t.cells().size()),
  _compatibilityMatrix(springs.size()*t.cells().size(), _FirstBoxDof + boxDofs.size()),
  _dynamicalMatrix(_FirstBoxDof + boxDofs.size()) {

  // prepare maps
  for(unsigned int i=0; i<t.cells().size(); ++i) {
    _cellToCellIndex[t.cells()[i]] = i;
  }
  for(unsigned int i=0; i<boxDofs.size(); ++i) {
    _boxDofToMatrixIndex[boxDofs[i]] = _FirstBoxDof + i;
  }
  int currentSpringIndex = 0;
  for(Cell *c : t.cells()) {
    for(SpringType st : springs) {
      _cellSpringToMatrixIndex[c][st] = currentSpringIndex;
      ++currentSpringIndex;
    }
  }
}

template<class Backend>
void Matrices<Backend>::computeCompatibilityMatrix() {
  for(Cell *c : _Tessellation.cells()) {
    for(auto springTypeAndIndex : _cellSpringToMatrixIndex.at(c)) {
      SpringType springType = springTypeAndIndex.first;
      int currentSpringIndex = springTypeAndIndex.second;
      if(springType==SpringType::Surface) {
        Vector3D derSurfaceWrtOwnPosition(0,0,0);
        for(DirectedFace *f : c->faces()) {
          addToCompatibilityMatrixCell(currentSpringIndex, f->otherCell(), f->derCellSurfaceWrtCellCenterConnection());
          derSurfaceWrtOwnPosition -= f->derCellSurfaceWrtCellCenterConnection();
          if(_boxDofToMatrixIndex.count(PeriodicBoxDof::PureShearXy)>0) {
            addToCompatibilityMatrixBoxDof(currentSpringIndex, dofIndexForBoxDof(PeriodicBoxDof::PureShearXy), 
                      _Box.x*f->periodicity().x*f->derCellSurfaceWrtCellCenterConnection().x() 
                    - _Box.y*f->periodicity().y*f->derCellSurfaceWrtCellCenterConnection().y());
          }
          if(_boxDofToMatrixIndex.count(PeriodicBoxDof::SimpleShearYx)>0) {
            addToCompatibilityMatrixBoxDof(currentSpringIndex, dofIndexForBoxDof(PeriodicBoxDof::SimpleShearYx), 
                    _Box.y*f->periodicity().y*f->derCellSurfaceWrtCellCenterConnection().x());
          }
        }
        addToCompatibilityMatrixCell(currentSpringIndex, c, derSurfaceWrtOwnPosition);
      } else if(springType==SpringType::Volume) {
        Vector3D derVolumeWrtOwnPosition(0,0,0);
        for(DirectedFace *f : c->faces()) {
          addToCompatibilityMatrixCell(currentSpringIndex, f->otherCell(), f->derCellVolumeWrtCellCenterConnection());
          derVolumeWrtOwnPosition -= f->derCellVolumeWrtCellCenterConnection();
          if(_boxDofToMatrixIndex.count(PeriodicBoxDof::PureShearXy)>0) {
            addToCompatibilityMatrixBoxDof(currentSpringIndex, dofIndexForBoxDof(PeriodicBoxDof::PureShearXy), 
                      _Box.x*f->periodicity().x*f->derCellVolumeWrtCellCenterConnection().x() 
                    - _Box.y*f->periodicity().y*f->derCellVolumeWrtCellCenterConnection().y());
          }
          if(_boxDofToMatrixIndex.count(PeriodicBoxDof::SimpleShearYx)>0) {
            addToCompatibilityMatrixBoxDof(currentSpringIndex, dofIndexForBoxDof(PeriodicBoxDof::SimpleShearYx), 
                    _Box.y*f->periodicity().y*f->derCellVolumeWrtCellCenterConnection().x());
          }
        }
        addToCompatibilityMatrixCell(currentSpringIndex, c, derVolumeWrtOwnPosition);
      }
    }
  }
  
  _compatibilityMatrix.performSingularValueDecomposition();
}

template<class Backend>
void Matrices<Backend>::computeDynamicalMatrix() {
  // compute most of the contributions
  for(Cell *c : _Tessellation.cells()) {
    computeCellContributionsToDynamicalMatrix(c);
  }

  // part pure shear derivatives
  if(_boxDofToMatrixIndex.count(PeriodicBoxDof::PureShearXy)>0) {
    int pureShearIndex = dofIndexForBoxDof(PeriodicBoxDof::PureShearXy);
    for(Cell *c : _Tessellation.cells()) {
      for(DirectedFace *f : c->faces()) {
        addToDynamicalMatrixBoxDofBoxDof(pureShearIndex, pureShearIndex, 
                _Box.x*f->periodicity().x*f->derCellEnergyWrtCellCenterConnection().x()
              + _Box.y*f->periodicity().y*f->derCellEnergyWrtCellCenterConnection().y());
      }
    }
  }

  _dynamicalMatrix.buildAndDiagonalize();
}

template<class Backend>
void Matrices<Backend>::computeCellContributionsToDynamicalMatrix(const Cell *c) {
  // this needs to come before the vertices to prepare treatment of the results fed by the vertices
  for(DirectedFace *f : c->faces()) {
    f->createDoublyLinkedEdgeList();
  }
  
  for(VertexOfCell *v : c->vertices()) {
    computeVertexContributionsToDynamicalMatrix(v);
  }
  for(DirectedEdgeOfCell *e : c->edges()) {
    if(e<e->conjugated()) {
      computeEdgeContributionsToDynamicalMatrix(e);
    }
  }
  for(DirectedFace *f : c->faces()) {
    computeFaceContributionsToDynamicalMatrix(f);
  }
  
  // terms from first area derivative
  for(DirectedFace *f : c->faces()) {
    const double AreaNorm = f->area().norm();
    const Vector3D NormalizedArea = f->area()/AreaNorm;
    
    Matrix3x3 derivativeFF(f->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(f->derFaceAreaWrtCellCenterConnection()));
    derivativeFF -= dyadicProduct(f->derFaceAreaWrtCellCenterConnection()*NormalizedArea, f->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
    derivativeFF *= c->derOwnEnergyWrtSurface()/AreaNorm;
    derivativeFF += (c->derOwnEnergyWrtVolume()/6.0)*(f->derFaceAreaWrtCellCenterConnection()+f->derFaceAreaWrtCellCenterConnection().transposed());
    addCellEnergyDerivativeToDynamicalMatrix(f, f, derivativeFF);
    
    // so far, surface terms only (not even accounting for the first derivatives of the surface)
    DirectedEdgeOfCell *e1 = f->firstEdge();
    do {
      Matrix3x3 derivativeE1F(e1->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(f->derFaceAreaWrtCellCenterConnection()));
      derivativeE1F -= dyadicProduct(e1->derFaceAreaWrtCellCenterConnection()*NormalizedArea, f->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
      derivativeE1F *= c->derOwnEnergyWrtSurface()/AreaNorm;
      derivativeE1F += (c->derOwnEnergyWrtVolume()/6.0)*e1->derFaceAreaWrtCellCenterConnection();
      addCellEnergyDerivativeToDynamicalMatrix(e1->conjugated()->face(), f, derivativeE1F);

      Matrix3x3 derivativeFE1(f->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(e1->derFaceAreaWrtCellCenterConnection()));
      derivativeFE1 -= dyadicProduct(f->derFaceAreaWrtCellCenterConnection()*NormalizedArea, e1->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
      derivativeFE1 *= c->derOwnEnergyWrtSurface()/AreaNorm;
      derivativeFE1 += (c->derOwnEnergyWrtVolume()/6.0)*e1->derFaceAreaWrtCellCenterConnection().transposed();
      addCellEnergyDerivativeToDynamicalMatrix(f, e1->conjugated()->face(), derivativeFE1);
      
      DirectedEdgeOfCell *e2 = f->firstEdge();
      do {
        Matrix3x3 derivativeE1E2(e1->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(e2->derFaceAreaWrtCellCenterConnection()));
        derivativeE1E2 -= dyadicProduct(e1->derFaceAreaWrtCellCenterConnection()*NormalizedArea, e2->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
        derivativeE1E2 *= c->derOwnEnergyWrtSurface()/AreaNorm;
        addCellEnergyDerivativeToDynamicalMatrix(e1->conjugated()->face(), e2->conjugated()->face(), derivativeE1E2);
//        Matrix3x3 derivativeE1E2(e1->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(e2->derFaceAreaWrtCellCenterConnection()));
//        derivativeE1E2 -= dyadicProduct(e1->derFaceAreaWrtCellCenterConnection()*NormalizedArea, e2->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
//        derivativeE1E2 *= 1.0/AreaNorm;
//        _tessellation.addCellEnergyDerivativeToDynamicalMatrix(e1->conjugated()->face(), e2->conjugated()->face(), derivativeE1E2);

        e2 = e2->nextAroundFace();
      } while(e2!=f->firstEdge());
      e1 = e1->nextAroundFace();
    } while(e1!=f->firstEdge());
  }
  
  // terms from first surface and volume derivatives
  for(DirectedFace *f1 : c->faces()) {
    for(DirectedFace *f2 : c->faces()) {
      Matrix3x3 derivativeF1F2(dyadicProduct(f1->derCellVolumeWrtCellCenterConnection(), f2->derCellVolumeWrtCellCenterConnection()));
      derivativeF1F2 *= c->derOwnEnergyWrtVolume2nd();
      derivativeF1F2 += dyadicProduct(c->derOwnEnergyWrtSurface2nd()*f1->derCellSurfaceWrtCellCenterConnection(), f2->derCellSurfaceWrtCellCenterConnection());
      addCellEnergyDerivativeToDynamicalMatrix(f1, f2, derivativeF1F2);
    }
  }
}

template<class Backend>
void Matrices<Backend>::computeVertexContributionsToDynamicalMatrix(const VertexOfCell *v) {
  for(int i=0; i<3; ++i) {
    const DirectedEdgeOfCell *edgeB(v->edges()[i]); 
    const DirectedEdgeOfCell *edgeC(v->edges()[(i+1)%3]); 
    const DirectedEdgeOfCell *edgeD(v->edges()[(i+2)%3]); 
    const DirectedFace *faceB(edgeB->face()); 
    const DirectedFace *faceC(edgeC->face()); 
    const DirectedFace *faceD(edgeD->face()); 
    const Vector3D &rab(faceB->cellCenterConnectionVector()); 
    const Vector3D &rac(faceC->cellCenterConnectionVector()); 
    const Vector3D &rad(faceD->cellCenterConnectionVector()); 

    // second derivative BB
    Matrix3x3x3 secondDerivativeBB(dyadicProduct(Matrix3x3::Identity, crossProduct(2*v->positionRelativeToCellDenominator()*rac,rad)));
    secondDerivativeBB += dyadicProduct(edgeB->derVertexPositionDenominatorWrtCellCenterConnection(), edgeB->derVertexPositionNumeratorWrtCellCenterConnection());
    secondDerivativeBB += dyadicProductInBetween(edgeB->derVertexPositionDenominatorWrtCellCenterConnection(), edgeB->derVertexPositionNumeratorWrtCellCenterConnection());
    secondDerivativeBB += dyadicProduct(edgeB->derVertexPositionDenominatorWrtCellCenterConnection()/(0.5*v->positionRelativeToCellDenominator()),
                                        edgeB->derVertexPositionDenominatorWrtCellCenterConnection(),
                                        v->positionRelativeToCellNumerator());
    for(DirectedEdgeOfCell *e : v->edges()) {
      addVertexPositionDerivativeToDynamicalMatrix(faceB, faceB, e, secondDerivativeBB);
    }


    // second derivative BC
    Matrix3x3x3 secondDerivativeOfNumeratorBC(dyadicProduct(-2*rab, antisymmetricMatrixFromVector(rad)));
    secondDerivativeOfNumeratorBC += dyadicProductInBetween(2*rac, antisymmetricMatrixFromVector(rad));
    secondDerivativeOfNumeratorBC += rad.normSq()*Matrix3x3x3::Epsilon;
    Matrix3x3 secondDerivativeOfDenominatorBC(dyadicProduct(edgeB->derVertexPositionDenominatorWrtCellCenterConnection()/(0.5*v->positionRelativeToCellDenominator()),
                                                            edgeC->derVertexPositionDenominatorWrtCellCenterConnection()));
    secondDerivativeOfDenominatorBC -= antisymmetricMatrixFromVector(rad*(2*v->positionRelativeToCellDenominator()*v->positionRelativeToCellDenominator()));
    Matrix3x3x3 secondDerivativeBC(v->positionRelativeToCellDenominator()*secondDerivativeOfNumeratorBC);
    secondDerivativeBC += dyadicProduct(edgeB->derVertexPositionDenominatorWrtCellCenterConnection(), edgeC->derVertexPositionNumeratorWrtCellCenterConnection());
    secondDerivativeBC += dyadicProductInBetween(edgeC->derVertexPositionDenominatorWrtCellCenterConnection(), edgeB->derVertexPositionNumeratorWrtCellCenterConnection());
    secondDerivativeBC += dyadicProduct(secondDerivativeOfDenominatorBC, v->positionRelativeToCellNumerator());
    for(DirectedEdgeOfCell *e : v->edges()) {
      addVertexPositionDerivativeToDynamicalMatrix(faceB, faceC, e, secondDerivativeBC);
    }

    
    // second derivative BD
    Matrix3x3x3 secondDerivativeOfNumeratorBD(dyadicProduct(2*rab, antisymmetricMatrixFromVector(rac)));
    secondDerivativeOfNumeratorBD -= rac.normSq()*Matrix3x3x3::Epsilon;
    secondDerivativeOfNumeratorBD -= dyadicProductInBetween(2*rad, antisymmetricMatrixFromVector(rac));
    Matrix3x3 secondDerivativeOfDenominatorBD(dyadicProduct(edgeB->derVertexPositionDenominatorWrtCellCenterConnection()/(0.5*v->positionRelativeToCellDenominator()),
                                                            edgeD->derVertexPositionDenominatorWrtCellCenterConnection()));
    secondDerivativeOfDenominatorBD += antisymmetricMatrixFromVector(rac*(2*v->positionRelativeToCellDenominator()*v->positionRelativeToCellDenominator()));
    Matrix3x3x3 secondDerivativeBD(v->positionRelativeToCellDenominator()*secondDerivativeOfNumeratorBD);
    secondDerivativeBD += dyadicProduct(edgeB->derVertexPositionDenominatorWrtCellCenterConnection(), edgeD->derVertexPositionNumeratorWrtCellCenterConnection());
    secondDerivativeBD += dyadicProductInBetween(edgeD->derVertexPositionDenominatorWrtCellCenterConnection(), edgeB->derVertexPositionNumeratorWrtCellCenterConnection());
    secondDerivativeBD += dyadicProduct(secondDerivativeOfDenominatorBD, v->positionRelativeToCellNumerator());
    for(DirectedEdgeOfCell *e : v->edges()) {
      addVertexPositionDerivativeToDynamicalMatrix(faceB, faceD, e, secondDerivativeBD);
    }
  }
}

template<class Backend>
void Matrices<Backend>::addVertexPositionDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const DirectedEdgeOfCell *e, const Matrix3x3x3 &vertexPositionDerivative) {
  // edge energy
  Cell *c = f1->cell();
  if(c->type()->edgeSprings) {
    addCellEnergyDerivativeToDynamicalMatrix(f1, f2, vertexPositionDerivative * e->derEdgeEnergyWrtVertexPosition());
  }

  // surface and volume energy
  addFaceAreaVectorDerivativeToDynamicalMatrix(f1, f2, e->face(), vertexPositionDerivative
          *antisymmetricMatrixFromVector(0.5*(e->nextAroundFace()->vertex()->positionRelativeToCell()-e->previousAroundFace()->vertex()->positionRelativeToCell())));
}


template<class Backend>
void Matrices<Backend>::computeEdgeContributionsToDynamicalMatrix(const DirectedEdgeOfCell *e) {
  Cell *c = e->face()->cell();
  const CellType *parameters = c->type();
  if(parameters->edgeSprings) {
    Vector3D nEdgeVector(e->negativeLengthVector());
    double length = nEdgeVector.norm();
    double energyFactor = 2.0*parameters->edgeElasticity;
    Matrix3x3 derEnergyWrtEdgeVector(energyFactor*(1.0-e->restLength()/length)*Matrix3x3::Identity);
    derEnergyWrtEdgeVector += dyadicProduct(energyFactor*e->restLength()/(length*length*length) * nEdgeVector, nEdgeVector);

    // this vertex - this vertex:
    for(DirectedEdgeOfCell *e1 : e->vertex()->edges()) {
      for(DirectedEdgeOfCell *e2 : e->vertex()->edges()) {
        addCellEnergyDerivativeToDynamicalMatrix(e1->face(), e2->face(), 
                e1->derVertexPositionWrtCellCenterConnection() * derEnergyWrtEdgeVector.multiplyFromRightByMatrix22( e2->derVertexPositionWrtCellCenterConnection() )
        );
      }
    }

    // this vertex - other vertex:
    for(DirectedEdgeOfCell *e1 : e->vertex()->edges()) {
      for(DirectedEdgeOfCell *e2 : e->conjugated()->vertex()->edges()) {
        Matrix3x3 der(e1->derVertexPositionWrtCellCenterConnection() * derEnergyWrtEdgeVector.multiplyFromRightByMatrix22( e2->derVertexPositionWrtCellCenterConnection() ));
        der *= -1;
        addCellEnergyDerivativeToDynamicalMatrix(e1->face(), e2->face(), der);
        addCellEnergyDerivativeToDynamicalMatrix(e2->face(), e1->face(), der.transposed());
      }
    }

    // other vertex - other vertex:
    for(DirectedEdgeOfCell *e1 : e->conjugated()->vertex()->edges()) {
      for(DirectedEdgeOfCell *e2 : e->conjugated()->vertex()->edges()) {
        addCellEnergyDerivativeToDynamicalMatrix(e1->face(), e2->face(), 
                e1->derVertexPositionWrtCellCenterConnection() * derEnergyWrtEdgeVector.multiplyFromRightByMatrix22( e2->derVertexPositionWrtCellCenterConnection() )
        );
      }
    }
  }  
}

template<class Backend>
void Matrices<Backend>::computeFaceContributionsToDynamicalMatrix(const DirectedFace *f) {
  const Matrix3x3x3 HalfEpsilon(0.5*Matrix3x3x3::Epsilon);
  DirectedEdgeOfCell *edge = f->firstEdge();
  do {
    DirectedEdgeOfCell *nextEdge = edge->nextAroundFace();
    
    VertexOfCell *v1 = edge->vertex();
    VertexOfCell *v2 = nextEdge->vertex();
    for(DirectedEdgeOfCell *edgeV1 : v1->edges()) {
      for(DirectedEdgeOfCell *edgeV2 : v2->edges()) {
        addFaceAreaVectorDerivativeToDynamicalMatrix(edgeV1->face(), edgeV2->face(), f, 
              - HalfEpsilon
                  .multiplyFromLeftByMatrix12(edgeV1->derVertexPositionWrtCellCenterConnection())
                  .multiplyIntoMiddleByMatrix22(edgeV2->derVertexPositionWrtCellCenterConnection()));
        addFaceAreaVectorDerivativeToDynamicalMatrix(edgeV2->face(), edgeV1->face(), f, 
                HalfEpsilon
                  .multiplyFromLeftByMatrix12(edgeV2->derVertexPositionWrtCellCenterConnection())
                  .multiplyIntoMiddleByMatrix22(edgeV1->derVertexPositionWrtCellCenterConnection()));
      }
    }
    
    edge = nextEdge;
  } while(edge!=f->firstEdge());
}

template<class Backend>
void Matrices<Backend>::addFaceAreaVectorDerivativeToDynamicalMatrix(const DirectedFace* f1, const DirectedFace* f2, const DirectedFace* f, const Matrix3x3x3& areaVectorDerivative) {
  Matrix3x3 energyDerivative(areaVectorDerivative*((f->cell()->derOwnEnergyWrtSurface()/f->area().norm())*f->area()));
  energyDerivative += areaVectorDerivative*((f->cell()->derOwnEnergyWrtVolume()/6.0)*f->cellCenterConnectionVector());
  addCellEnergyDerivativeToDynamicalMatrix(f1, f2, energyDerivative);
}

template<class Backend>
void Matrices<Backend>::addCellEnergyDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const Matrix3x3 &cellEnergyDerivative) {
  // actual dynamical matrix
  addToDynamicalMatrixCellCell(f1->otherCell(), f2->otherCell(), cellEnergyDerivative);
  addToDynamicalMatrixCellCell(f1->cell(), f2->otherCell(), -cellEnergyDerivative);
  addToDynamicalMatrixCellCell(f1->otherCell(), f2->cell(), -cellEnergyDerivative);
  addToDynamicalMatrixCellCell(f1->cell(), f2->cell(), cellEnergyDerivative);
  
  if(_boxDofToMatrixIndex.count(PeriodicBoxDof::PureShearXy)>0) {
    int pureShearIndex = dofIndexForBoxDof(PeriodicBoxDof::PureShearXy);
    addToDynamicalMatrixBoxDofCell(pureShearIndex, f1->otherCell(), 
            _Box.x*f2->periodicity().x*cellEnergyDerivative.columnX() - _Box.y*f2->periodicity().y*cellEnergyDerivative.columnY());
    addToDynamicalMatrixBoxDofCell(pureShearIndex, f1->cell(), 
          - _Box.x*f2->periodicity().x*cellEnergyDerivative.columnX() + _Box.y*f2->periodicity().y*cellEnergyDerivative.columnY());
    addToDynamicalMatrixBoxDofBoxDof(pureShearIndex, pureShearIndex, 
            _Box.x*_Box.x*f1->periodicity().x*f2->periodicity().x*cellEnergyDerivative.xx()
         -2*_Box.x*_Box.y*f1->periodicity().x*f2->periodicity().y*cellEnergyDerivative.xy()
          + _Box.y*_Box.y*f1->periodicity().y*f2->periodicity().y*cellEnergyDerivative.yy());
  }

  if(_boxDofToMatrixIndex.count(PeriodicBoxDof::SimpleShearYx)>0) {
    int simpleShearIndex = dofIndexForBoxDof(PeriodicBoxDof::SimpleShearYx);
    addToDynamicalMatrixBoxDofCell(simpleShearIndex, f1->otherCell(), _Box.y*f2->periodicity().y*cellEnergyDerivative.columnX());
    addToDynamicalMatrixBoxDofCell(simpleShearIndex, f1->cell(), -_Box.y*f2->periodicity().y*cellEnergyDerivative.columnX());
    addToDynamicalMatrixBoxDofBoxDof(simpleShearIndex, simpleShearIndex, _Box.y*_Box.y*f1->periodicity().y*f2->periodicity().y*cellEnergyDerivative.xx());  
  }
}


template<class Backend>
void Matrices<Backend>::addToCompatibilityMatrixCell(int springIndex, const Cell *c, const Vector3D &derivative) {
  int i=dofIndexForCell(c);
  _compatibilityMatrix.addToElement(springIndex, i+0, derivative.x());
  _compatibilityMatrix.addToElement(springIndex, i+1, derivative.y());
  _compatibilityMatrix.addToElement(springIndex, i+2, derivative.z());  
}

template<class Backend>
void Matrices<Backend>::addToCompatibilityMatrixBoxDof(int springIndex, int boxDofIndex, double derivative) {
  _compatibilityMatrix.addToElement(springIndex, boxDofIndex, derivative);
}

template<class Backend>
void Matrices<Backend>::addToDynamicalMatrixCellCell(const Cell *c1, const Cell *c2, const Matrix3x3 &derivative) {
  int i1=dofIndexForCell(c1);
  int i2=dofIndexForCell(c2);
  _dynamicalMatrix.addToElement(i1+0, i2+0, derivative.xx());
  _dynamicalMatrix.addToElement(i1+1, i2+0, derivative.yx());
  _dynamicalMatrix.addToElement(i1+2, i2+0, derivative.zx());
  _dynamicalMatrix.addToElement(i1+0, i2+1, derivative.xy());
  _dynamicalMatrix.addToElement(i1+1, i2+1, derivative.yy());
  _dynamicalMatrix.addToElement(i1+2, i2+1, derivative.zy());
  _dynamicalMatrix.addToElement(i1+0, i2+2, derivative.xz());
  _dynamicalMatrix.addToElement(i1+1, i2+2, derivative.yz());
  _dynamicalMatrix.addToElement(i1+2, i2+2, derivative.zz());
}
  
template<class Backend>
void Matrices<Backend>::addToDynamicalMatrixBoxDofCell(int boxDofIndex, const Cell *c, const Vector3D &derivative) {
  int i=dofIndexForCell(c);
  _dynamicalMatrix.addToElement(boxDofIndex, i+0, derivative.x());
  _dynamicalMatrix.addToElement(boxDofIndex, i+1, derivative.y());
  _dynamicalMatrix.addToElement(boxDofIndex, i+2, derivative.z());  
  _dynamicalMatrix.addToElement(i+0, boxDofIndex, derivative.x());
  _dynamicalMatrix.addToElement(i+1, boxDofIndex, derivative.y());
  _dynamicalMatrix.addToElement(i+2, boxDofIndex, derivative.z());  
}

template<class Backend>
void Matrices<Backend>::addToDynamicalMatrixBoxDofBoxDof(int boxDofIndex1, int boxDofIndex2, double derivative) {
  _dynamicalMatrix.addToElement(boxDofIndex1, boxDofIndex2, derivative);
}


#endif /* MATRICES_INLINE_H */

