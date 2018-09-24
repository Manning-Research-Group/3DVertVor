#include "misc/geometry/Vector3D.h"
#include "misc/geometry/Matrix3x3.h"
#include "misc/geometry/Matrix3x3x3.h"

#include "VertexOfCell.h"
#include "DirectedEdgeOfCell.h"
#include "DirectedFace.h"
#include "Cell.h"
#include "Tessellation.h"

#include "Tessellation-inline.h"
#include "DirectedEdgeOfCell-inline.h"

void Tessellation::computeDynamicalMatrix() {
  if(_dynamicalMatrix) {
    delete _dynamicalMatrix;
  }
  
  // 2 "control parameters": 0=pure shear xy; 1=simple shear yx (normal,force); 2=volume
  _dynamicalMatrix = new DynamicalMatrixOld(3*_cells.size(), 3);
//  _dynamicalMatrix = new DynamicalMatrixTest(3*_cells.size(), 2);
  for(unsigned int i=0; i<_cells.size(); ++i) {
    _cellToDynamicalMatrixIndex[_cells[i]] = 3*i;
  }

  // for volume derivative
  const double Volume = _box.volume();
  _boxDimensionsVector.set(_box.x, _box.y, _box.z);
  _boxDimensionsVector /= 3*Volume;
  
  // compute most of the contributions
  for(Cell *c : _cells) {
    c->computeContributionsToDynamicalMatrix(*this);   // argument only for debugging!
  }

  // part of control parameter derivatives
  const double VolumePreFactor = -2.0/(3.0*Volume);
  for(Cell *c : _cells) {
    for(DirectedFace *f : c->faces()) {
      _dynamicalMatrix->add2ndDerivativeControlParameter(0, 
              _box.x*f->periodicity().x*f->derCellEnergyWrtCellCenterConnection().x()
            + _box.y*f->periodicity().y*f->derCellEnergyWrtCellCenterConnection().y());
      Vector3D periodicityAndBox(f->periodicity()*_boxDimensionsVector);
      _dynamicalMatrix->add2ndDerivativeControlParameter(2, 
            VolumePreFactor*periodicityAndBox*f->derCellEnergyWrtCellCenterConnection());
    }
  }
}

void Tessellation::computeAndSolveDynamicalMatrix() {
  computeDynamicalMatrix();
  // diagonalization
  _dynamicalMatrix->solve();
}

inline void Tessellation::addCellEnergyDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const Matrix3x3 &cellEnergyDerivative) {
  // actual dynamical matrix
  addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(*_dynamicalMatrix, f1->otherCell(), f2->otherCell(), cellEnergyDerivative);
  addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(*_dynamicalMatrix, f1->cell(), f2->otherCell(), -cellEnergyDerivative);
  addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(*_dynamicalMatrix, f1->otherCell(), f2->cell(), -cellEnergyDerivative);
  addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(*_dynamicalMatrix, f1->cell(), f2->cell(), cellEnergyDerivative);
  
  // 0 = pure shear xy
  addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(*_dynamicalMatrix, 0, f1->otherCell(), 
          _box.x*f2->periodicity().x*cellEnergyDerivative.columnX() - _box.y*f2->periodicity().y*cellEnergyDerivative.columnY());
  addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(*_dynamicalMatrix, 0, f1->cell(), 
        - _box.x*f2->periodicity().x*cellEnergyDerivative.columnX() + _box.y*f2->periodicity().y*cellEnergyDerivative.columnY());
  _dynamicalMatrix->add2ndDerivativeControlParameter(0, 
          _box.x*_box.x*f1->periodicity().x*f2->periodicity().x*cellEnergyDerivative.xx()
       -2*_box.x*_box.y*f1->periodicity().x*f2->periodicity().y*cellEnergyDerivative.xy()
        + _box.y*_box.y*f1->periodicity().y*f2->periodicity().y*cellEnergyDerivative.yy());

  // 1 = simple shear yx (normal,force)
  addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(*_dynamicalMatrix, 1, f1->otherCell(), _box.y*f2->periodicity().y*cellEnergyDerivative.columnX());
  addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(*_dynamicalMatrix, 1, f1->cell(), -_box.y*f2->periodicity().y*cellEnergyDerivative.columnX());
  _dynamicalMatrix->add2ndDerivativeControlParameter(1, _box.y*_box.y*f1->periodicity().y*f2->periodicity().y*cellEnergyDerivative.xx());

  // 2 = volume
  Vector3D periodicityAndBox1(f1->periodicity()*_boxDimensionsVector);
  Vector3D periodicityAndBox2(f2->periodicity()*_boxDimensionsVector);
  addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(*_dynamicalMatrix, 2, f1->otherCell(), cellEnergyDerivative*periodicityAndBox2);
  addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(*_dynamicalMatrix, 2, f1->cell(), -cellEnergyDerivative*periodicityAndBox2);
  _dynamicalMatrix->add2ndDerivativeControlParameter(2, periodicityAndBox1*cellEnergyDerivative*periodicityAndBox2);
}

inline void Cell::computeContributionsToDynamicalMatrix(Tessellation &t) /*argument only for debugging!*/ {
  // this needs to come before the vertices to prepare treatment of the results fed by the vertices
  for(DirectedFace *f : _faces) {
    f->createDoublyLinkedEdgeList();
  }
  for(VertexOfCell *v : _vertices) {
    v->computeContributionsToDynamicalMatrix(t);
  }
  for(DirectedFace *f : _faces) {
    f->computeContributionsToDynamicalMatrix(t);
  }
  for(DirectedEdgeOfCell *e : _edges) {
    if(e<e->conjugated()) {
      e->computeContributionsToDynamicalMatrix();
    }
  }
  
  // terms from first area derivative
  for(DirectedFace *f : _faces) {
    const double AreaNorm = f->area().norm();
    const Vector3D NormalizedArea = f->area()/AreaNorm;
    
    Matrix3x3 derivativeFF(f->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(f->derFaceAreaWrtCellCenterConnection()));
    derivativeFF -= dyadicProduct(f->derFaceAreaWrtCellCenterConnection()*NormalizedArea, f->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
    derivativeFF *= derOwnEnergyWrtSurface()/AreaNorm;
    derivativeFF += (derOwnEnergyWrtVolume()/6.0)*(f->derFaceAreaWrtCellCenterConnection()+f->derFaceAreaWrtCellCenterConnection().transposed());
//    Matrix3x3 derivativeFF(f->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(f->derFaceAreaWrtCellCenterConnection()));
//    derivativeFF -= dyadicProduct(f->derFaceAreaWrtCellCenterConnection()*NormalizedArea, f->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
//    derivativeFF *= 1.0/AreaNorm;
//    Matrix3x3 derivativeFF((1.0/6.0)*(f->derFaceAreaWrtCellCenterConnection()+f->derFaceAreaWrtCellCenterConnection().transposed()));
    _tessellation.addCellEnergyDerivativeToDynamicalMatrix(f, f, derivativeFF);
    
    // so far, surface terms only (not even accounting for the first derivatives of the surface)
    DirectedEdgeOfCell *e1 = f->firstEdge();
    do {
      Matrix3x3 derivativeE1F(e1->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(f->derFaceAreaWrtCellCenterConnection()));
      derivativeE1F -= dyadicProduct(e1->derFaceAreaWrtCellCenterConnection()*NormalizedArea, f->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
      derivativeE1F *= derOwnEnergyWrtSurface()/AreaNorm;
      derivativeE1F += (derOwnEnergyWrtVolume()/6.0)*e1->derFaceAreaWrtCellCenterConnection();
      _tessellation.addCellEnergyDerivativeToDynamicalMatrix(e1->conjugated()->face(), f, derivativeE1F);

      Matrix3x3 derivativeFE1(f->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(e1->derFaceAreaWrtCellCenterConnection()));
      derivativeFE1 -= dyadicProduct(f->derFaceAreaWrtCellCenterConnection()*NormalizedArea, e1->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
      derivativeFE1 *= derOwnEnergyWrtSurface()/AreaNorm;
      derivativeFE1 += (derOwnEnergyWrtVolume()/6.0)*e1->derFaceAreaWrtCellCenterConnection().transposed();
      _tessellation.addCellEnergyDerivativeToDynamicalMatrix(f, e1->conjugated()->face(), derivativeFE1);
      
      DirectedEdgeOfCell *e2 = f->firstEdge();
      do {
        Matrix3x3 derivativeE1E2(e1->derFaceAreaWrtCellCenterConnection().multiplyFromRightByMatrix22(e2->derFaceAreaWrtCellCenterConnection()));
        derivativeE1E2 -= dyadicProduct(e1->derFaceAreaWrtCellCenterConnection()*NormalizedArea, e2->derFaceAreaWrtCellCenterConnection()*NormalizedArea);
        derivativeE1E2 *= derOwnEnergyWrtSurface()/AreaNorm;
        _tessellation.addCellEnergyDerivativeToDynamicalMatrix(e1->conjugated()->face(), e2->conjugated()->face(), derivativeE1E2);
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
  for(DirectedFace *f1 : _faces) {
    for(DirectedFace *f2 : _faces) {
      Matrix3x3 derivativeF1F2(dyadicProduct(f1->_derCellVolumeWrtCellCenterConnection, f2->_derCellVolumeWrtCellCenterConnection));
      derivativeF1F2 *= derOwnEnergyWrtVolume2nd();
      derivativeF1F2 += dyadicProduct(derOwnEnergyWrtSurface2nd()*f1->_derCellSurfaceWrtCellCenterConnection, f2->_derCellSurfaceWrtCellCenterConnection);
      _tessellation.addCellEnergyDerivativeToDynamicalMatrix(f1, f2, derivativeF1F2);
    }
  }
}

inline void Cell::addFaceAreaVectorDerivativeToDynamicalMatrix(const DirectedFace* f1, const DirectedFace* f2, const DirectedFace* f, const Matrix3x3x3& areaVectorDerivative) {
  Matrix3x3 energyDerivative(areaVectorDerivative*((derOwnEnergyWrtSurface()/f->area().norm())*f->area()));
  energyDerivative += areaVectorDerivative*((derOwnEnergyWrtVolume()/6.0)*f->cellCenterConnectionVector());
//  Matrix3x3 energyDerivative(areaVectorDerivative*((1.0/f->area().norm())*f->area()));
//  Matrix3x3 energyDerivative(areaVectorDerivative*((1.0/6.0)*f->cellCenterConnectionVector()));
  _tessellation.addCellEnergyDerivativeToDynamicalMatrix(f1, f2, energyDerivative);
}

inline void DirectedFace::computeContributionsToDynamicalMatrix(Tessellation &t) /*argument only for debugging!*/ {
  const Matrix3x3x3 HalfEpsilon(0.5*Matrix3x3x3::Epsilon);
  DirectedEdgeOfCell *edge = _firstEdge;
  do {
    DirectedEdgeOfCell *nextEdge = edge->nextAroundFace();
    
    VertexOfCell *v1 = edge->vertex();
    VertexOfCell *v2 = nextEdge->vertex();
    for(DirectedEdgeOfCell *edgeV1 : v1->edges()) {
      for(DirectedEdgeOfCell *edgeV2 : v2->edges()) {
        _cell->addFaceAreaVectorDerivativeToDynamicalMatrix(edgeV1->face(), edgeV2->face(), this, 
              - HalfEpsilon
                  .multiplyFromLeftByMatrix12(edgeV1->_derVertexPositionWrtCellCenterConnection)
                  .multiplyIntoMiddleByMatrix22(edgeV2->_derVertexPositionWrtCellCenterConnection));
        _cell->addFaceAreaVectorDerivativeToDynamicalMatrix(edgeV2->face(), edgeV1->face(), this, 
                HalfEpsilon
                  .multiplyFromLeftByMatrix12(edgeV2->_derVertexPositionWrtCellCenterConnection)
                  .multiplyIntoMiddleByMatrix22(edgeV1->_derVertexPositionWrtCellCenterConnection));
      }
    }
    
    edge = nextEdge;
  } while(edge!=_firstEdge);
}

inline void DirectedFace::addVertexPositionDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const DirectedEdgeOfCell *e, const Matrix3x3x3 &vertexPositionDerivative) {
  _cell->addFaceAreaVectorDerivativeToDynamicalMatrix(f1, f2, this, 
          vertexPositionDerivative*antisymmetricMatrixFromVector(0.5*e->nextAroundFace()->vertex()->positionRelativeToCell()));
  _cell->addFaceAreaVectorDerivativeToDynamicalMatrix(f1, f2, this, 
          vertexPositionDerivative*antisymmetricMatrixFromVector(-0.5*e->_previousAroundFace->vertex()->positionRelativeToCell()));
}

inline void DirectedEdgeOfCell::computeContributionsToDynamicalMatrix() const {
  Cell *c = face()->cell();
  const CellType *parameters = c->type();
  if(parameters->edgeSprings) {
    Vector3D nEdgeVector(negativeLengthVector());
    double length = nEdgeVector.norm();
    double energyFactor = 2.0*parameters->edgeElasticity;
    Matrix3x3 derEnergyWrtEdgeVector(energyFactor*(1.0-_restLength/length)*Matrix3x3::Identity);
    derEnergyWrtEdgeVector += dyadicProduct(energyFactor*_restLength/(length*length*length) * nEdgeVector, nEdgeVector);

    // this vertex - this vertex:
    for(DirectedEdgeOfCell *e1 : vertex()->edges()) {
      for(DirectedEdgeOfCell *e2 : vertex()->edges()) {
        c->tesselation().addCellEnergyDerivativeToDynamicalMatrix(e1->face(), e2->face(), 
                e1->derVertexPositionWrtCellCenterConnection() * derEnergyWrtEdgeVector.multiplyFromRightByMatrix22( e2->derVertexPositionWrtCellCenterConnection() )
        );
      }
    }

    // this vertex - other vertex:
    for(DirectedEdgeOfCell *e1 : vertex()->edges()) {
      for(DirectedEdgeOfCell *e2 : conjugated()->vertex()->edges()) {
        Matrix3x3 der(e1->derVertexPositionWrtCellCenterConnection() * derEnergyWrtEdgeVector.multiplyFromRightByMatrix22( e2->derVertexPositionWrtCellCenterConnection() ));
        der *= -1;
        c->tesselation().addCellEnergyDerivativeToDynamicalMatrix(e1->face(), e2->face(), der);
        c->tesselation().addCellEnergyDerivativeToDynamicalMatrix(e2->face(), e1->face(), der.transposed());
      }
    }

    // other vertex - other vertex:
    for(DirectedEdgeOfCell *e1 : conjugated()->vertex()->edges()) {
      for(DirectedEdgeOfCell *e2 : conjugated()->vertex()->edges()) {
        c->tesselation().addCellEnergyDerivativeToDynamicalMatrix(e1->face(), e2->face(), 
                e1->derVertexPositionWrtCellCenterConnection() * derEnergyWrtEdgeVector.multiplyFromRightByMatrix22( e2->derVertexPositionWrtCellCenterConnection() )
        );
      }
    }
  }
}

inline void DirectedEdgeOfCell::addVertexPositionDerivativeToDynamicalMatrix(const DirectedFace* f1, const DirectedFace* f2, const Matrix3x3x3& vertexPositionDerivative) const {
  Cell *c = f1->cell();
  if(c->type()->edgeSprings) {
    c->tesselation().addCellEnergyDerivativeToDynamicalMatrix(f1, f2, vertexPositionDerivative * derEdgeEnergyWrtVertexPosition());
  }
}

inline void VertexOfCell::computeContributionsToDynamicalMatrix(Tessellation &t) /*argument only for debugging!*/ {
  for(int i=0; i<3; ++i) {
    const DirectedEdgeOfCell *edgeB(_edges[i]); 
    const DirectedEdgeOfCell *edgeC(_edges[(i+1)%3]); 
    const DirectedEdgeOfCell *edgeD(_edges[(i+2)%3]); 
    const DirectedFace *faceB(edgeB->face()); 
    const DirectedFace *faceC(edgeC->face()); 
    const DirectedFace *faceD(edgeD->face()); 
    const Vector3D &rab(faceB->cellCenterConnectionVector()); 
    const Vector3D &rac(faceC->cellCenterConnectionVector()); 
    const Vector3D &rad(faceD->cellCenterConnectionVector()); 

    // second derivative BB
    Matrix3x3x3 secondDerivativeBB(dyadicProduct(Matrix3x3::Identity, crossProduct(2*_positionRelativeToCellPositionDenominator*rac,rad)));
    secondDerivativeBB += dyadicProduct(edgeB->_derVertexPositionDenominatorWrtCellCenterConnection, edgeB->_derVertexPositionNumeratorWrtCellCenterConnection);
    secondDerivativeBB += dyadicProductInBetween(edgeB->_derVertexPositionDenominatorWrtCellCenterConnection, edgeB->_derVertexPositionNumeratorWrtCellCenterConnection);
    secondDerivativeBB += dyadicProduct(edgeB->_derVertexPositionDenominatorWrtCellCenterConnection/(0.5*_positionRelativeToCellPositionDenominator),
                                        edgeB->_derVertexPositionDenominatorWrtCellCenterConnection,
                                        _positionRelativeToCellPositionNumerator);
    for(DirectedEdgeOfCell *e : _edges) {
      e->face()->addVertexPositionDerivativeToDynamicalMatrix(faceB, faceB, e, secondDerivativeBB);
      e->addVertexPositionDerivativeToDynamicalMatrix(faceB, faceB, secondDerivativeBB);
    }


    // second derivative BC
    Matrix3x3x3 secondDerivativeOfNumeratorBC(dyadicProduct(-2*rab, antisymmetricMatrixFromVector(rad)));
    secondDerivativeOfNumeratorBC += dyadicProductInBetween(2*rac, antisymmetricMatrixFromVector(rad));
    secondDerivativeOfNumeratorBC += rad.normSq()*Matrix3x3x3::Epsilon;
    Matrix3x3 secondDerivativeOfDenominatorBC(dyadicProduct(edgeB->_derVertexPositionDenominatorWrtCellCenterConnection/(0.5*_positionRelativeToCellPositionDenominator),
                                                            edgeC->_derVertexPositionDenominatorWrtCellCenterConnection));
    secondDerivativeOfDenominatorBC -= antisymmetricMatrixFromVector(rad*(2*_positionRelativeToCellPositionDenominator*_positionRelativeToCellPositionDenominator));
    Matrix3x3x3 secondDerivativeBC(_positionRelativeToCellPositionDenominator*secondDerivativeOfNumeratorBC);
    secondDerivativeBC += dyadicProduct(edgeB->_derVertexPositionDenominatorWrtCellCenterConnection, edgeC->_derVertexPositionNumeratorWrtCellCenterConnection);
    secondDerivativeBC += dyadicProductInBetween(edgeC->_derVertexPositionDenominatorWrtCellCenterConnection, edgeB->_derVertexPositionNumeratorWrtCellCenterConnection);
    secondDerivativeBC += dyadicProduct(secondDerivativeOfDenominatorBC, _positionRelativeToCellPositionNumerator);
    for(DirectedEdgeOfCell *e : _edges) {
      e->face()->addVertexPositionDerivativeToDynamicalMatrix(faceB, faceC, e, secondDerivativeBC);
      e->addVertexPositionDerivativeToDynamicalMatrix(faceB, faceC, secondDerivativeBC);
    }

    
    // second derivative BD
    Matrix3x3x3 secondDerivativeOfNumeratorBD(dyadicProduct(2*rab, antisymmetricMatrixFromVector(rac)));
    secondDerivativeOfNumeratorBD -= rac.normSq()*Matrix3x3x3::Epsilon;
    secondDerivativeOfNumeratorBD -= dyadicProductInBetween(2*rad, antisymmetricMatrixFromVector(rac));
    Matrix3x3 secondDerivativeOfDenominatorBD(dyadicProduct(edgeB->_derVertexPositionDenominatorWrtCellCenterConnection/(0.5*_positionRelativeToCellPositionDenominator),
                                                            edgeD->_derVertexPositionDenominatorWrtCellCenterConnection));
    secondDerivativeOfDenominatorBD += antisymmetricMatrixFromVector(rac*(2*_positionRelativeToCellPositionDenominator*_positionRelativeToCellPositionDenominator));
    Matrix3x3x3 secondDerivativeBD(_positionRelativeToCellPositionDenominator*secondDerivativeOfNumeratorBD);
    secondDerivativeBD += dyadicProduct(edgeB->_derVertexPositionDenominatorWrtCellCenterConnection, edgeD->_derVertexPositionNumeratorWrtCellCenterConnection);
    secondDerivativeBD += dyadicProductInBetween(edgeD->_derVertexPositionDenominatorWrtCellCenterConnection, edgeB->_derVertexPositionNumeratorWrtCellCenterConnection);
    secondDerivativeBD += dyadicProduct(secondDerivativeOfDenominatorBD, _positionRelativeToCellPositionNumerator);
    for(DirectedEdgeOfCell *e : _edges) {
      e->face()->addVertexPositionDerivativeToDynamicalMatrix(faceB, faceD, e, secondDerivativeBD);
      e->addVertexPositionDerivativeToDynamicalMatrix(faceB, faceD, secondDerivativeBD);
    }
  }
}
