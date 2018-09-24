#include "misc/geometry/Matrix3x3.h"
#include "PeriodicBox.h"
#include "DirectedFace.h"
#include "DirectedEdgeOfCell.h"

#include "VertexOfCell-inline.h"

void DirectedFace::bufferCellCenterConnectionVector() {
  _bufferCenterConnectionVector = _otherCell->position()-_cell->position()+_box.periodicityOffset(_periodicity);
}

void DirectedFace::computeGeometry() {
  // the area vector, perpendicular to the face
  // and oriented such that it points from this cell towards the other cell
  _area.set(0,0,0);
  Vector3D lastRelativeVertexPosition(firstEdge()->vertex()->positionRelativeToCell());
  DirectedEdgeOfCell *edge=firstEdge();
  do {
    edge = edge->nextAroundFace();
    VertexOfCell *v=edge->vertex();
    Vector3D curRelativeVertexPosition(v->positionRelativeToCell());
    _area += crossProduct(lastRelativeVertexPosition, curRelativeVertexPosition);
    lastRelativeVertexPosition = curRelativeVertexPosition;
  } while(edge!=firstEdge());
  _area *= -0.5;
}

void DirectedFace::computeFaceAreaDerivatives() {
  // compute derivative of area with respect to this face's cell center connection vector (Ccv))
  // and with respect to neighboring face's cell center connection vector (Nccv)
  _derFaceAreaWrtCellCenterConnection = Matrix3x3::Zero;
  Vector3D previousVertexPosition(firstEdge()->vertex()->positionRelativeToCell());
  DirectedEdgeOfCell *edge, *nextEdge = firstEdge()->nextAroundFace();
  do {
    edge = nextEdge;
    nextEdge = edge->nextAroundFace();

    const Matrix3x3 &firstVertexPositionWrtNccv(edge->conjugated()->nextAroundFace()->derVertexPositionWrtCellCenterConnection());
    const Vector3D &firstVertexPosition(edge->vertex()->positionRelativeToCell());
    
    const Matrix3x3 &secondVertexPositionWrtNccv(edge->conjugated()->derVertexPositionWrtCellCenterConnection());
    const Vector3D &secondVertexPosition(nextEdge->vertex()->positionRelativeToCell());

    const Vector3D &followingVertexPosition(nextEdge->nextAroundFace()->vertex()->positionRelativeToCell());

    edge->_derFaceAreaWrtCellCenterConnection = 
              firstVertexPositionWrtNccv * antisymmetricMatrixFromVector(previousVertexPosition)
            - firstVertexPositionWrtNccv * antisymmetricMatrixFromVector(secondVertexPosition)
            + secondVertexPositionWrtNccv * antisymmetricMatrixFromVector(firstVertexPosition)
            - secondVertexPositionWrtNccv * antisymmetricMatrixFromVector(followingVertexPosition);
    edge->_derFaceAreaWrtCellCenterConnection *= -0.5;

    const Matrix3x3 &firstVertexPositionWrtCcv(edge->derVertexPositionWrtCellCenterConnection());
    _derFaceAreaWrtCellCenterConnection -= firstVertexPositionWrtCcv * antisymmetricMatrixFromVector(secondVertexPosition);
    _derFaceAreaWrtCellCenterConnection += firstVertexPositionWrtCcv * antisymmetricMatrixFromVector(previousVertexPosition);
    
    previousVertexPosition = firstVertexPosition;
  } while(edge!=firstEdge());
  _derFaceAreaWrtCellCenterConnection *= -0.5;
}
