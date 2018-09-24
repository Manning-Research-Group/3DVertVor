#ifndef DIRECTEDEDGEOFCELL_INLINE_H
#define DIRECTEDEDGEOFCELL_INLINE_H

#include "CellType.h"
#include "DirectedEdgeOfCell.h"
#include "VertexOfCell-inline.h"

inline void DirectedEdgeOfCell::createDoublyLinkEdgeList() {
  // this is to prepare treatment of the results feeded by the vertices
  DirectedEdgeOfCell *edge = this;
  do {
    DirectedEdgeOfCell *nextEdge = edge->nextAroundFace();
    nextEdge->_previousAroundFace = edge;
    edge = nextEdge;
  } while(edge!=this);
}

inline Vector3D DirectedEdgeOfCell::negativeLengthVector() const {
  Vector3D v(vertex()->positionRelativeToCell());
  v -= conjugated()->vertex()->positionRelativeToCell();
  return v;
}

inline double DirectedEdgeOfCell::length() const {
  return negativeLengthVector().norm();
}

inline void DirectedEdgeOfCell::computeEdgeEnergyDerivative() {
  _derEdgeEnergyTermVertexPosition = negativeLengthVector();
  const CellType *parameters = _face->cell()->type();
  _derEdgeEnergyTermVertexPosition *= 2.0*parameters->edgeElasticity * (1.0 - _restLength/_derEdgeEnergyTermVertexPosition.norm());
}


#endif /* DIRECTEDEDGEOFCELL_INLINE_H */

