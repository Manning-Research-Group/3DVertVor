#ifndef CELL_INLINE_H
#define	CELL_INLINE_H

#include "Cell.h"
#include "Tessellation.h"

inline DirectedFace *Cell::newFace(Cell *otherCell, DirectedEdgeOfCell *edge) { 
  DirectedFace *f = new DirectedFace(_tessellation.box(), this, otherCell, edge); 
  _faces.push_back(f); 
  return f; 
}

inline DirectedEdgeOfCell *Cell::newEdge(VertexOfCell *vertex) {
  DirectedEdgeOfCell *e = new DirectedEdgeOfCell(vertex);
  _edges.push_back(e);
  return e;
}

inline VertexOfCell *Cell::newVertex(const Vector3D &voroPosition) { 
  VertexOfCell *v = new VertexOfCell(this, voroPosition); 
  _vertices.push_back(v); 
  return v; 
}

inline void Cell::moveIntoBox() { 
  _periodicity += _tessellation.box().computePositionPeriodicityAndRest(_position);
}

inline const Vector3D Cell::positionWithoutBox() const { 
  return _position + _tessellation.box().periodicityOffset(_periodicity); 
}



#endif	/* CELL_INLINE_H */

