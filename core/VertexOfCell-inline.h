/* 
 * File:   VertexOfCell-inline.h
 * Author: arbeit
 *
 * Created on October 29, 2015, 7:49 PM
 */

#ifndef VERTEXOFCELL_INLINE_H
#define	VERTEXOFCELL_INLINE_H

#include "VertexOfCell.h"
#include "Cell.h"

inline const Vector3D VertexOfCell::voroPosition() const { 
  return _cell->position() + _voroPositionRelativeToCellCenter;
}

inline const Vector3D VertexOfCell::position() const { 
  return _cell->position() + _positionRelativeToCellPosition;
}

#endif	/* VERTEXOFCELL_INLINE_H */

