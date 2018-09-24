#ifndef VERTEX_H
#define	VERTEX_H

#include <vector>
#include <iostream>

#include "misc/geometry/Vector3D.h"
#include "misc/geometry/Matrix3x3.h"
#include "povray/PovrayGenerator.h"
#include "PeriodicBox.h"
#include "DirectedEdgeOfCell.h"

class Cell;
class DirectedFace;
class Tessellation;  // only for debugging!

class VertexOfCell {
public:
  VertexOfCell(Cell *c, const Vector3D &vp) : _cell(c), _voroPositionRelativeToCellCenter(vp) {}

  const Cell *cell() const { return _cell; }
  std::vector<DirectedEdgeOfCell*> edges() const { return _edges; }
  void addEdge(DirectedEdgeOfCell *e) { _edges.push_back(e); }
  const Vector3D voroPosition() const;

  void computeGeometry();
  void computeVertexPositionDerivatives();
  const Vector3D position() const;
  const Vector3D positionWithRespectTo(const Vector3D& referencePosition) const { return _positionRelativeToCellPosition + referencePosition; }
  const Vector3D &positionRelativeToCell() const {return _positionRelativeToCellPosition; }
  const Vector3D &positionRelativeToCellNumerator() const {return _positionRelativeToCellPositionNumerator; }
  double positionRelativeToCellDenominator() const {return _positionRelativeToCellPositionDenominator; }

  bool checkVoroPositionAndVoroCondition(const PeriodicBox &box, const std::vector<Cell*> &cells) const;
  
private:
  const static double MaxRadiusSqDeviationForVoroCheck;
  const static double MaxDistanceSqDeviationForComparisonToVoro;

  Cell *_cell;
  std::vector<DirectedEdgeOfCell*> _edges;

  // geometry
  Vector3D _voroPositionRelativeToCellCenter;
  double   _normalization; // own computation
  Vector3D _positionRelativeToCellPosition; // own computation
  Vector3D _positionRelativeToCellPositionNumerator; // own computation
  double   _positionRelativeToCellPositionDenominator; // own computation
  
  // dynamical matrix
  void computeContributionsToDynamicalMatrix(Tessellation &t);  // argument only for debugging!
  friend class Cell;
};

#endif	/* VERTEX_H */

