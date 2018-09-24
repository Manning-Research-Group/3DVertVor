#ifndef DIRECTEDEDGEOFCELL_H
#define DIRECTEDEDGEOFCELL_H

#include "misc/geometry/Matrix3x3.h"
#include "VertexOfCell.h"

class VertexOfCell;
class DirectedFace;

class DirectedEdgeOfCell {
public:
  DirectedEdgeOfCell(VertexOfCell *vertex) : _vertex(vertex), _conjugated(NULL), _nextAroundFace(NULL), _face(NULL), _restLength(1.0) {}; 
  void setConjugatedAndNext(DirectedEdgeOfCell *conjugated, DirectedEdgeOfCell *nextAroundFace) { _conjugated=conjugated; _nextAroundFace=nextAroundFace; };
  void setFace(DirectedFace *face) { _face=face; };
  
  VertexOfCell *vertex() const { return _vertex; }
  DirectedEdgeOfCell *conjugated() const { return _conjugated; }
  DirectedEdgeOfCell *nextAroundFace() const { return _nextAroundFace; }
  DirectedEdgeOfCell *previousAroundFace() const { return _previousAroundFace; }
  DirectedFace *face() const { return _face; }

  void createDoublyLinkEdgeList();

  Vector3D negativeLengthVector() const;
  double length() const;
  double differenceLengthRestLength() const { return length()-_restLength; }
  double restLength() const { return _restLength; }
  void setRestLengthTo(double l) { _restLength = l; }
  void setRestLengthToCurrentLength() { _restLength = length(); }
  
  void computeEdgeEnergyDerivative();
  const Vector3D &derEdgeEnergyWrtVertexPosition() const { return _derEdgeEnergyTermVertexPosition; }
  const Matrix3x3 &derVertexPositionWrtCellCenterConnection() const { return _derVertexPositionWrtCellCenterConnection; }
  const Matrix3x3 &derVertexPositionNumeratorWrtCellCenterConnection() const { return _derVertexPositionNumeratorWrtCellCenterConnection; }
  const Vector3D &derVertexPositionDenominatorWrtCellCenterConnection() const { return _derVertexPositionDenominatorWrtCellCenterConnection; }
  const Matrix3x3 &derFaceAreaWrtCellCenterConnection() const { return _derFaceAreaWrtCellCenterConnection; }
  
private:
  VertexOfCell *_vertex;
  DirectedEdgeOfCell *_conjugated;
  DirectedEdgeOfCell *_nextAroundFace; // in CW order!!!
  DirectedEdgeOfCell *_previousAroundFace; // in CW order!!! Only needed to compute dynamical matrix.
  DirectedFace *_face;
  
  // geometry / energy
  double _restLength;
  
  // first derivative of L^2 with respect to the position of _vertex
  Vector3D _derEdgeEnergyTermVertexPosition;

  // first (row) index:  cell center connection vector of face 
  // second (column) index:  position of vertex
  Matrix3x3 _derVertexPositionNumeratorWrtCellCenterConnection;
  Vector3D _derVertexPositionDenominatorWrtCellCenterConnection;
  Matrix3x3 _derVertexPositionWrtCellCenterConnection;

  // first (row) index:  cell center connection vector of conjugated edge's face 
  // second (column) index:  area vector of this face
  Matrix3x3 _derFaceAreaWrtCellCenterConnection;
  
  void computeContributionsToDynamicalMatrix() const;
  void addVertexPositionDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const Matrix3x3x3 &derivative) const;
  
  friend class VertexOfCell;
  friend class DirectedFace;
  friend class Cell;
};

#endif /* DIRECTEDEDGEOFCELL_H */

