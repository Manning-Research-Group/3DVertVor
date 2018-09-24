#ifndef DIRECTEDFACE_H
#define	DIRECTEDFACE_H

#include <stdlib.h>
#include <vector>
#include <iostream>

#include "misc/geometry/Vector3D.h"
#include "misc/geometry/Matrix3x3.h"
#include "povray/Color.h"
#include "povray/PovrayGenerator.h"
#include "VertexOfCell.h"
#include "DirectedEdgeOfCell.h"
#include "Cell.h"
#include "PeriodicBox.h"

class Cell;
class Tessellation;   // only for debugging!

/**
 * A Directed face belongs to a Cell and its area vector points outside of this cell. The vertices are sorted with respect to a left hand rule.
 * The conjugated face corresponds to the same face but it is associated with the cell  on the other side and correspondingly, its area normal points into the opposite direction.
 */
class DirectedFace {
public:
  DirectedFace(const PeriodicBox &box, Cell *c, Cell *oc, DirectedEdgeOfCell *firstEdge) : _box(box), _cell(c), _otherCell(oc), _conjugated(NULL), _firstEdge(firstEdge) {}

  Cell *cell() const { return _cell; }
  Cell *otherCell() const { return _otherCell; }
  void setConjugated(DirectedFace *c) { _conjugated = c; }
  DirectedFace *conjugated() const { return _conjugated; }
  DirectedEdgeOfCell *firstEdge() const { return _firstEdge; }
  
  // topology
  void createDoublyLinkedEdgeList() { firstEdge()->createDoublyLinkEdgeList(); };
  void updateAdditionalInterfacialTension();
  double additionalInterfacialTension() const { return _additionalInterfacialTension; }

  // helper to figure out the periodicity
  const Vector3D &averageVertexVoronoiPositions() const { return _bufferedAverageVoronoiPositions; }
  void bufferAverageVertexVoronoiPositions() {
    _bufferedAverageVoronoiPositions.set(0,0,0);
    int numEdges = 0;
    DirectedEdgeOfCell *edge = firstEdge();
    do {
      _bufferedAverageVoronoiPositions += edge->vertex()->voroPosition();
      ++numEdges;
      edge = edge->nextAroundFace();
    } while(edge!=firstEdge());
    _bufferedAverageVoronoiPositions /= numEdges;
  };
  
  void setPeriodicity(const PeriodicityVector3D &p) { _periodicity = p; }
  PeriodicityVector3D &periodicity() { return _periodicity; }
  const PeriodicityVector3D &periodicity() const { return _periodicity; }
  
  void bufferCellCenterConnectionVector();
  const Vector3D &cellCenterConnectionVector() const { return _bufferCenterConnectionVector; }
  
  // depends on vertex positions and cellCenterConnectionVectors
  void computeGeometry();
  void computeFaceAreaDerivatives();
  const Vector3D &area() const { return _area; }
  const Matrix3x3 &derFaceAreaWrtCellCenterConnection() const { return _derFaceAreaWrtCellCenterConnection; }
  const Vector3D &derCellEnergyWrtCellCenterConnection() const { return _derCellEnergyWrtCellCenterConnection; }
  const Vector3D &derCellSurfaceWrtCellCenterConnection() const { return _derCellSurfaceWrtCellCenterConnection; }
  const Vector3D &derCellVolumeWrtCellCenterConnection() const { return _derCellVolumeWrtCellCenterConnection; }
  
 
private:
  const PeriodicBox &_box;
  Cell *_cell, *_otherCell;
  DirectedFace *_conjugated;
  DirectedEdgeOfCell *_firstEdge;
  PeriodicityVector3D _periodicity;
  
  Vector3D _bufferedAverageVoronoiPositions;
  Vector3D _bufferCenterConnectionVector;
  
  // additional interfacial tension, buffered here for speed
  double _additionalInterfacialTension;
  
  // depends on vertex positions:
  Vector3D _area;
  
  // first (row) index:  cell center connection vector of this face 
  // second (column) index:  area vector of this face
  Matrix3x3 _derFaceAreaWrtCellCenterConnection;

  // derivatives of own cell's energy with respect to cell center connection vector
  Vector3D _derCellSurfaceWrtCellCenterConnection;
  Vector3D _derCellVolumeWrtCellCenterConnection;
  Vector3D _derCellEnergyWrtCellCenterConnection;
  
  friend class Cell;

  // dynamical matrix
  void computeContributionsToDynamicalMatrix(Tessellation &t);  // argument only for debugging!
  void addVertexPositionDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const DirectedEdgeOfCell *e, const Matrix3x3x3 &derivative);
  friend class VertexOfCell;
};

#endif	/* DIRECTEDFACE_H */

