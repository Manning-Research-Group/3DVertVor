#ifndef CELL_H
#define	CELL_H

#include <vector>
#include <iostream>

#include "misc/geometry/Vector3D.h"
#include "misc/geometry/Matrix3x3x3.h"
#include "misc/geometry/Ellipsoid.h"

#include "povray/PovrayMesh.h"

#include "config.h"
#include "PeriodicBox.h"
#include "VertexOfCell.h"
#include "DirectedEdgeOfCell.h"
#include "DirectedFace.h"
#include "CellType.h"

class Tessellation;

class Cell {
public:
  Cell(Tessellation &tessellation, const CellType &parameters, const Vector3D &initialPosition);
  virtual ~Cell();
  Tessellation &tesselation() const { return _tessellation; }
  
  // cell type
  void setType(const CellType &parameters) { _type = &parameters; }
  const CellType *type() const { return _type; }
  
  // topology
  void resetTopology();
  DirectedFace *newFace(Cell *otherCell, DirectedEdgeOfCell *edge);
  DirectedEdgeOfCell *newEdge(VertexOfCell *vertex);
  VertexOfCell *newVertex(const Vector3D &voroPosition);
  const std::vector<DirectedFace*> &faces() const { return _faces; }
  const std::vector<DirectedEdgeOfCell*> &edges() const { return _edges; }
  const std::vector<VertexOfCell*> &vertices() const { return _vertices; }
  
  // geometry
  void setPosition(const Vector3D &p) { _position = p; }
  void moveIntoBox();
  const Vector3D &position() const { return _position; }
  Vector3D &position() { return _position; }
  const PeriodicityVector3D &periodicity() const { return _periodicity; }
  PeriodicityVector3D &periodicity() { return _periodicity; }
  const Vector3D positionWithoutBox() const;
  void computeGeometry();
  double volume() const { return _volume; }
  double volumeDifference() const { return _volume-_type->preferredVolume; }
  double surface() const { return _surface; }
  double surfaceDifference() const { return _surface-_type->preferredSurface; }
  void setEdgeRestLengthsTo(double l);
  void setEdgeRestLengthsToCurrentLengths();
  Ellipsoid fitEllipsoid() const { Vector3D cm; return fitEllipsoid(cm); }
  Ellipsoid fitEllipsoid(Vector3D &centerOfMass) const;
#if defined(CONFIG_CREATE_ANGULAR_NOISE_WALKING_ALONG_SPHERE) || defined(CONFIG_CREATE_ANGULAR_NOISE_RAIBLE_WINKLER)
  void setDirectionOfSelfPropulsion(const Vector3D &n) { _directionOfSelfPropulsion = n/n.norm(); }
  const Vector3D &directionOfSelfPropulsion() const { return _directionOfSelfPropulsion; }
#endif

  // dynamics
  void timeStep(const double deltaT);
  const Vector3D &speed() const { return _speed; }

  // energy
  double edgeEnergy() const;
  double additionalSurfaceEnergy() const;
  double energy() const {
    const double volumeDiff = volumeDifference();
    const double surfaceDiff = surfaceDifference();
    const double energy0 = _type->volumeElasticity*volumeDiff*volumeDiff + _type->surfaceElasticity*surfaceDiff*surfaceDiff + additionalSurfaceEnergy();
    if(_type->edgeSprings) {
      return energy0 + edgeEnergy();
    } else {
      return energy0;
    }
  }

  // force
  void computeCellEnergyDerivatives();
  void computeTotalEnergyDerivative();
  Vector3D force() const { return -_derTotalEnergyWrtCellPosition; }
  const Vector3D &derTotalEnergyWrtCellPosition() const { return _derTotalEnergyWrtCellPosition; }
  const double derOwnEnergyWrtVolume() const { return 2*_type->volumeElasticity*(_volume-_type->preferredVolume); }
  const double derOwnEnergyWrtVolume2nd() const { return 2*_type->volumeElasticity; }
  const double derOwnEnergyWrtSurface() const { return 2*_type->surfaceElasticity*(_surface-_type->preferredSurface); }
  const double derOwnEnergyWrtSurface2nd() const { return 2*_type->surfaceElasticity; }

  // consistency checks
  bool checkVertexPositionsAndVoroCondition(const std::vector<Cell*> &cells) const;
  
  // for drawing
  PovrayMesh<const VertexOfCell*> createMesh() const;
 
protected:
  // force
  Vector3D _derTotalEnergyWrtCellPosition;

private:
  Tessellation &_tessellation;
  const CellType *_type;

  // topology
  std::vector<DirectedFace*> _faces;
  std::vector<DirectedEdgeOfCell*> _edges;
  std::vector<VertexOfCell*> _vertices;
  
  // geometry
  Vector3D _position, _speed;
  PeriodicityVector3D _periodicity;
#if defined(CONFIG_CREATE_ANGULAR_NOISE_WALKING_ALONG_SPHERE) || defined(CONFIG_CREATE_ANGULAR_NOISE_RAIBLE_WINKLER)
  Vector3D _directionOfSelfPropulsion;
#endif
#ifdef CONFIG_CREATE_ANGULAR_NOISE_ANGLE_DYNAMICS
  double _selfPropulsionTheta, _selfPropulsionPhi;
#endif
  double _volume, _surface;

  // dynamical matrix
  void computeContributionsToDynamicalMatrix(Tessellation &t);  // argument only for debugging!
  void addFaceAreaVectorDerivativeToDynamicalMatrix(const DirectedFace *f1, const DirectedFace *f2, const DirectedFace *f, const Matrix3x3x3 &derivative);
  friend class Tessellation;
  friend class DirectedFace;
};

#endif	/* CELL_H */

