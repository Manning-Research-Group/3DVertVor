#include <set>

#include "Cell.h"
#include "Tessellation.h"

#include "DirectedEdgeOfCell-inline.h"

Cell::Cell(Tessellation &tessellation, const CellType &type, const Vector3D &initialPosition) 
  : _tessellation(tessellation), _type(&type), _position(initialPosition), _periodicity(0,0,0) {
  // initialize the direction of self-propulsion
#if defined(CONFIG_CREATE_ANGULAR_NOISE_WALKING_ALONG_SPHERE) || defined(CONFIG_CREATE_ANGULAR_NOISE_RAIBLE_WINKLER)
  double norm = 0.0;
  do {
    _directionOfSelfPropulsion = Random::gaussianVector3D();
    norm = _directionOfSelfPropulsion.norm();
  } while(norm<1e-8); // avoid division by zero
  _directionOfSelfPropulsion /= norm;
#endif
#ifdef CONFIG_CREATE_ANGULAR_NOISE_ANGLE_DYNAMICS
  Vector3D directionOfSelfPropulsion;
  do { directionOfSelfPropulsion = Random::gaussianVector3D(); } while(directionOfSelfPropulsion.norm()<1e-8);
  _selfPropulsionTheta = directionOfSelfPropulsion.theta();
  _selfPropulsionPhi = directionOfSelfPropulsion.phi();
#endif
}

Cell::~Cell() {
  resetTopology();
}

void Cell::resetTopology() {
  for(DirectedFace *f : _faces) {
    delete f;
  }
  _faces.clear();
  for(DirectedEdgeOfCell *e : _edges) {
    delete e;
  }
  _edges.clear();
  for(VertexOfCell *v : _vertices) {
    delete v;
  }
  _vertices.clear();
}


void Cell::computeGeometry() {
  _volume = 0;
  _surface = 0;
  for(unsigned int i=0; i<_faces.size(); ++i) {
    DirectedFace *f=_faces[i];
    _volume += f->cellCenterConnectionVector() * f->area();
    _surface += f->area().norm();
  }
  _volume /= 6;
}

Ellipsoid Cell::fitEllipsoid(Vector3D &centerOfMass) const {
  const Matrix3x3 InverseRegularTetrahedron(Vector3D(1,1,1)/sqrt(6), Vector3D(2,-1,-1)/sqrt(3), Vector3D(0, -1, 1));
  const double VolumeRegularTetrahedron = sqrt(2)/12.0;
  const Vector3D CenterRegularTetrahedron(0.25*sqrt(6), 0, 0);
  const Matrix3x3 IntegratedDyadicProductRegularTetrahedron(Matrix3x3::createDiagonalMatrix(0.4, 0.025, 0.025));
  
  double totalRelativeVolume = 0;
  centerOfMass.set(0,0,0);
  Matrix3x3 integratedDyadicProduct(Matrix3x3::Zero);
  
  for(DirectedFace *face : faces()) {
    Vector3D centerOfFace(0.5*face->cellCenterConnectionVector());
    
    DirectedEdgeOfCell *edge = face->firstEdge();
    VertexOfCell *lastVertex = edge->vertex();
    do {
      edge = edge->nextAroundFace();
      VertexOfCell *vertex = edge->vertex();
      
      const Matrix3x3 TetrahedralVectors(centerOfFace, lastVertex->positionRelativeToCell(), vertex->positionRelativeToCell());
      const Matrix3x3 TrafoFromRegularTetrahedron(TetrahedralVectors*InverseRegularTetrahedron);
      const double RelativeVolume = TrafoFromRegularTetrahedron.determinant();
      totalRelativeVolume += RelativeVolume;
      centerOfMass += RelativeVolume * (TrafoFromRegularTetrahedron * CenterRegularTetrahedron);
      integratedDyadicProduct += RelativeVolume * (TrafoFromRegularTetrahedron * IntegratedDyadicProductRegularTetrahedron * TrafoFromRegularTetrahedron.transposed());
      
      lastVertex = vertex;
    } while(edge!=face->firstEdge());
  }
  centerOfMass /= totalRelativeVolume;
//  integratedDyadicProduct /= totalRelativeVolume;
  integratedDyadicProduct *= VolumeRegularTetrahedron;
  
//  std::cout << "volume (old formula): " << volume() << std::endl;
//  std::cout << "volume (tetrahedra):  " << VolumeRegularTetrahedron*totalRelativeVolume << std::endl;
//  std::cout << "center of mass: " << centerOfMass << std::endl;
//  std::cout << "matrix: " << integratedDyadicProduct - dyadicProduct(VolumeRegularTetrahedron*totalRelativeVolume*centerOfMass, centerOfMass) << std::endl;
  
  return Ellipsoid(integratedDyadicProduct - dyadicProduct(VolumeRegularTetrahedron*totalRelativeVolume*centerOfMass, centerOfMass));
}

/*
 * Here, the moment of inertia considers a point mass at each vertex. Other formulations of
 * moment of inertia exist, but here we use the simplest which seems to work very well for our
 * needs of defining an ellipsoid for each cell. For more details, see Dobrovolskis (1996) - Inertia
 * of any polyhedron.
*/
EllipsoidByUnitPointMassPolyhedron Cell::fitEllipsoidByUnitPointMassPolyhedron() {
   
  // six terms of symmetric moment of intertia tensor
  double Ixx = 0;
  double Ixy = 0;
  double Ixz = 0;
  double Iyy = 0;
  double Iyz = 0;
  double Izz = 0;

  int vertexCounter = 0;
  for(DirectedFace *face : faces()) {
    
    DirectedEdgeOfCell *edge = face->firstEdge();
    do {
      edge = edge->nextAroundFace();
      VertexOfCell *vertex = edge->vertex();

      double x = vertex->positionRelativeToCell().x();
      double y = vertex->positionRelativeToCell().y();
      double z = vertex->positionRelativeToCell().z();
      
      Ixx += y*y + z*z;
      Iyy += x*x + z*z;
      Izz += x*x + y*y;
      Ixy += x*y;
      Ixz += x*z;
      Iyz += y*z;

      vertexCounter++;
      
    } while(edge!=face->firstEdge());
  } //end face loop
  
  Ixx /= vertexCounter;
  Iyy /= vertexCounter;
  Izz /= vertexCounter;
  Ixy /= vertexCounter;
  Ixz /= vertexCounter;
  Iyz /= vertexCounter;
  
  _unitPointMassMomentOfInertiaTensor = Matrix3x3(Vector3D(Ixx,-Ixy,-Ixz), Vector3D(-Ixy,Iyy,-Iyz), 
                                     Vector3D(-Ixz,-Iyz,Izz)); 
  
  return EllipsoidByUnitPointMassPolyhedron(_unitPointMassMomentOfInertiaTensor);
}

void Cell::setEdgeRestLengthsTo(double l) {
  for(DirectedEdgeOfCell *e : _edges) {
    e->setRestLengthTo(l);
  }
}

void Cell::setEdgeRestLengthsToCurrentLengths() {
  for(DirectedEdgeOfCell *e : _edges) {
    e->setRestLengthToCurrentLength();
  }
}


double Cell::edgeEnergy() const {
  double sumEdgeEnergies = 0.0;
  for(DirectedEdgeOfCell *e : edges()) {
    if(e<e->conjugated()) {
      double diff = e->differenceLengthRestLength();
      sumEdgeEnergies += diff*diff;
    }
  }
  return _type->edgeElasticity * sumEdgeEnergies;
}

double Cell::additionalSurfaceEnergy() const {
  double sum = 0.0;
  for(DirectedFace *f : _faces) {
    sum += f->additionalInterfacialTension() * f->area().norm();
  }
  return sum;
}

void Cell::computeCellEnergyDerivatives() {
  // compute geometrical derivatives of vertex positions
  for(VertexOfCell *v : vertices()) {
    v->computeVertexPositionDerivatives();
  }
  
  // compute derivatives of face area vectors
  for(DirectedFace *f : faces()) {
    f->computeFaceAreaDerivatives();
  }
  
  // compute derivative of edge energy
  for(DirectedEdgeOfCell *e : edges()) {
    e->computeEdgeEnergyDerivative();
  }
  
  // reset all derivatives
  for(DirectedFace *face : faces()) {
    face->_derCellSurfaceWrtCellCenterConnection.set(0,0,0);
    face->_derCellVolumeWrtCellCenterConnection.set(0,0,0);
    face->_derCellEnergyWrtCellCenterConnection.set(0,0,0);  //  first only contains the cell-cell interfacial energy derivative, then the rest
  }
  
  for(DirectedFace *face : faces()) {
    double areaNorm = face->area().norm();
    areaNorm = (areaNorm>0)?areaNorm:1.0;
    Vector3D areaDirection(face->area()/areaNorm);
    Vector3D areaDerivative(face->derFaceAreaWrtCellCenterConnection()*areaDirection);
    face->_derCellSurfaceWrtCellCenterConnection += areaDerivative;
    face->_derCellEnergyWrtCellCenterConnection += face->additionalInterfacialTension() * areaDerivative;
    face->_derCellVolumeWrtCellCenterConnection += face->area() + face->derFaceAreaWrtCellCenterConnection()*face->cellCenterConnectionVector();

    DirectedEdgeOfCell *edge = face->firstEdge();
    do {
      edge = edge->nextAroundFace();
      DirectedFace *otherFace = edge->conjugated()->face();
      Vector3D areaDerivative(edge->derFaceAreaWrtCellCenterConnection()*areaDirection);
      otherFace->_derCellSurfaceWrtCellCenterConnection += areaDerivative;
      otherFace->_derCellEnergyWrtCellCenterConnection += face->additionalInterfacialTension() * areaDerivative;
      otherFace->_derCellVolumeWrtCellCenterConnection += edge->derFaceAreaWrtCellCenterConnection()*face->cellCenterConnectionVector();
    } while(edge!=face->firstEdge());
  }

  for(DirectedFace *face : faces()) {
    face->_derCellVolumeWrtCellCenterConnection /= 6.0;
    face->_derCellEnergyWrtCellCenterConnection +=   derOwnEnergyWrtVolume()*face->_derCellVolumeWrtCellCenterConnection
                                                   + derOwnEnergyWrtSurface()*face->_derCellSurfaceWrtCellCenterConnection;
    if(face->_derCellEnergyWrtCellCenterConnection.isnan()) {
      std::cerr << "Cell::computeEnergyDerivatives: energy derivative is nan: " << face->_derCellEnergyWrtCellCenterConnection 
                << "; volume derivative: " << face->_derCellVolumeWrtCellCenterConnection
                << "; surface derivative: " << face->_derCellSurfaceWrtCellCenterConnection 
                << "; volume: " << _volume
                << "; surface: " << _surface
                << std::endl;
      std::cerr << "Cell::computeEnergyDerivatives: "
                << "formula: " << derOwnEnergyWrtVolume()*face->_derCellVolumeWrtCellCenterConnection + derOwnEnergyWrtSurface()*face->_derCellSurfaceWrtCellCenterConnection
                << "; derOwnEnergyWrtVolume(): " << derOwnEnergyWrtVolume()
                << "; derOwnEnergyWrtSurface(): " << derOwnEnergyWrtSurface()
                << "; kV: " << _type->volumeElasticity
                << "; V0: " << _type->preferredVolume
                << "; kS: " << _type->surfaceElasticity
                << "; S0: " << _type->preferredSurface
                << std::endl;
      exit(1);
    }
  }
  
  // add edge derivatives
  if(_type->edgeSprings) {
    for(DirectedFace *face : faces()) {
      Vector3D derTotalEdgeEnergyWrtCellCenterConnection(0,0,0);
      DirectedEdgeOfCell *fEdge = face->firstEdge();
      do {
        Vector3D derTotalEdgeEnergyWrtVertexPosition(0,0,0);
        for(DirectedEdgeOfCell *e : fEdge->vertex()->edges()) {
          derTotalEdgeEnergyWrtVertexPosition += e->derEdgeEnergyWrtVertexPosition();
        }
        derTotalEdgeEnergyWrtCellCenterConnection += fEdge->derVertexPositionWrtCellCenterConnection() * derTotalEdgeEnergyWrtVertexPosition;
        fEdge = fEdge->nextAroundFace();
      } while(fEdge!=face->firstEdge());
      face->_derCellEnergyWrtCellCenterConnection += derTotalEdgeEnergyWrtCellCenterConnection;
    }
  }
}

void Cell::computeTotalEnergyDerivative() {
  // derivative of own energy with respect to own position
  _derTotalEnergyWrtCellPosition.set(0,0,0);
  for(DirectedFace *face : faces()) {
    _derTotalEnergyWrtCellPosition -= face->derCellEnergyWrtCellCenterConnection();
  }
  // derivative of neighbor's energies with respect to my position
  for(DirectedFace *face : faces()) {
    _derTotalEnergyWrtCellPosition += face->conjugated()->derCellEnergyWrtCellCenterConnection();
  }
}

void Cell::timeStep(const double deltaT) {
  // *** update direction of self-propulsion ***
  
#ifdef CONFIG_CREATE_ANGULAR_NOISE_WALKING_ALONG_SPHERE
  // get noise vector
  Vector3D noise(Random::gaussianVector3D(sqrt(2*_type->angularDiffusion*deltaT)));
  // project perpendicular to direction of self-propulsion
  Vector3D normalVector(noise - _directionOfSelfPropulsion*(_directionOfSelfPropulsion*noise));
  // get its norm
  double beta = normalVector.norm();
  // avoid division by zero
  if(beta>1e-10) { 
    // change direction of self-propulsion, taking non-linear effects into account
    _directionOfSelfPropulsion = cos(beta)*_directionOfSelfPropulsion + sin(beta)/beta*normalVector;
    // re-normalize to avoid numerical drift
    _directionOfSelfPropulsion /= _directionOfSelfPropulsion.norm();
  }
#endif          

#ifdef CONFIG_CREATE_ANGULAR_NOISE_RAIBLE_WINKLER
  // get noise vector
  Vector3D noise(Random::gaussianVector3D(sqrt(2*_type->angularDiffusion*deltaT)));
  // project perpendicular to direction of self-propulsion
  Vector3D normalVector(noise - _directionOfSelfPropulsion*(_directionOfSelfPropulsion*noise));
  // change direction of self-propulsion
  _directionOfSelfPropulsion += normalVector;
  // re-normalize
  _directionOfSelfPropulsion /= _directionOfSelfPropulsion.norm();
#endif          
  
#ifdef CONFIG_CREATE_ANGULAR_NOISE_ANGLE_DYNAMICS
  // set noise standard deviation
  double NoiseStdDevition = sqrt(2*_type->angularDiffusion*deltaT);
  // add noise; TODO: deal with divisions by zero, and situations close to it
  _selfPropulsionTheta += _type->angularDiffusion/tan(_selfPropulsionTheta) + Random::normal(NoiseStdDevition);
  _selfPropulsionPhi += Random::normal(NoiseStdDevition)/sin(_selfPropulsionTheta);
#endif  
  
  // *** update position ***
//  Vector3D deltaPosition = deltaT * (
//#if defined(CONFIG_CREATE_ANGULAR_NOISE_WALKING_ALONG_SPHERE) || defined(CONFIG_CREATE_ANGULAR_NOISE_RAIBLE_WINKLER)
//          + _parameters.Speed * _directionOfSelfPropulsion
//#endif          
//#ifdef CONFIG_CREATE_ANGULAR_NOISE_ANGLE_DYNAMICS
//          + _parameters.Speed * Vector3D::fromNormThetaPhi(1.0, _selfPropulsionTheta, _selfPropulsionPhi)
//#endif  
//          + force());
//  _position += deltaPosition;
//  _positionWithoutBox += deltaPosition;
//  _speed = deltaPosition/deltaT;

  //DEBUG
//  if(_directionOfSelfPropulsion.x()!=_directionOfSelfPropulsion.x()) {
//    std::cerr << "Cell::timeStep:  direction of self propulsion is nan!" << std::endl;
//  }

//  Vector3D forceValue(force());
//  if(forceValue.x()!=forceValue.x()) {
//    std::cerr << "Cell::timeStep:  force is nan!" << std::endl;
//  }
  
#if defined(CONFIG_CREATE_ANGULAR_NOISE_WALKING_ALONG_SPHERE) || defined(CONFIG_CREATE_ANGULAR_NOISE_RAIBLE_WINKLER)
  Vector3D deltaPosition = _type->speed * _directionOfSelfPropulsion;
#endif          
#ifdef CONFIG_CREATE_ANGULAR_NOISE_ANGLE_DYNAMICS
  Vector3D deltaPosition = _type->speed * Vector3D::fromNormThetaPhi(1.0, _selfPropulsionTheta, _selfPropulsionPhi);
#endif  
  if(_type->experiencesForce) {
    deltaPosition += force();
  }
  deltaPosition *= deltaT;
  CellType type3;
  if(_cellstype==2){
    //_position += 1e-2*Vector3D(0, 0, 1)*deltaT;
    _position += 1e-4*Vector3D(0, 0, 1);
  }
  else{_position += deltaPosition;}
  
  _speed = deltaPosition/deltaT;
//  if(_position.x()!=_position.x()) {
//    std::cerr << "Cell::timeStep:  position is nan!" << std::endl;
//  }
}

bool Cell::checkVertexPositionsAndVoroCondition(const std::vector<Cell*> &cells) const {
  for(unsigned int i=0; i<_vertices.size(); ++i) {
    VertexOfCell *v=_vertices[i];
    if(!v->checkVoroPositionAndVoroCondition(_tessellation.box(), cells)) {
      return false;
    }
  }
  return true;
}

PovrayMesh<const VertexOfCell*> Cell::createMesh() const {
  PovrayMesh<const VertexOfCell*> mesh;
  for(const VertexOfCell *v : _vertices) mesh.addPoint(v, v->position());
  for(const DirectedFace *f : _faces) {
    DirectedEdgeOfCell *e1 = f->firstEdge();
    DirectedEdgeOfCell *e2 = e1->nextAroundFace();
    DirectedEdgeOfCell *e3 = e2->nextAroundFace();
    do {
      mesh.addTriangle(e1->vertex(), e2->vertex(), e3->vertex());
      e2 = e3;
      e3 = e3->nextAroundFace();
    } while(e3!=f->firstEdge());
  }
  return mesh;
}
