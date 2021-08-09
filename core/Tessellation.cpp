#include <stdlib.h>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <math.h>

#include "misc/other/misc-math.h"
#include "misc/other/fileIo.h"

#include "PeriodicBox.h"
#include "Cell.h"
#include "DirectedFace.h"
#include "VertexOfCell.h"
#include "DirectedEdgeOfCell.h"
#include "Tessellation.h"

#include "VertexOfCell-inline.h"
#include "Cell-inline.h"
#include "DirectedFace-inline.h"

const double Tessellation::NumericalDerivativeDifference = 1e-8;

Tessellation::Tessellation(PeriodicBox &box) : _dynamicalMatrix(NULL), _box(box), _time(0.0), _topologicalElementsPresent(false) {
}

Tessellation::~Tessellation() {
  cleanup();
}

void Tessellation::cleanup() {
  for(Cell *c : _cells) delete c;
  _cells.clear();
  if(_dynamicalMatrix) {
    delete _dynamicalMatrix;
  }
}

// ***** creation ***** //

void Tessellation::addCell(const CellType &cellType, const Vector3D &position) {
  _cells.push_back(new Cell(*this, cellType, position));
}

void Tessellation::addCellsAtRandomPositions(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells
  for(int i=0; i<NumberOfCells; ++i) {
    addCell(cellType, Vector3D(Random::uniform(_box.x), Random::uniform(_box.y), Random::uniform(_box.z)));
//    std::cout << _cells[i].position() << std::endl;
  }
}

void Tessellation::addCellsTopHalf(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells
  for(int i=0; i<NumberOfCells; ++i) {
    addCell(cellType, Vector3D(Random::uniform(_box.x), Random::uniform(_box.y), Random::uniform(0.5*_box.z)));
//    std::cout << _cells[i].position() << std::endl;
  }
}

void Tessellation::addCellsBottomHalf(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells
  for(int i=0; i<NumberOfCells; ++i) {
    addCell(cellType, Vector3D(Random::uniform(_box.x), Random::uniform(_box.y), 0.5*_box.z+Random::uniform(0.5*_box.z)));
//    std::cout << _cells[i].position() << std::endl;
  }
}

void Tessellation::squareLattice(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells

  int sidelength = int(cbrt(2*NumberOfCells));
  for(int i=0; i<NumberOfCells; ++i) {
    //  if(i%) {
    //}
    double xpos = 0.0;
    double ypos = 0.0;
    double zpos = 0.0;
    int ii,jj,kk;
    int sidedist = sidelength/2;

    ii = i%(sidedist);
    jj = int(i/sidedist)%sidelength;
    kk = int(i/(sidedist*sidelength));

    //std::cout << i << std::endl;
    //std::cout << ii << std::endl;
    //std::cout << jj << std::endl;
    //std::cout << kk << std::endl;
    //=2*B3+IF(MOD(B4-1,2)=0,1,0)+0.5

    xpos = 2*ii+0.5+((jj-1)%2==0)*(kk%2==0)+(jj%2==0)*((kk-1)%2==0)+0.1*Random::uniform(_box.x)/sidelength;
    ypos = jj+0.5+0.1*Random::uniform(_box.y)/sidelength;
    zpos = kk+0.5+0.1*Random::uniform(_box.z)/sidelength;

    //std::cout << i << std::endl;
    //std::cout << xpos << std::endl;
    //std::cout << ypos << std::endl;
    //std::cout << zpos << std::endl;

    addCell(cellType, Vector3D(xpos, ypos, zpos));
    //addCell(cellType, Vector3D(Random::uniform(_box.x), Random::uniform(_box.y), 0.5*_box.z+Random::uniform(0.5*_box.z)));
//    std::cout << _cells[i].position() << std::endl;
  }
}

void Tessellation::squareLattice2(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells

  int sidelength = int(cbrt(2*NumberOfCells));
  for(int i=0; i<NumberOfCells; ++i) {
    double xpos = 0.0;
    double ypos = 0.0;
    double zpos = 0.0;
    int ii,jj,kk;
    int sidedist = sidelength/2;

    ii = i%(sidedist);
    jj = int(i/sidedist)%sidelength;
    kk = int(i/(sidedist*sidelength));

    xpos = 2*ii+0.5+((jj-1)%2==0)*((kk-1)%2==0)+(jj%2==0)*((kk)%2==0)+0.1*Random::uniform(_box.x)/sidelength;
    ypos = jj+0.5+0.1*Random::uniform(_box.y)/sidelength;
    zpos = kk+0.5+0.1*Random::uniform(_box.z)/sidelength;

    addCell(cellType, Vector3D(xpos, ypos, zpos));
  }
}

void Tessellation::basementmembrane(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells
  for(int i=0; i<NumberOfCells; ++i) {
    addCell(cellType, Vector3D(Random::uniform(_box.x), Random::uniform(_box.y), Random::uniform(0.25*_box.z)));
//    std::cout << _cells[i].position() << std::endl;
  }
}

void Tessellation::basal(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells
  for(int i=0; i<NumberOfCells; ++i) {
    addCell(cellType, Vector3D(Random::uniform(_box.x), Random::uniform(_box.y), 0.25*_box.z+Random::uniform(0.25*_box.z)));
//    std::cout << _cells[i].position() << std::endl;
  }
}

void Tessellation::suprabasal(const CellType &cellType, const int NumberOfCells) {
  // create randomly placed cells
  for(int i=0; i<NumberOfCells; ++i) {
    addCell(cellType, Vector3D(Random::uniform(_box.x), Random::uniform(_box.y), 0.5*_box.z+Random::uniform(0.5*_box.z)));
//    std::cout << _cells[i].position() << std::endl;
  }
}


// ***** minimization ***** //

void Tessellation::setDofs(MinimizerWithDerivative &minimizer) {
  minimizer.clearDofs();
  for(Cell *c : _cells) {
    minimizer.addDof(c->position().x(), c->derTotalEnergyWrtCellPosition().x());
    minimizer.addDof(c->position().y(), c->derTotalEnergyWrtCellPosition().y());
    minimizer.addDof(c->position().z(), c->derTotalEnergyWrtCellPosition().z());
  }
}

bool Tessellation::relax(MinimizerWithDerivative &minimizer) {
  return minimizer.minimize([this]()->double { 
    this->computeTopologyAndGeometry();
    return this->energy();
  }, [this]()->void {
    this->computeTopologyAndGeometry();
    this->computeEnergyDerivatives();
  }, [this]()->double { 
    this->computeTopologyAndGeometry();
    this->computeEnergyDerivatives();
    return this->energy();
  });
}


// ***** dynamics ***** //

void Tessellation::timeStep(const double deltaT) {
  for(unsigned int ci=0; ci<_cells.size(); ++ci) {
    Cell *c=_cells[ci];
    c->timeStep(deltaT);
  }
  _time += deltaT;
}

// ***** geometry ***** //

void Tessellation::computeVertexPositions() {
  for(Cell *c : _cells) {
    for(unsigned int fi=0; fi<c->faces().size(); ++fi) {
      DirectedFace *f=c->faces()[fi];
      // NOTE: this has to come *after* fixing the periodicity for this face!
      f->bufferCellCenterConnectionVector();
    }
    // NOTE: this has to come *after* buffering the cell center connection vectors!
    for(unsigned int i=0; i<c->vertices().size(); ++i) {
      VertexOfCell *v=c->vertices()[i];
      v->computeGeometry();
    }
  }
}

void Tessellation::computeGeometry() {
  computeVertexPositions();

  for(Cell *c : _cells) {
    // NOTE: this has to come *after* having computed the vertex positions and *after* having the cell center connection vectors!
    for(unsigned int i=0; i<c->faces().size(); ++i) {
      DirectedFace *f=c->faces()[i];
      f->computeGeometry();
    }
    // NOTE: this has to come *after* having computed the directed bond area and the cell center connection vectors!
    c->computeGeometry();
  }
}

bool Tessellation::topologyChanged() {
  if(!_topologicalElementsPresent) {
    return true;
  }
    
  for(Cell *c : _cells) {
    c->moveIntoBox();
  }
  
  computeVertexPositions();
          
  // refill mailbox
  _mailbox.refill(_box, _cells);

  // check
  double maxDimension = (_box.x>_box.y)?_box.x:_box.y;
  maxDimension = (maxDimension>_box.z)?maxDimension:_box.z;
  for(Cell *c : _cells) {
    for(VertexOfCell *v : c->vertices()) {
      // compute radius of sphere around vertex according to old tessellation
      const double SqRadius = v->positionRelativeToCell().normSq();
      double radius = sqrt(SqRadius);
      if(radius>maxDimension) radius = maxDimension;
      
      // check if cell position lies within this radius
      const BoxIndex MaxIndices(_mailbox.maxIndicesForRadius(radius));
      const Vector3D VPosition(v->position());
      const BoxIndex VIndex(_mailbox.boxIndexForPosition(VPosition));
      BoxIndex cIndex;
      for(int iOffset=-MaxIndices.i; iOffset<=MaxIndices.i; ++iOffset) {
        cIndex.i = modulo(VIndex.i + iOffset, _mailbox.nx());
        for(int jOffset=-MaxIndices.j; jOffset<=MaxIndices.j; ++jOffset) {
          cIndex.j = modulo(VIndex.j + jOffset, _mailbox.ny());
          for(int kOffset=-MaxIndices.k; kOffset<=MaxIndices.k; ++kOffset) {
            cIndex.k = modulo(VIndex.k + kOffset, _mailbox.nz());
            
            // loop over all cells at this index
            for(Cell *oc : _mailbox.box(cIndex)) {
              Vector3D distanceVector(VPosition-oc->position());
              _box.minimalDistanceVector(distanceVector);
              if(distanceVector.normSq()<SqRadius) {
                // check whether it's one of the cells already touching
                if(oc==c) {
                  // already touching -> ignore
                  continue;
                }
                bool isAmongstDefiningCells = false;
                for(auto *e : v->edges()) {
                  if(oc==e->face()->otherCell()) {
                    // already touching -> ignore
                    isAmongstDefiningCells = true;
                    break;
                  }
                }
                if(!isAmongstDefiningCells) {
                  // Then it's another cell inside of the sphere -> topology changed
                  return true;
                }
              }
            }
          }
        }
      }
    }
  }
  
  return false;
}

bool Tessellation::topologyChangedBruteForce() {
  if(!_topologicalElementsPresent) {
    return true;
  }

  for(Cell *c : _cells) {
    c->moveIntoBox();
  }
  
  computeVertexPositions();

  // check
  for(Cell *c : _cells) {
    for(VertexOfCell *v : c->vertices()) {
      // compute radius of sphere around vertex according to old tessellation
      const double SqRadius = v->positionRelativeToCell().normSq();
      
      // check if any cell position lies within this radius
      const Vector3D VPosition(v->position());
      for(Cell *oc : _cells) {
        Vector3D distanceVector(VPosition-oc->position());
        _box.minimalDistanceVector(distanceVector);
        if(distanceVector.normSq()<SqRadius) {
          // check whether it's one of the cells already touching
          if(oc==c) {
            // already touching -> ignore
            continue;
          }
          bool isAmongstDefiningCells = false;
          for(auto *e : v->edges()) {
            if(oc==e->face()->otherCell()) {
              // already touching -> ignore
              isAmongstDefiningCells = true;
              break;
            }
          }
          if(!isAmongstDefiningCells) {
            // Then it's another cell inside of the sphere -> topology changed
            return true;
          }
        }
      }
    }
  }
  
  return false;
}

bool Tessellation::computeTopologyAndGeometry() {
  // check if topology changed
  bool tc = topologyChanged();
  if(tc) {
    // DEBUG
//    std::cout << "Topology changed!" << std::endl;
    // first, remove everything up to the cell positions
    resetTopology();
    // voronoi tesselation
    setTopologyByVoronoiTesselation();
    // update interfacial tensions
    for(Cell *c : cells()) {
      for(DirectedFace *f : c->faces()) {
        f->updateAdditionalInterfacialTension();
      }
    }
  }

  // geometry
  computeGeometry();
  
  return tc;
//  return true;
}

bool Tessellation::checkVoronoiTesselationFromVoroLibrary() const {
  for(unsigned int i=0; i<_cells.size(); ++i) {
    Cell *c=_cells[i];
    if(!c->checkVertexPositionsAndVoroCondition(_cells)) {
      return false;
    }
  }
  return true;
}

void Tessellation::affinelyRescaleBox(double factorX, double factorY, double factorZ) {
  _box.x *= factorX;  _box.y *= factorY;  _box.z *= factorZ;
  for(Cell *c : _cells) {
    Vector3D p(c->position());
    c->position().set(factorX*p.x(), factorY*p.y(), factorZ*p.z());
  }
}

void Tessellation::setEdgeRestLengthsTo(double l) {
  for(Cell *c : _cells) {
    c->setEdgeRestLengthsTo(l);
  }
}

void Tessellation::setEdgeRestLengthsToCurrentLengths() {
  for(Cell *c : _cells) {
    c->setEdgeRestLengthsToCurrentLengths();
  }
}


// ***** energy ***** //

double Tessellation::energy() const {
  double e=0.0;
  for(const Cell *c : cells()) e += c->energy();
  return e;
}
  
  
// ***** forces ***** //

void Tessellation::computeEnergyDerivatives() {
  // set cell energy derivatives
  for(Cell *c : _cells) {
    c->computeCellEnergyDerivatives();
  }
  for(Cell *c : _cells) {
    c->computeTotalEnergyDerivative();
  }
  
  // set periodic box derivatives
  _box.derEnergyWrtShearYx = 0.0;
  for(Cell *c : _cells) {
    for(DirectedFace *f : c->faces()) {
      _box.derEnergyWrtShearYx += f->periodicity().y*_box.y*f->derCellEnergyWrtCellCenterConnection().x();
    }
  }
}

#define STRESS_COMPONENT_COMPUTATION(normalAxis,forceAxis) { \
    double tension=0.0; \
    for(Cell *c : _cells) { \
      for(DirectedFace *f : c->faces()) { \
        tension += f->periodicity().normalAxis*_box.normalAxis*f->derCellEnergyWrtCellCenterConnection().forceAxis(); \
      } \
    } \
    return tension/_box.volume(); \
  }

double Tessellation::stressXx() const STRESS_COMPONENT_COMPUTATION(x,x)
double Tessellation::stressYy() const STRESS_COMPONENT_COMPUTATION(y,y)
double Tessellation::stressZz() const STRESS_COMPONENT_COMPUTATION(z,z)
double Tessellation::stressYx() const { return _box.derEnergyWrtShearYx/_box.volume(); }
double Tessellation::pressure() const { return -(stressXx() + stressYy() + stressZz())/3.0; }
double Tessellation::pressureAlt() const { 
  Vector3D boxDimensionsVector(_box.x, _box.y, _box.z);
  double pressure = 0;
  for(Cell *c : _cells) {
    for(DirectedFace *f : c->faces()) {
      Vector3D periodicityAndBox(f->periodicity()*boxDimensionsVector);
      pressure -= periodicityAndBox*f->derCellEnergyWrtCellCenterConnection();
    }
  }
  return pressure/3.0/_box.volume(); 
}

Matrix3x3 Tessellation::numericalDerVertexPositionWrtCellPosition(VertexOfCell *v, Cell *c) {
  Vector3D cp0(c->position());
  Vector3D vp0(v->positionRelativeToCell());

  c->setPosition(cp0+Vector3D(NumericalDerivativeDifference,0,0));
  computeGeometry();
  Vector3D dvpx(v->positionRelativeToCell()-vp0);
  
  c->setPosition(cp0+Vector3D(0,NumericalDerivativeDifference,0));
  computeGeometry();
  Vector3D dvpy(v->positionRelativeToCell()-vp0);
  
  c->setPosition(cp0+Vector3D(0,0,NumericalDerivativeDifference));
  computeGeometry();
  Vector3D dvpz(v->positionRelativeToCell()-vp0);
  
  c->setPosition(cp0);
  computeGeometry();

  return Matrix3x3(Vector3D(dvpx.x(), dvpy.x(), dvpz.x()), Vector3D(dvpx.y(), dvpy.y(), dvpz.y()), Vector3D(dvpx.z(), dvpy.z(), dvpz.z()))/NumericalDerivativeDifference;
}

Matrix3x3 Tessellation::numericalDerFaceAreaWrtCellPosition(DirectedFace *f, Cell *c) {
  Vector3D cp0(c->position());
  Vector3D area0(f->area());

  c->setPosition(cp0+Vector3D(NumericalDerivativeDifference,0,0));
  computeGeometry();
  Vector3D dAreax(f->area()-area0);
  
  c->setPosition(cp0+Vector3D(0,NumericalDerivativeDifference,0));
  computeGeometry();
  Vector3D dAreay(f->area()-area0);
  
  c->setPosition(cp0+Vector3D(0,0,NumericalDerivativeDifference));
  computeGeometry();
  Vector3D dAreaz(f->area()-area0);
  
  c->setPosition(cp0);
  computeGeometry();

  return Matrix3x3(Vector3D(dAreax.x(), dAreay.x(), dAreaz.x()), Vector3D(dAreax.y(), dAreay.y(), dAreaz.y()), Vector3D(dAreax.z(), dAreay.z(), dAreaz.z()))/NumericalDerivativeDifference;
}

Vector3D Tessellation::numericalDerCellEnergyWrtCellPosition(Cell *ce, Cell *cp) {
  Vector3D p0(cp->position());
  double e0 = ce->energy();

  cp->setPosition(p0+Vector3D(NumericalDerivativeDifference,0,0));
  computeGeometry();
  double ex=ce->energy()-e0;
  
  cp->setPosition(p0+Vector3D(0,NumericalDerivativeDifference,0));
  computeGeometry();
  double ey=ce->energy()-e0;
  
  cp->setPosition(p0+Vector3D(0,0,NumericalDerivativeDifference));
  computeGeometry();
  double ez=ce->energy()-e0;
  
  cp->setPosition(p0);
  computeGeometry();
  
  return Vector3D(ex,ey,ez)/NumericalDerivativeDifference;
}

Vector3D Tessellation::numericalForceOnCell(Cell *c) {
  Vector3D p0(c->position());
  double e0 = energy();

  c->setPosition(p0+Vector3D(NumericalDerivativeDifference,0,0));
  computeGeometry();
  double ex=energy()-e0;
  
  c->setPosition(p0+Vector3D(0,NumericalDerivativeDifference,0));
  computeGeometry();
  double ey=energy()-e0;
  
  c->setPosition(p0+Vector3D(0,0,NumericalDerivativeDifference));
  computeGeometry();
  double ez=energy()-e0;
  
  c->setPosition(p0);
  computeGeometry();
  
  return -Vector3D(ex,ey,ez)/NumericalDerivativeDifference;
}

double Tessellation::totalCellForceNormSq() const {
  double sumNormSq = 0.0;
  for(Cell *c : cells()) {
    sumNormSq += c->force().normSq();
  }
  return sumNormSq;
}

  
// ***** dynamical matrix ***** //
  
double Tessellation::cellPairModulus(const Cell *c1, const Cell *c2, const double cutoffEVal, double cutoffEVecComponentSq) const {
  if((cutoffEVal>=0) && (cutoffEVecComponentSq<0)) {
    cutoffEVecComponentSq = cutoffEVal;
  }
  Vector3D dist(_box.minimalDistanceVector(c2->position()-c1->position()));
  int i1 = _cellToDynamicalMatrixIndex.at(c1);
  int i2 = _cellToDynamicalMatrixIndex.at(c2);
  double sum = 0.0;
  for(int q=0; q<_dynamicalMatrix->eigenValues().rows(); ++q) {
    const double omegaSq = _dynamicalMatrix->eigenValues()(q);
    const DynamicalMatrixOld::Vector &ev = _dynamicalMatrix->eigenVectors().col(q);
    const double element = (ev(i2)-ev(i1))*dist.x() + (ev(i2+1)-ev(i1+1))*dist.y() + (ev(i2+2)-ev(i1+2))*dist.z();
    if(omegaSq>=cutoffEVal) {
      sum += element*element/omegaSq;
    } else if(element*element>=cutoffEVecComponentSq) {
      return 0.0;
    }
  }
  return dist.normSq()/sum;
}

//double Tessellation::cellPairModulusNumerical(const Cell *c1const, const Cell *c2const) {
//  // this is BAD, but it's just a numerical check and we restore the positions
//  Cell *c1=(Cell*)c1const, *c2=(Cell*)c2const;
//  
//  // set up minimizer
//  GslMinimizer cg;
//  cg.setInitialStepSize(0.01);
//  cg.setLineMinimizationTolerance(0);
//  cg.setForceTolerancePerDof(1e-8);
//  cg.setMaxIterationsPerDof(100);
////  cg.setLogInterval(1);
//  
//  // save current cell positions
//  std::vector<Vector3D> oldPositions;
//  for(Cell *c : _cells) oldPositions.push_back(c->position());
//
//  // compute periodicity and distance
//  Vector3D distVec(c2->position()-c1->position());
//  PeriodicityVector3D periodicity(_box.computeDistancePeriodicityAndRest(distVec));
//
//  // define rotated coordinate system and initial distance
//  double dist0 = distVec.norm();
//  double dist = dist0;
//  double phi = atan2(distVec.y(), distVec.x());
//  double theta = atan2(distVec.z(), sqrt(distVec.x()*distVec.x() + distVec.y()*distVec.y()));
//
//  // forces on rotated coordinates
//  double derEnergyWrtDist, derEnergyWrtPhi, derEnergyWrtTheta;
//  Vector3D derEnergyWrtPositionOfCell1;
//  
//  // update functions
//  std::function<void()> updateCellPositions = [this, &dist, &phi, &theta, periodicity, c1, c2]() -> void {
//    double cosTheta = cos(theta);
//    c2->position() = c1->position() + dist*Vector3D(cosTheta*cos(phi), cosTheta*sin(phi), sin(theta)) + this->_box.periodicityOffset(periodicity);
//  };
//  std::function<void()> recomputeForces = [&derEnergyWrtDist, &derEnergyWrtPhi, &derEnergyWrtTheta, &derEnergyWrtPositionOfCell1, &dist, &phi, &theta, periodicity, c1, c2]() -> void {
//    double cosTheta = cos(theta);
//    double sinTheta = sin(theta);
//    double cosPhi = cos(phi);
//    double sinPhi = sin(phi);
//    derEnergyWrtPositionOfCell1 = c1->derTotalEnergyWrtCellPosition() + c2->derTotalEnergyWrtCellPosition();
//    derEnergyWrtDist = c2->derTotalEnergyWrtCellPosition()*Vector3D(cosTheta*cosPhi, cosTheta*sinPhi, sinTheta);
//    derEnergyWrtPhi = dist*c2->derTotalEnergyWrtCellPosition()*Vector3D(-cosTheta*sinPhi, cosTheta*cosPhi, 0);
//    derEnergyWrtTheta = dist*c2->derTotalEnergyWrtCellPosition()*Vector3D(-sinTheta*cosPhi, -sinTheta*sinPhi, cosTheta);
//  };
//  
//  // get force
//  dist = dist0;
//  updateCellPositions();
//  computeGeometry();
//  computeEnergyDerivatives();
//  recomputeForces();
//  const double derEnergyWrtDist0 = derEnergyWrtDist;
//
//  // change length
//  const double Delta = 1e-4;
//  dist = dist0 + Delta;
//
//  // set up degrees of freedom
//  cg.clearDofs();
//  for(Cell *c : _cells) {
//    if(((c!=c1) && (c!=c2))) {
//      cg.addDof(c->position().x(), c->derTotalEnergyWrtCellPosition().x());
//      cg.addDof(c->position().y(), c->derTotalEnergyWrtCellPosition().y());
//      cg.addDof(c->position().z(), c->derTotalEnergyWrtCellPosition().z());
//    }
//  }
//  cg.addDof(c1->position().x(), derEnergyWrtPositionOfCell1.x());
//  cg.addDof(c1->position().y(), derEnergyWrtPositionOfCell1.y());
//  cg.addDof(c1->position().z(), derEnergyWrtPositionOfCell1.z());
//  cg.addDof(phi, derEnergyWrtPhi);
//  cg.addDof(theta, derEnergyWrtTheta);
//  
//  // minimize
//  cg.minimize([this,updateCellPositions,recomputeForces]()->double { 
//      updateCellPositions();
//      this->computeGeometry();
//      return this->energy();
//    }, [this,updateCellPositions,recomputeForces]()->void {
//      updateCellPositions();
//      this->computeGeometry();
//      this->computeEnergyDerivatives();
//      recomputeForces();
//    }, [this,updateCellPositions,recomputeForces]()->double { 
//      updateCellPositions();
//      this->computeGeometry();
//      this->computeEnergyDerivatives();
//      recomputeForces();
//      return this->energy();
//    });
//    
//  // compute new force
//  updateCellPositions();
//  computeGeometry();
//  computeEnergyDerivatives();
//  recomputeForces();
//  const double derEnergyWrtDist1 = derEnergyWrtDist;
//
//  // reset to old positions
//  for(unsigned int i=0; i<_cells.size(); ++i) _cells[i]->setPosition(oldPositions[i]);
//  
//  // return second energy derivative
//  return (derEnergyWrtDist1-derEnergyWrtDist0)/Delta;
//}

double Tessellation::cellPairModulusX(const Cell *c1, const Cell *c2, const double cutoffEVal, double cutoffEVecComponentSq) const {
  if((cutoffEVal>=0) && (cutoffEVecComponentSq<0)) {
    cutoffEVecComponentSq = cutoffEVal;
  }  
  int i1 = this->_cellToDynamicalMatrixIndex.at(c1);
  int i2 = this->_cellToDynamicalMatrixIndex.at(c2);
  double sum = 0.0;
  for(int q=0; q<_dynamicalMatrix->eigenValues().rows(); ++q) {
    const double omegaSq = _dynamicalMatrix->eigenValues()(q);
    const DynamicalMatrixOld::Vector &ev = _dynamicalMatrix->eigenVectors().col(q);
    const double element = ev(i2)-ev(i1);
    if(omegaSq>=cutoffEVal) {
      sum += element*element/omegaSq;
    } else if(element*element>=cutoffEVecComponentSq) {
      return 0.0;
    }
  }
  return 1.0/sum;  
}

double Tessellation::cellPairModulusY(const Cell *c1, const Cell *c2, const double cutoffEVal, double cutoffEVecComponentSq) const {
  if((cutoffEVal>=0) && (cutoffEVecComponentSq<0)) {
    cutoffEVecComponentSq = cutoffEVal;
  }  
  int i1 = this->_cellToDynamicalMatrixIndex.at(c1);
  int i2 = this->_cellToDynamicalMatrixIndex.at(c2);
  double sum = 0.0;
  for(int q=0; q<_dynamicalMatrix->eigenValues().rows(); ++q) {
    const double omegaSq = _dynamicalMatrix->eigenValues()(q);
    const DynamicalMatrixOld::Vector &ev = _dynamicalMatrix->eigenVectors().col(q);
    const double element = ev(i2+1)-ev(i1+1);
    if(omegaSq>=cutoffEVal) {
      sum += element*element/omegaSq;
    } else if(element*element>=cutoffEVecComponentSq) {
      return 0.0;
    }
  }
  return 1.0/sum;  
}

double Tessellation::cellPairModulusZ(const Cell *c1, const Cell *c2, const double cutoffEVal, double cutoffEVecComponentSq) const {
  if((cutoffEVal>=0) && (cutoffEVecComponentSq<0)) {
    cutoffEVecComponentSq = cutoffEVal;
  }  
  int i1 = this->_cellToDynamicalMatrixIndex.at(c1);
  int i2 = this->_cellToDynamicalMatrixIndex.at(c2);
  double sum = 0.0;
  for(int q=0; q<_dynamicalMatrix->eigenValues().rows(); ++q) {
    const double omegaSq = _dynamicalMatrix->eigenValues()(q);
    const DynamicalMatrixOld::Vector &ev = _dynamicalMatrix->eigenVectors().col(q);
    const double element = ev(i2+2)-ev(i1+2);
    if(omegaSq>=cutoffEVal) {
      sum += element*element/omegaSq;
    } else if(element*element>=cutoffEVecComponentSq) {
      return 0.0;
    }
  }
  return 1.0/sum;  
}

// ***** drawing ***** //

void Tessellation::saveAsImageDefault(const std::string &filename, int width, int height, const Vector3D &CameraPositionOffset, const Vector3D &CameraLooksAt,
                                    const Color &Background) const {
  const double Scaling = cubicRoot(_box.volume());
  saveAsImage(filename, [this,Scaling](PovrayGenerator &g){
    g.addLight(Vector3D(10*Scaling,0,0), Color(0.75,0.70,0.70)); 
    g.addLight(Vector3D(0,-10*Scaling,0), Color(0.30,0.30,0.30)); 
    g.addLight(Vector3D(0,0,10*Scaling), Color(0.95,0.90,0.90)); 
    this->drawEdgesDefault(g);
    this->drawVerticesDefault(g);
    this->drawCellPositionsDefault(g);
  }, width, height, CameraPositionOffset, CameraLooksAt);
}

void Tessellation::saveAsImage(const std::string &filename, std::function<void(PovrayGenerator &g)> drawElements, int width, int height, const Vector3D &CameraPositionOffset, const Vector3D &CameraLooksAt,
                                    const Color &Background) const {
  const double Scaling = cubicRoot(_box.volume());
  PovrayGenerator g;
  g.start(Background);
  Vector3D ref(_box.x*CameraLooksAt.x(), _box.y*CameraLooksAt.y(), _box.z*CameraLooksAt.z());
  g.setCamera(Scaling*CameraPositionOffset+ref, ref, double(width)/height);
  drawElements(g);
  g.create(filename, width);
}

void Tessellation::drawCellPositionsDefault(PovrayGenerator &g, const Color &color, double radius) const {
  drawCells(g, color, [radius](PovrayGenerator &g, const Cell *c){ g.drawSphere(c->position(), radius); });
}

void Tessellation::drawCells(PovrayGenerator &g, const Color &color, std::function<void(PovrayGenerator &g, const Cell *c)> drawCell) const {
  g.startUnion();
  for(const Cell *c : cells()) drawCell(g, c);
  g.endUnionWithColor(color);
}

void Tessellation::drawEdgesDefault(PovrayGenerator &g, const Color &color, double radius) const {
  drawEdges(g, color, [radius](PovrayGenerator &g, const Cell *c, const DirectedEdgeOfCell *e){ 
      VertexOfCell *v1 = e->vertex();
      VertexOfCell *v2 = e->nextAroundFace()->vertex();
      g.drawCylinder(v1->positionWithRespectTo(c->position()), v2->positionWithRespectTo(c->position()), radius);
  });
}

void Tessellation::drawEdges(PovrayGenerator &g, const Color &color, std::function<void(PovrayGenerator &g, const Cell *c, const DirectedEdgeOfCell *e)> drawEdge) const {
  g.startUnion();
  for(const Cell *c : cells()) {
    for(const DirectedFace *f : c->faces()) {
      DirectedEdgeOfCell *edge=f->firstEdge();
      do {
        drawEdge(g, c, edge);
        edge = edge->nextAroundFace();
      } while(edge!=f->firstEdge());
    }
  }
  g.endUnionWithColor(color);
}

void Tessellation::drawVerticesDefault(PovrayGenerator &g, const Color &color, double radius) const {
  drawVertices(g, color, [radius](PovrayGenerator &g, const Cell *c, const VertexOfCell *v){ 
    g.drawSphere(v->positionWithRespectTo(c->position()), radius);
  });
}

void Tessellation::drawVertices(PovrayGenerator &g, const Color &color, std::function<void(PovrayGenerator &g, const Cell *c, const VertexOfCell *v)> drawVertex) const {
  g.startUnion();
  for(const Cell *c : cells()) {
    for(const VertexOfCell *v : c->vertices()) drawVertex(g, c, v);
  }
  g.endUnionWithColor(color);
}

// ***** creation by loading box and positions from dump ***** //
bool Tessellation::loadFromDump(const CellType &cellParameters, const std::string &filename) {
  return loadFromDumpHelper(filename, [this,&cellParameters](const Vector3D &initialPosition) { return new Cell(*this, cellParameters, initialPosition); });
}

// (Default LAMMPS dump format, see lammps.sandia.gov/doc/dump.html )
/*
ITEM: TIMESTEP
1530100
ITEM: NUMBER OF ATOMS
54
ITEM: BOX BOUNDS xy xz yz pp pp pp
0.0 8.0 0.0
0.0 8.0 0.0
0.0 8.0 0.0
ITEM: ATOMS id x y z
2 .00250 .00250 .00250
...
*/
bool Tessellation::loadFromDumpHelper(const std::string &filename, std::function<Cell*(const Vector3D &initialPosition)> createCell) {
  cleanup();
  
  //Open the input file
  std::ifstream Input_File(filename.c_str());	
  if(!(Input_File.is_open() && Input_File.good())) {
      std::cerr << "Was not able to open " << filename << std::endl;
    return false;
  }

  // Parameters for the dump style to be read
  const int nheader=9, lineofnspheres=4;
  std::string style1cols("ITEM: ATOMS id x y z "), tmpstr;
  char line[256];
  int ncells, index;
  double minx,miny,minz, maxx,maxy,maxz, tiltxy,tiltxz,tiltyz;
  double x,y,z;

  // In this dump file style, skip the top lines (TIMESTEP, VALUE, NUMBER OF ATOMS)
  for (int iread = 0; iread < lineofnspheres-1; iread++)
    Input_File.getline(line,256);

  // Read number of cells
  Input_File >> ncells >> tmpstr;
  std::cout << ncells << " is the number of particles in the top of dump file " << std::endl;

  //Skip line that describes "ITEM: BOX BOUNDS pp pp pp"
  Input_File.getline(line,256);

  // Read box dimensions
  Input_File >> minx >> maxx >> tiltxy;
  Input_File >> miny >> maxy >> tiltxz;
  Input_File >> minz >> maxz >> tiltyz;
  if ((tiltxz != 0) || (tiltyz != 0)) {
    std::cerr << "No box tilt in xz or yz directions currently allowed." << std::endl;
    return -1;
  }
  _box.z = maxz - minz;
  _box.y = maxy - miny;
  _box.x = maxx - minx - std::abs(tiltxy);

  // And skip the end line character
  Input_File.getline(line,256);
  //std::cout << line << std::endl;

  // Skip the rest of the header down to the column description
  for (int iread = 0; iread < nheader-lineofnspheres-4; iread++)
    Input_File.getline(line,256);
  std::string linestring(line);
  if (0 != style1cols.compare(linestring)) {
    std::cerr << "Expected cols: " << std::endl << style1cols << std::endl <<
                " but found " << std::endl << linestring << std::endl;
    std::cerr << "(Difference could be white space at end of line.)" << std::endl;
  }

  for(int iparticle = 0 ; iparticle < ncells ; iparticle++) {

    Input_File >> index >> x >> y >> z;
    if ((_box.x < x) || (_box.y < y) || (_box.z < z) || (x < 0) || (y < 0) || (z < 0))
      std::cerr << "Read in particles that are outside of the box." << std::endl;

    // IDs begin at 1 whereas array indexing begins at 0
    _cells.push_back(createCell(Vector3D(x, y, z)));
  }

  Input_File.close();
  std::cout << "Finished reading dump file" << filename << std::endl;

  return true;
}

bool Tessellation::writeToDump(const std::string &filename) {
  
  // Open the input file
  std::ofstream Output_File(filename.c_str());	
  if(!(Output_File.is_open() && Output_File.good())) {
      std::cerr << "Was not able to open " << filename << std::endl;
    return false;
  }

  // (LAMMPS dump format, see lammps.sandia.gov/doc/dump.html )
  Output_File << "ITEM: TIMESTEP" << std::endl;
  Output_File << "0" << std::endl;
  Output_File << "ITEM: NUMBER OF ATOMS" << std::endl;
  Output_File << _cells.size() << std::endl;

  // Because box is tri-clinic (sheared), BOX BOUNDS style is:
  Output_File << "ITEM: BOX BOUNDS xy xz yz pp pp pp" << std::endl;
  // pp referes to fully periodic in the x, y, or z directions.
  // Give bounding-box limits, then tilt factors.  
  // If PeriodicBox class includes tiltxz or tiltyz, will include here.
  double tiltxy = _box.shearYx*_box.y;
  double minxbound = std::min(0.0,             tiltxy);
  double maxxbound = std::max(_box.x, _box.x + tiltxy);
  Output_File << minxbound << " " << maxxbound << " " << tiltxy << std::endl;
  Output_File << 0.0 << " " << _box.y << " " << 0.0 << std::endl;
  Output_File << 0.0 << " " << _box.z << " " << 0.0 << std::endl;

  Output_File << "ITEM: ATOMS id x y z" << std::endl;
  int cellid = 1;
  for(Cell *c : _cells) {
    Vector3D position = c->position();
    Output_File << cellid++ << " " << 
            position.x() << " " << position.y() << " " << position.z() << std::endl;
  }

  Output_File.close();
  std::cout << "Closed dump file " << filename << std::endl;
  return true;
}

#if USE_NETCDF

bool Tessellation::loadFromNetCdfFileWithCellParameters(const std::string &Path, CellType &parameters, int record) {
  cleanup();

  const int Dim = 3;

  std::cout << "Loading from NetCDF file " << Path << "..." << std::endl;

  NcError err(NcError::verbose_fatal);

  // open file:
  NcFile file(Path.c_str(), NcFile::ReadOnly);
  if(!file.is_valid()) {
    std::cerr << "Problem loading file!" << std::endl;
    return false;
  }
  
  // get dimensions
  NcDim *recDim, *dimDim, *dm2Dim, *nPDim, *dofDim, *strDim;
  if(!(recDim = file.get_dim("rec"))) return false;
  if(!(dimDim = file.get_dim("dim"))) return false;
  if(!(dm2Dim = file.get_dim("dim2"))) return false;
  if(!(nPDim  = file.get_dim("NP"))) return false;
  if(!(dofDim = file.get_dim("dof"))) return false;
  if(!(strDim = file.get_dim("StringSize"))) return false;

  // check record size
  const int NumRecords = recDim->size();
  if(record>=NumRecords) {
    std::cerr << "Record " << record << " not found! Only " << NumRecords << " present in file." << std::endl;
    return false;
  }

  // check dimensions
  const int NumDims = dimDim->size();
  if(Dim!=NumDims) {
    std::cerr << "Found " << NumDims << "D != " << Dim << "D data!" << std::endl;
    return false;
  }
  if(Dim*Dim!=dm2Dim->size()) {
    std::cerr << "Wrong dim2 dimension: " << dm2Dim->size() << std::endl;
    return false;
  }

  // check number of cells
  const int NumCells = nPDim->size();
  if(NumCells<1)  {
    std::cerr << "No cells there to be loaded!" << std::endl;
    return false;
  }
  if(Dim*NumCells!=dofDim->size()) {
    std::cerr << "Wrong number of degrees of freedom for " << NumCells << " cells: " << dofDim->size() << std::endl;
    return false;
  }

  
  // get variables
  NcVar *posVar, *boxMatrixVar, *boxStringVar;
  if(!(posVar = file.get_var("pos"))) return false;
  if(!(boxMatrixVar = file.get_var("BoxMatrix"))) return false;
  if(!(boxStringVar = file.get_var("BoxString"))) return false;

  // check periodic box string
  char *boxStringC = new char[strDim->size()];
  boxStringVar->get(boxStringC, 1, strDim->size());
  std::istringstream boxStringS(boxStringC);
  std::string boxName;
  std::getline(boxStringS, boxName, ':');
  delete[] boxStringC;
  if(0!=boxName.compare(PeriodicBox::NetCdfName)) {
    std::cerr << "Invalid box: " << boxName << "! Should be: " << PeriodicBox::NetCdfName << "!" << std::endl;
    return false;
  }

  // set box dimensions
  double transformationColumnMajor[Dim*Dim];
  boxMatrixVar->get(transformationColumnMajor, 1, Dim*Dim);
  box().x = transformationColumnMajor[0];
  box().y = transformationColumnMajor[4];
  box().z = transformationColumnMajor[8];
  box().shearYx = transformationColumnMajor[3]/box().y;
  
  // load cell positions and create cells
  double *positions = new double[Dim*NumCells];
  posVar->get(positions, 1, Dim*NumCells);
  for(int i=0; i<NumCells; ++i) {
    double *pos = positions + Dim*i;
    _cells.push_back(new Cell(*this, parameters, box().fromNormalizedPosition(Vector3D(pos[0], pos[1], pos[2]))));
  }
  delete[] positions;
          
  // close file
  if(!file.close()) return false;
  std::cout << "Done." << std::endl;
  return true;
}

bool Tessellation::saveAsNetCdfFile(const std::string &Path) const {
  const int Dim = 3;
  const int NP = _cells.size();
  const int STRING_SIZE = 128;
//  const char STRING_DELIMITER = ':';
  const char STRING_TERMINATOR = '\000';

  std::cout << "Saving as NetCDF file " << Path << "..." << std::endl;
  
	NcError err(NcError::verbose_fatal);
  
	// create file:
  NcFile file(Path.c_str(), fileExists(Path)?NcFile::Replace:NcFile::New);
	if(!file.is_valid()) {
    std::cerr << "Problem creating file!" << std::endl;
    return false;
  }

	// set the dimensions:
	NcDim *recDim, *dimDim, *dm2Dim, *NPDim, *dofDim, *strDim;
	if(!(recDim = file.add_dim("rec"))) return false;
	if(!(dimDim = file.add_dim("dim", Dim))) return false;
	if(!(dm2Dim = file.add_dim("dim2", Dim*Dim))) return false;
	if(!(NPDim  = file.add_dim("NP",  NP))) return false;
	if(!(dofDim = file.add_dim("dof", NP*Dim))) return false;
	if(!(strDim = file.add_dim("StringSize", STRING_SIZE))) return false;

	// set the variables:
	NcVar *posVar, *boxMatrixVar, *boxStringVar;
	if(!(posVar = file.add_var("pos", ncDouble, recDim, dofDim))) return false;
	if(!(boxMatrixVar = file.add_var("BoxMatrix", ncDouble, recDim, dm2Dim))) return false;
	if(!(boxStringVar = file.add_var("BoxString", ncChar, recDim, strDim))) return false;
  
  // write positions:
  double *normalizedPositions = new double[Dim*NP];
  for(int i=0; i<NP; ++i) {
    Cell *c = _cells[i];
    c->moveIntoBox();
    const Vector3D p = box().toNormalizedPosition(_cells[i]->position());
    double *np = normalizedPositions + Dim*i;
    np[0] = p.x();
    np[1] = p.y();
    np[2] = p.z();
  }
  if(!posVar->put_rec(normalizedPositions, 0)) return false;
  delete[] normalizedPositions;

  // write periodic box transformation:
  double transformationColumnMajor[9] = { box().x, 0.0, 0.0, box().shearYx*box().y, box().y, 0.0, 0.0, 0.0, box().z };
	if(!boxMatrixVar->put_rec(transformationColumnMajor, 0)) return false;
  
  // write periodic box identifier:
	std::string boxString(PeriodicBox::NetCdfName);
  boxString += ":1";
  boxString.append(STRING_SIZE-boxString.size(), STRING_TERMINATOR);
	if(!boxStringVar->put_rec(boxString.c_str(),	0)) return false;
  
  // close file:
  if(!file.sync()) return false;
  if(!file.close()) return false;
  std::cout << "Done." << std::endl;
  return true;
}

#endif /* USE_NETCDF */
