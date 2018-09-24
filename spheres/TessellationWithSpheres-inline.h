#ifndef TESSELLATIONWITHSPHERES_INLINE_H
#define TESSELLATIONWITHSPHERES_INLINE_H

#include <set>
#include <vector>
#include <string>
#include "misc/other/DerivedContainer.h"
#include "misc/other/fileIo.h"
#include "core/Tessellation.h"
#include "core/Tessellation-inline.h"
#include "core/Cell-inline.h"

#include "CellWithSphere.h"
#include "DirectedSphereOverlap.h"
#include "TessellationWithSpheres.h"

template<typename Potential> 
TessellationWithSpheres<Potential>::TessellationWithSpheres(PeriodicBox &box, Potential &p, double alpha) 
  : Tessellation(box), _potential(p), _alpha(alpha) {
}

template<typename Potential> 
TessellationWithSpheres<Potential>::~TessellationWithSpheres() { removeSphereOverlaps(); }

// ***** creation ***** //

template<typename Potential> 
void TessellationWithSpheres<Potential>::addCellsAtRandomPositions(const CellWithSphereParameters &cellParameters, const int NumberOfCells) {
  for(int i=0; i<NumberOfCells; ++i) {
    _cells.push_back(new CellWithSphere(*this, cellParameters, Vector3D(Random::uniform(box().x), Random::uniform(box().y), Random::uniform(box().z))));
  }
}

template<typename Potential> 
int TessellationWithSpheres<Potential>::bruteForceCheckForMissingOverlaps() {
  int missing = 0;
  int total = 0, ok = 0;
  for(CellWithSphere *c1 : cellsWithSpheres()) {
    for(CellWithSphere *c2 : cellsWithSpheres()) {
      if(c1<c2) {
        ++total;
        double distance = box().minimalDistanceVector(c2->position()-c1->position()).norm();
        if(distance<=(c1->parameters()->radius+c2->parameters()->radius)) {
          bool found=false;
          for(DirectedSphereOverlap *so : _sphereOverlaps) {
            if((c1==so->From) && (c2==so->To)) {
              found = true;
              break;
            }
          }
          if(!found) {
            ++missing;
          } else {
            ++ok;
          }
        }
      }
    }
  }
  std::cout << "ok: " << ok << std::endl;
  std::cout << "total: " << total << std::endl;
  return missing;
}

// ***** geometry ***** //

template<typename Potential> 
void TessellationWithSpheres<Potential>::computeGeometry() {
  // do tissue stuff
  Tessellation::computeGeometry();

  // remove sphere overlaps
  removeSphereOverlaps();
  // create sphere overlaps and compute stuff
  std::set<CellPointerWithPeriodicity> visitedCells;
  for(CellWithSphere *c1 : cellsWithSpheres()) {
    Vector3D totalDistance;
    PeriodicityVector3D totalPeriodicity;
    recursivelyCreateSphereOverlaps(c1, c1, visitedCells, totalDistance, totalPeriodicity);
    visitedCells.clear();
  }
}

template<typename Potential> 
void TessellationWithSpheres<Potential>::removeSphereOverlaps() {
  for(DirectedSphereOverlap *so : _sphereOverlaps) delete so;
  _sphereOverlaps.clear();
}

template<typename Potential> 
void TessellationWithSpheres<Potential>::recursivelyCreateSphereOverlaps(CellWithSphere *c1, CellWithSphere *c2, std::set<CellPointerWithPeriodicity> &visitedCells, const Vector3D &totalDistance, const PeriodicityVector3D &totalPeriodicity) {
  CellPointerWithPeriodicity index(c2, totalPeriodicity);
  if(visitedCells.count(index)==0) {
    visitedCells.insert(index);
    // check distance
    double radiusSum = c1->parameters()->radius + c2->parameters()->radius;
    if(totalDistance.norm() <= radiusSum) {
      // add overlap, but only one per cell pair (c1<c2)
      // also, (c1<c2) excludes a cell being connected to itself
      if(c1<c2) {
        _sphereOverlaps.push_back(new DirectedSphereOverlap(c1, c2, totalDistance, totalPeriodicity, radiusSum));
      }
      // loop over all neighbors of cell c2
      for(DirectedFace *f : c2->faces()) {
        recursivelyCreateSphereOverlaps(c1, (CellWithSphere *)f->otherCell(), visitedCells, totalDistance+f->cellCenterConnectionVector(), totalPeriodicity+f->periodicity());
      }
    }
  }
}


// ***** energy ***** //

template<typename Potential> 
double TessellationWithSpheres<Potential>::energy() const {
  double eTissue = Tessellation::energy();
  // add overlap energies
  double eSpheres = 0.0;
  for(DirectedSphereOverlap *so : _sphereOverlaps) {
    eSpheres += _potential.energy(so->RelativeOverlap);
  }
  return _alpha*eTissue + (1.0-_alpha)*eSpheres;
}


// ***** forces ***** //

template<typename Potential> 
void TessellationWithSpheres<Potential>::computeEnergyDerivatives() {
  Tessellation::computeEnergyDerivatives();
  
  // multiply all cell energy derivatives by alpha
  for(CellWithSphere *c : cellsWithSpheres()) {
    c->_derTotalEnergyWrtCellPosition *= _alpha;
  }
  // add overlap energy derivatives to cell forces (and something else?)
  for(DirectedSphereOverlap *so : _sphereOverlaps) {
    const Vector3D Derivative(-(1.0-_alpha)*_potential.firstDerivative(so->RelativeOverlap)/so->Sigma*so->DistanceVector/so->DistanceVector.norm());
    so->From->_derTotalEnergyWrtCellPosition -= Derivative;
    so->To->_derTotalEnergyWrtCellPosition += Derivative;
  }
  
  // update periodic box derivatives
  box().derEnergyWrtShearYx *= _alpha;
  double derSpheres=0.0;
  for(DirectedSphereOverlap *so : _sphereOverlaps) {
    derSpheres -= _potential.firstDerivative(so->RelativeOverlap)/so->Sigma  *so->DistanceVector.x()/so->DistanceVector.norm() *so->Periodicity.y*box().y;
  }
  box().derEnergyWrtShearYx += (1.0-_alpha)*derSpheres;
}

#define OVERLOAD_STRESS_COMPONENT_COMPUTATION_WITH_OVERLAPS(function,normalAxis,forceAxis) \
  template<typename Potential> \
  double TessellationWithSpheres<Potential>::function() const { \
    double tensionTissue=Tessellation::function(); \
    double tensionSpheres=0.0; \
    for(DirectedSphereOverlap *so : _sphereOverlaps) { \
      tensionSpheres -= _potential.firstDerivative(so->RelativeOverlap)/so->Sigma \
                        *so->DistanceVector.forceAxis()/so->DistanceVector.norm() \
                        *so->Periodicity.normalAxis*box().normalAxis; \
    } \
    tensionSpheres /= box().volume(); \
    return _alpha*tensionTissue + (1.0-_alpha)*tensionSpheres; \
  }
OVERLOAD_STRESS_COMPONENT_COMPUTATION_WITH_OVERLAPS(stressXx, x, x)
OVERLOAD_STRESS_COMPONENT_COMPUTATION_WITH_OVERLAPS(stressYy, y, y)
OVERLOAD_STRESS_COMPONENT_COMPUTATION_WITH_OVERLAPS(stressZz, z, z)
OVERLOAD_STRESS_COMPONENT_COMPUTATION_WITH_OVERLAPS(stressYx, y, x) // first: normal vector, second: force vector


// ***** dynamical matrix ***** //

template<typename Potential> 
void TessellationWithSpheres<Potential>::computeDynamicalMatrix() {
  Tessellation::computeDynamicalMatrix();
  *_dynamicalMatrix *= _alpha;
  
  // 2 "control parameters": 0=pure shear xy; 1=simple shear yx (normal,force)
  DynamicalMatrixOld dmSpheres(3*_cells.size(), 2);

  // compute most of the contributions
  for(DirectedSphereOverlap *so : _sphereOverlaps) {
    const double LengthSq = so->DistanceVector.normSq();
    const double Length = sqrt(LengthSq);
    Vector3D unitVector(so->DistanceVector/Length);
    
    Matrix3x3 cellCellTerm(dyadicProduct((_potential.secondDerivative(so->RelativeOverlap)/(so->Sigma*so->Sigma))*unitVector, unitVector));
    Matrix3x3 secondCellCellTerm(Matrix3x3::Identity);
    secondCellCellTerm -= dyadicProduct(unitVector, unitVector);
    secondCellCellTerm *= -_potential.firstDerivative(so->RelativeOverlap)/(Length*so->Sigma);
    cellCellTerm += secondCellCellTerm;
    addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(dmSpheres, so->From, so->From, cellCellTerm);
    addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(dmSpheres, so->To, so->To, cellCellTerm);
    Matrix3x3 minusCellCellTerm(-cellCellTerm);
    addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(dmSpheres, so->From, so->To, minusCellCellTerm);
    addTotalEnergyDerivativeToDynamicalMatrixPositionPosition(dmSpheres, so->To, so->From, minusCellCellTerm);
    
    // 0 = pure shear xy
    addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(dmSpheres, 0, so->To,
            so->Periodicity.x*box().x*cellCellTerm.columnX() - so->Periodicity.y*box().y*cellCellTerm.columnY());
    addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(dmSpheres, 0, so->From,
          - so->Periodicity.x*box().x*cellCellTerm.columnX() + so->Periodicity.y*box().y*cellCellTerm.columnY());
    dmSpheres.add2ndDerivativeControlParameter(0, 
            so->Periodicity.x *so->Periodicity.x *box().x *box().x *cellCellTerm.xx()
        -2* so->Periodicity.x *so->Periodicity.y *box().x *box().y *cellCellTerm.xy()
          + so->Periodicity.y *so->Periodicity.y *box().y *box().y *cellCellTerm.yy());
    // need to add "pressure" part to second pure shear derivative
    dmSpheres.add2ndDerivativeControlParameter(0, -_potential.firstDerivative(so->RelativeOverlap)/so->Sigma 
            *(so->Periodicity.x *box().x *unitVector.x() + so->Periodicity.y *box().y *unitVector.y()) );

    // 1 = simple shear yx (normal,force)
    addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(dmSpheres, 1, so->To, box().y*so->Periodicity.y*cellCellTerm.columnX());
    addTotalEnergyDerivativeToDynamicalMatrixControlParameterPosition(dmSpheres, 1, so->From, -box().y*so->Periodicity.y*cellCellTerm.columnX());
    dmSpheres.add2ndDerivativeControlParameter(1, box().y*box().y*so->Periodicity.y*so->Periodicity.y*cellCellTerm.xx());    
  }
  
  // add to dynamical matrix
  dmSpheres *= 1.0-_alpha;
  *_dynamicalMatrix += dmSpheres;
}


// ***** drawing ***** //

template<typename Potential> 
void TessellationWithSpheres<Potential>::saveAsImageDefault(const std::string &filename, int width, int height, const Vector3D &CameraPositionOffset, const Vector3D &CameraLooksAt) const {
  saveAsImage(filename, [this](PovrayGenerator &g){
    this->drawEdgesDefault(g);
    this->drawVerticesDefault(g);
    this->drawCellSpheresDefault(g);
  }, width, height, CameraPositionOffset, CameraLooksAt);
}

template<typename Potential> 
void TessellationWithSpheres<Potential>::drawCellSpheresDefault(PovrayGenerator &g, const Color &color) const {
  drawCells(g, color, [](PovrayGenerator &g, const Cell *c){ g.drawSphere(c->position(), ((CellWithSphere*)c)->parameters()->radius); });
}


// ***** loading and saving ***** //
template<typename Potential> 
bool TessellationWithSpheres<Potential>::loadFromDump(const CellWithSphereParameters &cellParameters, const std::string &filename) {
  return loadFromDumpHelper(filename, [this,&cellParameters](const Vector3D &initialPosition) { return new CellWithSphere(*this, cellParameters, initialPosition); });
}

#if USE_NETCDF

template<typename Potential> 
bool TessellationWithSpheres<Potential>::loadFromNetCdfFileWithCellParameters(const std::string &Path, CellWithSphereParameters &parameters, int record) {
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
	NcDim *recDim, *dimDim, *dm2Dim, *nPDim, *dofDim, *strDim, *potStrSizeDim;
	if(!(recDim = file.get_dim("rec"))) return false;
	if(!(dimDim = file.get_dim("dim"))) return false;
	if(!(dm2Dim = file.get_dim("dim2"))) return false;
	if(!(nPDim  = file.get_dim("NP"))) return false;
	if(!(dofDim = file.get_dim("dof"))) return false;
	if(!(strDim = file.get_dim("StringSize"))) return false;
	if(!(potStrSizeDim = file.get_dim("PotStrSize"))) return false;

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
  NcVar *posVar, *radVar, *boxMatrixVar, *boxStringVar, *potStringVar;
	if(!(posVar = file.get_var("pos"))) return false;
	if(!(radVar = file.get_var("rad"))) return false;
	if(!(boxMatrixVar = file.get_var("BoxMatrix"))) return false;
	if(!(boxStringVar = file.get_var("BoxString"))) return false;
  if(!(potStringVar = file.get_var("PotString"))) return false;

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

  // check periodic box string
  char *potStringC = new char[potStrSizeDim->size()];
  potStringVar->get(potStringC, 1, potStrSizeDim->size());
  std::istringstream potStringS(potStringC);
  std::string potName;
  std::getline(potStringS, potName, ':');
  delete[] potStringC;
  if(0!=potName.compare(Potential::NetCdfName)) {
    std::cout << "Warning: Different potentials! In file: " << potName << ".  I have: " << Potential::NetCdfName << "." << std::endl;
  } else {
    std::string potParameterString;
    std::getline(potStringS, potParameterString, '\000');
    _potential.setParametersFromNetCdfString(potParameterString);
  }

  // check whether all radii are the same
  const double MaxRelDeviation=1e-12;
  double *radii = new double[NumCells];
  radVar->get(radii, 1, NumCells);
  const double Radius = radii[0];
  const double MaxRadiusDeviation = Radius*MaxRelDeviation;
  for(int i=1; i<NumCells; ++i) {
    if(fabs(radii[i]-Radius)>MaxRadiusDeviation) {
      std::cerr << "Not all radii are the same, which is currently required by the loading routine!" << std::endl;
      return false;
    }
  }
  delete[] radii;
  std::cout << "Using radius = " << Radius << std::endl;

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
  parameters.radius = Radius;
  for(int i=0; i<NumCells; ++i) {
    double *pos = positions + Dim*i;
    _cells.push_back(new CellWithSphere(*this, parameters, box().fromNormalizedPosition(Vector3D(pos[0], pos[1], pos[2]))));
  }
  delete[] positions;
          
  // close file
  if(!file.close()) return false;
  std::cout << "Done." << std::endl;
  return true;
}

template<typename Potential> 
bool TessellationWithSpheres<Potential>::saveAsNetCdfFile(const std::string &Path) const {
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
	NcDim *recDim, *dimDim, *dm2Dim, *NPDim, *dofDim, *strDim, *potStrSizeDim;
	if(!(recDim = file.add_dim("rec"))) return false;
	if(!(dimDim = file.add_dim("dim", Dim))) return false;
	if(!(dm2Dim = file.add_dim("dim2", Dim*Dim))) return false;
	if(!(NPDim  = file.add_dim("NP",  NP))) return false;
	if(!(dofDim = file.add_dim("dof", NP*Dim))) return false;
	if(!(strDim = file.add_dim("StringSize", STRING_SIZE))) return false;
	if(!(potStrSizeDim = file.add_dim("PotStrSize", STRING_SIZE))) return false;

	// set the variables:
	NcVar *posVar, *radVar, *boxMatrixVar, *boxStringVar, *potStringVar;
	if(!(posVar = file.add_var("pos", ncDouble, recDim, dofDim))) return false;
	if(!(radVar = file.add_var("rad", ncDouble, recDim, NPDim))) return false;
	if(!(boxMatrixVar = file.add_var("BoxMatrix", ncDouble, recDim, dm2Dim))) return false;
	if(!(boxStringVar = file.add_var("BoxString", ncChar, recDim, strDim))) return false;
  if(!(potStringVar = file.add_var("PotString", ncChar, recDim, potStrSizeDim))) return false;
  
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

  // write radii:
  double *radii = new double[NP];
  for(int i=0; i<NP; ++i) {
    radii[i] = cellWithSphere(i)->parameters()->radius;
  }
  if(!radVar->put_rec(radii, 0)) return false;
  delete[] radii;

  // write periodic box transformation:
  double transformationColumnMajor[9] = { box().x, 0.0, 0.0, box().shearYx*box().y, box().y, 0.0, 0.0, 0.0, box().z };
	if(!boxMatrixVar->put_rec(transformationColumnMajor, 0)) return false;
  
  // write periodic box identifier:
	std::string boxString(PeriodicBox::NetCdfName);
  boxString += ":1";
  boxString.append(STRING_SIZE-boxString.size(), STRING_TERMINATOR);
	if(!boxStringVar->put_rec(boxString.c_str(),	0)) return false;
  
  // write potential identifier:
	std::string potentialString(_potential.getNetCdfNameAndParameters());
  potentialString.append(STRING_SIZE-potentialString.size(), STRING_TERMINATOR);
  if(!potStringVar->put_rec(potentialString.c_str(), 0)) return false;
  
  // close file:
  if(!file.sync()) return false;
  if(!file.close()) return false;
  std::cout << "Done." << std::endl;
  return true;
}

#endif /* USE_NETCDF */


#endif /* TESSELLATIONWITHSPHERES_INLINE_H */

