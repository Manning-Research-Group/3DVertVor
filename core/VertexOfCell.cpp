#include "VertexOfCell.h"
#include "Cell.h"

#include "VertexOfCell-inline.h"

const double VertexOfCell::MaxRadiusSqDeviationForVoroCheck = 1e-14;
const double VertexOfCell::MaxDistanceSqDeviationForComparisonToVoro = 1e-14;

void VertexOfCell::computeGeometry() {
  if(_edges.size()!=3) {
    std::cerr << "VertexOfCell::computeGeometry:  have " << _edges.size() << " edges!" << std::endl;
    exit(1);
  }
  
  // relative vectors between cell centers, *including* the right periodicities
  const Vector3D &rab(_edges[0]->face()->cellCenterConnectionVector()); 
  const Vector3D &rac(_edges[1]->face()->cellCenterConnectionVector()); 
  const Vector3D &rad(_edges[2]->face()->cellCenterConnectionVector()); 
  
  // this is in the denominator for several computations, the order matters!
  _normalization = rab*crossProduct(rac, rad);
  if(_normalization==0) {
    std::cerr << "VertexOfCell::computeGeometry:  normalization == 0!" << std::endl;
  }
  _positionRelativeToCellPositionDenominator = 0.5/_normalization;
  
  // position
  _positionRelativeToCellPositionNumerator =  rab.normSq()*crossProduct(rac, rad) 
                                            + rac.normSq()*crossProduct(rad, rab) 
                                            + rad.normSq()*crossProduct(rab, rac);
  _positionRelativeToCellPosition = _positionRelativeToCellPositionNumerator*_positionRelativeToCellPositionDenominator;
  if(_positionRelativeToCellPosition.isnan()) {
    std::cerr << "VertexOfCell::computeGeometry:  vertex position is nan:" << _positionRelativeToCellPosition << std::endl;
    exit(1);
  }
}

void VertexOfCell::computeVertexPositionDerivatives() {
  const Vector3D &rab(_edges[0]->face()->cellCenterConnectionVector()); 
  const Vector3D &rac(_edges[1]->face()->cellCenterConnectionVector()); 
  const Vector3D &rad(_edges[2]->face()->cellCenterConnectionVector()); 

  // derivatives
  double denominatorOfDenominatorDerivative(_positionRelativeToCellPositionDenominator/_normalization);

  _edges[0]->_derVertexPositionNumeratorWrtCellCenterConnection  = dyadicProduct(2*rab, crossProduct(rac, rad));
  _edges[0]->_derVertexPositionNumeratorWrtCellCenterConnection += rac.normSq()*antisymmetricMatrixFromVector(rad);
  _edges[0]->_derVertexPositionNumeratorWrtCellCenterConnection -= rad.normSq()*antisymmetricMatrixFromVector(rac);
  _edges[0]->_derVertexPositionDenominatorWrtCellCenterConnection = -crossProduct(rac, rad)*denominatorOfDenominatorDerivative;
  _edges[0]->_derVertexPositionWrtCellCenterConnection  = _positionRelativeToCellPositionDenominator*_edges[0]->_derVertexPositionNumeratorWrtCellCenterConnection;
  _edges[0]->_derVertexPositionWrtCellCenterConnection += dyadicProduct(_edges[0]->_derVertexPositionDenominatorWrtCellCenterConnection, _positionRelativeToCellPositionNumerator);

  _edges[1]->_derVertexPositionNumeratorWrtCellCenterConnection  = dyadicProduct(2*rac, crossProduct(rad, rab));
  _edges[1]->_derVertexPositionNumeratorWrtCellCenterConnection += rad.normSq()*antisymmetricMatrixFromVector(rab);
  _edges[1]->_derVertexPositionNumeratorWrtCellCenterConnection -= rab.normSq()*antisymmetricMatrixFromVector(rad);
  _edges[1]->_derVertexPositionDenominatorWrtCellCenterConnection = -crossProduct(rad, rab)*denominatorOfDenominatorDerivative;
  _edges[1]->_derVertexPositionWrtCellCenterConnection  = _positionRelativeToCellPositionDenominator*_edges[1]->_derVertexPositionNumeratorWrtCellCenterConnection;
  _edges[1]->_derVertexPositionWrtCellCenterConnection += dyadicProduct(_edges[1]->_derVertexPositionDenominatorWrtCellCenterConnection, _positionRelativeToCellPositionNumerator);

  _edges[2]->_derVertexPositionNumeratorWrtCellCenterConnection  = dyadicProduct(2*rad, crossProduct(rab, rac));
  _edges[2]->_derVertexPositionNumeratorWrtCellCenterConnection += rab.normSq()*antisymmetricMatrixFromVector(rac);
  _edges[2]->_derVertexPositionNumeratorWrtCellCenterConnection -= rac.normSq()*antisymmetricMatrixFromVector(rab);
  _edges[2]->_derVertexPositionDenominatorWrtCellCenterConnection = -crossProduct(rab, rac)*denominatorOfDenominatorDerivative;
  _edges[2]->_derVertexPositionWrtCellCenterConnection  = _positionRelativeToCellPositionDenominator*_edges[2]->_derVertexPositionNumeratorWrtCellCenterConnection;
  _edges[2]->_derVertexPositionWrtCellCenterConnection += dyadicProduct(_edges[2]->_derVertexPositionDenominatorWrtCellCenterConnection, _positionRelativeToCellPositionNumerator);
}

bool VertexOfCell::checkVoroPositionAndVoroCondition(const PeriodicBox &box, const std::vector<Cell*>& cells) const {
  // compare voro position to computed position
  if((_positionRelativeToCellPosition-_voroPositionRelativeToCellCenter).normSq()>MaxDistanceSqDeviationForComparisonToVoro) {
    std::cout << "Computed vertex position is different from that of the voro++ code!" << std::endl;
    return false;
  }
  
  // check radius
  Vector3D v(voroPosition()-_cell->position());
  double radiusSq = v.normSq();
  for(DirectedEdgeOfCell *e : _edges) {
    DirectedFace *f = e->face();
    double distanceSq = (f->cellCenterConnectionVector()-v).normSq();
    if(fabs(distanceSq-radiusSq)>MaxRadiusSqDeviationForVoroCheck) {
      std::cout << "Voro Position is not at the center of the tetrahedron!" << std::endl;
//      std::cout << "radiusSq: " << radiusSq << "; distanceSq: " << distanceSq << "; cell center connection vector norm Sq: " << f->cellCenterConnectionVector().normSq() << std::endl;
//      std::cout << "cell center connection vector: " << f->cellCenterConnectionVector() << std::endl;
//      std::cout << "cell center connection vector: " << f->periodicity() << std::endl;
      return false;
    }
  }
  
  // check all cell centers
  for(unsigned int ci=0; ci<cells.size(); ++ci) {
    Cell *c=cells[ci];
    Vector3D diffVector(c->position()-voroPosition());
    // take the shortest of all possible distance vectors:
    box.computeDistancePeriodicityAndRest(diffVector);
    double distanceSq = diffVector.normSq();
    if(distanceSq<radiusSq-MaxRadiusSqDeviationForVoroCheck) {
      std::cout << "Other cell is closer than sphere radius!" << std::endl;
      return false;
    }
  }
  
  return true;
}
