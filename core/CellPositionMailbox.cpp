#include "CellPositionMailbox.h"


void CellPositionMailbox::refill(const PeriodicBox &pBox, const std::vector<Cell*> &cellContainer) {
  // *** check if we have to recreate boxes ***
  const double AverageNumberOfCellsPerBox = 1;
  const double BoxVolume = AverageNumberOfCellsPerBox*pBox.volume()/cellContainer.size();
  const double BoxSide = exp(log(BoxVolume)/3);
  const int NewNx = round(pBox.x/BoxSide);
  const int NewNy = round(pBox.y/BoxSide);
  const int NewNz = round(pBox.z/BoxSide);
  if((NewNx!=_nx) || (NewNy!=_ny) || (NewNz!=_nz)) {
    // recreate boxes
    _nx = NewNx; _ny=NewNy; _nz=NewNz;
    _boxes.assign(_nx*_ny*_nz, std::vector<Cell*>());
    for(auto &b : _boxes) {
      b.reserve(10*AverageNumberOfCellsPerBox);
    }
  } else {
    for(auto &b : _boxes) {
      b.clear();
      b.reserve(10*AverageNumberOfCellsPerBox);
    }
  }
  _bx=pBox.x/_nx, _by=pBox.y/_ny, _bz=pBox.z/_nz;
  
//  std::cout << "size: " <<  _nx << ", " << _ny << ", " << _nz << "; " << _boxes.size() << std::endl;
  // *** fill boxes *** //
  for(Cell *c : cellContainer) {
//    BoxIndex bi(boxIndexForPosition(c->position()));
//    std::cout << bi.i << ", " << bi.j << ", " << bi.k << "; " << bi.i+_nx*(bi.j+_ny*bi.k) << std::endl;
    boxRef(boxIndexForPosition(c->position())).push_back(c);
  }
}
  