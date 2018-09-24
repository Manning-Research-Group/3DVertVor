#ifndef CELLPOSITIONMAILBOX_H
#define CELLPOSITIONMAILBOX_H

#include <vector>
#include "misc/other/misc-math.h"
#include "PeriodicBox.h"
#include "Cell.h"

struct BoxIndex {
  BoxIndex() {}
  BoxIndex(int i_, int j_, int k_) : i(i_), j(j_), k(k_) {}
  int i, j, k;
};

class CellPositionMailbox {
public:
  CellPositionMailbox() : _nx(-1), _ny(-1), _nz(-1) {}
  void refill(const PeriodicBox &box, const std::vector<Cell*> &cellContainer);
  const std::vector<Cell*> &box(const BoxIndex &bi) const { return _boxes[bi.i+_nx*(bi.j+_ny*bi.k)]; }
  BoxIndex boxIndexForPosition(const Vector3D &p) { 
    return BoxIndex(modulo(p.x()/_bx, _nx), modulo(p.y()/_by, _ny), modulo(p.z()/_bz, _nz)); 
  }
  BoxIndex maxIndicesForRadius(double radius) { return BoxIndex(ceil(radius/_bx), ceil(radius/_by), ceil(radius/_bz)); }
  int nx() const { return _nx; }
  int ny() const { return _ny; }
  int nz() const { return _nz; }
  
private:  
  std::vector<Cell*> &boxRef(const BoxIndex &bi) { return _boxes[bi.i+_nx*(bi.j+_ny*bi.k)]; }
  
  int _nx, _ny, _nz;
  std::vector< std::vector<Cell*> > _boxes;
  double _bx, _by, _bz;
};

#endif /* CELLPOSITIONMAILBOX_H */

