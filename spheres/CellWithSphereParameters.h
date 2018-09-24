#ifndef CELLWITHSPHEREPARAMETERS_H
#define CELLWITHSPHEREPARAMETERS_H

#include "core/CellType.h"

struct CellWithSphereParameters : public CellType {
public:
  CellWithSphereParameters() : CellType(), radius(1.0) {}
  double radius;
};

#endif /* CELLWITHSPHEREPARAMETERS_H */

