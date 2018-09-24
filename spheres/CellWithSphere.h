#ifndef CELLWITHSPHERE_H
#define CELLWITHSPHERE_H

#include "core/Cell.h"

#include "CellWithSphereParameters.h"

template<typename Potential> class TessellationWithSpheres;

class CellWithSphere : public Cell {
public:
  template<typename P> CellWithSphere(TessellationWithSpheres<P> &tessellation, const CellWithSphereParameters &parameters, const Vector3D &initialPosition)
        : Cell(tessellation, parameters, initialPosition) {}
  virtual ~CellWithSphere() {}
  
  void setParameters(const CellWithSphereParameters &parameters) { Cell::setType(parameters); }
  const CellWithSphereParameters *parameters() const { return (const CellWithSphereParameters *)Cell::type(); }

private:
  template<typename Potential> friend class TessellationWithSpheres;
};

#endif /* CELLWITHSPHERE_H */

