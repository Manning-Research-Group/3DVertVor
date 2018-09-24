#ifndef ELLIPSOID_H
#define ELLIPSOID_H

#include "Matrix3x3.h"

class Ellipsoid {
public:
  Ellipsoid(const Matrix3x3 &integratedDyadicProductWrtCenterOfMass);
  double ev0() const { return _ev0; }
  double ev1() const { return _ev1; }
  double ev2() const { return _ev2; }
  double a() const { return _a; }
  double b() const { return _b; }
  double c() const { return _c; }
  double volume() const { return _volume; }
//  double surface() const;
  double surfaceApproximation() const;
  
private:
  double _ev0, _ev1, _ev2; // eigen values
  double _volume;
  double _a, _b, _c; // half-axes, sorted in descending order
};

#endif /* ELLIPSOID_H */

