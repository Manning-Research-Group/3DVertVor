#ifndef ELLIPSOIDBYUNITPOINTMASSPOLYHEDRON_H
#define ELLIPSOIDBYUNITPOINTMASSPOLYHEDRON_H

#include "Matrix3x3.h"

class EllipsoidByUnitPointMassPolyhedron {
public:
  EllipsoidByUnitPointMassPolyhedron(const Matrix3x3 &unitPointMassMomentOfInertiaTensor);
  double ev0() const { return _ev0; }
  double ev1() const { return _ev1; }
  double ev2() const { return _ev2; }
  double a() const { return _a; }
  double b() const { return _b; }
  double c() const { return _c; }
  Vector3D aAxis() const { return _aAxis; }
  Vector3D bAxis() const { return _bAxis; }
  Vector3D cAxis() const { return _cAxis; }
  
private:
  double _ev0, _ev1, _ev2;         // eigen values in increasing order
  double _a, _b, _c;               // ellipsoid principal radii in descending order
  Vector3D _aAxis, _bAxis, _cAxis; // direction of principal radii _a, _b, _c
};

#endif /* ELLIPSOIDBYUNITPOINTMASSPOLYHEDRON_H */