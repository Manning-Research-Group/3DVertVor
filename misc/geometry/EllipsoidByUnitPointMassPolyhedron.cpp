#include <iostream>
#include <Eigen/Dense>

#include "EllipsoidByUnitPointMassPolyhedron.h"

/* 
 * Fit ellipsoid based on moment of inertia tensor of a polyhedron. For details see
 * Dobrovolskis (1996) - Inertia of any polyhedron
*/
EllipsoidByUnitPointMassPolyhedron::EllipsoidByUnitPointMassPolyhedron(const Matrix3x3 &j) {
  Eigen::Matrix3d matrix;
  matrix << j.xx(), j.xy(), j.xz(), j.yx(), j.yy(), j.yz(), j.zx(), j.zy(), j.zz();
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigensolver(matrix);
  auto &ev = eigensolver.eigenvalues();
  auto &evec = eigensolver.eigenvectors();

  // from Dobrovolskis, A<=B<=C and Eigen3 stores eigenvalues in increasing order
  _ev0=ev(0);  // eigenvalue A
  _ev1=ev(1);  // eigenvalue B
  _ev2=ev(2);  // eigenvalue C

  // principal radii of fitted ellipsoid
  _a = sqrt(2.5*(_ev1 + _ev2 - _ev0));   // major
  _b = sqrt(2.5*(_ev0 + _ev2 - _ev1));   // intermediate
  _c = sqrt(2.5*(_ev0 + _ev1 - _ev2));   // minor
  
  // pricipal axes of fitted ellipsoid
  _aAxis = Vector3D(evec(0), evec(1), evec(2));
  _bAxis = Vector3D(evec(3), evec(4), evec(5));
  _cAxis = Vector3D(evec(6), evec(7), evec(8));

}
