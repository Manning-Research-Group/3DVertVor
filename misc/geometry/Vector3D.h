#ifndef VECTOR_3D_H
#define VECTOR_3D_H

#include <ostream>
#include <sstream>
#include <complex>
#include <cmath>
#include <math.h>
#include <stdio.h>

class Matrix3x3;

/**
 * A vector in three dimensions with implementations for the basic algebraic operations.
 */
class Vector3D {
public:
  Vector3D() : _x(0.0), _y(0.0), _z(0.0) {}
  Vector3D(const double x, const double y, const double z) : _x(x), _y(y), _z(z) {}
  Vector3D(const Vector3D &other) : _x(other._x), _y(other._y), _z(other._z) {}
  
  // read
  const double &x() const { return _x; }
  const double &y() const { return _y; }
  const double &z() const { return _z; }
  double &x() { return _x; }
  double &y() { return _y; }
  double &z() { return _z; }
  double normSq() const { return _x*_x + _y*_y + _z*_z; }
  double norm() const { return sqrt(normSq()); }
  double theta() const { return atan2(sqrt(_x*_x+_y*_y), _z); }
  double phi() const { return atan2(_y, _x); }
  bool isnan() const { return std::isnan(_x) || std::isnan(_y) || std::isnan(_z); }
  
//  std::complex<double> sphericalHarmonic(int l, int m) const {
//    double normBuffered = norm();
//    double P = std::tr1::assoc_legendre(l, m, _z/normBuffered);
//    return 0;
//  }

  // manipulate
  void set(const double x, const double y, const double z) { _x = x; _y = y; _z = z; }
  Vector3D operator-() { Vector3D tmp(*this); tmp*=-1; return tmp; }
  void operator+=(const Vector3D &o) { _x += o._x; _y += o._y; _z += o._z; }
  void operator-=(const Vector3D &o) { _x -= o._x; _y -= o._y; _z -= o._z; }
  void operator*=(const double d) { _x *= d; _y *= d; _z *= d; }
  void operator*=(const Matrix3x3 &o);
  void operator/=(const double d) { _x /= d; _y /= d; _z /= d; }

  friend Vector3D operator+(const Vector3D &l, const Vector3D &r) { return Vector3D(l._x+r._x, l._y+r._y, l._z+r._z); }
  friend Vector3D operator-(const Vector3D &l, const Vector3D &r) { return Vector3D(l._x-r._x, l._y-r._y, l._z-r._z); }
  friend Vector3D operator-(const Vector3D &r) { return Vector3D(-r._x, -r._y, -r._z); }
  friend Vector3D operator*(const Vector3D &l, const double d) { return Vector3D(d*l._x, d*l._y, d*l._z); }
  friend Vector3D operator/(const Vector3D &l, const double d) { return Vector3D(l._x/d, l._y/d, l._z/d); }
  friend Vector3D operator*(const double d, const Vector3D &r) { return Vector3D(d*r._x, d*r._y, d*r._z); }
  friend double operator*(const Vector3D &l, const Vector3D &r) { return l._x*r._x + l._y*r._y + l._z*r._z; }
  friend Vector3D crossProduct(const Vector3D &l, const Vector3D &r) { return Vector3D(l._y*r._z-l._z*r._y,  l._z*r._x-l._x*r._z, l._x*r._y-l._y*r._x); }
  friend std::ostream& operator<<(std::ostream& os, const Vector3D &r) { return os << r._x << ", " << r._y << ", " << r._z; }

  static Vector3D fromNormThetaPhi(double norm, double theta, double phi) { 
    double rho = norm*sin(theta);
    return Vector3D(rho*cos(phi), rho*sin(phi), norm*cos(theta));
  }
private:
  double _x, _y, _z;
};


 
#endif
