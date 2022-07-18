#ifndef MATRIX3X3_H
#define	MATRIX3X3_H

#include "Vector3D.h"

class Matrix3x3x3;

/**
 * A 3x3 matrix.
 */
class Matrix3x3 {
public:
  // creation
  Matrix3x3() : _x(0.0, 0.0, 0.0), _y(0.0, 0.0, 0.0), _z(0.0, 0.0, 0.0) {}
  Matrix3x3(const Vector3D &x, const Vector3D &y, const Vector3D &z) : _x(x), _y(y), _z(z) {}
  Matrix3x3(const Matrix3x3 &other) : _x(other._x), _y(other._y), _z(other._z) {}
  static Matrix3x3 createDiagonalMatrix(double xx, double yy, double zz) { return Matrix3x3(Vector3D(xx,0,0), Vector3D(0,yy,0), Vector3D(0,0,zz)); }

  // element access
  const Vector3D &columnX() const { return _x; }
  const Vector3D &columnY() const { return _y; }
  const Vector3D &columnZ() const { return _z; }
  double xx() const { return _x.x(); }
  double xy() const { return _y.x(); }
  double xz() const { return _z.x(); }
  double yx() const { return _x.y(); }
  double yy() const { return _y.y(); }
  double yz() const { return _z.y(); }
  double zx() const { return _x.z(); }
  double zy() const { return _y.z(); }
  double zz() const { return _z.z(); }
  
  double determinant() const {
    return _x.x()*(_y.y()*_z.z() - _y.z()*_z.y()) + _y.x()*(_z.y()*_x.z() - _z.z()*_x.y()) + _z.x()*(_x.y()*_y.z() - _x.z()*_y.y());
  }
  Matrix3x3 transposed() const { return Matrix3x3(Vector3D(_x.x(),_y.x(),_z.x()), Vector3D(_x.y(),_y.y(),_z.y()), Vector3D(_x.z(),_y.z(),_z.z())); }
  
  // addition
  Matrix3x3 operator-() { Matrix3x3 tmp(*this); tmp*=-1; return tmp; }
  void operator+=(const Matrix3x3 &o) { _x += o._x; _y += o._y; _z += o._z; }
  void operator-=(const Matrix3x3 &o) { _x -= o._x; _y -= o._y; _z -= o._z; }
  friend Matrix3x3 operator-(const Matrix3x3 &r) { return Matrix3x3(-r._x, -r._y, -r._z); }
  friend Matrix3x3 operator+(const Matrix3x3 &l, const Matrix3x3 &r) { return Matrix3x3(l._x+r._x, l._y+r._y, l._z+r._z); }
  friend Matrix3x3 operator-(const Matrix3x3 &l, const Matrix3x3 &r) { return Matrix3x3(l._x-r._x, l._y-r._y, l._z-r._z); }

  // multiplication with scalar
  void operator*=(const double d) { _x *= d; _y *= d; _z *= d; }
  void operator/=(const double d) { _x /= d; _y /= d; _z /= d; }
  friend Matrix3x3 operator*(const Matrix3x3 &l, const double d) { return Matrix3x3(d*l._x, d*l._y, d*l._z); }
  friend Matrix3x3 operator*(const double d, const Matrix3x3 &r) { return Matrix3x3(d*r._x, d*r._y, d*r._z); }
  friend Matrix3x3 operator/(const Matrix3x3 &l, const double d) { return Matrix3x3(l._x/d, l._y/d, l._z/d); }
  
  // multiplication with vector
  friend Vector3D operator*(const Vector3D &l, const Matrix3x3 &r) { return Vector3D(l*r._x, l*r._y, l*r._z); }
  friend Vector3D operator*(const Matrix3x3 &l, const Vector3D &r) { Vector3D v(l._x*r.x()); v+=l._y*r.y(); v+=l._z*r.z(); return v; }

  // multiplication with matrix
  void operator*=(const Matrix3x3 &o) { _x *= o; _y *= o; _z *= o; }
  Matrix3x3 multiplyFromRightByMatrix21(const Matrix3x3 &m) const { return Matrix3x3(*this*m._x, *this*m._y, *this*m._z); }
//  Matrix3x3 multiplyFromRightByMatrix22(const Matrix3x3 &m) const { Matrix3x3 r(dyadicProduct(m._x,_x)); r+=dyadicProduct(m._y,_y); r+=dyadicProduct(m._z,_z); return r; }
  Matrix3x3 multiplyFromRightByMatrix22(const Matrix3x3 &m) const { return Matrix3x3(_x*m._x.x()+_y*m._y.x()+_z*m._z.x(), _x*m._x.y()+_y*m._y.y()+_z*m._z.y(), _x*m._x.z()+_y*m._y.z()+_z*m._z.z()); }
  friend Matrix3x3 operator*(const Matrix3x3 &l, const Matrix3x3 &r) { return l.multiplyFromRightByMatrix21(r); }
  double traceOfDotProductWith(const Matrix3x3 &m) const { return _x*m._x + _y*m._y + _z*m._z; }
  
  friend std::ostream& operator<<(std::ostream& os, const Matrix3x3 &r) { return os << r._x << "; " << r._y << "; " << r._z; }

  const static Matrix3x3 Zero, Identity;
private:
  Vector3D _x, _y, _z;  // these are column vectors
  
  friend class Matrix3x3x3;
  friend Matrix3x3x3 dyadicProduct(const Vector3D &l, const Matrix3x3 &r);
  friend Matrix3x3x3 dyadicProductInBetween(const Vector3D &i, const Matrix3x3 &o);
  friend Matrix3x3x3 operator*(const Matrix3x3x3 &l, const Matrix3x3 &r);
};

inline void Vector3D::operator*=(const Matrix3x3 &o) { *this = (*this)*o; }

// this yields \epsilon_{ijk}v_k, where i is the row index and j the column index of the resulting matrix
inline Matrix3x3 antisymmetricMatrixFromVector(const Vector3D &v) { 
  return Matrix3x3(Vector3D(0, -v.z(), v.y()), Vector3D(v.z(), 0, -v.x()), Vector3D(-v.y(), v.x(), 0));
}

// this yields l_ir_j, where i is the row index and j the column index of the resulting matrix
inline Matrix3x3 dyadicProduct(const Vector3D &l, const Vector3D &r) { 
  return Matrix3x3(r.x()*l, r.y()*l, r.z()*l);
}

#endif	/* MATRIX3X3_H */

