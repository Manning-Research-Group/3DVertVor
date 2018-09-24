#ifndef MATRIX3X3X3_H
#define	MATRIX3X3X3_H

#include "Vector3D.h"
#include "Matrix3x3.h"

class Matrix3x3x3 {
public:
  // creation
  Matrix3x3x3() : _x(Matrix3x3::Zero), _y(Matrix3x3::Zero), _z(Matrix3x3::Zero) {}
  Matrix3x3x3(const Matrix3x3 &x, const Matrix3x3 &y, const Matrix3x3 &z) : _x(x), _y(y), _z(z) {}
  Matrix3x3x3(const Matrix3x3x3 &other) : _x(other._x), _y(other._y), _z(other._z) {}

  // addition
//  Matrix3x3x3 operator-() { Matrix3x3x3 tmp(*this); tmp*=-1; return tmp; }
  void operator+=(const Matrix3x3x3 &o) { _x += o._x; _y += o._y; _z += o._z; }
  void operator-=(const Matrix3x3x3 &o) { _x -= o._x; _y -= o._y; _z -= o._z; }
  friend Matrix3x3x3 operator-(const Matrix3x3x3 &r) { return Matrix3x3x3(-r._x, -r._y, -r._z); }
  friend Matrix3x3x3 operator+(const Matrix3x3x3 &l, const Matrix3x3x3 &r) { return Matrix3x3x3(l._x+r._x, l._y+r._y, l._z+r._z); }
  friend Matrix3x3x3 operator-(const Matrix3x3x3 &l, const Matrix3x3x3 &r) { return Matrix3x3x3(l._x-r._x, l._y-r._y, l._z-r._z); }

  // multiplication with scalar
  void operator*=(const double d) { _x *= d; _y *= d; _z *= d; }
  void operator/=(const double d) { _x /= d; _y /= d; _z /= d; }
  friend Matrix3x3x3 operator*(const Matrix3x3x3 &l, const double d) { return Matrix3x3x3(d*l._x, d*l._y, d*l._z); }
  friend Matrix3x3x3 operator*(const double d, const Matrix3x3x3 &r) { return Matrix3x3x3(d*r._x, d*r._y, d*r._z); }
  friend Matrix3x3x3 operator/(const Matrix3x3x3 &l, const double d) { return Matrix3x3x3(l._x/d, l._y/d, l._z/d); }
  
  // multiplication with vector
  Matrix3x3 multiplyByVector1stIndex(const Vector3D &v) const { return Matrix3x3(v*_x, v*_y, v*_z); }
  Matrix3x3 multiplyByVector2ndIndex(const Vector3D &v) const { return Matrix3x3(_x*v, _y*v, _z*v); }
  Matrix3x3 multiplyByVector3rdIndex(const Vector3D &v) const { Matrix3x3 m(_x*v.x()); m+=_y*v.y(); m+=_z*v.z(); return m; }
  friend Matrix3x3 operator*(const Vector3D &l, const Matrix3x3x3 &r) { return r.multiplyByVector1stIndex(l); }
  friend Matrix3x3 operator*(const Matrix3x3x3 &l, const Vector3D &r) { return l.multiplyByVector3rdIndex(r); }

  // multiplication with matrix
  Matrix3x3x3 multiplyFromLeftByMatrix12(const Matrix3x3 &m) const { return Matrix3x3x3(m*_x, m*_y, m*_z); }
  Matrix3x3x3 multiplyIntoMiddleByMatrix22(const Matrix3x3 &m) const { 
    return Matrix3x3x3(_x.multiplyFromRightByMatrix22(m), _y.multiplyFromRightByMatrix22(m), _z.multiplyFromRightByMatrix22(m));
  }
//  Matrix3x3x3 multiplyFromRightByMatrix22(const Matrix3x3 &m) const { 
//    return Matrix3x3x3(Matrix3x3(_x._x*m._x.x()+_x._y*m._y.x()+_x._z*m._z.x(), _y._x*m._x.x()+_y._y*m._y.x()+_y._z*m._z.x(), _z._x*m._x.x()+_z._y*m._y.x()+_z._z*m._z.x()),
//                       Matrix3x3(_x._x*m._x.y()+_x._y*m._y.y()+_x._z*m._z.y(), _y._x*m._x.y()+_y._y*m._y.y()+_y._z*m._z.y(), _z._x*m._x.y()+_z._y*m._y.y()+_z._z*m._z.y()),
//                       Matrix3x3(_x._x*m._x.z()+_x._y*m._y.z()+_x._z*m._z.z(), _y._x*m._x.z()+_y._y*m._y.z()+_y._z*m._z.z(), _z._x*m._x.z()+_z._y*m._y.z()+_z._z*m._z.z()));
//  }
//  Matrix3x3x3 multiplyFromRightByMatrix32(const Matrix3x3 &m) const { 
//    return Matrix3x3x3(_x*m._x.x()+_y*m._y.x()+_z*m._z.x(), _x*m._x.y()+_y*m._y.y()+_z*m._z.y(), _x*m._x.z()+_y*m._y.z()+_z*m._z.z());
//  }
  Matrix3x3x3 multiplyByMatrix31(const Matrix3x3 &m) const { return Matrix3x3x3(*this*m._x, *this*m._y, *this*m._z); }
  friend Matrix3x3x3 operator*(const Matrix3x3 &l, const Matrix3x3x3 &r) { return r.multiplyFromLeftByMatrix12(l); }
  friend Matrix3x3x3 operator*(const Matrix3x3x3 &l, const Matrix3x3 &r) { return l.multiplyByMatrix31(r); }

  friend std::ostream& operator<<(std::ostream& os, const Matrix3x3x3 &r) { return os << r._x << std::endl << r._y << std::endl << r._z; }

  const static Matrix3x3x3 Zero, Epsilon;
private:
  Matrix3x3 _x, _y, _z;  // each 3x3 matrix covers the first two indices; 
};

inline Matrix3x3x3 dyadicProduct(const Vector3D &l, const Matrix3x3 &r) { 
  return Matrix3x3x3(dyadicProduct(l,r._x), dyadicProduct(l,r._y), dyadicProduct(l,r._z));
}

inline Matrix3x3x3 dyadicProductInBetween(const Vector3D &i, const Matrix3x3 &o) { 
  return Matrix3x3x3(dyadicProduct(o._x,i), dyadicProduct(o._y,i), dyadicProduct(o._z,i));
}

inline Matrix3x3x3 dyadicProduct(const Matrix3x3 &l, const Vector3D &r) { 
  return Matrix3x3x3(l*r.x(), l*r.y(), l*r.z());
}

inline Matrix3x3x3 dyadicProduct(const Vector3D &l, const Vector3D &m, const Vector3D &r) { 
  return dyadicProduct(dyadicProduct(l,m),r);
}

#endif	/* MATRIX3X3X3_H */

