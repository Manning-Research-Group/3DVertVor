/* 
 * File:   Periodicity.h
 * Author: arbeit
 *
 * Created on October 9, 2015, 7:55 PM
 */

#ifndef PERIODICITY_H
#define	PERIODICITY_H

#include <ostream>
#include <sstream>

struct PeriodicityVector3D {
  PeriodicityVector3D() : x(0), y(0), z(0) {}
  PeriodicityVector3D(int px, int py, int pz) : x(px), y(py), z(pz) {}
  int x, y, z;

  bool operator==(const PeriodicityVector3D& other) { return (x==other.x) && (y==other.y) && (z==other.z); }
  bool operator!=(const PeriodicityVector3D& other) { return !(*this==other); }
  PeriodicityVector3D operator+=(const PeriodicityVector3D& other) { x+=other.x; y+=other.y; z+=other.z; return *this; }
  friend bool operator<(const PeriodicityVector3D& l, const PeriodicityVector3D &r) { 
    return (l.x<r.x) || (
              (l.x==r.x) && (
                  (l.y<r.y) || ( (l.y==r.y) && (l.z<r.z) )
              )
           );
  }
  
  friend std::ostream& operator<<(std::ostream& os, const PeriodicityVector3D &r) { return os << "(" << r.withoutBraces() << ")"; }
  std::string withoutBraces() const { 
    std::stringstream ss;
    ss << x << ", " << y << ", " << z;
    return ss.str();
  }

};

inline PeriodicityVector3D operator+(const PeriodicityVector3D& l, const PeriodicityVector3D& r) { PeriodicityVector3D tmp(l); tmp+=r; return tmp; }
inline Vector3D operator*(const PeriodicityVector3D& l, const Vector3D& r) { return Vector3D(l.x*r.x(), l.y*r.y(), l.z*r.z()); }

#endif	/* PERIODICITY_H */

