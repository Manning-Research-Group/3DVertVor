#ifndef PERIODICBOX_H
#define	PERIODICBOX_H

#include <math.h>

#include "misc/geometry/Vector3D.h"

#include "PeriodicityVector3D.h"

struct PeriodicBox {
  PeriodicBox(double s) : x(s), y(s), z(s), shearYx(0.0) {}
  PeriodicBox(double sX, double sY, double sZ) : x(sX), y(sY), z(sZ), shearYx(0.0) {}
  double volume() const { return x*y*z; }
  Vector3D periodicityOffset(const PeriodicityVector3D &p) const { return Vector3D(x*p.x + shearYx*y*p.y, y*p.y, z*p.z); }
  PeriodicityVector3D computePositionPeriodicityAndRest(Vector3D &position) const {
    int py = floor(position.y()/y);
    PeriodicityVector3D pv(floor((position.x()-shearYx*y*py)/x), py, floor(position.z()/z));
    position -= periodicityOffset(pv);
    return pv;
  }
  PeriodicityVector3D computeDistancePeriodicityAndRest(Vector3D &distance) const {
    int py = round(distance.y()/y);
    PeriodicityVector3D pv(round((distance.x()-shearYx*y*py)/x), py, round(distance.z()/z));
    distance -= periodicityOffset(pv);
    return pv;
  }
  Vector3D minimalDistanceVector(const Vector3D &distance) const {
    Vector3D md(distance);
    computeDistancePeriodicityAndRest(md);
    return md;
  }
  Vector3D toNormalizedPosition(const Vector3D &p) const { return Vector3D((p.x()-shearYx*p.y())/x, p.y()/y, p.z()/z);  }
  Vector3D fromNormalizedPosition(const Vector3D &p) const { return Vector3D(x*p.x()+y*shearYx*p.y(), y*p.y(), z*p.z());  }
  double x, y, z;
  double shearYx; 
  double derEnergyWrtShearYx; 
#if USE_NETCDF
  static const std::string NetCdfName;
#endif
};

#endif	/* PERIODICBOX_H */

