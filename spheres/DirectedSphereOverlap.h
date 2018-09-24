#ifndef DIRECTEDSPHEREOVERLAP_H
#define DIRECTEDSPHEREOVERLAP_H

#include "core/PeriodicityVector3D.h"
#include "CellWithSphere.h"

struct DirectedSphereOverlap {
  DirectedSphereOverlap(CellWithSphere *from, CellWithSphere *to, Vector3D distanceVector, PeriodicityVector3D periodicity, double sigma) 
      : From(from), To(to), DistanceVector(distanceVector), Periodicity(periodicity), Sigma(sigma), RelativeOverlap(1.0 - DistanceVector.norm()/Sigma) {}
  CellWithSphere * const From, * const To;
  const Vector3D DistanceVector;
  const PeriodicityVector3D Periodicity;
  const double Sigma;
  const double RelativeOverlap;
};

#endif /* SPHEREOVERLAP_H */

