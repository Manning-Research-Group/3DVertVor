#ifndef TESSELLATIONWITHSPHERES_H
#define TESSELLATIONWITHSPHERES_H

#include <set>
#include <vector>
#if USE_NETCDF
#include <netcdfcpp.h>
#endif


#include "misc/other/DerivedContainer.h"
#include "core/Tessellation.h"
#include "CellWithSphere.h"
#include "DirectedSphereOverlap.h"

template<typename Potential> class TessellationWithSpheres : public Tessellation {
public:
  TessellationWithSpheres(PeriodicBox &box, Potential &p, double alpha);
  virtual ~TessellationWithSpheres();
  
  ConstDerivedContainer<const std::vector<Cell*>,CellWithSphere*> cellsWithSpheres() const { return ConstDerivedContainer<const std::vector<Cell*>,CellWithSphere*>(_cells); }
  DerivedContainer<std::vector<Cell*>,CellWithSphere*> cellsWithSpheres() { return DerivedContainer<std::vector<Cell*>,CellWithSphere*>(_cells); }
  const CellWithSphere *cellWithSphere(int i) const { return (const CellWithSphere*)cells()[i]; }
  CellWithSphere *cellWithSphere(int i) { return (CellWithSphere*)cells()[i]; }

  // creation
  void addCellsAtRandomPositions(const CellWithSphereParameters &cellParameters, const int NumberOfCells);
  int bruteForceCheckForMissingOverlaps();
  unsigned int numberOfSphereOverlaps() const { return _sphereOverlaps.size(); }
  
  // geometry
  virtual void computeGeometry() override;

  // energy
  virtual double energy() const override;
  
  // forces
  virtual void computeEnergyDerivatives() override;
  virtual double stressXx() const override;
  virtual double stressYy() const override;
  virtual double stressZz() const override;
  virtual double stressYx() const override; // first: normal vector, second: force vector

  // drawing
  virtual void saveAsImageDefault(const std::string &filename, int width=800, int height=800, const Vector3D &CameraPositionOffset=Vector3D(1.6, -0.8, 1.37), const Vector3D &CameraLooksAt=Vector3D(0.5, 0.5, 0.375)) const override;
  void drawCellSpheresDefault(PovrayGenerator &g, const Color &color=Color(0.0, 0.0, 0.8)) const;

  // loading and saving; compatible with the LiuJamming code
  bool loadFromDump(const CellWithSphereParameters &cellParameters, const std::string &filename);
#if USE_NETCDF
  /** this uses the format of the CStaticDatabase class of the LiuJamming code */
  bool loadFromNetCdfFileWithCellParameters(const std::string &Path, CellWithSphereParameters &parameters, int record=0);
  /** this uses the format of the CStaticDatabase class of the LiuJamming code */
  bool saveAsNetCdfFile(const std::string &Path) const;
#endif
  
private:
  Potential &_potential;
  double _alpha; // alpha==0 -> spheres part only;  alpha==1 -> tissue part only
  std::vector<DirectedSphereOverlap*> _sphereOverlaps;

  class CellPointerWithPeriodicity {
  public:
    CellPointerWithPeriodicity(const CellWithSphere *c_, const PeriodicityVector3D &pv_) : c(c_), pv(pv_) {}
    friend bool operator<(const CellPointerWithPeriodicity& l, const CellPointerWithPeriodicity &r) { return (l.c<r.c) || ((l.c==r.c) && (l.pv<r.pv)); }
  private:
    const CellWithSphere *c;
    const PeriodicityVector3D pv;
  };
  
  void removeSphereOverlaps();
  void recursivelyCreateSphereOverlaps(CellWithSphere *c1, CellWithSphere *c2, std::set<CellPointerWithPeriodicity> &visitedCells, const Vector3D &totalDistance, const PeriodicityVector3D &totalPeriodicity);
  
  // dynamical matrix
  virtual void computeDynamicalMatrix() override;
};


#endif /* TESSELLATIONWITHSPHERES_H */

