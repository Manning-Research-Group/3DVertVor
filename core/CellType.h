#ifndef CELLTYPE_H
#define	CELLTYPE_H

#include <unordered_map>

struct CellType {
public:
  CellType() : experiencesForce(true), 
          speed(0.0), angularDiffusion(1.0), 
          volumeElasticity(1.0), surfaceElasticity(1.0), preferredSurface(5.4), preferredVolume(1.0), 
          edgeSprings(false), edgeElasticity(1.0) {}

  bool experiencesForce;
  double speed;
  double angularDiffusion;
  double volumeElasticity;
  double surfaceElasticity;
  double preferredSurface;
  double preferredVolume;

  void setAdditionalInterfacialTensionWith(const CellType &other, double at) { _additionalInterfacialTension[&other] = at; }
  double additionalInterfacialTensionWith(const CellType *other) const { return (_additionalInterfacialTension.count(other))?_additionalInterfacialTension.at(other):0.0; }
  
  bool edgeSprings;
  double edgeElasticity;

private:
  std::unordered_map<const CellType*,double> _additionalInterfacialTension;
};

#endif	/* CELLPARAMETERS_H */

