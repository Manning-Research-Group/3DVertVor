#include "misc/minimization/LmDerivativeNewton.h"
#include "misc/minimization/CgMinimizer.h"
#include "misc/minimization/CgMinimizerSearchDirectionUpdates.h"
#include "misc/minimization/CgMinimizer-inline.h"

#include "core/PeriodicBox.h"
#include "core/Tessellation.h"
#include "core/CellType.h"

#include "core/Tessellation-inline.h"

int main(int argc, char** argv) {
  const int NumberOfCells = 128;
  
  // create cell type
  CellType type;
  type.volumeElasticity = 1;  //  K_V
  type.preferredVolume = 1;  //  V_0
  type.surfaceElasticity = 1;  //  K_S
  type.preferredSurface = 5.0;  //  S_0 
  
  // create box, tesselation, and add cells at random positions:
  const double Side = cubicRoot(NumberOfCells);
  PeriodicBox box(Side);
  Tessellation t(box);
  t.addCellsAtRandomPositions(type, NumberOfCells);
 
  
  // initialize conjugated gradient minimizer;  first the line minimization
  LmDerivativeNewton lmDn;
  lmDn.setC2(1e-12);
  lmDn.setMinMaxAbsSlope(1e-12);
  // By changing the template parameters in the line below, you can change the line minimization (currently only LmDerivativeNewton) 
  // and the kind of conjugated gradient method (alternatives:  PolakRibiere [often best], FletcherReeves, MemorylessBfgs).
  CgMinimizer<LmDerivativeNewton, PolakRibiere> minimizer(lmDn);
  minimizer.setGradientTolerancePerDof(1e-12);
  minimizer.setMaxIterationsPerDof(100);
  
  // alternatively, you can also use a simple gradient minimizer with step length adaptation (SimpleGradient) or the FIRE minimizer (FireMinimizer)
  
  // tell minimizer which degrees of freedom to use for minimization
  t.setDofs(minimizer);

  // minimize
  t.relax(minimizer);

  // now minimize also with respect to simple shear degree of freedom
  // comment out the following two lines if you don't want this
  minimizer.addDof(box.shearYx, box.derEnergyWrtShearYx);
  t.relax(minimizer);
  
  // after any change, always run this function before reading anything new from the configuration:
  t.computeTopologyAndGeometry();
  
  // if you are interested in any energy derivatives, also call this if you want to know moduli
  t.computeEnergyDerivatives();
   
  // if you want to know the dynamical matrix and / or moduli
  t.computeAndSolveDynamicalMatrix();
  
  // draw final configuration:
  t.saveAsImageDefault("minimized.png");

  // write data:
  std::cout << "Energy: " << t.energy() << std::endl;

  long sumNeighborNumbers=0;
  for(Cell *c : t.cells()) sumNeighborNumbers += c->faces().size();
  std::cout << "Average neighbor number: " << ((double)sumNeighborNumbers)/t.cells().size() << std::endl;

  double sum=0;
  for(Cell *c : t.cells()) sum += c->surface();
  std::cout << "Average cell surface: " << sum/t.cells().size() << std::endl;
  
  sum=0;
  for(Cell *c : t.cells()) sum += c->volume();
  std::cout << "Average cell volume: " << sum/t.cells().size() << std::endl;
  
  std::cout << "bond orientational order Q_6:  " << t.bondOrientationalOrder<6>() << std::endl;
  std::cout << "bond orientational order Q_8:  " << t.bondOrientationalOrder<8>() << std::endl;
  std::cout << "bond orientational order Q_10:  " << t.bondOrientationalOrder<10>() << std::endl;
  
  std::cout << "smallest eigen value of Hessian:  " << t.dynamicalMatrix()->eigenValues()(0) << std::endl;
  std::cout << "smallest nontrivial eigen value of Hessian:  " << t.dynamicalMatrix()->eigenValues()(3) << std::endl;
  
  std::cout << "residual force per dof without shear dof:  " << sqrt(t.totalCellForceNormSq()/(3*t.cells().size())) << std::endl;
  std::cout << "residual force per dof with shear dof:  " << sqrt((t.totalCellForceNormSq() + box.derEnergyWrtShearYx*box.derEnergyWrtShearYx)/(3*t.cells().size()+1)) << std::endl;
  
  std::cout << "Pure shear modulus: " << t.pureShearModulusXy(1e-14) << std::endl;  // computation according to Eq. (A48) in Merkel, Manning, New Journal of Physics, 2018;  the value passed is the cutoff for declaring an eigen value of the dynamical matrix zero
  std::cout << "Simple shear modulus: " << t.simpleShearModulusYx(1e-14) << std::endl;  // computation according to Eq. (A48) in Merkel, Manning, New Journal of Physics, 2018;  the value passed is the cutoff for declaring an eigen value of the dynamical matrix zero
  
  return 0;
}

