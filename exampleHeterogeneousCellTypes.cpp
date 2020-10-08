#include "core/PeriodicBox.h"
#include "core/Tessellation.h"
#include "core/CellType.h"

#include "core/Cell-inline.h"

int main(int argc, char** argv) {
  Random::seed(0x508A6472);  // sets the random seed
  
  const int NumberOfCellsType1 = 128;
  const int NumberOfCellsType2 = 128;
  const int SimulationsSteps = 2000;
  const int MovieFrameEachNSteps = 10;
  const double DeltaT = 0.01;  
  
  // create cell types
  CellType type1;
  type1.volumeElasticity = 1;  //  K_V
  type1.preferredVolume = 1;  //  V_0
  type1.surfaceElasticity = 1;  //  K_S
  type1.preferredSurface = 5.5;  //  S_0 
  type1.angularDiffusion = 1;  // D_r
  type1.speed = 0.1;  // v_0
  CellType type2;
  type2.volumeElasticity = 1;  //  K_V
  type2.preferredVolume = 1;  //  V_0
  type2.surfaceElasticity = 1;  //  K_S
  type2.preferredSurface = 5.5;  //  S_0 
  type2.angularDiffusion = 1;  // D_r
  type2.speed = 0.1;  // v_0
  
  // set additional interfacial tension between type1 and type2 cells
  type1.setAdditionalInterfacialTensionWith(type2, 2);
  // ^^^^^ For most purposes, this does the same as:
  //  type1.setAdditionalInterfacialTensionWith(type2, 1);
  //  type2.setAdditionalInterfacialTensionWith(type1, 1);

  // create box, tesselation, and add cells at random positions:
  const double Side = cubicRoot(NumberOfCellsType1 + NumberOfCellsType2);
  PeriodicBox box(Side, Side, Side);
  Tessellation t(box);
  t.addCellsAtRandomPositions(type1, NumberOfCellsType1);
  t.addCellsAtRandomPositions(type2, NumberOfCellsType2);

  // alternatively, this creates for instance a spherical arrangement of type-1 cells within type-2 cells as an initial state:
//  const Vector3D MidPoint(0.5*Side, 0.5*Side, 0.5*Side);
//  const double Radius = 0.3*Side;
//  for(int i=0; i<NumberOfCellsType1; ++i) {
//    Vector3D position;
//    do {
//      position.set(Random::uniform(box.x), Random::uniform(box.y), Random::uniform(box.z));
//    } while((position-MidPoint).norm()>Radius);
//    t.addCell(type1, position);
//  }
//  for(int i=0; i<NumberOfCellsType2; ++i) {
//    Vector3D position;
//    do {
//      position.set(Random::uniform(box.x), Random::uniform(box.y), Random::uniform(box.z));
//    } while((position-MidPoint).norm()<Radius);
//    t.addCell(type2, position);
//  }
  
  for(int i=0; i<=SimulationsSteps; ++i) {
    // updates of the Voronoi tesselation
    t.computeTopologyAndGeometry();
    t.computeEnergyDerivatives();
    
    if(i%MovieFrameEachNSteps==0) {
      std::cout << "t = " << t.time() << " / " << DeltaT*SimulationsSteps << std::endl;
 
      // this is how one can loop over all cells:
//      for(const Cell *c : t.cells()) {
        // and this is how one can access the cell position,
        // more precisely the position that is not put back into the periodic box, 
        // which is ideal to compute the MSD
//        c->positionWithoutBox()
        // the x component of this position is for instance:
//      c->positionWithoutBox().x()
//      }

      // this is to create the sum of all cell-cell interfaces in the system (totalSurface)
      // and the sum of all interfaces between cells of type 1 and type 2 (surface12)
      double totalSurface=0.0, surface12=0.0;
      for(const Cell *c : t.cells()) {
        totalSurface += c->surface();
        if(c->type()==&type1) {
          for(const DirectedFace *f : c->faces()) {
            if(f->otherCell()->type()==&type2) {
              surface12 += f->area().norm();
            }
          }
        }
      }
      totalSurface *= 0.5;  // counted twice:  sum of all cell surfaces = 2* total cell-cell interface
      std::cout << "12 interfaces: " << surface12 << " / " << totalSurface << " (" << surface12/totalSurface << ")" << std::endl;
     
    }
    // compute time step
    t.timeStep(DeltaT);
  }

  return 0;
}

