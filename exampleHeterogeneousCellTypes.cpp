#include "core/PeriodicBox.h"
#include "core/Tessellation.h"
#include "core/CellType.h"
#include "misc/other/MovieCreator.h"

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
  
  MovieCreator mc;
  mc.start();
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
     
      
      // this is to create the movie
      t.saveAsImage(mc.nextFrameFileName(), [&t, &type1, Side](PovrayGenerator &g) {
        g.addLight(Vector3D(10*Side,0,0), Color(0.75,0.70,0.70)); 
        g.addLight(Vector3D(0,-10*Side,0), Color(0.30,0.30,0.30)); 
        g.addLight(Vector3D(0,0,10*Side), Color(0.95,0.90,0.90)); 

        const Vector3D BoxPosition1(-2*Side, -2*Side, -2*Side), BoxPosition2(2*Side, 2*Side, 0.5*Side);
        PovrayMesh<int> meshForCutting;
        meshForCutting.addPoint(0, Vector3D(0, 0, 0));
        meshForCutting.addPoint(1, Vector3D(Side, 0, 0));
        meshForCutting.addPoint(2, Vector3D(0, Side, 0));
        meshForCutting.addPoint(3, Vector3D(0, 0, Side));
        meshForCutting.addTriangle(0, 1, 2);
        meshForCutting.addTriangle(0, 1, 3);
        meshForCutting.addTriangle(0, 2, 3);
        meshForCutting.addTriangle(1, 2, 3);
        
        for(const Cell *c : t.cells()) {
//          if(c->type()==&type1) {
//            g.startUnion();
//            g.drawMesh(c->createMesh());
//            g.endUnionWithColor((c->type()==&type1)?Color(0,0,1):Color(1,0,0));

//            g.startUnionWithBoxIntersection(BoxPosition1, BoxPosition2);
            g.startUnionWithMeshIntersection(meshForCutting);
            g.drawMesh(c->createMesh());
            g.endUnionWithColor((c->type()==&type1)?Color(0,0,1):Color(1,0,0));
          //          }
        }
        
//        g.startUnion();
////        g.startUnionWithBoxIntersection(BoxPosition1, BoxPosition2);
//        for(const Cell *c : t.cells()) {
////          if(c->type()==&type1) {
//            for(const DirectedFace *f : c->faces()) {
//              DirectedEdgeOfCell *edge=f->firstEdge();
//              do {
//                VertexOfCell *v1 = edge->vertex();
//                VertexOfCell *v2 = edge->nextAroundFace()->vertex();
//                g.drawCylinder(v1->positionWithRespectTo(c->position()), v2->positionWithRespectTo(c->position()), 0.05);
//                edge = edge->nextAroundFace();
//              } while(edge!=f->firstEdge());
//            }
////          }
//        }
//        g.endUnionWithColor(Color(0.5, 0.5, 0.5));
////        t.drawVerticesDefault(g, Color(0.5, 0.5, 0.5));
//        g.startUnion();
////        g.startUnionWithBoxIntersection(BoxPosition1, BoxPosition2);
//        for(const Cell *c : t.cells()) {
//          for(const VertexOfCell *v : c->vertices()) g.drawSphere(v->positionWithRespectTo(c->position()), 0.04);
//        }
//        g.endUnionWithColor(Color(0.5, 0.5, 0.5));
//  
//        for(const Cell *c : t.cells()) {
//          g.startUnion();
////          g.startUnionWithBoxIntersection(BoxPosition1, BoxPosition2);
//          g.drawSphere(c->position(), 0.4);
//          g.endUnionWithColor((c->type()==&type1)?Color(0,0,1):Color(1,0,0));
//        }
      }, 800, 800, Vector3D(1.6, -0.8, 1.37));
    }
    // compute time step
    t.timeStep(DeltaT);
  }
  mc.createMovie("movie-cut.mp4", 20);

  return 0;
}

