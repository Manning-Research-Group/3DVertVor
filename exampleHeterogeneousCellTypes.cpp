#include "core/PeriodicBox.h"
#include "core/Tessellation.h"
#include "core/CellType.h"

#include "core/Cell-inline.h"
#include "core/VertexOfCell-inline.h"

#include "misc/other/MovieCreator.h"
#include "misc/geometry/Ellipsoid.h"
#include "misc/geometry/EllipsoidByUnitPointMassPolyhedron.h"

#include "voro++/cell.hh"
#include "voro++/voro++.hh"
#include <math.h>

#include <experimental/filesystem>
#include <exception>
#include <unordered_map>

#include <vtkCellArray.h>
#include <vtkCellData.h>
#include <vtkDataArray.h>
#include <vtkDataSetMapper.h>
#include <vtkDoubleArray.h>
#include <vtkIntArray.h>
#include <vtkIdList.h>
#include <vtkNamedColors.h>
#include <vtkPointData.h>
#include <vtkPoints.h>
#include <vtkPolyhedron.h>
#include <vtkProperty.h>
#include <vtkSmartPointer.h>
#include <vtkUnstructuredGrid.h>
#include <vtkVersion.h>
#include <vtkXMLUnstructuredGridWriter.h>

int main(int argc, char** argv) {
  // comment this out for random results
  // sets the random seed for reproducible results
  //Random::seed(0x508A6472);

  // ------------------------------------------
  // Paraview output stuff
  // ------------------------------------------

  // choose whether to output to Paraview
  bool outputVtkVis = true;

  // check if paraview output directory exists, if not, create it
  if( outputVtkVis ) {
    if( !std::experimental::filesystem::exists( "paraview" ) ) {
      std::experimental::filesystem::create_directory("paraview");
    }
    else {
      // check if is indeed a directory (it could be a file)
      if( !std::experimental::filesystem::is_directory("paraview")) {
        std::cerr << "\"paraview\" is a file, not a directory" << std::endl;
        std::cerr << "Remove \"paraview\" file, and re-run the program." << std::endl;
        throw std::exception();
      }
    }
  }

  // header of timeseries file
  ofstream vtkTimeseries ("timeseries.pvd");
  if( outputVtkVis ) {
    vtkTimeseries << "<?xml version=\"1.0\"?>" << endl;
    vtkTimeseries << "<VTKFile type=\"Collection\" version=\"0.1\"" << endl;
    vtkTimeseries << "         byte_order=\"LittleEndian\"" << endl;
    vtkTimeseries << "         compressor=\"vtkZLibDataCompressor\">" << endl;
    vtkTimeseries << "  <Collection>" << endl;
  }

  // ------------------------------------------
  // user input
  // ------------------------------------------

  const int TotalCells = 1728; //1728
  const int NumberOfCellsType1 = 864;
  const int NumberOfCellsType2 = 864; //864
  //const int SimulationsSteps = 5e6;
  //const int MovieFrameEachNSteps = 100000;
  //const int SimulationsSteps = 5100;
  //const int MovieFrameEachNSteps = 1;
  //const int SimulationsSteps = 500100;
  //const int MovieFrameEachNSteps = 100;

  //From Preeti Paper
  const int SimulationsSteps = 20015; //50001
  const int MovieFrameEachNSteps = 2e3;

  //const int SimulationsSteps = 500000;
  //const int MovieFrameEachNSteps = 50000;
  const double DeltaT = 0.01;  
  const int maxNeighbors = 60;  // max neighbor each cell has
  
  // create cell types
  CellType type1;
  type1.volumeElasticity = 1;  //  K_V
  type1.preferredVolume = 1.0;  //  V_0
  type1.surfaceElasticity = 1;  //  K_S
  type1.preferredSurface = 5.8;  //  S_0 
  type1.angularDiffusion = 1;  // D_rz
  type1.speed = 0.0;  // v_0
  CellType type2;
  type2.volumeElasticity = 1;  //  K_V
  type2.preferredVolume = 1.0;  //  V_0
  type2.surfaceElasticity = 1;  //  K_S
  type2.preferredSurface = 5.8;  //  S_0 
  type2.angularDiffusion = 1;  // D_r
  type2.speed = 0.0;  // v_0
  CellType type3;
  type3.volumeElasticity = 1;  //  K_V
  type3.preferredVolume = 1.0;  //  V_0
  type3.surfaceElasticity = 1;  //  K_S
  type3.preferredSurface = 5.8;  //  S_0 
  type3.angularDiffusion = 1;  // D_r
  type3.speed = 0.00;  // v_0  
  
  // set additional interfacial tension between type1 and type2 cells
  type1.setAdditionalInterfacialTensionWith(type2, 0.64);
  type3.setAdditionalInterfacialTensionWith(type2, 0.64);
  // ^^^^^ For most purposes, this does the same as:
  //  type1.setAdditionalInterfacialTensionWith(type2, 1);
  //  type2.setAdditionalInterfacialTensionWith(type1, 1);

  // ------------------------------------------
  // create tessellation
  // ------------------------------------------

  // create box, tesselation, and add cells at random positions:
  const double Side = cubicRoot(NumberOfCellsType1 + NumberOfCellsType2);
  PeriodicBox box(Side, Side, Side);
  Tessellation t(box);
  //t.addCellsAtRandomPositions(type1, NumberOfCellsType1);
  //t.addCellsAtRandomPositions(type2, NumberOfCellsType2);
  //t.squareLattice(type1, NumberOfCellsType1);
  //t.squareLattice2(type2, NumberOfCellsType2);

  //Top half type 1 bot half type 2
  t.addCellsTopHalf(type1, NumberOfCellsType1);
  t.addCellsBottomHalf(type2, NumberOfCellsType2);  

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

  // ------------------------------------------
  // loop in time
  // ------------------------------------------

  // get number of cells
  int nCells = t.cells().size();

  // loop in time to stabilize system before v_0 is added to any of the cells
   for(int i=0; i<SimulationsSteps; ++i) {
      if(i==7500){
        type1.speed = 0.00;
        type2.speed = 0.00;
      }
      if(i==20000){
        int findcell = 0;
         for(Cell *c : t.cells()) {
            if(c->type()==&type1 && findcell==0 && c->position().z()>5) {
              findcell = 1;
              c->setType(type3);
              c->setcellType(2);
            }
            else if(c->type()==&type1){
              c->setcellType(0);
            }    
            else{c->setcellType(1);}
      }}


      // updates of the Voronoi tesselation
      t.computeTopologyAndGeometry();
      t.computeEnergyDerivatives();
      
      //if( (i%MovieFrameEachNSteps==0) && outputVtkVis ) {
      if( (i>2e4 || i%MovieFrameEachNSteps==0) && outputVtkVis ) {
         std::cout << "t = " << t.time() << " / " << DeltaT*SimulationsSteps << " | " << sqrt(t.totalCellForceNormSq()/(3*t.cells().size()))  << " | " <<  sqrt(t.totalCellForceNormSq()) << std::endl;
         
         // Paraview file name
         char filename[256]; 
         sprintf(filename,"paraview/cells_t0_%06f.vtu", t.time());

         int nCells = t.cells().size();
         std::cout << "Number of cells " << nCells << std::endl;
         
         // unstructured grid and points (i.e. vertices)
         vtkSmartPointer<vtkUnstructuredGrid> uGrid = vtkSmartPointer<vtkUnstructuredGrid>::New();
         vtkSmartPointer<vtkPoints> points = vtkSmartPointer<vtkPoints>::New();
         
         // time stamp
         vtkSmartPointer<vtkDoubleArray> vtkTimeArray = vtkSmartPointer<vtkDoubleArray>::New();
         vtkTimeArray->SetNumberOfComponents(1);
         vtkTimeArray->SetNumberOfTuples(1);
         vtkTimeArray->SetName("TIME");

         // cell ID
         vtkSmartPointer<vtkIntArray> cellID = vtkSmartPointer<vtkIntArray>::New();
         cellID->SetNumberOfComponents(1);
         cellID->SetNumberOfTuples(nCells);
         cellID->SetName("cellID");
         
         // cell type
         vtkSmartPointer<vtkIntArray> cellType = vtkSmartPointer<vtkIntArray>::New();
         cellType->SetNumberOfComponents(1);
         cellType->SetNumberOfTuples(nCells);
         cellType->SetName("cellType");

         // cell position (center)
         vtkSmartPointer<vtkDoubleArray> cellPosition = vtkSmartPointer<vtkDoubleArray>::New();
         cellPosition->SetNumberOfComponents(3);
         cellPosition->SetNumberOfTuples(nCells);
         cellPosition->SetName("cellPosition");

          // cell volume
         vtkSmartPointer<vtkDoubleArray> cellVolume = vtkSmartPointer<vtkDoubleArray>::New();
         cellVolume->SetNumberOfComponents(1);
         cellVolume->SetNumberOfTuples(nCells);
         cellVolume->SetName("cellVolume");

         // cell surface area
         vtkSmartPointer<vtkDoubleArray> cellSurfArea = vtkSmartPointer<vtkDoubleArray>::New();
         cellSurfArea->SetNumberOfComponents(1);
         cellSurfArea->SetNumberOfTuples(nCells);
         cellSurfArea->SetName("cellSurfArea");

         // cell energy
         vtkSmartPointer<vtkDoubleArray> cellEnergy = vtkSmartPointer<vtkDoubleArray>::New();
         cellEnergy->SetNumberOfComponents(1);
         cellEnergy->SetNumberOfTuples(nCells);
         cellEnergy->SetName("cellEnergy");

         // cell periodicity
        // +1 indicates over positive side once; -1 indicates over negative side once
        vtkSmartPointer<vtkDoubleArray> cellPeriodicity = vtkSmartPointer<vtkDoubleArray>::New();
        cellPeriodicity->SetNumberOfComponents(3);
        cellPeriodicity->SetNumberOfTuples(nCells);
        cellPeriodicity->SetName("cellPeriodicity");

         // cell velocity
        // +1 indicates over positive side once; -1 indicates over negative side once
        vtkSmartPointer<vtkDoubleArray> cellSP = vtkSmartPointer<vtkDoubleArray>::New();
        cellSP->SetNumberOfComponents(3);
        cellSP->SetNumberOfTuples(nCells);
        cellSP->SetName("cellSP");

        // cell neighbors
        // When numberOfFacesOfCell exceeds maxNeighbors, the cellNeighbor vtkDataArray
        // becomes invalid.
        vtkSmartPointer<vtkIntArray> cellNeighbor = vtkSmartPointer<vtkIntArray>::New();
        cellNeighbor->SetNumberOfComponents(maxNeighbors);
        cellNeighbor->SetNumberOfTuples(nCells);
        cellNeighbor->SetName("cellNeighbor");

         // cell ellipsoid principal radii
         vtkSmartPointer<vtkDoubleArray> cellEllipsoidRadii = vtkSmartPointer<vtkDoubleArray>::New();
         cellEllipsoidRadii->SetNumberOfComponents(3);
         cellEllipsoidRadii->SetNumberOfTuples(nCells);
         cellEllipsoidRadii->SetName("cellEllipsoidRadii");
         
         // cell ellipsoid principal a-axis (major)
         vtkSmartPointer<vtkDoubleArray> cellEllipsoidAaxis = vtkSmartPointer<vtkDoubleArray>::New();
         cellEllipsoidAaxis->SetNumberOfComponents(3);
         cellEllipsoidAaxis->SetNumberOfTuples(nCells);
         cellEllipsoidAaxis->SetName("cellEllipsoidAaxis");
         
         // cell ellipsoid principal a-axis (major)
         vtkSmartPointer<vtkDoubleArray> cellEllipsoidBaxis = vtkSmartPointer<vtkDoubleArray>::New();
         cellEllipsoidBaxis->SetNumberOfComponents(3);
         cellEllipsoidBaxis->SetNumberOfTuples(nCells);
         cellEllipsoidBaxis->SetName("cellEllipsoidBaxis");
        
         // cell ellipsoid principal c-axis (minor)
         vtkSmartPointer<vtkDoubleArray> cellEllipsoidCaxis = vtkSmartPointer<vtkDoubleArray>::New();
         cellEllipsoidCaxis->SetNumberOfComponents(3);
         cellEllipsoidCaxis->SetNumberOfTuples(nCells);
         cellEllipsoidCaxis->SetName("cellEllipsoidCaxis");

         // cell ellipsoid principal a-axis orientation
         vtkSmartPointer<vtkDoubleArray> cellorient = vtkSmartPointer<vtkDoubleArray>::New();
         cellorient->SetNumberOfComponents(1);
         cellorient->SetNumberOfTuples(nCells);
         cellorient->SetName("cellOrientation");

         //Interfacial area between basal and basement
         vtkSmartPointer<vtkDoubleArray> InterfacialCellArea12 = vtkSmartPointer<vtkDoubleArray>::New();
         InterfacialCellArea12->SetNumberOfComponents(maxNeighbors);
         InterfacialCellArea12->SetNumberOfTuples(nCells);
         InterfacialCellArea12->SetName("InterfacialCellArea12");

         // cell registration
         vtkSmartPointer<vtkDoubleArray> cellRegister12 = vtkSmartPointer<vtkDoubleArray>::New();
         cellRegister12->SetNumberOfComponents(maxNeighbors);
         cellRegister12->SetNumberOfTuples(nCells);
         cellRegister12->SetName("cellRegister12");

        std::unordered_map<Cell*, unsigned int> cellToIndex;
        for(unsigned int i=0; i<t.cells().size(); ++i) {
          cellToIndex[t.cells()[i]] = i;
        }

         int cellCounter = 0;
         int numberTotalVerticesWithDups = 0;
         int cellTypeInt = 0;
         

        for(Cell *c : t.cells()) {
          if(c->type()==&type1) {
            cellTypeInt = 0;
          }
          else if( c->type()==&type2 ) {
            cellTypeInt = 1;
          }
          else if( c->type()==&type3 ) {
            cellTypeInt = 2;
          }

           cellCounter++;
            int numberOfFacesOfCell = c->faces().size();
            const int numberOfVerticesOfCell = c->vertices().size();

            const double tmpVol = c->volume(); 
            
//            std::cout << "Cell " << cellCounter << std::endl;
//            std::cout << "   Number of faces: " << numberOfFacesOfCell << std::endl;
//            std::cout << "   Number of vert:  " << numberOfVerticesOfCell << std::endl;
            
            vtkSmartPointer<vtkIdList> vtkFaces = vtkSmartPointer<vtkIdList>::New(); // vtk faces
            vtkFaces->InsertNextId(numberOfFacesOfCell);
            
            int faceCounter = 0;
            int neighborIndices[maxNeighbors];
          	std::fill_n(neighborIndices,maxNeighbors,-1); // -1 means no neighbor
            // loop through faces and edges to store vertices
            for(const DirectedFace *f : c->faces()) {
               int numberOfVerticesOfFace = 0;

            	Cell *neighbor=f->otherCell();
           		neighborIndices[faceCounter] = cellToIndex[neighbor];

               // loop through edges to count face vertices and store vertices' positions
               DirectedEdgeOfCell *edge=f->firstEdge();
               do {
                  VertexOfCell *v1 = edge->vertex();
                  points->InsertNextPoint(v1->position().x(),v1->position().y(),v1->position().z());
                  edge = edge->nextAroundFace();
                  ++numberOfVerticesOfFace;
                  ++numberTotalVerticesWithDups;
               } while(edge!=f->firstEdge());
                
               // loop through edges again to store vertices' indexes
               // insert vertices' indexes in reverse order 
               // 3D-Voronoi uses left-hand rule; Paraview/vtk uses right-hand rule
               vtkFaces->InsertNextId(numberOfVerticesOfFace);
               for( int jj=0 ; jj < numberOfVerticesOfFace ; ++jj ) { 
                  vtkFaces->InsertNextId(numberTotalVerticesWithDups-1-jj);
               }
               ++faceCounter;
            } // end face loop
 
            double interfacialAreaArray12[maxNeighbors];
            double RegisterArray12[maxNeighbors];

            std::fill_n(interfacialAreaArray12,maxNeighbors,-1); // -1 means no interface
            std::fill_n(RegisterArray12,maxNeighbors,-1); // -1 means no interface

            double interfacialArea;
            int count12=0;          
            int othertype;
            double cellvols;
            double tempVec1[3];
            tempVec1[0] = c->position().x();
            tempVec1[1] = c->position().y();
            tempVec1[2] = c->position().z();
            double tempVec2[3];
            double tempx;
            double tempy;
            // loop through faces of basal cells
            //std::cout <<  cellTypeInt << std::endl;
            if(cellTypeInt==1 or cellTypeInt==0){
              for(const DirectedFace *f : c->faces()) {
                    if(f->otherCell()->type()==&type1) {
                      othertype = 0;
                    }
                    else if(f->otherCell()->type()==&type2 ) {
                      othertype = 1;
                    }
                  
                  //std::cout <<  cellTypeInt << " "<< othertype << std::endl;
                  
				  interfacialArea = f->area().norm();
				  interfacialAreaArray12[count12]=interfacialArea;
                  if(othertype!=cellTypeInt) {                     
                      cellvols = c->volume();
                      //std::cout <<  interfacialArea << std::endl;
                      tempVec2[0] = f->otherCell()->position().x();
                      tempVec2[1] = f->otherCell()->position().y();
                      tempVec2[2] = f->otherCell()->position().z();                
                      if(othertype==0)
                      {
                        
                        
                        if(interfacialArea>0.5){
                        tempy=tempVec2[1]-tempVec1[1];
                        tempx=tempVec2[0]-tempVec1[0];
                        if(abs(tempx)>pow(TotalCells,1/3)-2){
                          tempx=abs(tempx)-pow(TotalCells,1/3);
                        }
                        if(abs(tempy)>pow(TotalCells,1/3)-2){
                          tempy=abs(tempy)-pow(TotalCells,1/3);
                        }
                        RegisterArray12[count12]=1-pow(pow(tempx,2)+pow(tempy,2),0.5)/(pow(cellvols,1/3));
                        }
                        
                      }
                  }
              count12 +=1;}       
            }


            // add cell to unstructure grid
            uGrid->InsertNextCell(VTK_POLYHEDRON,vtkFaces);

            // add attributes to cell
            cellID->InsertValue(cellCounter-1, cellCounter);
            cellType->InsertValue(cellCounter-1, cellTypeInt);

            double tempVec[3];
            double temporient;
            double tempSPP[3];
            tempVec[0] = c->position().x();
            tempVec[1] = c->position().y();
            tempVec[2] = c->position().z();
            cellPosition->InsertTuple(cellCounter-1, tempVec); 
            cellNeighbor->InsertTupleValue(cellCounter-1,neighborIndices);
            cellVolume->InsertValue(cellCounter-1, c->volume());
            cellSurfArea->InsertValue(cellCounter-1, c->surface());
            cellEnergy->InsertValue(cellCounter-1, c->energy());
            InterfacialCellArea12->InsertTupleValue(cellCounter-1,interfacialAreaArray12);
            cellRegister12->InsertTupleValue(cellCounter-1,RegisterArray12);
           
            tempSPP[0] = c->force().x();
            tempSPP[1] = c->force().y();
            tempSPP[2] = c->force().z();
            cellSP->InsertTuple(cellCounter-1, tempSPP);

            EllipsoidByUnitPointMassPolyhedron fittedCell = 
               c->fitEllipsoidByUnitPointMassPolyhedron();
            tempVec[0] = fittedCell.a();
            tempVec[1] = fittedCell.b();
            tempVec[2] = fittedCell.c();
            cellEllipsoidRadii->InsertTuple(cellCounter-1, tempVec);
            tempVec[0] = fittedCell.aAxis().x();
            tempVec[1] = fittedCell.aAxis().y();
            tempVec[2] = fittedCell.aAxis().z();

            temporient = acos(abs(tempVec[2]));
            cellorient->InsertValue(cellCounter-1, temporient); 
            cellEllipsoidAaxis->InsertTuple(cellCounter-1,tempVec);
            tempVec[0] = fittedCell.bAxis().x();
            tempVec[1] = fittedCell.bAxis().y();
            tempVec[2] = fittedCell.bAxis().z();
            cellEllipsoidBaxis->InsertTuple(cellCounter-1,tempVec);
            tempVec[0] = fittedCell.cAxis().x();
            tempVec[1] = fittedCell.cAxis().y();
            tempVec[2] = fittedCell.cAxis().z();
            cellEllipsoidCaxis->InsertTuple(cellCounter-1,tempVec);

              
            uGrid->GetCellData()->AddArray(cellSP); 
            uGrid->GetCellData()->AddArray(cellorient);              
            uGrid->GetCellData()->AddArray(cellID);
            uGrid->GetCellData()->AddArray(cellType);
            uGrid->GetCellData()->AddArray(cellPosition);
            uGrid->GetCellData()->AddArray(cellVolume);
            uGrid->GetCellData()->AddArray(cellSurfArea);
            uGrid->GetCellData()->AddArray(cellEnergy);
            uGrid->GetCellData()->AddArray(cellEllipsoidRadii);
            uGrid->GetCellData()->AddArray(cellEllipsoidAaxis);
            uGrid->GetCellData()->AddArray(cellEllipsoidBaxis);
            uGrid->GetCellData()->AddArray(cellEllipsoidCaxis);
            uGrid->GetCellData()->AddArray(cellNeighbor);
            uGrid->GetCellData()->AddArray(InterfacialCellArea12);
            uGrid->GetCellData()->AddArray(cellRegister12);
            
         } // end cell loop
         
         //std::cout << "Total number of vertices: " << numberTotalVerticesWithDups << std::endl;

         // add time stamp to vtk data set
         vtkTimeArray->InsertValue(0,t.time());
         uGrid->GetFieldData()->AddArray(vtkTimeArray);

         uGrid->SetPoints(points);

         vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter =
            vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
         vtkWriter->SetInputData(uGrid);
         vtkWriter->SetFileName(filename);
         vtkWriter->SetDataModeToBinary();
         vtkWriter->Update();
 
      } // end outputVtkVis if

      // compute time step
      t.timeStep(DeltaT);
   
   } // end time step loop

   return 0;
}