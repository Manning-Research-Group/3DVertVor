#include "core/PeriodicBox.h"
#include "core/Tessellation.h"
#include "core/CellType.h"

#include "core/Cell-inline.h"
#include "core/VertexOfCell-inline.h"

#include "voro++/cell.hh"
#include "voro++/voro++.hh"

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
  Random::seed(0x508A6472);

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

  const int NumberOfCellsType1 = 128;
  const int NumberOfCellsType2 = 128;
  const int SimulationsSteps = 200;
  const int MovieFrameEachNSteps = 10;
  const double DeltaT = 0.01;  
  const int maxNeighbors = 60;  // max neighbor each cell has
  
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

  // ------------------------------------------
  // create tessellation
  // ------------------------------------------

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

  // ------------------------------------------
  // loop in time
  // ------------------------------------------

  // get number of cells
  int nCells = t.cells().size();

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

      if( outputVtkVis ) {

        // Paraview file name
        char filename[256];
        sprintf(filename,"paraview/cells_t0_%06d.vtu", i);

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

        // create map from Cell* to integer
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

          int numberOfFacesOfCell = c->faces().size();
          // output warning message to error file in case maxNeighbors is too small
          if( numberOfFacesOfCell > maxNeighbors ) {
            std::clog << "************************************************";
            std::clog << "***************************************************" << std::endl;
            std::clog << "WARNING: number of faces exceeded maxNeighbors. ";
            std::clog << "cellNeighbor vtkDataArray in *.vtu files is invalid" << std::endl;
            std::clog << "   cellCounter  = " << cellCounter << std::endl;
            std::clog << "   maxNeighbors = " << maxNeighbors << std::endl;
            std::clog << "   numberOfFacesOfCell = " << numberOfFacesOfCell << std::endl;
            std::clog << "************************************************";
            std::clog << "***************************************************" << std::endl;
          }

          // vtk faces
          // [numberOfCellFaces, (numberOfPointsOfFace0, pointId0, pointId1, … ),
          //                     (numberOfPointsOfFace1, pointId0, pointId1, …), … ].
          vtkSmartPointer<vtkIdList> vtkFaces = vtkSmartPointer<vtkIdList>::New();
          vtkFaces->InsertNextId(numberOfFacesOfCell);

          // loop through faces and edges to store neighbors and vertices
          int faceCounter = 0;
          int neighborIndices[maxNeighbors];
          std::fill_n(neighborIndices,maxNeighbors,-1); // -1 means no neighbor
          for(const DirectedFace *f : c->faces()) {
            int numberOfVerticesOfFace = 0;

            // store neighbors
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

          // add cell to unstructure grid
          uGrid->InsertNextCell(VTK_POLYHEDRON,vtkFaces);

          // add attributes to cell
          cellID->InsertValue(cellCounter, cellCounter);
          cellType->InsertValue(cellCounter, cellTypeInt);
          cellVolume->InsertValue(cellCounter, c->volume());
          cellSurfArea->InsertValue(cellCounter, c->surface());
          cellEnergy->InsertValue(cellCounter, c->energy());

          cellCounter++;
        }// end cell loop

        // add cell data to vtk unstructured grid
        uGrid->GetCellData()->AddArray(cellID);
        uGrid->GetCellData()->AddArray(cellType);
        uGrid->GetCellData()->AddArray(cellVolume);
        uGrid->GetCellData()->AddArray(cellSurfArea);
        uGrid->GetCellData()->AddArray(cellEnergy);

        // add time stamp to vtk unstructured grid
        vtkTimeArray->InsertValue(0,t.time());
        uGrid->GetFieldData()->AddArray(vtkTimeArray);

        // add vertices to vtk unstructured grid
        uGrid->SetPoints(points);

        vtkSmartPointer<vtkXMLUnstructuredGridWriter> vtkWriter =
           vtkSmartPointer<vtkXMLUnstructuredGridWriter>::New();
        vtkWriter->SetInputData(uGrid);
        vtkWriter->SetFileName(filename);
        vtkWriter->SetDataModeToBinary();
        //vtkWriter->SetDataModeToAscii(); // if you want to see numeric files, but it's 10x larger
        vtkWriter->Update();

        // add .vtu file to timeseries file
        vtkTimeseries << "    <DataSet timestep=\"" << t.time() <<
                             "\" group=\"\" part=\"0\" file=\"" << filename << "\"/>" << endl;

      }// end outputVtkVis if
    }
    // compute time step
    t.timeStep(DeltaT);
  }

  // finish timeseries files
  vtkTimeseries << "  </Collection>" << endl;
  vtkTimeseries << "</VTKFile>" << endl;
  vtkTimeseries.close();

  return 0;
}

