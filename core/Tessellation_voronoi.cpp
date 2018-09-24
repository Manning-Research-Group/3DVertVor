#include <voro++/voro++.hh>

#include "misc/other/misc-math.h"
#include "misc/geometry/Vector3D.h"

#include "Tessellation.h"

#include "VertexOfCell-inline.h"
#include "Cell-inline.h"

//const double Tessellation::DistanceSqCutoffForPeriodicityInference = 2e-10;
const double Tessellation::DistanceSqCutoffForPeriodicityInference = 1e-8;

void Tessellation::resetTopology() {
  for(unsigned int i=0; i<_cells.size(); ++i) {
    Cell *c=_cells[i];
    c->resetTopology();
  }
  _topologicalElementsPresent = false;
}

void Tessellation::setTopologyByVoronoiTesselation() {
//  std::cout << "Creating and filling voro container..." << std::endl;
  // it is recommended that the number of cells per computational bock should be on the order of 5, maybe a bit higher
  const double ComputationalBlockSize = cubicRoot(5.0*_box.volume()/_cells.size());
  const int NumComputationalBlocksX = ceil(_box.x/ComputationalBlockSize);
  const int NumComputationalBlocksY = ceil(_box.y/ComputationalBlockSize);
  const int NumComputationalBlocksZ = ceil(_box.z/ComputationalBlockSize);

  // create voronoi container with periodic boundary conditions with skew
  voro::container_periodic voronoiContainer(_box.x, 
                                            _box.shearYx*_box.y, _box.y, 
                                            0.0, 0.0, _box.z, 
                                            NumComputationalBlocksX, NumComputationalBlocksY, NumComputationalBlocksZ,
                                            16);
  
  // fill with cells
  for(unsigned int i=0; i<_cells.size(); ++i) {
    Cell *c = _cells[i];
//    c->moveIntoBox();
    const Vector3D &p = c->position();
    // check for nans:
    if(std::isnan(p.x()) || std::isnan(p.y()) || std::isnan(p.z())) {
      std::cerr << "Tessellation::setTopologyByVoronoiTesselation: Particle " << i << " has position containing nans: " << p << std::endl;
      exit(1);
    }
    // otherwise add to container
    voronoiContainer.put(i, p.x(), p.y(), p.z());
  }

  // Here, I rely on the assumption that the indices in the voro module coincide with the indices in the _cells array!
  // loop over all cells in the voro module
//  std::cout << "Translating topology..." << std::endl;
  voro::c_loop_all_periodic voroParticleLoop(voronoiContainer);
  if(voroParticleLoop.start()) {
    do {
      const int CurrentCellIndex = voroParticleLoop.pid();
      Cell *cell = _cells[CurrentCellIndex];

      voro::voronoicell_neighbor voroCell;
      if(voronoiContainer.compute_cell(voroCell, voroParticleLoop)) {
        // to map the combined vertex + bond index to an edge pointer
        DirectedEdgeOfCell **cellEdges = new DirectedEdgeOfCell*[3*voroCell.p];
        
        // create the vertices and edges
        for(int vertexIndex=0; vertexIndex<voroCell.p; ++vertexIndex) {
          // create vertex within cell
          Vector3D twoTimesRelativeVertexPosition(voroCell.pts[3*vertexIndex], voroCell.pts[3*vertexIndex+1], voroCell.pts[3*vertexIndex+2]);
          VertexOfCell *vertex = cell->newVertex(0.5*twoTimesRelativeVertexPosition);

          // check number of bonds
          if(voroCell.nu[vertexIndex]!=3) {
            std::cerr << "Tessellation::setTopologyByVoronoiTesselation: Wrong number of bonds per vertex (is " << voroCell.nu[vertexIndex] << " but should be 3)!" << std::endl;
            exit(1);
          }
          for(int relEdgeIndex=0; relEdgeIndex<3; ++relEdgeIndex) {
            DirectedEdgeOfCell *edge = cell->newEdge(vertex);
            cellEdges[3*vertexIndex + relEdgeIndex] = edge;
            vertex->addEdge(edge);
          }
        }
        
        // wire bonds
        for(int vertexIndex=0; vertexIndex<voroCell.p; ++vertexIndex) {
          for(int edgeIndex=0; edgeIndex<3; ++edgeIndex) {
            int nextVertexIndex = voroCell.ed[vertexIndex][edgeIndex];
            int conjugatedEdgeIndex = voroCell.ed[vertexIndex][3+edgeIndex];
            DirectedEdgeOfCell *conjEdge = cellEdges[3*nextVertexIndex + conjugatedEdgeIndex];
                        
            int nextEdgeIndex = (conjugatedEdgeIndex+1)%3;
            DirectedEdgeOfCell *nextEdge = cellEdges[3*nextVertexIndex + nextEdgeIndex];

            cellEdges[3*vertexIndex + edgeIndex]->setConjugatedAndNext(conjEdge, nextEdge);
          }
        }
        
        // create faces
        for(int vertexIndex=0; vertexIndex<voroCell.p; ++vertexIndex) {
          // loop over all bonds
          for(int edgeIndex = 0; edgeIndex < voroCell.nu[vertexIndex]; ++edgeIndex) {
            DirectedEdgeOfCell *firstEdge = cellEdges[3*vertexIndex + edgeIndex];
            if(!firstEdge->face()) {
              // first, create face
              int otherCellIndex = voroCell.ne[vertexIndex][edgeIndex];
              DirectedFace *face = cell->newFace(_cells[otherCellIndex], firstEdge);
              
              // loop around face and wire
              DirectedEdgeOfCell *edge = firstEdge;
              do {
                edge->setFace(face);
                edge = edge->nextAroundFace();
              } while(edge!=firstEdge);
            }
          }
        }
        
        // free allocated memory
        delete[] cellEdges;
      } else {
        Vector3D voroPos(voroParticleLoop.x(), voroParticleLoop.y(), voroParticleLoop.z());
        std::cerr << "Tessellation::setTopologyByVoronoiTesselation: Error computing cell for particle " << voroParticleLoop.pid() << " at voro pos: " << voroPos << ", my pos: " << _cells[voroParticleLoop.pid()]->position() << "!" << std::endl;
        exit(1);
      }
    } while(voroParticleLoop.inc());
  }
  
//  // get cells out of the box
//  for(Cell *c : _cells) {
////    c->outOfBox();
//    c->restorePosition();
//  }
  
//  std::cout << "Figuring out conjugated face pairs and periodicities..." << std::endl;
  // figure out conjugated face pairs and compute periodicity vectors
  // buffer average positions
  for(unsigned int i=0; i<_cells.size(); ++i) {
    Cell *c=_cells[i];
    for(unsigned int j=0; j<c->faces().size(); ++j) {
      DirectedFace *f=c->faces()[j];
      f->bufferAverageVertexVoronoiPositions();
    }
  }
  // now match directed faces and set periodicities
  for(Cell *c : _cells) {
    for(DirectedFace *f : c->faces()) {
      Cell *oc = f->otherCell();
      DirectedFace *firstFaceFound = NULL;
      int facesFound = 0;
      for(DirectedFace *f2: oc->faces()) {
        if(f2->otherCell()==c) {
          firstFaceFound = f2;
          ++facesFound;
        }
      }
      if(facesFound==0) {
        std::cerr << "Tessellation::setTopologyByVoronoiTesselation: Could not find conjugated face!" << std::endl;
        exit(1);
      } else if(facesFound==1) {
        Vector3D dist = f->averageVertexVoronoiPositions() - firstFaceFound->averageVertexVoronoiPositions();
        PeriodicityVector3D p = _box.computeDistancePeriodicityAndRest(dist);
        double normSq = dist.normSq();
        if(normSq>DistanceSqCutoffForPeriodicityInference) {
          std::cerr << "Tessellation::setTopologyByVoronoiTesselation: Corresponding faces are too far away! DistanceSq: " << normSq << std::endl;
          exit(1);
        } else {
          f->setPeriodicity(p);
          f->setConjugated(firstFaceFound);
        }
      } else { // In this case we have to do it the ugly way.
        Vector3D centerOfF = f->averageVertexVoronoiPositions();
        bool found = false;
        for(DirectedFace *f2: oc->faces()) {
          Vector3D dist = centerOfF-f2->averageVertexVoronoiPositions();
          PeriodicityVector3D p = _box.computeDistancePeriodicityAndRest(dist);
          double normSq = dist.normSq();
          if(normSq<=DistanceSqCutoffForPeriodicityInference) {
            if(found) {
              std::cerr << "Tessellation::setTopologyByVoronoiTesselation: Two faces are too close! DistanceSq: " << normSq << std::endl;
              exit(1);
            } else {
              found = true;
              f->setPeriodicity(p);
              f->setConjugated(f2);
            }
          }
        }
        if(!found) {
          std::cerr << "Tessellation::setTopologyByVoronoiTesselation: Could not find conjugated face!" << std::endl;
          exit(1);
        }
      }
    }
  }

  _topologicalElementsPresent = true;
//  std::cout << "Done." << std::endl;
}