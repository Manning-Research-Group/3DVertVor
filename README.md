This code was originally developed by Matthias Merkel and described below. There have been minor alterations to the 3D Voronoi code but it primarily remains unchanged from its original state.

With minor changes to the following:  <br />
Core: Cell.h / Cell.cpp 

# 3D-Voronoi
If this software is being used for a scientific publication, the user is encouraged to cite the following publication:
> Matthias Merkel, M. Lisa Manning, "A geometrically controlled rigidity transition in a model for confluent 3D tissues", _New Journal of Physics_ **20**, p. 022002 (2018).

Please also see the `LICENSE` file.

## Dependencies
This software has been tested on Linux and MacOs.  You might have to do some work in order to use it on Windows.

- an at least partly C++11-compatible compiler (code was tested with gcc versions 4.6.3 and 5.4.0)
- cmake
- the Eigen3 library (http://eigen.tuxfamily.org/index.php?title=Main_Page)
- to render 3D images:  povray
- to create movies:  ffmpeg
- _[optional] netcdf in order to exchange configurations with the LiuJamming code_

## Compilation
1. In some parent folder (`parent/`), create the folders `parent/src/` and `parent/build/`.
2. Git clone (or unpack) this code into the `parent/src/` folder.
3. Cd into the `parent/build/` folder and run `cmake ../src`.
4. In the `parent/build/` folder, run `make`.

### Recompilation
Unless any of the `CMakeLists.txt` files is changed, step 4 above is sufficient for recompilation.  Otherwise, in folder `parent/build/`, run `rm -r`, and then steps 3 & 4 above.

### If there is a problem with the Eigen3 library
1. In `parent/src/CMakeLists.txt`, comment out the line `find_package(Eigen3 REQUIRED)` and uncomment the line below, pointing it to your Eigen3 folder.
2. In `parent/build/` run `rm -r`.
3. Repeat steps 3 & 4 above.

## Directory structure
- `parent/src/`:  all source code files, including cmake files.  This folder always remains free of any build-related files (object files, binaries, etc).
  - `cmake/`:  files needed by cmake.
  - `computation/`:  code for dynamical matrix and compatibility matrix-related computations
  - `core/`:  main code
  - `dynamicalMatrix/`:  used for old way of computing dynamical matrix (turns out to be better)
  - `spheres/`:  version of code including spheres around the Voronoi centers
  - `misc/`:  miscellaneous general routines and classes used by the code
  - `povray/`:  an little wrapper around povray
  - `voro++`:  the slightly modified voro++ code by Christ Rycroft (see `LICENSE` file within folder, and http://math.lbl.gov/voro++/)
  - `CMakeFile.txt`:  controls the code compilation
  - `exampleEnergyMinimization.cpp`:  example illustrating simple energy minimization
  - `exampleHeterogeneousCellTypes.cpp`:  example file for how to use the code for a time-dependent simulation with two different cell types, and more functions for the 3D rendering
- `parent/build/`:  all compiled binaries will be here

## Quickstart
The quickest way to get started is to inspect the example main files.

### `exampleEnergyMinimization.cpp`
This file shows how to:
- initialize a random configuration of cells with specific homogeneous properties
- minimize this configuration,
- and read information about the final configuration

More specifically, it shows how to:
- define cell mechanical properties
- define the periodic box
- initialize a tessellation
- fill it with the cells
- initializes a conjugated gradient minimizer and set its numerical parameters
- tell the minimizer which degrees of freedom to minimize
- carry out shear-stabilized minimizations
- recompute the topology, geometry, and first and second derivative
- simply draw the configuration without much ado
- compute quantities like total system energy, elastic moduli, etc
- loop over all cells in a configuration and access their surface and volume

### `exampleHeterogeneousCellTypes.cpp`
This file shows how to:
- use several cell types with an additional interfacial tension between the
- carry out a time-dependent simulation
- create a movie from that

More specifically, it shows how to:
- define different cell types
- create an initial configuration out of randomly placed cells of different types
- create an initial configuration out of specifically placed cells of different types
- evolve the system by a time step
- access cell-cell interfacial areas
- drawing:  more specifically define the drawing, draw cells as filled polyhedra, and create cuts when drawing for better visibility of interesting features in the 3D rendering

### Pitfalls
- Always call the functions `Tessellation::computeTopologyAndGeometry()` after the configuration changed (e.g. via `Tessellation::timeStep()`, through a minimizer, or via directly changing cell positions).
- Always call `Tessellation::computeEnergyDerivatives()` before computing forces, stresses, or moduli.
- Always call `Tessellation::computeAndSolveDynamicalMatrix()` before accessing the dynamical matrix or computing elastic moduli.

## Structure of the code
- The central class is **`Tessellation`**.  It contains all **`Cell`** s, and all the topological elements needed for the bookkeeping.  It also contains a **_reference_** to the **`PeriodicBox`**.  Because of this, any changes you do to the box after the initialization **_still affects the Tessellation_**.
- The **`PeriodicBox`** contains the box dimensions and is reponsible for turing cell positions into distances between cell positions.  (For that, it uses co called "periodicity vectors".  For more on "periodicity vectors", see Ph.D. thesis of Matthias Merkel, appendix C.1 and/or appendix A.2 in Merkel & Manning, New Journal of Physics, 2018.)
- The **`Cell`** class represents a cell and the **`CellType`** struct contains the cell properties.
- The **`DirectedFace`** class represents "one half" of the interface between two cells.  It point's to it's own cell and to its "conjugated" face.
- The **`DirectedEdgeOfCell`** class represents a directed edge of a cell.  The edges of a cell are stored in a half edge structure.  To loop around the edges of a face, first call `DirectedFace::firstEdge` and then iterate calling `DirectedEdgeOfCell::nextAroundFace` on the current edge until the first edge is encountered again.
- The **`VertexOfCell`** class represents a vertex of a cell.  Note that for performance reasons, the same vertex appears four times, once for each of the four adjacent cells.

### I want to know more
The best way to learn more about the code is, e.g. if you want to access some information not contained in the examples, just check the header files of the corresponding class in `src/core/` and see whether you find an appropriate getter function.
