# Heterotypic interfacial tension inside the 3D vertex model
This work adds updated visualization and multiple cell types to the open-source 3D vertex model code created by Zhang and Schwarz. The original open-source code can be found at https://github.com/ZhangTao-SJTU/tvm and the original paper describing the code and it's applications is found at https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.4.043148 or on arXiv at https://arxiv.org/abs/2204.07081. 

We add the capibility for multiple cell types as well as allow for heterotypic interfacial tension between different cell types. In addition, we change the visualization to be cell focused instead of polygon focused and use .vtu files instead of .vtk files. 

## Additions and alterations to the tvm code
The primary additions to the code are the following:

Run: Run.h / Run.cpp  <br />
Energy: Interface.h / Interface.cpp  <br />
Polygon: Polygon.h / Polygon.cpp <br />
Vertex: Vertex.h / Vertex.cpp
                
<br />                
With minor changes to the following:  <br />
main.cpp  <br />
tvm.cpp  <br />
Energy: Volume.cpp
