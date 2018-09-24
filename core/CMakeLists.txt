### library
add_library(core Cell.cpp DirectedFace.cpp VertexOfCell.cpp Tessellation.cpp Tessellation_voronoi.cpp 
        PeriodicBox.cpp CellPositionMailbox.cpp 2ndEnergyDerivatives.cpp)
target_link_libraries(core povray miscGeometry miscOther miscMinimization voro++ ${NETCDF_LIBRARIES})
