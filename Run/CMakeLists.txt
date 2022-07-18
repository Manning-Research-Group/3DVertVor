set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

aux_source_directory(. DIR_LIB_SRCS)
add_library (run ${DIR_LIB_SRCS})
target_link_libraries(run miscGeometry miscOther miscMinimization ${NETCDF_LIBRARIES})
