cmake_minimum_required (VERSION 3.1)
project (tvm)

set (CMAKE_CXX_STANDARD 11)
set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# for make: cmake -G "Unix Makefiles" ../src
#set(CMAKE_BUILD_TYPE Debug)
set(CMAKE_BUILD_TYPE Release)
#set(CMAKE_BUILD_TYPE RelWithDebInfo)
# or just: cmake -DCMAKE_BUILD_TYPE=Release -G "Unix Makefiles" ../src

set(CMAKE_CONFIGURATION_TYPES "RelWithDebInfo;Debug;Release;MinSizeRel")
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake/Modules)

### packages
find_package(Eigen3 REQUIRED)
#set(EIGEN3_INCLUDE_DIR "~/Code/lib/")

set(NETCDF_CXX "YES")
find_package(NetCDF)
if(NETCDF_FOUND)
  message(STATUS "NetCDF library found!")
  add_definitions(-DUSE_NETCDF=1)
else()
  message(STATUS "NetCDF library not found!")
  add_definitions(-DUSE_NETCDF=0)
endif()

find_package(VTK COMPONENTS
  vtkCommonColor
  vtkCommonCore
  vtkCommonDataModel
  vtkIOXML
  vtkInteractionStyle
  vtkRenderingCore
  vtkRenderingFreeType)
if (NOT VTK_FOUND)
  message("VTK library not found")
  return ()
endif()
message (STATUS "VTK_VERSION: ${VTK_VERSION}")
include(${VTK_USE_FILE})

### includes
include_directories(${CMAKE_CURRENT_SOURCE_DIR} ${EIGEN3_INCLUDE_DIR} ${NETCDF_INCLUDE_DIRS})


aux_source_directory(. DIR_SRCS)
add_subdirectory(Run)
add_subdirectory(Vertex)
add_subdirectory(Edge)
add_subdirectory(Polygon)
add_subdirectory(Cell)
add_subdirectory(Energy)
add_subdirectory(Reconnection)
add_subdirectory(misc)
add_executable(tvm tvm.cpp)
target_link_libraries(tvm run vertex edge polygon cell energy reconnection ${VTK_LIBRARIES} -lstdc++fs)

