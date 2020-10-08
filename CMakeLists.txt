cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project("Implementation of the 3D Voronoi model")


### compiler & flags
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED on)
set(CMAKE_CC_FLAGS "${CMAKE_CC_FLAGS} -Wall -pedantic")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -pedantic -std=c++0x -Wno-deprecated-declarations")

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


### subdirectories
add_subdirectory(voro++)
add_subdirectory(misc)
add_subdirectory(povray)
add_subdirectory(core)
add_subdirectory(spheres)


### executables
add_executable(exampleHeterogeneousCellTypes.exe exampleHeterogeneousCellTypes.cpp)
target_link_libraries(exampleHeterogeneousCellTypes.exe core povray ${VTK_LIBRARIES} -lstdc++fs)

add_executable(exampleEnergyMinimization.exe exampleEnergyMinimization.cpp)
target_link_libraries(exampleEnergyMinimization.exe core miscMinimization)
