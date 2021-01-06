cmake_minimum_required(VERSION 3.15.0)

project(SimpleMPM3)

# Pathes of external dependencies (user defined part)
# Qt5
set(Qt5_DIR "D:/Softwares/Qt/5.15.0/msvc2019_64/lib/cmake/Qt5")
# hdf5
set(HDF5_DIR "D:/Vendors/HDF5-1.12.0-win64/cmake/hdf5/")
# free type
set(FREETYPE_DIR "D:/Vendors/freetype/")
# tbb
set(TBB_DIR "D:/Vendors/tbb/cmake/")
# add external dependencies
list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Vendors/") # cmake file path
include(Vendors REQUIRED) # call Vendors.cmake

# Utilities
add_subdirectory(Utilities/)

# Material Models
add_subdirectory(MaterialModels/)

# Geometry
add_subdirectory(Geometry/)

# Simulations
add_subdirectory(Simulations/)

# Simulations OMP
add_subdirectory(SimulationsOMP/)

# ModelViewer
add_subdirectory(ModelViewer/)

# Tests
add_subdirectory(Tests/)

# TestsParallel
add_subdirectory(TestsParallel/)
