# need to include bin folder in system path !!
cmake_minimum_required(VERSION 3.15.0)

project(SimpleMPM3)

# External dependencies
# Qt5
set(Qt5_DIR D:/Softwares/Qt/5.14.1/msvc2015_64/)
# hdf5
set(HDF5_DIR D:/Vendors/HDF5-1.12.0-win64/cmake/hdf5/)
# add external dependencies
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Vendors/) # cmake file path
include(Vendors REQUIRED)

# Utilities
add_subdirectory(Utilities/)

# Material Models
add_subdirectory(MaterialModels/)

# Geometry
add_subdirectory(Geometry/)

# Simulations
add_subdirectory(Simulations/)

# ModelViewer
add_subdirectory(ModelViewer/)

# Tests
add_subdirectory(Tests/)
