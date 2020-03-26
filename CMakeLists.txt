# need to include bin folder in system path !!
cmake_minimum_required(VERSION 3.15.0)

project(SimpleMPM3)

# External dependencies
# packages path
set(Qt5_DIR D:/Softwares/Qt/5.14.1/msvc2015_64/) # Qt5
set(HDF5_DIR D:/Vendors/HDF5-1.10.6-win64/cmake/hdf5/) # Hdf5
# add external dependencies
list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/Vendors/)
include(Vendors REQUIRED)

# MemoryUtils
set(MEMORYUTILS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/MemoryUtils/)
add_subdirectory(MemoryUtils/)

# File Operations
set(FILEOPS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/FileOps/)
add_subdirectory(FileOps/)

# Material Models
set(MATERIALMODELS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/MaterialModels/)
add_subdirectory(MaterialModels/)

# Simulations
set(SIMULATIONS_INCLUDE_DIR ${CMAKE_SOURCE_DIR}/Simulations/)
add_subdirectory(Simulations/)

# Tests

# ResultViewerQt
set(RESULTVIEWERQT_INCLUDE_DIR
    ${CMAKE_SOURCE_DIR}/ResultViewerQt/
    )
add_subdirectory(ResultViewerQt/)
