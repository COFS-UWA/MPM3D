cmake_minimum_required(VERSION 3.15.0)

project(SimpleMPM3)

# External dependencies
include(${PROJECT_SOURCE_DIR}/Vendors/VendorsConfig.cmake)

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
