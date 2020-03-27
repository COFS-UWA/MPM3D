#===================================
#
# External Dependencies
#
#===================================
set(VENDORS_DIR ${CMAKE_SOURCE_DIR}/Vendors/)

# Qt5
find_package(Qt5 REQUIRED COMPONENTS Widgets OpenGL)

#set(CMAKE_AUTOMOC ON) # Qt5 c++ extension
#set(CMAKE_AUTORCC ON)
#set(CMAKE_AUTOUIC ON)

# hdf5
set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED)

# Eigen
set(EIGEN_INCLUDE_DIR ${VENDORS_DIR}/Eigen/include/)

# OpenGL
# find_library(OPENGL_LIB
    # NAMES opengl32.lib
    # )
# set(OPENGL_LIBRARIES ${OPENGL_LIB})