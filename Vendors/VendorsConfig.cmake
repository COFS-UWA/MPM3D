#===================================
#
# External Dependencies
# 
#===================================

# vcpkg config
include($ENV{VCPKG_ROOT}/scripts/buildsystems/vcpkg.cmake)
# tbb
find_package(TBB CONFIG REQUIRED)
# hdf5
find_package(hdf5 CONFIG REQUIRED)
# freetype
find_package(freetype CONFIG REQUIRED)

# Qt
set(Qt5_DIR $ENV{Qt5_ROOT}/5.15.2/msvc2019_64/lib/cmake/Qt5)
find_package(Qt5 REQUIRED COMPONENTS Widgets OpenGL)

# Vendors
# Eigen
set(EIGEN_DIR $ENV{VENDORS_ROOT}/Eigen)
set(EIGEN_INCLUDE_DIR ${EIGEN_DIR}/include/)

# Deprecated
# # FreeType
# set(FREETYPE_DIR $ENV{VENDORS_ROOT}/freetype/)
# set(FREETYPE_INCLUDE_DIR "${FREETYPE_DIR}include/freetype2/") # must to freetype2 dir
# # debug
# set(FREETYPE_DEBUG_LIBRARIES_DIR "${FREETYPE_DIR}lib_debug/")
# find_library(FREETYPE_DEBUG_LIBRARIES
    # NAMES freetyped
    # PATHS ${FREETYPE_DEBUG_LIBRARIES_DIR}
    # )
# unset(FREETYPE_DEBUG_LIBRARIES_DIR CACHE)
# # release
# set(FREETYPE_RELEASE_LIBRARIES_DIR "${FREETYPE_DIR}lib/")
# find_library(FREETYPE_RELEASE_LIBRARIES
    # NAMES freetype
    # PATHS ${FREETYPE_RELEASE_LIBRARIES_DIR}
    # )
# unset(FREETYPE_RELEASE_LIBRARIES_DIR CACHE)

# # vcpkg
# set(VCPKG_ROOT "C:/Softwares/vcpkg/scripts/buildsystems/vcpkg.cmake" CACHE STRING "")
# set(CMAKE_TOOLCHAIN_FILE ${VCPKG_ROOT} CACHE STRING "")

# # Pathes of external dependencies (user defined part)
# # Qt5
# set(Qt5_DIR "C:/Softwares/Qt/5.15.0/msvc2019_64/lib/cmake/Qt5")
# # hdf5
# set(HDF5_DIR "C:/Vendors/HDF5-1.12.0-win64/cmake/hdf5/")
# # free type
# set(FREETYPE_DIR "C:/Vendors/freetype/")
# # Eigen
# set(EIGEN_INCLUDE_DIR "C:/Vendors/Eigen/include/")
# # tbb
# find_package(TBB CONFIG REQUIRED)
# # add external dependencies
# list(APPEND CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/Vendors/") # cmake file path
# include(Vendors REQUIRED) # call Vendors.cmake

# set(VENDORS_DIR "${CMAKE_SOURCE_DIR}/Vendors/")

# # Qt5
# find_package(Qt5 REQUIRED COMPONENTS Widgets OpenGL)

# #set(CMAKE_AUTOMOC ON) # Qt5 c++ extension
# #set(CMAKE_AUTORCC ON)
# #set(CMAKE_AUTOUIC ON)

# # hdf5
# set(HDF5_USE_STATIC_LIBRARIES ON)
# find_package(HDF5 REQUIRED)

# # OpenGL
# # find_library(OPENGL_LIB
    # # NAMES opengl32.lib
    # # )
# # set(OPENGL_LIBRARIES ${OPENGL_LIB})

# # FreeType
# set(FREETYPE_INCLUDE_DIR "${FREETYPE_DIR}include/freetype2/") # must to freetype2 dir
# # debug
# set(FREETYPE_DEBUG_LIBRARIES_DIR "${FREETYPE_DIR}lib_debug/")
# find_library(FREETYPE_DEBUG_LIBRARIES
    # NAMES freetyped
    # PATHS ${FREETYPE_DEBUG_LIBRARIES_DIR}
    # )
# unset(FREETYPE_DEBUG_LIBRARIES_DIR CACHE)
# # release
# set(FREETYPE_RELEASE_LIBRARIES_DIR "${FREETYPE_DIR}lib/")
# find_library(FREETYPE_RELEASE_LIBRARIES
    # NAMES freetype
    # PATHS ${FREETYPE_RELEASE_LIBRARIES_DIR}
    # )
# unset(FREETYPE_RELEASE_LIBRARIES_DIR CACHE)
