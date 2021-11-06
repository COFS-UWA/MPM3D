#===================================
#
# External Dependencies
#
#===================================
set(VENDORS_DIR "${CMAKE_SOURCE_DIR}/Vendors/")

# Qt5
find_package(Qt5 REQUIRED COMPONENTS Widgets OpenGL)

#set(CMAKE_AUTOMOC ON) # Qt5 c++ extension
#set(CMAKE_AUTORCC ON)
#set(CMAKE_AUTOUIC ON)

# hdf5
set(HDF5_USE_STATIC_LIBRARIES ON)
find_package(HDF5 REQUIRED)

# OpenGL
# find_library(OPENGL_LIB
    # NAMES opengl32.lib
    # )
# set(OPENGL_LIBRARIES ${OPENGL_LIB})

# FreeType
set(FREETYPE_INCLUDE_DIR "${FREETYPE_DIR}include/freetype2/") # must to freetype2 dir
# debug
set(FREETYPE_DEBUG_LIBRARIES_DIR "${FREETYPE_DIR}lib_debug/")
find_library(FREETYPE_DEBUG_LIBRARIES
    NAMES freetyped
    PATHS ${FREETYPE_DEBUG_LIBRARIES_DIR}
    )
unset(FREETYPE_DEBUG_LIBRARIES_DIR CACHE)
# release
set(FREETYPE_RELEASE_LIBRARIES_DIR "${FREETYPE_DIR}lib/")
find_library(FREETYPE_RELEASE_LIBRARIES
    NAMES freetype
    PATHS ${FREETYPE_RELEASE_LIBRARIES_DIR}
    )
unset(FREETYPE_RELEASE_LIBRARIES_DIR CACHE)

# Tbb
# set(TBB_INCLUDE_DIR "${TBB_DIR}../include/")
# set(TBB_LIBRARIES_DIR "${TBB_DIR}../lib/intel64/vc14/")
# find_package(TBB REQUIRED)
