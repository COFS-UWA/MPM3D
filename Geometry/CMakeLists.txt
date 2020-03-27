#===================================
#
# Geometry Project
#
#===================================
add_library(
    Geometry STATIC
    # precomplied header
    Geometry_pcp.h Geometry_pcp.cpp
    #
    Geometry.h Geometry.cpp
    TetrahedronMesh.h TetrahedronMesh.cpp
    )

set(GEOMETRY_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    Geometry PUBLIC
    # Internal
    ${MEMORYUTILS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    )
