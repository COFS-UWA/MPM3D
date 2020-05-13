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
    NumPairHashTable.hpp
    TetrahedronMeshTemplate.hpp
    TetrahedronMesh.h #TetrahedronMesh.cpp
    ParticleGenerator3D.hpp
    SearchingGrid3D.hpp
    )

set(GEOMETRY_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    Geometry PUBLIC
    # Internal
    ${COMMONUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    # External
    ${HDF5_INCLUDE_DIR}
    )
