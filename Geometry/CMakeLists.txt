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
    GeometryUtils.h GeometryUtils.cpp
    Geometry2D.h
    Geometry3D.h
    NumPairHashTable.hpp
    # 2D geometry
    TriangleUtils.h
    TriangleUtils_f.h
    TriangleMeshTemplate.hpp
    SearchingGrid2D.hpp
    PDSBgGrid2D.hpp
    TriangleMesh.h TriangleMesh.cpp
    #RandomPointQueue.h RandomPointQueue.cpp
    ParticleGenerator2D.hpp
    # 3D geometry
    TetrahedronUtils.h
    TetrahedronMeshTemplate.hpp
    TetrahedronMesh.h
    SearchingGrid3D.hpp
    ParticleGenerator3D.hpp
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
