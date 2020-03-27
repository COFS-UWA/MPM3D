#===================================
#
# Simulations Project
#
#===================================
add_library(
    Simulations STATIC
    # precomplied header
    Simulations_pcp.h Simulations_pcp.cpp
    # model
    Model.h
    Step.h Step.cpp
    TimeHistory.h TimeHistory.cpp
    ResultFile.h
    # boundary conditions
    BCs.h
    # utilities
    ParticleGenerator3D.h ParticleGenerator3D.cpp
    )

set(SIMULATIONS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    Simulations PUBLIC
    # Internal dependencies
    ${COMMONUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${FILEUTILS_INCLUDE_DIR}
    ${MATERIALMODELS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    ${SIMULATIONS_INCLUDE_DIR}
    # External dependencies
    ${EIGEN_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    )

target_link_libraries(
    Simulations PUBLIC
    # Internal dependencies
    MemoryUtils
    MaterialModels
    Geometry
    # External dependencies
    ${HDF5_LIBRARIES}
    )

if (MSVC)
    target_compile_definitions(Simulations PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(Simulations PRIVATE -YcSimulations_pcp.h)
endif()
