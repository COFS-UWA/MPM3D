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
    BC.h
    # utilities
    ParticleGenerator.h ParticleGenerator.cpp
    )

target_include_directories(
    Simulations PUBLIC
    # Internal dependencies
    ${SIMULATIONS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${CONSTITUTIVEMODEL_INCLUDE_DIR}
    # External dependencies
    ${EIGEN_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    )

target_link_libraries(
    Simulations PUBLIC
    # Internal dependencies
    MemoryUtils
    MaterialModels
    # External dependencies
    ${HDF5_LIBRARIES}
    )

if (MSVC)
    target_compile_definitions(Simulations PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(Simulations PRIVATE -YcSimulations_pcp.h)
endif()
