#===================================
#
# TestsParallel Project
#
#===================================
add_executable(
    TestsParallel
    # precompiled header
    TestsParallel_pcp.h TestsParallel_pcp.cpp
    # entry
    TestsParallel_main.cpp
    )

target_include_directories(
    TestsParallel PUBLIC
    # Internal dependencies
    ${COMMONUTILS_INCLUDE_DIR}
    ${FILEUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${MATERIALMODELS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    ${SIMULATIONS_INCLUDE_DIR}
    ${SIMULATIONS_OMP_INCLUDE_DIR}
    ${MODELVIEWER_INCLUDE_DIR}
    # External dependencies
    ${EIGEN_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    )

target_link_libraries(
    TestsParallel PUBLIC
    # Internal dependencies
    CommonUtils
    MemoryUtils
    MaterialModels
    Geometry
    Simulations
    SimulationsOMP
    ModelViewer
    # External dependencies
    Qt5::Widgets
    Qt5::OpenGL
    ${HDF5_LIBRARIES}
    )

if (WIN32)
    # supress warning for some C functions in visual studio
    target_compile_definitions(TestsParallel PRIVATE _CRT_SECURE_NO_WARNINGS)
    # use precompiled header
    target_compile_options(TestsParallel PRIVATE -YcTestsParallel_pcp.h)
endif()
