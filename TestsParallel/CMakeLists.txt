#===================================
#
# TestsParallel Project
#
#===================================
add_executable(
    TestsParallel
    # precompiled header
    TestsParallel_pcp.h TestsParallel_pcp.cpp
    test_parallel_utils.h test_parallel_utils2.h
    #test_parallel_utils2.h
    test_simulations_omp.h
    test_model_view_omp.h
    # entry
    TestsParallel_main.cpp
    test_rigid_object.cpp
    test_contact_model_3d.cpp
    test_omp_barrier_time.cpp
    test_von_mises.cpp
    test_Step_OMP.cpp
    test_t2d_me_mt_test1.cpp
    test_t2d_me_mt_test2.cpp
    test_t2d_me_s_test2.cpp
    test_t2d_me_mt_cap_compression.cpp
    test_t2d_me_mt_strip_footing_smaller.cpp
    test_t2d_me_mt_strip_footing.cpp
    test_t3d_me_s_1d_compression.cpp
    test_t3d_me_mt_1d_compression.cpp
    test_t3d_me_mt_cap_compression.cpp
    test_t3d_me_mt_cylinder_foundation.cpp
    test_t3d_me_mt_cylinder_foundation_den.cpp
    test_t3d_me_mt_block_sliding.cpp
    test_t3d_me_mt_test_rigid_mesh.cpp
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
