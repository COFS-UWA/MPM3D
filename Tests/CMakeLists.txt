#===================================
#
# Tests Project
#
#===================================
add_executable(
    Tests 
    # precompiled header
    Tests_pcp.h Tests_pcp.cpp
    # entry
    Tests_main.cpp
    # header
    utils.h utils.cpp
    test_utilities.h
    test_geometry.h
    test_material_models.h
    test_simulations.h
    test_model_view.h
    # unit tests
    # utilities
    #test_stack_and_link_list.cpp
    #test_solve_functions.cpp
    # material models
    #test_model_container.cpp
    #test_mcc_get_Su.cpp
    #test_mcc_compression.cpp
    #test_undrained_mcc.cpp
    # 2D
    # back ground searching grid
    #test_searching_grid2d.cpp
    # random point generation
    #test_random_point_for_Model_T2D_ME_s.cpp
    # display model
    #test_Model_T2D_ME_s_display.cpp
    # geostatic calculation
    test_t2d_me_s_geostatic.cpp
    #test_t2d_me_s_geostatic_restart.cpp
    test_t2d_chm_s_geostatic_mcc.cpp
    # deformation calculation
    test_t2d_me_s_1d_compression.cpp
    #test_t2d_me_s_1d_compression_horizontal.cpp
    #test_t2d_me_s_plate_with_hole.cpp
    test_t2d_me_s_cap_compression.cpp
    test_t2d_me_s_t_bar_smaller_soil.cpp
    test_t2d_me_s_shallow_foundation_smaller_soil.cpp
    test_t2d_me_s_shallow_foundation.cpp
    test_t2d_chm_s_1d_consolidation.cpp
    #test_t2d_chm_s_undrained_compression.cpp
    test_t2d_chm_s_t_bar_smaller_soil.cpp
    test_t2d_chm_s_ud_t_bar_smaller_soil.cpp
    # restart utilities
    #test_t2d_me_s_test_restart.cpp
    #test_t2d_chm_s_test_restart.cpp
    # rigid circle
    #test_t2d_me_s_test_rigid_circle.cpp
    #test_t2d_chm_s_test_rigid_circle.cpp
    # parallel version
    #test_t2d_me_p_geostatic.cpp
    #test_t2d_me_p_test.cpp
    #test_t2d_me_p_1d_compression.cpp
    #test_t2d_me_p_rigid_circle_penetration.cpp
    # shallow pipe embedment for conference
    # coupled hydro-mechanics
    test_t2d_chm_s_pipe_conference_geo.cpp
    test_t2d_chm_s_pipe_conference_restart1.cpp
    test_t2d_chm_s_pipe_conference_restart2.cpp
    # completely drained
    #test_t2d_me_s_pipe_conference_geo.cpp
    #test_t2d_me_s_pipe_conference_restart.cpp
    #test_t2d_me_s_pipe_conference_restart2.cpp
    #test_t2d_me_p_pipe_conference_geo.cpp
    #test_t2d_me_p_pipe_conference_restart.cpp
    # no geostatic stress
    test_t2d_me_s_pipe_conference_drained_no_geostress.cpp
    test_t2d_chm_s_pipe_conference_undrained_no_geostress.cpp
    # undrained
    #test_t2d_me_s_pipe_conference_geo_undrained.cpp
    # 3D
    test_tetrahedron_mesh.cpp
    test_Model_T3D_ME_s_display.cpp
    test_RigidTetrahedronMesh.cpp
    test_RigidTetrahedronMesh_display.cpp
    test_RigidTetrahedronMesh_rotation.cpp
    test_searching_grid3d.cpp
    test_t3d_me_s_1d_compression.cpp
    test_t3d_me_s_cap_compression.cpp
    test_t3d_chm_s_1d_consolidation.cpp
    test_fem_t3d_me_s_1d_compression.cpp
    test_fem_t3d_me_s_1d_compression2.cpp
    test_mat_model_utils.cpp
    test_von_mises.cpp
    test_tresca.cpp
    # sand hypoplasticity
    test_hypoplasticity_match_Wichtmann.cpp
    )

target_include_directories(
    Tests PUBLIC
    # Internal dependencies
    ${COMMONUTILS_INCLUDE_DIR}
    ${FILEUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${MATERIALMODELS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    ${SIMULATIONS_INCLUDE_DIR}
    ${MODELVIEWER_INCLUDE_DIR}
    # External dependencies
    ${EIGEN_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${TBB_INCLUDE_DIR}
    )

target_link_directories(
    Tests PUBLIC
    ${TBB_LIBRARIES_DIR}
    )

target_link_libraries(
    Tests PUBLIC
    # Internal dependencies
    CommonUtils
    MemoryUtils
    MaterialModels
    Geometry
    Simulations
    ModelViewer
    # External dependencies
    Qt5::Widgets
    Qt5::OpenGL
    ${HDF5_LIBRARIES}
    TBB
    )

if (WIN32)
    # supress warning for some C functions in visual studio
    target_compile_definitions(Tests PRIVATE _CRT_SECURE_NO_WARNINGS)
    # use precompiled header
    target_compile_options(Tests PRIVATE -YcTests_pcp.h)
endif()
