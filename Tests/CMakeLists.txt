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
    test_stack_and_link_list.cpp
    test_model_container.cpp
    test_mcc_get_Su.cpp
    test_solve_functions.cpp
    test_tetrahedron_mesh.cpp
    test_Model_T3D_ME_s_display.cpp
    test_searching_grid3d.cpp
    test_PropsMPM3DApp.cpp
    test_PospSingleFrame_display.cpp
    test_t3d_me_s_1d_compression.cpp
    test_t3d_chm_s_1d_consolidation.cpp
    test_fem_t3d_me_s_1d_compression.cpp
    test_fem_t3d_me_s_1d_compression2.cpp
    )

target_include_directories(
    Tests PUBLIC
    ${COMMONUTILS_INCLUDE_DIR}
    ${FILEUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${MATERIALMODELS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    ${SIMULATIONS_INCLUDE_DIR}
    ${MODELVIEWER_INCLUDE_DIR}
    )

#link_directories(${GLFW_LIBRARIES_DIR})
target_link_libraries(
    Tests PUBLIC
    CommonUtils
    MemoryUtils
    MaterialModels
    Geometry
    Simulations
    ModelViewer
    )

if (WIN32)
    # supress warning for certain functions in visual studio
    target_compile_definitions(Tests PRIVATE _CRT_SECURE_NO_WARNINGS)
    # use precompiled header
    target_compile_options(Tests PRIVATE -YcTests_pcp.h)
endif()
