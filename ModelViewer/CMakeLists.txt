#===================================
#
# ModelViewer Project
#
#===================================
# This project display model before and after calculation
set(CMAKE_AUTOMOC ON)

add_library(
    ModelViewer STATIC
    # Precomplied header
    ModelViewer_pcp.h ModelViewer_pcp.cpp
    # architecture
    QtGLView.h QtGLView.cpp
    QtController.h
    # prep processing
    QtSceneFromModel.h
    QtController_Prep.h QtController_Prep.cpp
    # post processing
    QtSceneFromHdf5.h
    QtController_Posp_Static.h QtController_Posp_Static.cpp
    QtController_Posp_Animation.h QtController_Posp_Animation.cpp
    # utilities
    # color map
    UniformColorMap.h UniformColorMap.cpp
    UniformColorMap_Abaqus.h UniformColorMap_Abaqus.cpp
    # display chars
    QtCharBitmapLoader.h QtCharBitmapLoader.cpp
    QtUniformColorMapObject.h QtUniformColorMapObject.cpp
    # field data extraction from hdf5 file
    Hdf5DataLoader.h Hdf5DataLoader.cpp
    Hdf5Field.h Hdf5Field.cpp
    Hdf5FieldExtraction.h
    Hdf5FieldExtraction_x.h Hdf5FieldExtraction_x.cpp
    Hdf5FieldExtraction_y.h Hdf5FieldExtraction_y.cpp
    Hdf5FieldExtraction_z.h Hdf5FieldExtraction_z.cpp
    Hdf5FieldExtraction_vol.h Hdf5FieldExtraction_vol.cpp
    Hdf5FieldExtraction_m.h Hdf5FieldExtraction_m.cpp
    Hdf5FieldExtraction_density.h Hdf5FieldExtraction_density.cpp
    Hdf5FieldExtraction_vx.h Hdf5FieldExtraction_vx.cpp
    Hdf5FieldExtraction_vy.h Hdf5FieldExtraction_vy.cpp
    Hdf5FieldExtraction_vz.h Hdf5FieldExtraction_vz.cpp
    Hdf5FieldExtraction_s11.h Hdf5FieldExtraction_s11.cpp
    Hdf5FieldExtraction_s22.h Hdf5FieldExtraction_s22.cpp
    Hdf5FieldExtraction_s33.h Hdf5FieldExtraction_s33.cpp
    Hdf5FieldExtraction_s12.h Hdf5FieldExtraction_s12.cpp
    Hdf5FieldExtraction_s23.h Hdf5FieldExtraction_s23.cpp
    Hdf5FieldExtraction_s31.h Hdf5FieldExtraction_s31.cpp
    Hdf5FieldExtraction_e11.h Hdf5FieldExtraction_e11.cpp
    Hdf5FieldExtraction_e22.h Hdf5FieldExtraction_e22.cpp
    Hdf5FieldExtraction_e33.h Hdf5FieldExtraction_e33.cpp
    Hdf5FieldExtraction_e12.h Hdf5FieldExtraction_e12.cpp
    Hdf5FieldExtraction_e23.h Hdf5FieldExtraction_e23.cpp
    Hdf5FieldExtraction_e31.h Hdf5FieldExtraction_e31.cpp
    Hdf5FieldExtraction_mises_strain_2d.h Hdf5FieldExtraction_mises_strain_2d.cpp
    Hdf5FieldExtraction_ee11.h Hdf5FieldExtraction_ee11.cpp
    Hdf5FieldExtraction_ee22.h Hdf5FieldExtraction_ee22.cpp
    Hdf5FieldExtraction_ee33.h Hdf5FieldExtraction_ee33.cpp
    Hdf5FieldExtraction_ee12.h Hdf5FieldExtraction_ee12.cpp
    Hdf5FieldExtraction_ee23.h Hdf5FieldExtraction_ee23.cpp
    Hdf5FieldExtraction_ee31.h Hdf5FieldExtraction_ee31.cpp
    Hdf5FieldExtraction_pe11.h Hdf5FieldExtraction_pe11.cpp
    Hdf5FieldExtraction_pe22.h Hdf5FieldExtraction_pe22.cpp
    Hdf5FieldExtraction_pe33.h Hdf5FieldExtraction_pe33.cpp
    Hdf5FieldExtraction_pe12.h Hdf5FieldExtraction_pe12.cpp
    Hdf5FieldExtraction_pe23.h Hdf5FieldExtraction_pe23.cpp
    Hdf5FieldExtraction_pe31.h Hdf5FieldExtraction_pe31.cpp
    Hdf5FieldExtraction_plastic_mises_strain_2d.h Hdf5FieldExtraction_plastic_mises_strain_2d.cpp
    #Hdf5FieldExtraction_mat_mean_stress.h Hdf5FieldExtraction_mat_mean_stress.cpp
    #Hdf5FieldExtraction_mat_mises_stress.h Hdf5FieldExtraction_mat_mises_stress.cpp
    # CHM MPM
    Hdf5FieldExtraction_n.h Hdf5FieldExtraction_n.cpp
    Hdf5FieldExtraction_m_s.h Hdf5FieldExtraction_m_s.cpp
    Hdf5FieldExtraction_m_f.h Hdf5FieldExtraction_m_f.cpp
    Hdf5FieldExtraction_density_s.h Hdf5FieldExtraction_density_s.cpp
    Hdf5FieldExtraction_density_f.h Hdf5FieldExtraction_density_f.cpp
    Hdf5FieldExtraction_vol_m.h Hdf5FieldExtraction_vol_m.cpp
    Hdf5FieldExtraction_vx_s.h Hdf5FieldExtraction_vx_s.cpp
    Hdf5FieldExtraction_vy_s.h Hdf5FieldExtraction_vy_s.cpp
    Hdf5FieldExtraction_vz_s.h Hdf5FieldExtraction_vz_s.cpp
    Hdf5FieldExtraction_vx_f.h Hdf5FieldExtraction_vx_f.cpp
    Hdf5FieldExtraction_vy_f.h Hdf5FieldExtraction_vy_f.cpp
    Hdf5FieldExtraction_vz_f.h Hdf5FieldExtraction_vz_f.cpp
    Hdf5FieldExtraction_p.h Hdf5FieldExtraction_p.cpp
    ################## 2D MPM model display ######################
    # 2D OpenGL Buffer Objects
    QtTriangleMeshGLObject.h QtTriangleMeshGLObject.cpp
    QtMonoColorCircleGLObject.h QtMonoColorCircleGLObject.cpp
    QtMultiColorCircleGLObject.h QtMultiColorCircleGLObject.cpp
    QtRigidCircleObject.h QtRigidCircleObject.cpp
    QtRigidRectObject.h QtRigidRectObject.cpp
    # Model_T2D_ME_s
    # prep process module
    QtSceneFromModel_T2D_ME_s.h QtSceneFromModel_T2D_ME_s.cpp
    QtApp_Prep_T2D_ME_s.h QtApp_Prep_T2D_ME_s.cpp
    # post process module
    QtSceneFromHdf5_T2D_ME_s.h QtSceneFromHdf5_T2D_ME_s.cpp
    QtApp_Posp_T2D_ME_s.h QtApp_Posp_T2D_ME_s.cpp
    # Model_T2D_ME_mt
    # prep process module
    QtSceneFromModel_T2D_ME_mt.h QtSceneFromModel_T2D_ME_mt.cpp
    QtApp_Prep_T2D_ME_mt.h QtApp_Prep_T2D_ME_mt.cpp
    # post process module
    QtSceneFromHdf5_T2D_ME_mt.h QtSceneFromHdf5_T2D_ME_mt.cpp
    QtApp_Posp_T2D_ME_mt.h QtApp_Posp_T2D_ME_mt.cpp
    # Model_T2D_ME_p
    # prep process module
    #QtSceneFromModel_T2D_ME_p.h QtSceneFromModel_T2D_ME_p.cpp
    #QtApp_Prep_T2D_ME_p.h QtApp_Prep_T2D_ME_p.cpp
    # post process module
    #QtSceneFromHdf5_T2D_ME_p.h QtSceneFromHdf5_T2D_ME_p.cpp
    #QtApp_Posp_T2D_ME_p.h QtApp_Posp_T2D_ME_p.cpp
    # Model_T2D_CHM_s
    # prep process module
    QtSceneFromModel_T2D_CHM_s.h QtSceneFromModel_T2D_CHM_s.cpp
    QtApp_Prep_T2D_CHM_s.h QtApp_Prep_T2D_CHM_s.cpp
    # post process module
    QtSceneFromHdf5_T2D_CHM_s.h QtSceneFromHdf5_T2D_CHM_s.cpp
    QtApp_Posp_T2D_CHM_s.h QtApp_Posp_T2D_CHM_s.cpp
    ################## 3D MPM model display ####################
    # 3D OpenGL Buffer Objects
    QtTetrahedronMeshGLObject.h QtTetrahedronMeshGLObject.cpp
    QtMonoColorBallGLObject.h QtMonoColorBallGLObject.cpp
    QtMultiColorBallGLObject.h QtMultiColorBallGLObject.cpp
    QtRigidTetrahedronMeshGLObject.h QtRigidTetrahedronMeshGLObject.cpp
    # Model_T3D_ME_s
    # prep process module
    QtSceneFromModel_T3D_ME_s.h QtSceneFromModel_T3D_ME_s.cpp
    QtApp_Prep_T3D_ME_s.h QtApp_Prep_T3D_ME_s.cpp
    # post process module
    QtSceneFromHdf5_T3D_ME_s.h QtSceneFromHdf5_T3D_ME_s.cpp
    QtApp_Posp_T3D_ME_s.h QtApp_Posp_T3D_ME_s.cpp
    # Model_T3D_CHM_s
    # prep process module
    QtSceneFromModel_T3D_CHM_s.h QtSceneFromModel_T3D_CHM_s.cpp
    QtApp_Prep_T3D_CHM_s.h QtApp_Prep_T3D_CHM_s.cpp
    # post process module
    QtSceneFromHdf5_T3D_CHM_s.h QtSceneFromHdf5_T3D_CHM_s.cpp
    QtApp_Posp_T3D_CHM_s.h QtApp_Posp_T3D_CHM_s.cpp
    )

set(MODELVIEWER_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    ModelViewer PUBLIC
    # Internal
    ${COMMONUTILS_INCLUDE_DIR}
    ${FILEUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${MATERIALMODELS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    ${SIMULATIONS_INCLUDE_DIR}
    ${SIMULATIONS_OMP_INCLUDE_DIR}
    ${MODELVIEWER_INCLUDE_DIR}
    # External
    ${EIGEN_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    ${FREETYPE_INCLUDE_DIR}
    )

target_link_libraries(
    ModelViewer PUBLIC
    # Internal dependencies
    FileUtils
    # External dependencies
    Qt5::Widgets
    Qt5::OpenGL
    ${HDF5_LIBRARIES}
    )

target_link_libraries(
    ModelViewer PUBLIC
    debug ${FREETYPE_DEBUG_LIBRARIES}
    optimized ${FREETYPE_RELEASE_LIBRARIES}
    )

if (MSVC)
    # get rid of console window
    # set_target_properties(main
    # PROPERTIES
        # WIN32_EXECUTABLE YES
        # LINK_FLAGS "/ENTRY:mainCRTStartup"
	# )
    # suppressed warning for C functions
    target_compile_definitions(ModelViewer PRIVATE _CRT_SECURE_NO_WARNINGS)
    # precompiled header
    target_compile_options(ModelViewer PRIVATE -YcModelViewer_pcp.h)
    # exclude moc files from precompilation
    set_source_files_properties(${CMAKE_BINARY_DIR}\\ModelViewer\\ModelViewer_autogen\\mocs_compilation.cpp PROPERTIES COMPILE_FLAGS /Y-)
endif()

set(CMAKE_AUTOMOC OFF)