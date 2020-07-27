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
    ################## 2D MPM model display ######################
    # 2D OpenGL Buffer Objects
    QtTriangleMeshGLObject.h QtTriangleMeshGLObject.cpp
    QtMonoColorCircleGLObject.h QtMonoColorCircleGLObject.cpp
    QtMultiColorCircleGLObject.h QtMultiColorCircleGLObject.cpp
    QtRigidCircleObject.h QtRigidCircleObject.cpp
    # display chars
    QtCharBitmapLoader.h QtCharBitmapLoader.cpp
    QtUniformColorMapObject.h QtUniformColorMapObject.cpp
    # Model_T2D_ME_s
    # prep process module
    QtSceneFromModel_T2D_ME_s.h QtSceneFromModel_T2D_ME_s.cpp
    QtApp_Prep_T2D_ME_s.h QtApp_Prep_T2D_ME_s.cpp
    # post process module
    QtSceneFromHdf5_T2D_ME_s.h QtSceneFromHdf5_T2D_ME_s.cpp
    QtApp_Posp_T2D_ME_s.h QtApp_Posp_T2D_ME_s.cpp
    # Model_T2D_ME_p
    # prep process module
    QtSceneFromModel_T2D_ME_p.h QtSceneFromModel_T2D_ME_p.cpp
    QtApp_Prep_T2D_ME_p.h QtApp_Prep_T2D_ME_p.cpp
    # post process module
    QtSceneFromHdf5_T2D_ME_p.h QtSceneFromHdf5_T2D_ME_p.cpp
    QtApp_Posp_T2D_ME_p.h QtApp_Posp_T2D_ME_p.cpp
    # Model_T2D_CHM_s
    # prep process module
    QtSceneFromModel_T2D_CHM_s.h QtSceneFromModel_T2D_CHM_s.cpp
    QtApp_Prep_T2D_CHM_s.h QtApp_Prep_T2D_CHM_s.cpp
    # post process module
    QtSceneFromHdf5_T2D_CHM_s.h QtSceneFromHdf5_T2D_CHM_s.cpp
    QtApp_Posp_T2D_CHM_s.h QtApp_Posp_T2D_CHM_s.cpp
    ################## 3D MPM model display ######################
    MPM3DModelView.h MPM3DModelView.cpp
    # Mesh gl buffer
    TetrahedronMeshBuffer.h TetrahedronMeshBuffer.cpp
    # Particle gl buffer - unicolor
    MonoColorCubeParticleSystem.h
    MonoColorParticleBuffer.h MonoColorParticleBuffer.cpp
    # Particle gl buffer with Phong model
    ValueToColor.h
    PhongCubeParticleSystem.h PhongCubeParticleSystem.cpp
    PhongParticleBuffer.h PhongParticleBuffer.cpp
    # Particle as balls
    ball_mesh_data.h
    BallParticleBuffer.h BallParticleBuffer.cpp
    # Preprocess single frame controller
    PrepSingleFrameControllerBase.h
    PrepSingleFrameControllerTemplate.hpp
    # Preprocess window
    PrepMPM3DApp.h PrepMPM3DApp.cpp
    # Post-process single frame controller
    PospSingleFrameController.h PospSingleFrameController.cpp
    # Animation controller
    AnimationGenerationController.h AnimationGenerationController.cpp
    # Post-process window
    PospMPM3DApp.h PospMPM3DApp.cpp
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
    ${MODELVIEWER_INCLUDE_DIR}
    # External
    ${HDF5_INCLUDE_DIR}
    ${FREETYPE_INCLUDE_DIR}
    )

target_link_libraries(
    ModelViewer
    FileUtils
    Qt5::Widgets
    Qt5::OpenGL
    ${HDF5_LIBRARIES}
    )

target_link_libraries(
    ModelViewer
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