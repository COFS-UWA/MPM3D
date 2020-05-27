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
    # utils
    hdf5_utils.h hdf5_utils.cpp
    # 3D MPM model display
    # Mesh gl buffer
    TetrahedronMeshBuffer.h TetrahedronMeshBuffer.cpp
    # Particle gl buffer
    MonoColorCubeParticleSystem.h
    MonoColorParticleBuffer.h MonoColorParticleBuffer.cpp
    MultiColorCubeParticleSystem.h ValueToColor.h
    MultiColorParticleBuffer.h MultiColorParticleBuffer.cpp
    # Point gl buffer
    PointBuffer.h PointBuffer.cpp
    MPM3DModelView.h MPM3DModelView.cpp
    PrepMPM3DApp.h PrepMPM3DApp.cpp
    TimeHistory_ModelView.h TimeHistory_ModelView.cpp
    # Helper of MPM3DModelView: input data to ModelView
    TimeHistory_ModelView_T3D_ME_s.h TimeHistory_ModelView_T3D_ME_s.cpp
    # Single frame
    PospSingleFrameController.h PospSingleFrameController.cpp
    # Animation
    AnimationGenerationController.h AnimationGenerationController.cpp
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
    )

target_link_libraries(
    ModelViewer
    FileUtils
    Qt5::Widgets
    Qt5::OpenGL
    ${HDF5_LIBRARIES}
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