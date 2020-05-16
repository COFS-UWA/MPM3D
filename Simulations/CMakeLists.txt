#===================================
#
# Simulations Project
#
#===================================
set(CMAKE_AUTOMOC ON)

add_library(
    Simulations STATIC
    # precomplied header
    Simulations_pcp.h Simulations_pcp.cpp
    macro_utils.h
    # model
    Model.h
    Step.h Step.cpp
    ModelData.h ModelData.cpp
    TimeHistory.h TimeHistory.cpp
    ResultFile.h
    ResultFile_hdf5.h ResultFile_hdf5.cpp
    ResultFile_XML.h
    Model_hdf5_utilities.h
    # boundary conditions
    BCs.h
    TimeHistory_ConsoleProgressBar.h TimeHistory_ConsoleProgressBar.cpp
    # 3D ME MPM
    Model_T3D_ME_s.h Model_T3D_ME_s.cpp
    Step_T3D_ME_s.h Step_T3D_ME_s.cpp
    ModelData_T3D_ME_s.h ModelData_T3D_ME_s.cpp
    Model_T3D_ME_s_hdf5_utilities.h Model_T3D_ME_s_hdf5_utilities.cpp
    TimeHistory_T3D_ME_s_complete.h TimeHistory_T3D_ME_s_complete.cpp
    # 3D CHM MPM
    Model_T3D_CHM_s.h Model_T3D_CHM_s.cpp
    Step_T3D_CHM_s.h Step_T3D_CHM_s.cpp
    Model_T3D_CHM_s_hdf5_utilities.h Model_T3D_CHM_s_hdf5_utilities.cpp
    # display model
    ModelWindow.h ModelWindow.cpp
    GLModelWindow.h GLModelWindow.cpp
    ModelViewer3D.h ModelViewer3D.cpp
    ModelToViewerBase.h # transfer model data to viewer class
    ModelToViewer3D.hpp
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
    Qt5::Widgets
    Qt5::OpenGL
    ${HDF5_LIBRARIES}
    )

if (MSVC)
    target_compile_definitions(Simulations PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(Simulations PRIVATE -YcSimulations_pcp.h)
    # exclude moc files from precompilation
    set_source_files_properties(${CMAKE_BINARY_DIR}\\Simulations\\Simulations_autogen\\mocs_compilation.cpp PROPERTIES COMPILE_FLAGS /Y-)
endif()

set(CMAKE_AUTOMOC OFF)
