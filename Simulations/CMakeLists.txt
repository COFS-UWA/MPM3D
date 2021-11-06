#===================================
#
# Simulations Project
#
#===================================
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
    TimeHistory_ConsoleProgressBar.h TimeHistory_ConsoleProgressBar.cpp
    ResultFile.h
    ResultFile_hdf5.h ResultFile_hdf5.cpp
    ResultFile_XML.h
    Model_hdf5_utilities.h Model_hdf5_utilities.cpp
    # boundary conditions
    BCs.h
    # 2D ME MPM serial version
    Model_T2D_ME_s.h Model_T2D_ME_s.cpp
    Model_T2D_ME_s_hdf5_utilities.h Model_T2D_ME_s_hdf5_utilities.cpp
    ModelData_T2D_ME_s.h ModelData_T2D_ME_s.cpp
    ## deformation calculation
    Step_T2D_ME_s.h Step_T2D_ME_s.cpp Step_T2D_ME_s_avg.cpp
    TimeHistory_T2D_ME_s_complete.h TimeHistory_T2D_ME_s_complete.cpp
    Step_T2D_CHM_s_ud.h Step_T2D_CHM_s_ud.cpp Step_T2D_CHM_s_ud_avg.cpp
    TimeHistory_T2D_CHM_s_ud_complete.h TimeHistory_T2D_CHM_s_ud_complete.cpp
    ## geostatic calculation
    Step_T2D_ME_s_Geo.h Step_T2D_ME_s_Geo.cpp Step_T2D_ME_s_Geo_avg.cpp
    TimeHistory_T2D_ME_s_Geo_complete.h TimeHistory_T2D_ME_s_Geo_complete.cpp
    #TimeHistory_T2D_ME_p_complete.h TimeHistory_T2D_ME_p_complete.cpp
    # 2D CHM MPM serial version
    Model_T2D_CHM_s.h Model_T2D_CHM_s.cpp
    Model_T2D_CHM_s_hdf5_utilities.h Model_T2D_CHM_s_hdf5_utilities.cpp
    ModelData_T2D_CHM_s.h ModelData_T2D_CHM_s.cpp
    ## deformation calculation
    Step_T2D_CHM_s.h Step_T2D_CHM_s.cpp Step_T2D_CHM_s_avg.cpp
    TimeHistory_T2D_CHM_s_complete.h TimeHistory_T2D_CHM_s_complete.cpp
    ## geostatic calculation
    Step_T2D_CHM_s_Geo.h Step_T2D_CHM_s_Geo.cpp Step_T2D_CHM_s_Geo_avg.cpp
    TimeHistory_T2D_CHM_s_Geo_complete.h TimeHistory_T2D_CHM_s_Geo_complete.cpp
    # 2D CHM MPM double particle sets
    Model_T2D_CHM_d.h Model_T2D_CHM_d.cpp
    Step_T2D_CHM_d.h Step_T2D_CHM_d.cpp
    # 3D ME MPM
    Model_T3D_ME_s.h Model_T3D_ME_s.cpp
    Step_T3D_ME_s.h Step_T3D_ME_s.cpp Step_T3D_ME_s_rb.cpp
    Model_T3D_ME_s_hdf5_utilities.h Model_T3D_ME_s_hdf5_utilities.cpp
    ModelData_T3D_ME_s.h ModelData_T3D_ME_s.cpp
    TimeHistory_T3D_ME_s_complete.h TimeHistory_T3D_ME_s_complete.cpp
    # 3D CHM MPM
    Model_T3D_CHM_s.h Model_T3D_CHM_s.cpp
    Step_T3D_CHM_s.h Step_T3D_CHM_s.cpp
    Model_T3D_CHM_s_hdf5_utilities.h Model_T3D_CHM_s_hdf5_utilities.cpp
    ModelData_T3D_CHM_s.h ModelData_T3D_CHM_s.cpp
    TimeHistory_T3D_CHM_s_complete.h TimeHistory_T3D_CHM_s_complete.cpp
    # 3D ME FEM
    Model_FEM_T3D_ME_s.h Model_FEM_T3D_ME_s.cpp
    Step_FEM_T3D_ME_s.h Step_FEM_T3D_ME_s.cpp
    Model_FEM_T3D_ME_s_hdf5_utilities.h Model_FEM_T3D_ME_s_hdf5_utilities.cpp
    ModelData_FEM_T3D_ME_s.h ModelData_FEM_T3D_ME_s.cpp
    TimeHistory_FEM_T3D_ME_s_complete.h TimeHistory_FEM_T3D_ME_s_complete.cpp
    )

add_subdirectory(RigidBody/)

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
	#TBB::tbb
    ${HDF5_LIBRARIES}
    )

if (MSVC)
    target_compile_definitions(Simulations PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(Simulations PRIVATE -YcSimulations_pcp.h)
endif()
