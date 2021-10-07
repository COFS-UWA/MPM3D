#===================================
#
# SimulationsOMP Project
#
#===================================
add_library(
    SimulationsOMP STATIC
    # precomplied header
    SimulationsOMP_pcp.h SimulationsOMP_pcp.cpp
    SimulationsOMPUtils.h
    ParticleVariablesGetter.h ParticleVariablesGetter.cpp
    # OpenMP multi-thread version
    Step_OMP.h Step_OMP.cpp
    # 2D ME MPM
    Model_T2D_ME_mt.h Model_T2D_ME_mt.cpp
    Model_T2D_ME_mt_hdf5_utilities.h Model_T2D_ME_mt_hdf5_utilities.cpp
    ModelData_T2D_ME_mt.h ModelData_T2D_ME_mt.cpp
    Step_T2D_ME_mt.h Step_T2D_ME_mt.cpp
    TimeHistory_T2D_ME_mt_complete.h TimeHistory_T2D_ME_mt_complete.cpp
    ## Rectangle elements
    Model_R2D_ME_mt.h Model_R2D_ME_mt.cpp
    Step_R2D_ME_mt.h Step_R2D_ME_mt.cpp
    # 3D ME MPM
    Model_T3D_ME_mt.h Model_T3D_ME_mt.cpp
    Model_T3D_ME_mt_hdf5_utilities.h Model_T3D_ME_mt_hdf5_utilities.cpp
    ModelData_T3D_ME_mt.h ModelData_T3D_ME_mt.cpp
    Step_T3D_ME_mt.h Step_T3D_ME_mt.cpp
    TimeHistory_T3D_ME_mt_complete.h TimeHistory_T3D_ME_mt_complete.cpp
    Step_T3D_ME_mt_Geo.h Step_T3D_ME_mt_Geo.cpp
    TimeHistory_T3D_ME_mt_Geo_complete.h TimeHistory_T3D_ME_mt_Geo_complete.cpp
    TimeHistory_T3D_ME_mt_Geo_ratio.h TimeHistory_T3D_ME_mt_Geo_ratio.cpp
    # 2D CHM MPM
    Model_T2D_CHM_mt.h Model_T2D_CHM_mt.cpp
    Step_T2D_CHM_mt.h Step_T2D_CHM_mt.cpp Step_T2D_CHM_mt2.cpp
    Model_T2D_CHM_mt_hdf5_utilities.h Model_T2D_CHM_mt_hdf5_utilities.cpp
    ModelData_T2D_CHM_mt.h ModelData_T2D_CHM_mt.cpp
    TimeHistory_T2D_CHM_mt_complete.h TimeHistory_T2D_CHM_mt_complete.cpp
    # 3D CHM MPM
    Model_T3D_CHM_mt.h Model_T3D_CHM_mt.cpp
    Model_T3D_CHM_mt_hdf5_utilities.h Model_T3D_CHM_mt_hdf5_utilities.cpp
    ModelData_T3D_CHM_mt.h ModelData_T3D_CHM_mt.cpp
    Step_T3D_CHM_mt.h Step_T3D_CHM_mt.cpp Step_T3D_CHM_mt2.cpp
    Step_T3D_CHM_ud_mt.h Step_T3D_CHM_ud_mt.cpp
    TimeHistory_T3D_CHM_mt_complete.h TimeHistory_T3D_CHM_mt_complete.cpp
    Step_T3D_CHM_mt_Geo.h Step_T3D_CHM_mt_Geo.cpp
    TimeHistory_T3D_CHM_mt_Geo_complete.h TimeHistory_T3D_CHM_mt_Geo_complete.cpp
    TimeHistory_T3D_CHM_mt_Geo_ratio.h TimeHistory_T3D_CHM_mt_Geo_ratio.cpp
    Step_T3D_CHM_ud_mt_subiter.h Step_T3D_CHM_ud_mt_subiter.cpp
    Step_T3D_CHM_ud_mt_subiter_utils.cpp
    TimeHistory_T3D_CHM_ud_mt_subiter_complete.h TimeHistory_T3D_CHM_ud_mt_subiter_complete.cpp
    TimeHistory_T3D_CHM_ud_mt_subiter_ratio.h TimeHistory_T3D_CHM_ud_mt_subiter_ratio.cpp
    # Parallelism using TBB
    DivideTask.hpp
    MergeTask.hpp
    ParallelUtils.h
    MSDRadixSortUtils.h MSDRadixSortUtils.cpp
    MSDRadixSortTask.h MSDRadixSortTask.cpp
    SortTriMeshNodeTask.h SortTriMeshNodeTask.cpp
    SortTehMeshNodeTask.h SortTehMeshNodeTask.cpp
    SortParticleTask.h SortParticleTask.cpp
    Step_TBB.h Step_TBB.cpp
    #
    ParaUtil.h ParaUtil.cpp
    ParallelForTask.hpp ParallelReduceTask.hpp
    DataMem.h DataMem.cpp
    ParaUtilSerialSort.h ParaUtilSerialSort.cpp
    ParaUtilMSDRadix16Sort.h ParaUtilMSDRadix16Sort.cpp
    NodeElemSortInternal.hpp NodeElemSort.hpp
    PclSort.h PclSort.cpp
    # 2D ME MPM
    Step_T2D_ME_Task.h Step_T2D_ME_Task.cpp
    Step_T2D_ME_TBB.h Step_T2D_ME_TBB.cpp
    TimeHistory_T2D_ME_TBB_complete.h TimeHistory_T2D_ME_TBB_complete.cpp
    # 2D CHM MPM
    Step_T2D_CHM_Task.h Step_T2D_CHM_Task.cpp
    Step_T2D_CHM_TBB.h Step_T2D_CHM_TBB.cpp
    TimeHistory_T2D_CHM_TBB_complete.h TimeHistory_T2D_CHM_TBB_complete.cpp
    # 3D ME MPM
    Step_T3D_ME_TBB_Task.h Step_T3D_ME_TBB_Task.cpp
    Step_T3D_ME_TBB.h Step_T3D_ME_TBB.cpp
    TimeHistory_T3D_ME_TBB_complete.h TimeHistory_T3D_ME_TBB_complete.cpp
    # 3D CHM MPM
    Step_T3D_CHM_TBB_Task.h Step_T3D_CHM_TBB_Task.cpp
    Step_T3D_CHM_TBB.h Step_T3D_CHM_TBB.cpp
    TimeHistory_T3D_CHM_TBB_complete.h TimeHistory_T3D_CHM_TBB_complete.cpp
    )

add_subdirectory(RigidObject/)

set(SIMULATIONS_OMP_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    SimulationsOMP PUBLIC
    # Internal dependencies
    ${COMMONUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    ${FILEUTILS_INCLUDE_DIR}
    ${MATERIALMODELS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    ${SIMULATIONS_INCLUDE_DIR}
    ${SIMULATIONS_OMP_INCLUDE_DIR}
    # External dependencies
    ${EIGEN_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
<<<<<<< HEAD
	${TBB_INCLUDE_DIR}
=======
    ${TBB_INCLUDE_DIR}
>>>>>>> 012e995fcd67451d8607327eb2ee9b3d588a1610
    )

find_package(TBB CONFIG REQUIRED)

target_link_libraries(
    SimulationsOMP PUBLIC
    # Internal dependencies
    MemoryUtils
    MaterialModels
    Geometry
    Simulations
    # External dependencies
    ${HDF5_LIBRARIES}
    TBB::tbb
    )

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    target_link_libraries(SimulationsOMP PUBLIC OpenMP::OpenMP_CXX)
endif()

if (MSVC)
    target_compile_definitions(SimulationsOMP PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(SimulationsOMP PRIVATE -YcSimulationsOMP_pcp.h)
endif()
