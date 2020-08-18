#===================================
#
# MaterialModels Project
#
#===================================
add_library(
    MaterialModels STATIC
    # precompiled header
    MaterialModels_pcp.h MaterialModels_pcp.cpp
    # Base class
    MaterialModel.h
    # Linear Elasticity
    LinearElasticity.h LinearElasticity.cpp
    # Modified Cam Clay
    ModifiedCamClay.h ModifiedCamClay.cpp
    # Modified Cam Clay for completely undrained deformation
    UndrainedModifiedCamClay.h UndrainedModifiedCamClay.cpp
    # ModelContainer
    MatModelContainer.h MatModelContainer.cpp
    # Map from material model id to pointer
    MatModelIdToPointerMap.h MatModelIdToPointerMap.cpp
    )

set(MATERIALMODELS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    MaterialModels PUBLIC
    ${MATERIALMODELS_INCLUDE_DIR}
    ${COMMONUTILS_INCLUDE_DIR}
    ${MEMORYUTILS_INCLUDE_DIR}
    )
    
if (MSVC)
    target_compile_definitions(MaterialModels PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(MaterialModels PRIVATE -YcMaterialModels_pcp.h)
endif()
