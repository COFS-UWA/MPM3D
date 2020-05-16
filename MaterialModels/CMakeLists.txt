#===================================
#
# MaterialModels Project
#
#===================================
add_library(
    MaterialModels STATIC
    # p
    MaterialModels_pcp.h MaterialModels_pcp.cpp
    # Base class
    MaterialModel.h
    # ModelContainer
    MatModelContainer.h MatModelContainer.cpp
    # Linear Elasticity
    LinearElasticity.h LinearElasticity.cpp
    # Modified Cam Clay
    ModifiedCamClay.h ModifiedCamClay.cpp
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
