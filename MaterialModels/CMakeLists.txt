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
    MaterialModelType.h
    # ModelContainer
    ModelContainer.h ModelContainer.cpp
    # Linear Elasticity
    LinearElasticity.h LinearElasticity.cpp
    # Modified Cam Clay
    ModifiedCamClay.h ModifiedCamClay.cpp
    )

target_include_directories(
    MaterialModels PUBLIC
    ${MEMORYUTILS_INCLUDE_DIR}
    )
    
if (MSVC)
    target_compile_definitions(MaterialModels PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(MaterialModels PRIVATE -YcMaterialModels_pcp.h)
endif()
