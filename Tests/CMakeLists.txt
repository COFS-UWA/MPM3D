#===================================
#
# Tests Project
#
#===================================
add_executable(
    Tests Tests_main.cpp
    # precompiled header
    Tests_pcp.h Tests_pcp.cpp
    # unit tests file
    )

target_include_directories(
    Tests
    PRIVATE
    ${SIMULATIONCORE_INCLUDE_DIR}
    ${CONSTITUTIVEMODEL_INCLUDE_DIR}
    ${CMAKE_SOURCE_DIR}/PostProcessor
    ${UTILITIES_INCLUDE_DIR}
    ${EIGEN_INCLUDE_DIR}
    ${OPENGL_INCLUDE_DIR}
    ${HDF5_INCLUDE_DIR}
    )

link_directories(${GLFW_LIBRARIES_DIR})
target_link_libraries(
    Tests PRIVATE
    SimulationCore
    ConstitutiveModel
    PostProcessor
    Utilities
    ${GLFW_LIBRARIES}
    )

if (WIN32)
    # supress warning for certain functions in visual studio
    target_compile_definitions(Tests PRIVATE _CRT_SECURE_NO_WARNINGS)
    # use precompiled header
    target_compile_options(Tests PRIVATE -YcTests_pcp.h)
endif()
