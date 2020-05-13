#===================================
#
# ResultViewerQt Project
#
#===================================
set(CMAKE_AUTOMOC ON)

add_executable(
    ResultViewerQt
    # Entry
    ResultViewerQt_main.cpp
    # Precomplied header
    ResultViewerQt_pcp.h ResultViewerQt_pcp.cpp
    hdf5_utils.h hdf5_utils.cpp
    MainWindow.h MainWindow.cpp
    GLBufferObject.h GLBufferObject.cpp
    GLWindow.h GLWindow.cpp
    )

set(RESULTVIEWERQT_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    ResultViewerQt PUBLIC
    # Internal
    ${MEMORYUTILS_INCLUDE_DIR}
    ${MATERIALMODELS_INCLUDE_DIR}
    ${GEOMETRY_INCLUDE_DIR}
    ${RESULTVIEWERQT_INCLUDE_DIR}
    # External
    ${HDF5_INCLUDE_DIR}
    )

target_link_libraries(
    ResultViewerQt
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
    target_compile_definitions(ResultViewerQt PRIVATE _CRT_SECURE_NO_WARNINGS)
    # precompiled header
    target_compile_options(ResultViewerQt PRIVATE -YcResultViewerQt_pcp.h)
    # exclude moc files from precompilation
    set_source_files_properties(${CMAKE_BINARY_DIR}\\ResultViewerQt\\ResultViewerQt_autogen\\mocs_compilation.cpp PROPERTIES COMPILE_FLAGS /Y-)
endif()

set(CMAKE_AUTOMOC OFF)