#===================================
#
# File Operations
#
#===================================
add_library(
    FileUtils STATIC
    # precomplied header
    FileUtils_pcp.h FileUtils_pcp.cpp
    #
    file_utils.h file_utils.cpp
    GifCreator.h GifCreator.cpp
    )

set(FILEUTILS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    FileUtils PUBLIC
    ${FILEUTILS_INCLUDE_DIR}
    )

if (MSVC)
    target_compile_definitions(FileUtils PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(FileUtils PRIVATE -YcFileUtils_pcp.h)
endif()
