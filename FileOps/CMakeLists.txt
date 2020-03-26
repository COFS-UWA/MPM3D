#===================================
#
# File Operations
#
#===================================
add_library(
    FileOps STATIC
    # precomplied header
    FileOps_pcp.h FileOps_pcp.cpp
    # common utilties
    common_utilities.h common_utilities.cpp
    )

if (MSVC)
    target_compile_definitions(FileOps PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(FileOps PRIVATE -YcFileOps_pcp.h)
endif()
