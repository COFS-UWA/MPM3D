#===================================
#
# MemoryUtils
#
#===================================
add_library(
    CommonUtils STATIC
    # precompiled header
    CommonUtils_pcp.h CommonUtils_pcp.cpp
    #
    LinkList.hpp
    StackList.hpp
    RandNumGenerator.h RandNumGenerator.cpp
    ThreadBarrierFixedNum.h
    )

set(COMMONUTILS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    CommonUtils PUBLIC
    ${COMMONUTILS_INCLUDE_DIR}
    )

if (MSVC)
    target_compile_definitions(CommonUtils PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(CommonUtils PRIVATE -YcCommonUtils_pcp.h)
endif()
