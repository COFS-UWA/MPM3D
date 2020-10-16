#===================================
#
# MemoryUtils
#
#===================================
add_library(
    MemoryUtils STATIC
    # precompiled header
    MemoryUtils_pcp.h MemoryUtils_pcp.cpp
    # memory
    ItemArray.hpp
    ItemArrayFast.hpp
    ItemBuffer.hpp
    ItemStack.hpp
    MemoryBuffer.h
    StringBuffer.h
    PreAllocStringBuffer.hpp
    CacheAlignedMem.h
    )

set(MEMORYUTILS_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING INTERNAL FORCE)

target_include_directories(
    MemoryUtils PUBLIC
    ${MEMORYUTILS_INCLUDE_DIR}
    )

if (MSVC)
    target_compile_definitions(MemoryUtils PRIVATE _CRT_SECURE_NO_WARNINGS)
    target_compile_options(MemoryUtils PRIVATE -YcMemoryUtils_pcp.h)
endif()
