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
    )
