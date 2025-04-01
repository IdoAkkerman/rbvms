# This file is part of the RBVMS application.
#
# RBVMSis free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# Add Files
#-------------------------------------------------------------------------------
set(SOURCES
    coefficients.cpp
    solver.cpp
    printInfo.cpp)

set(HEADERS
    coefficients.hpp
    solver.hpp)

#-------------------------------------------------------------------------------
# Build library
#-------------------------------------------------------------------------------
add_library(util ${SOURCES} ${HEADERS})
target_link_libraries(util PRIVATE ${MFEM_LIBRARIES} ${CMAKE_DL_LIBS})

#-------------------------------------------------------------------------------
# Convert headers
#-------------------------------------------------------------------------------
unset(tmp_names)
foreach(name ${HEADERS})
   list(APPEND tmp_names ${CMAKE_CURRENT_SOURCE_DIR}/${name})
endforeach()
set(GENERIC_HEADERS ${tmp_names} PARENT_SCOPE)
