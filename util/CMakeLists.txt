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
    solver.hpp
    buildInfo.hpp)

#-------------------------------------------------------------------------------
# Get buildInfo
#-------------------------------------------------------------------------------
get_filename_component(tmp1 ${MFEM_DIR} REALPATH)
get_filename_component(tmp2 ${PROJECT_BINARY_DIR}/.. REALPATH)

if (${tmp1} STREQUAL ${tmp2})
   message(STATUS "Report RBVMS, MFEM and hypre in build info.")
   add_custom_command(OUTPUT buildInfo.hpp
                      PRE_BUILD
                      COMMAND /bin/bash
                              "${PROJECT_SOURCE_DIR}/config/buildInfo.sh"
                              "${PROJECT_SOURCE_DIR}/util" "mfem" "hypre"
                      COMMENT "Getting git info")
else()
   message(STATUS "External MFEM and hypre libraries used.")
   message(STATUS "Only report RBVMS in build info.")
   add_custom_command(OUTPUT buildInfo.hpp
                      PRE_BUILD
                      COMMAND /bin/bash
                              "${PROJECT_SOURCE_DIR}/config/buildInfo.sh"
                              "${PROJECT_SOURCE_DIR}/util"
                      COMMENT "Getting git info")
endif()

#-------------------------------------------------------------------------------
# Build library
#-------------------------------------------------------------------------------
add_library(util ${SOURCES} ${HEADERS})
target_link_libraries(util PRIVATE ${LIBRARIES} ${CMAKE_DL_LIBS})
