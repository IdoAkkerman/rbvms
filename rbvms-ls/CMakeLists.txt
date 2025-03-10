# This file is part of the RBVMS application.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#-------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.8)
project(RBVMS-LS LANGUAGES CXX C)

#-------------------------------------------------------------------------------
# Check MFEM_DIR
#-------------------------------------------------------------------------------
if(NOT DEFINED MFEM_DIR)
  message(FATAL_ERROR "MFEM_DIR not defined")
endif()

message(STATUS "MFEM_DIR = ${MFEM_DIR}")
if(NOT IS_DIRECTORY ${MFEM_DIR})
  message(FATAL_ERROR "MFEM_DIR is not a directory")
endif()

#-------------------------------------------------------------------------------
# Add Files
#-------------------------------------------------------------------------------
set(SOURCES
    coefficients.cpp
    evolution.cpp
    monitor.cpp
    precon.cpp
    weakform.cpp
    printInfo.cpp
    dist_solver.cpp)

set(HEADERS
    coefficients.hpp
    evolution.hpp
    monitor.hpp
    precon.hpp
    weakform.hpp
    buildInfo.hpp
    dist_solver.hpp)

#-------------------------------------------------------------------------------
# Create header
#-------------------------------------------------------------------------------
get_filename_component(tmp1 ${MFEM_DIR} REALPATH)
get_filename_component(tmp2 ${PROJECT_BINARY_DIR}/../mfem REALPATH)

if (${tmp1} STREQUAL ${tmp2})
  message(STATUS "Report RBVMS-LS, MFEM and hypre in build info.")
  add_custom_command(OUTPUT buildInfo.hpp
                     PRE_BUILD
                     COMMAND /bin/bash
                             "${PROJECT_SOURCE_DIR}/../config/buildInfo.sh"
                             "${PROJECT_SOURCE_DIR}" "mfem" "hypre"
                     COMMENT "Getting git info")
else()
  message(STATUS "External MFEM and hypre libraries used.")
  message(STATUS "Only report RBVMS-Ls in build info.")
  add_custom_command(OUTPUT buildInfo.hpp
                     PRE_BUILD
                     COMMAND /bin/bash
                             "${PROJECT_SOURCE_DIR}/../config/buildInfo.sh"
                             "${PROJECT_SOURCE_DIR}"
                     COMMENT "Getting git info")
endif()

#-------------------------------------------------------------------------------
# Add mfem
#-------------------------------------------------------------------------------
set(mfem_DIR ${MFEM_DIR})
find_package(mfem NAMES MFEM)

message(STATUS "Found mfem: in ${MFEM_DIR} (version ${MFEM_VERSION})")

# Add mfem common
add_library(mfem-common STATIC IMPORTED GLOBAL)
set_property(TARGET mfem-common PROPERTY
             IMPORTED_LOCATION "${MFEM_DIR}/miniapps/common/libmfem-common.a")

set_property(TARGET mfem-common PROPERTY
             INCLUDE_DIRECTORIES "${MFEM_DIR}/../../mfem/miniapps/common")

list(PREPEND  MFEM_LIBRARIES mfem-common)

# Add mfem common
message("MFEM_CXX_FLAGS = ${MFEM_CXX_FLAGS}")
message("MFEM_INCLUDE_DIRS = ${MFEM_INCLUDE_DIRS}")
message("MFEM_LIBRARIES = ${MFEM_LIBRARIES}")
message("MFEM_COMMON_LIBRARIES = ${MFEM_DIR}/mfem/miniapps/common")
message("MFEM_COMMON_INCLUDES= ${MFEM_DIR}/../mfem/miniapps/common")



#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
add_executable(rbvms-ls rbvms-ls.cpp ${SOURCES} ${HEADERS})

target_link_libraries(rbvms-ls PRIVATE ${MFEM_LIBRARIES}  ${CMAKE_DL_LIBS})

install(TARGETS rbvms-ls DESTINATION bin)

