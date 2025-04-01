#-------------------------------------------------------------------------------
# This file is part of the RBVMS application.
#
# This is the superbuild file.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#-------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.8)
include(ExternalProject)
project("RBVMS Superbuild" LANGUAGES CXX C)

set(INSTALL_DIR ${PROJECT_BINARY_DIR})

#----------------------------------------------------------------------------
# IF MFEM_DIR is not given build a local version here
#----------------------------------------------------------------------------
if (NOT DEFINED MFEM_DIR)
   #----------------------------------------------------------------------------
   # Get git submodules
   #----------------------------------------------------------------------------
   find_package(Git REQUIRED)
   if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
      option(GIT_SUBMODULE "Check submodules during build" ON)
      if(GIT_SUBMODULE)
         message(STATUS "Submodule update")
         execute_process(COMMAND ${GIT_EXECUTABLE}
                         submodule update --init --recursive --single-branch
                         WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                         RESULT_VARIABLE GIT_SUBMOD_RESULT)
         if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive "
                    "failed with ${GIT_SUBMOD_RESULT}, "
                    "please checkout submodules")
         endif()
      endif()
   endif()

   #----------------------------------------------------------------------------
   # Build hypre
   #----------------------------------------------------------------------------
   ExternalProject_Add(hypre
                       PREFIX ${INSTALL_DIR}
                       CONFIGURE_COMMAND ${CMAKE_COMMAND}
                                         "${PROJECT_SOURCE_DIR}/hypre/src"
                                         "-DHYPRE_INSTALL_PREFIX=${INSTALL_DIR}"
                                         "-DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}"
                       SOURCE_DIR ${PROJECT_SOURCE_DIR}/hypre
                       BINARY_DIR ${PROJECT_BINARY_DIR}/hypre)

   #----------------------------------------------------------------------------
   # Build mfem
   #----------------------------------------------------------------------------
   ExternalProject_Add(mfem
                       PREFIX ${INSTALL_DIR}
                       DEPENDS hypre
                       CMAKE_ARGS -DMFEM_USE_MPI=YES
                                  -DMFEM_USE_LAPACK=YES
                                  -DHYPRE_DIR=${INSTALL_DIR}
                                  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
                                  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                  -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                       SOURCE_DIR ${PROJECT_SOURCE_DIR}/mfem
                       BINARY_DIR ${PROJECT_BINARY_DIR}/mfem)

   #----------------------------------------------------------------------------
   # Set mfem directory
   #----------------------------------------------------------------------------
   set(mfem_DIR ${MFEM_DIR}/mfem)
else()
   message(STATUS "Use provided MFEM_DIR")
   set(mfem_DIR ${MFEM_DIR}/mfem)
endif()

#-------------------------------------------------------------------------------
# Add mfem
#-------------------------------------------------------------------------------
find_package(mfem NAMES MFEM)
message(STATUS "Found mfem: in ${MFEM_DIR} (version ${MFEM_VERSION})")

#-------------------------------------------------------------------------------
# Define Function BuildInfo 
#-------------------------------------------------------------------------------
function(buildInfo)
   get_filename_component(tmp1 ${MFEM_DIR} REALPATH)
   get_filename_component(tmp2 ${PROJECT_BINARY_DIR}/.. REALPATH)

   if (${tmp1} STREQUAL ${tmp2})
      message(STATUS "Report RBVMS, MFEM and hypre in build info.")
      add_custom_command(OUTPUT buildInfo.hpp
                         PRE_BUILD
                         COMMAND /bin/bash
                                 "${PROJECT_SOURCE_DIR}/config/buildInfo.sh"
                                 "${PROJECT_SOURCE_DIR}" "mfem" "hypre"
                         COMMENT "Getting git info")
   else()
      message(STATUS "External MFEM and hypre libraries used.")
      message(STATUS "Only report RBVMS in build info.")
      add_custom_command(OUTPUT buildInfo.hpp
                         PRE_BUILD
                         COMMAND /bin/bash
                                 "${PROJECT_SOURCE_DIR}/config/buildInfo.sh"
                                 "${PROJECT_SOURCE_DIR}"
                         COMMENT "Getting git info")
   endif()
endfunction()

#-------------------------------------------------------------------------------
# Add directories
#-------------------------------------------------------------------------------
add_subdirectory(util)
add_subdirectory(rbvms)
add_subdirectory(rbvms-ls)

