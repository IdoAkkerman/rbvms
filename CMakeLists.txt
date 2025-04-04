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

#-------------------------------------------------------------------------------
# Build executables when MFEM_DIR is given
#-------------------------------------------------------------------------------
if (DEFINED MFEM_DIR)
   #----------------------------------------------------------------------------
   # Check MFEM_DIR
   #----------------------------------------------------------------------------
   if(NOT IS_DIRECTORY ${MFEM_DIR})
      message(FATAL_ERROR "MFEM_DIR is not a directory")
   endif()

   #----------------------------------------------------------------------------
   # Add mfem
   #----------------------------------------------------------------------------
   message(STATUS "Use provided MFEM_DIR: ${MFEM_DIR} ")
   set(mfem_DIR ${MFEM_DIR}/mfem)
   find_package(mfem NAMES MFEM)
   message(STATUS "Found mfem: in ${MFEM_DIR} (version ${MFEM_VERSION})")
   set(LIBRARIES ${MFEM_LIBRARIES})

   #----------------------------------------------------------------------------
   # Add
   #----------------------------------------------------------------------------
   add_subdirectory(util)
   add_subdirectory(rbvms)
   add_subdirectory(rbvms-ls)
   add_subdirectory(rbvms-ls2)

#-------------------------------------------------------------------------------
# If no MFEM_DIR is given --> build a local version here
#-------------------------------------------------------------------------------
else()
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
   # Build rbvms with local MFEM
   #----------------------------------------------------------------------------
   ExternalProject_Add(rbvms
                       PREFIX ${INSTALL_DIR}
                       DEPENDS mfem
                       CMAKE_ARGS -DMFEM_DIR=${INSTALL_DIR}
                                  -DCMAKE_INSTALL_PREFIX=${INSTALL_DIR}
                                  -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
                                  -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}
                       SOURCE_DIR ${PROJECT_SOURCE_DIR}
                       BINARY_DIR ${PROJECT_BINARY_DIR})
endif()
