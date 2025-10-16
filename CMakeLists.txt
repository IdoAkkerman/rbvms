#-------------------------------------------------------------------------------
# This file is part of the RBVMS application.
#
# This is the superbuild file.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.
#-------------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.8...4.99)
include(ExternalProject)
project("RBVMS Superbuild" LANGUAGES CXX C)

# Prohibit in-source builds
if ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are prohibited.")
endif ("${CMAKE_SOURCE_DIR}" STREQUAL "${CMAKE_BINARY_DIR}")

# Defaults
set(INSTALL_DIR ${PROJECT_BINARY_DIR})
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
include("${CMAKE_CURRENT_SOURCE_DIR}/config/defaults.cmake")

#-------------------------------------------------------------------------------
# Build executables when MFEM_DIR is given
#-------------------------------------------------------------------------------
if (DEFINED MFEM_DIR)

   #----------------------------------------------------------------------------
   # Load mfem
   #----------------------------------------------------------------------------
   message(STATUS "Looking for mfem ...")
   set(MFEM_DIR "" CACHE PATH "Path to the MFEM build or install prefix.")
   if (MFEM_DIR)
      find_package(mfem REQUIRED NAMES MFEM HINTS "${MFEM_DIR}"
                  "${MFEM_DIR}/lib/cmake/mfem" NO_DEFAULT_PATH)
   else()
      find_package(mfem REQUIRED NAMES MFEM)
   endif()
   message(STATUS "Found mfem config in: ${mfem_DIR} (version ${MFEM_VERSION})")
   message(STATUS "Mfem compiler is: ${MFEM_CXX_COMPILER}")
   if (NOT CMAKE_CXX_COMPILER AND MFEM_CXX_COMPILER)
     set(CMAKE_CXX_COMPILER "${MFEM_CXX_COMPILER}")
   endif()
   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MFEM_CXX_FLAGS}")
   set(LIBRARIES ${MFEM_LIBRARIES})
   set(INCLUDE_DIRS "${MFEM_INCLUDE_DIRS}")

   #----------------------------------------------------------------------------
   # Load CUDA
   #----------------------------------------------------------------------------
   if (MFEM_USE_CUDA)
      # Check for conflicts
      if (MFEM_USE_HIP)
         message(FATAL_ERROR " *** MFEM_USE_HIP cannot be combined with MFEM_USE_CUDA.")
      endif()

      # Use ${CMAKE_CXX_COMPILER} as the cuda host compiler.
      if (NOT CMAKE_CUDA_HOST_COMPILER)
         set(CMAKE_CUDA_HOST_COMPILER ${CMAKE_CXX_COMPILER})
      endif()

      # Obtain target CUDA architecture
      if (NOT CMAKE_CUDA_ARCHITECTURES)
         string(REGEX REPLACE "^sm_" "" ARCH_NUMBER "${CUDA_ARCH}")
         if ("${CUDA_ARCH}" STREQUAL "sm_${ARCH_NUMBER}")
            set(CMAKE_CUDA_ARCHITECTURES "${ARCH_NUMBER}")
         else()
            message(FATAL_ERROR "Unknown CUDA_ARCH: ${CUDA_ARCH}")
         endif()
      else()
         set(CUDA_ARCH "CMAKE_CUDA_ARCHITECTURES: ${CMAKE_CUDA_ARCHITECTURES}")
      endif()
      message(STATUS "Using CUDA architecture: ${CUDA_ARCH}")
      enable_language(CUDA)

      # Set CMAKE variables
      if (CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
         set(MFEM_CUDA_COMPILER_IS_NVCC ON)
         set(CUDA_FLAGS "${CUDA_FLAGS} --expt-extended-lambda --expt-relaxed-constexpr")
      endif()
      set(CMAKE_CUDA_STANDARD ${CMAKE_CXX_STANDARD} CACHE STRING
         "CUDA standard to use.")
      set(CMAKE_CUDA_STANDARD_REQUIRED ON CACHE BOOL
         "Force the use of the chosen CUDA standard.")
      set(CMAKE_CUDA_EXTENSIONS OFF CACHE BOOL "Enable CUDA standard extensions.")
      set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CUDA_FLAGS}")
      find_package(CUDAToolkit REQUIRED)

      # Set cusparse library
      set(CUSPARSE_LIBRARIES ${CUDA_cusparse_LIBRARY})
      if (NOT CUSPARSE_LIBRARIES)
         message(FATAL_ERROR
         "Error setting CUSPARSE_LIBRARIES: ${CUSPARSE_LIBRARIES}")
      endif()
      set(CUSPARSE_FOUND TRUE)

      # Set cublas library
      set(CUBLAS_LIBRARIES ${CUDA_cublas_LIBRARY})
      if (NOT CUBLAS_LIBRARIES)
         message(FATAL_ERROR
         "Error setting CUBLAS_LIBRARIES: ${CUBLAS_LIBRARIES}")
      endif() 
      set(CUBLAS_FOUND TRUE)
   endif()

   #----------------------------------------------------------------------------
   # Add directories
   #----------------------------------------------------------------------------
   add_subdirectory(util)
   add_subdirectory(rbvms)
   add_subdirectory(rbvms-ls)
   add_subdirectory(supg)

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
