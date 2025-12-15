# This file is part of the RBVMS application.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

#-------------------------------------------------------------------------------
# Add Files
#-------------------------------------------------------------------------------
set(SOURCES
    supg.cpp
    integrator.cpp)

set(HEADERS
    integrator.hpp)

#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
include_directories(../util)

if (MFEM_USE_CUDA)
   enable_language(CUDA)
   set_source_files_properties(${SOURCES} PROPERTIES LANGUAGE CUDA)
endif()

add_executable(supg ${SOURCES} ${HEADERS})
target_link_libraries(supg PRIVATE ${LIBRARIES} util ${CMAKE_DL_LIBS})
install(TARGETS supg DESTINATION bin)

# Small test utility to check tau variance under vertex-order permutations
add_executable(tau_variance tau_variance.cpp)
target_link_libraries(tau_variance PRIVATE ${LIBRARIES} util ${CMAKE_DL_LIBS})
install(TARGETS tau_variance DESTINATION bin)

# Alternate variant with direction sweeps and separate convective/diffusive checks
add_executable(tau_variance_alt tau_variance_alt.cpp)
target_link_libraries(tau_variance_alt PRIVATE ${LIBRARIES} util ${CMAKE_DL_LIBS})
install(TARGETS tau_variance_alt DESTINATION bin)
