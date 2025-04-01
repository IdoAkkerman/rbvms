# This file is part of the RBVMS application.
#
# RBVMSis free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

#-------------------------------------------------------------------------------
# Add Files
#-------------------------------------------------------------------------------
set(SOURCES
    formulation.cpp
    integrator.cpp)

set(HEADERS
    formulation.hpp
    integrator.hpp
    buildInfo.hpp)

#-------------------------------------------------------------------------------
# Create info header
#-------------------------------------------------------------------------------
buildInfo()

#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
add_executable(rbvms rbvms.cpp
               ${SOURCES}
               ${HEADERS} ${GENERIC_HEADERS})

target_link_libraries(rbvms PRIVATE ${MFEM_LIBRARIES} util ${CMAKE_DL_LIBS})
install(TARGETS rbvms DESTINATION bin)
