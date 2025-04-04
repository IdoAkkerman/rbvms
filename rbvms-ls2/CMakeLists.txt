# This file is part of the RBVMS application.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.

#-------------------------------------------------------------------------------
# Add Files
#-------------------------------------------------------------------------------
set(SOURCES
    formulation.cpp
    integrator.cpp
    dist_solver.cpp)

set(HEADERS
    formulation.hpp
    integrator.hpp
    dist_solver.hpp)

#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
include_directories(../util)
add_executable(rbvms-ls2 rbvms-ls2.cpp ${SOURCES} ${HEADERS})
target_link_libraries(rbvms-ls2 PRIVATE ${LIBRARIES} util ${CMAKE_DL_LIBS})
install(TARGETS rbvms-ls2 DESTINATION bin)
