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


message("${SOURCES}")
message("${GENERIC_SOURCES}")
message("${HEADERS}")
message("${GENERIC_HEADERS}")

#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
add_executable(rbvms-ls rbvms-ls.cpp
               ${SOURCES} ${GENERIC_SOURCES}
               ${HEADERS} ${GENERIC_HEADERS})

target_link_libraries(rbvms-ls PRIVATE ${MFEM_LIBRARIES} ${CMAKE_DL_LIBS})
install(TARGETS rbvms-ls DESTINATION bin)
