# This file is part of the RBVMS application.
#
# RBVMS is free software; you can redistribute it and/or modify it under the
# terms of the BSD-3 license. We welcome feedback and contributions, see file
# CONTRIBUTING.md for details.


#-------------------------------------------------------------------------------
# Build executable
#-------------------------------------------------------------------------------
include_directories(../util)

add_executable(vvp_NS_mm vvp_NS_mm.cpp)
target_link_libraries(vvp_NS_mm PRIVATE ${LIBRARIES} util ${CMAKE_DL_LIBS})
install(TARGETS vvp_NS_mm DESTINATION bin)
