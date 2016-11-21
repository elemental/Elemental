#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#

# Include the header directory
# ----------------------------
message(STATUS 
  "Prepending ${CMAKE_CURRENT_SOURCE_DIR}/include for ElSuiteSparse headers")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Define the header files installation rules
# ------------------------------------------
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
  FILES_MATCHING PATTERN "*.h" PATTERN "*.hpp")
    
set(CMAKE_C_FLAGS_${UPPER_BUILD_TYPE} "${C_FLAGS}")

# Define the main library and its link libraries
# ----------------------------------------------
file(GLOB_RECURSE EL_SUITESPARSE_SRC RELATIVE 
  ${CMAKE_CURRENT_SOURCE_DIR} "*.c" "*.cpp" "*.h" "*.hpp")
add_library(ElSuiteSparse ${LIBRARY_TYPE} ${EL_SUITESPARSE_SRC})
if (LINUX)
   target_link_libraries(ElSuiteSparse m)
endif()

if(EL_LINK_FLAGS)
  set_target_properties(ElSuiteSparse PROPERTIES LINK_FLAGS ${EL_LINK_FLAGS})
  set_target_properties(ElSuiteSparse PROPERTIES VERSION ${EL_VERSION_MINOR} SOVERSION ${EL_VERSION_MAJOR})
endif()
install(TARGETS ElSuiteSparse EXPORT ElementalTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
