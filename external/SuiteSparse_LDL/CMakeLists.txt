#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#

# TODO: Add support for 'long' version of SuiteSparse LDL?

# Include the header directory
# ----------------------------
message(STATUS 
  "Appending ${CMAKE_CURRENT_SOURCE_DIR} for SuiteSparse_LDL's headers")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}")

# Define the header files installation rules
# ------------------------------------------
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
  FILES_MATCHING PATTERN "*.h")
    
set(CMAKE_C_FLAGS_${UPPER_BUILD_TYPE} "${C_FLAGS}")

# Define the main library and its link libraries
# ----------------------------------------------
file(GLOB_RECURSE SUITESPARSE_LDL_SRC RELATIVE 
  ${CMAKE_CURRENT_SOURCE_DIR} "*.c" "*.h")
add_library(SuiteSparse_LDL ${LIBRARY_TYPE} ${SUITESPARSE_LDL_SRC})
if(EL_LINK_FLAGS)
  set_target_properties(SuiteSparse_LDL PROPERTIES LINK_FLAGS ${EL_LINK_FLAGS})
endif()
install(TARGETS SuiteSparse_LDL EXPORT ElTargets DESTINATION lib)
