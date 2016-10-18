#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  Copyright (c) 2013, Jed Brown 
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#

option(HAVE_SPINLOCKS "Enable if pthread lib supports spinlocks" OFF)
MARK_AS_ADVANCED(HAVE_SPINLOCKS)
if(NOT HAVE_SPINLOCKS)
  add_definitions(-DNOSPINLOCKS)
endif()

# Include the header directory
# ----------------------------
message(STATUS 
  "Prepending ${CMAKE_CURRENT_SOURCE_DIR}/include for PMRRR's headers")
include_directories(BEFORE "${CMAKE_CURRENT_SOURCE_DIR}/include")

# Define the header files installation rules
# ------------------------------------------
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
	DESTINATION ${CMAKE_INSTALL_INCLUDEDIR} 
  FILES_MATCHING PATTERN "*.h")
if(WIN32 OR APPLE OR NOT EL_HYBRID)
  # Disable for Windows because of lack of POSIX support;
  # disable for Mac OS X because of deprecation of unnamed semaphores
  # (and therefore sem_init and sem_destroy)
  option(DISABLE_PTHREADS "Disable pthreads?" ON)
else()
  option(DISABLE_PTHREADS "Disable pthreads?" OFF)
endif()
if(DISABLE_PTHREADS)
  add_definitions(-DDISABLE_PTHREADS)
endif()

# Ensure that an MPI C compiler was found
if(NOT MPI_C_FOUND)
  message(FATAL_ERROR "No MPI C compiler was found, so PMRRR cannot be built")
endif()

if(NOT DISABLE_PTHREADS)
  include(ElCheckFunctionExists)
  set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS}")
  set(CMAKE_REQUIRED_LINKER_FLAGS
    "${MPI_C_LINK_FLAGS} ${CMAKE_EXE_LINKER_FLAGS}")
  set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
  set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
  El_check_function_exists(MPI_Query_thread HAVE_MPI_QUERY_THREAD)
  if(NOT HAVE_MPI_QUERY_THREAD)
    message(FATAL_ERROR "Do not have MPI_Query_thread, so cannot build PMRRR")
  endif()
  set(CMAKE_REQUIRED_FLAGS)
  set(CMAKE_REQUIRED_LINKER_FLAGS)
  set(CMAKE_REQUIRED_INCLUDES)
  set(CMAKE_REQUIRED_LIBRARIES)
endif()
    
set(CMAKE_C_FLAGS_${UPPER_BUILD_TYPE} "${C_FLAGS}")

if(NOT DISABLE_PTHREADS AND NOT CMAKE_THREAD_LIBS_INIT)
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads)
  if(NOT CMAKE_USE_PTHREADS_INIT)
    message(FATAL_ERROR "Could not find a pthreads library, cannot build PMRRR")
  endif()
endif()
set(CMAKE_THREAD_LIBS_INIT ${CMAKE_THREAD_LIBS_INIT} PARENT_SCOPE)

# Define the main library and its link libraries
# ----------------------------------------------
file(GLOB_RECURSE PMRRR_SRC RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.c" "*.h")
add_library(pmrrr ${LIBRARY_TYPE} ${PMRRR_SRC})
if(DISABLE_PTHREADS)
  target_link_libraries(pmrrr ${MPI_C_LIBRARIES} ${MATH_LIBS})
else()
  target_link_libraries(pmrrr 
    ${MPI_C_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MATH_LIBS})
endif()
if(EL_LINK_FLAGS)
  set_target_properties(pmrrr PROPERTIES LINK_FLAGS ${EL_LINK_FLAGS})
endif()
install(TARGETS pmrrr EXPORT ElementalTargets DESTINATION ${CMAKE_INSTALL_LIBDIR})
