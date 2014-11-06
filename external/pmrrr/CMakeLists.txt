#
#  Copyright (c) 2009-2013, Jack Poulson
#                     2013, Jed Brown 
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#

option(HAVE_SPINLOCKS "Enable if pthread lib supports spinlocks" OFF)
MARK_AS_ADVANCED(HAVE_SPINLOCKS)

set(CMAKE_REQUIRED_FLAGS "${MPI_C_COMPILE_FLAGS} ${MPI_LINK_FLAGS}")
set(CMAKE_REQUIRED_INCLUDES ${MPI_C_INCLUDE_PATH})
set(CMAKE_REQUIRED_LIBRARIES ${MPI_C_LIBRARIES})
check_function_exists(MPI_Query_thread HAVE_MPI_QUERY_THREAD)
if(NOT HAVE_MPI_QUERY_THREAD)
  message(FATAL_ERROR "Do not have MPI_Query_thread, so cannot build PMRRR")
endif()

# Add in the preprocessor definitions and fortran name-mangling
set(C_BASE_FLAGS "${ELEM_BLAS_DEFS} ${ELEM_LAPACK_DEFS}")
if(NOT HAVE_SPINLOCKS)
  set(C_BASE_FLAGS "${C_BASE_FLAGS} -DNOSPINLOCKS")
endif()

# Ensure that an MPI C compiler was found
if(NOT MPI_C_FOUND)
  message(FATAL_ERROR "No MPI C compiler was found, so PMRRR cannot be built")
endif()
include_directories(${MPI_C_INCLUDE_PATH})
    
set(CMAKE_C_FLAGS_${UPPER_BUILD_TYPE} "${C_BASE_FLAGS} ${C_FLAGS}")

# Look for a PTHREADS library unless one was already specified
if(NOT CMAKE_THREAD_LIBS_INIT)
  set(CMAKE_THREAD_PREFER_PTHREAD ON)
  find_package(Threads)
  if(NOT CMAKE_USE_PTHREADS_INIT)
    message(WARNING "Could not find a pthreads library, cannot build PMRRR")
    return()
  endif()
endif()

# Grab all of the .c and .h files
file(GLOB_RECURSE PMRRR_C RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.c")
file(GLOB_RECURSE PMRRR_H RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.h")
set(PMRRR_SRC "${PMRRR_C};${PMRRR_H}")

# Define the main library and its link libraries
add_library(pmrrr ${LIBRARY_TYPE} ${PMRRR_SRC})
target_link_libraries(pmrrr 
  ${MPI_C_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT} ${MATH_LIBS})
if(MPI_LINK_FLAGS)
  set_target_properties(pmrrr PROPERTIES LINK_FLAGS ${MPI_LINK_FLAGS})
endif()
install(TARGETS pmrrr DESTINATION lib)
set(CMAKE_THREAD_LIBS_INIT ${CMAKE_THREAD_LIBS_INIT} PARENT_SCOPE)

# Define the header-preparation rules
set(PREPARED_HEADERS)
foreach(HEADER ${PMRRR_H})
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
    COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}
            ${CMAKE_CURRENT_BINARY_DIR}/${HEADER}
    DEPENDS "${CMAKE_CURRENT_SOURCE_DIR}/${HEADER}")
  list(APPEND PREPARED_HEADERS ${CMAKE_CURRENT_BINARY_DIR}/${HEADER})

  get_filename_component(HEADER_PATH ${HEADER} PATH)
  install(FILES ${CMAKE_CURRENT_BINARY_DIR}/${HEADER} 
          DESTINATION ${HEADER_PATH})
endforeach()
add_custom_target(prepare_pmrrr_headers DEPENDS ${PREPARED_HEADERS})
add_dependencies(pmrrr prepare_pmrrr_headers)

# Make sure that the PMRRR headers can be found
include_directories("${CMAKE_CURRENT_BINARY_DIR}/include")
