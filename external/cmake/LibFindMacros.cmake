# Version 1.0 (2013-04-12)
# Public Domain, originally written by Lasse Karkkainen <tronic@zi.fi>
# Published at http://www.cmake.org/Wiki/CMake:How_To_Find_Libraries

# If you improve the script, please modify the forementioned wiki page because
# I no longer maintain my scripts (hosted as static files at zi.fi). Feel free
# to remove this entire header if you use real version control instead.

# Changelog:
# 2013-04-12  Added version number (1.0) and this header, no other changes
# 2009-10-08  Originally published

# Works the same as find_package, but forwards the "REQUIRED" and "QUIET" arguments
# used for the current package. For this to work, the first parameter must be the
# prefix of the current package, then the prefix of the new package etc, which are
# passed to find_package.
macro (libfind_package PREFIX)
  set (LIBFIND_PACKAGE_ARGS ${ARGN})
  if (${PREFIX}_FIND_QUIETLY)
    set (LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} QUIET)
  endif (${PREFIX}_FIND_QUIETLY)
  if (${PREFIX}_FIND_REQUIRED)
    set (LIBFIND_PACKAGE_ARGS ${LIBFIND_PACKAGE_ARGS} REQUIRED)
  endif (${PREFIX}_FIND_REQUIRED)
  find_package(${LIBFIND_PACKAGE_ARGS})
endmacro (libfind_package)

# CMake developers made the UsePkgConfig system deprecated in the same release (2.6)
# where they added pkg_check_modules. Consequently I need to support both in my scripts
# to avoid those deprecated warnings. Here's a helper that does just that.
# Works identically to pkg_check_modules, except that no checks are needed prior to use.
macro (libfind_pkg_check_modules PREFIX PKGNAME)
  if (${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} EQUAL 4)
    include(UsePkgConfig)
    pkgconfig(${PKGNAME} ${PREFIX}_INCLUDE_DIRS ${PREFIX}_LIBRARY_DIRS ${PREFIX}_LDFLAGS ${PREFIX}_CFLAGS)
  else (${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} EQUAL 4)
    find_package(PkgConfig)
    if (PKG_CONFIG_FOUND)
      pkg_check_modules(${PREFIX} ${PKGNAME})
    endif (PKG_CONFIG_FOUND)
  endif (${CMAKE_MAJOR_VERSION} EQUAL 2 AND ${CMAKE_MINOR_VERSION} EQUAL 4)
endmacro (libfind_pkg_check_modules)

# This macro searchs the _lib_list list for _lib, and, if found, adds any 
# missing dependencies in the _deps list. Also, the ${PREFIX}_FIND_REQUIRED_${_dep}
# variable is set to true if the ${PREFIX}_FIND_REQUIRED_${_lib} is true
macro(libfind_add_dep PREFIX _lib_list _lib _deps)
  
  list(FIND ${_lib_list} ${_lib} _lib_find)
  if(NOT _lib_find EQUAL -1)
    foreach(_dep ${_deps})
      # Add missing dependencies to the list.
      list(FIND ${_lib_list} ${_dep} _dep_find)
      if(_dep_find EQUAL -1)
        list(APPEND ${_lib_list} ${_dep})
      endif()

      # Set the find required flag for the dependency 
      if(${PREFIX}_FIND_REQUIRED_${_lib})
        set(${PREFIX}_FIND_REQUIRED_${_dep} TRUE)
      else()
        set(${PREFIX}_FIND_REQUIRED_${_dep} FALSE)
      endif()
    endforeach()
  endif()

endmacro()

# Do the final processing once the paths have been detected.
# If include dirs are needed, ${PREFIX}_PROCESS_INCLUDES should be set to contain
# all the variables, each of which contain one include directory.
# Ditto for ${PREFIX}_PROCESS_LIBS and library files.
# Will set ${PREFIX}_FOUND, ${PREFIX}_INCLUDE_DIRS and ${PREFIX}_LIBRARIES.
# Also handles errors in case library detection was required, etc.
macro (libfind_process PREFIX)
  # Skip processing if already processed during this run
  if (NOT ${PREFIX}_FOUND)
    # Start with the assumption that the library was found
    set (${PREFIX}_FOUND TRUE)

    # Process all includes and set _FOUND to false if any are missing
    foreach (i ${${PREFIX}_PROCESS_INCLUDES})
      if (${i})
        set (${PREFIX}_INCLUDE_DIRS ${${PREFIX}_INCLUDE_DIRS} ${${i}})
        mark_as_advanced(${i})
      else (${i})
        set (${PREFIX}_FOUND FALSE)
      endif (${i})
    endforeach (i)

    # Process all libraries and set _FOUND to false if any are missing
    foreach (i ${${PREFIX}_PROCESS_LIBS})
      if (${i})
        set (${PREFIX}_LIBRARIES ${${PREFIX}_LIBRARIES} ${${i}})
        mark_as_advanced(${i})
      else (${i})
        set (${PREFIX}_FOUND FALSE)
      endif (${i})
    endforeach (i)

    # Print message and/or exit on fatal error
    if (${PREFIX}_FOUND)
      if (NOT ${PREFIX}_FIND_QUIETLY)
        message (STATUS "Found ${PREFIX} ${${PREFIX}_VERSION}")
      endif (NOT ${PREFIX}_FIND_QUIETLY)
    else (${PREFIX}_FOUND)
      if (${PREFIX}_FIND_REQUIRED)
        foreach (i ${${PREFIX}_PROCESS_INCLUDES} ${${PREFIX}_PROCESS_LIBS})
          message("${i}=${${i}}")
        endforeach (i)
        message (FATAL_ERROR "Required library ${PREFIX} NOT FOUND.\n"
                             "Install the library (dev version) and try again. "
                             "If the library is already installed, "
                             "use ccmake to set the missing variables manually.")
      endif (${PREFIX}_FIND_REQUIRED)
    endif (${PREFIX}_FOUND)
  endif (NOT ${PREFIX}_FOUND)
endmacro (libfind_process)

macro(libfind_header PREFIX _var _header)
  set(${PREFIX}_INCLUDE_SEARCH_DIR)
  if(${PREFIX}_ROOT_DIR)
    set(${PREFIX}_INCLUDE_SEARCH_DIR "${${PREFIX}_ROOT_DIR}/include")
  endif()

  find_path(${_var} ${_header}
      HINTS ${${PREFIX}_INCLUDE_DIR} ${${PREFIX}_INCLUDE_SEARCH_DIR}
      NO_CMAKE_SYSTEM_PATH)
endmacro()

macro(libfind_library PREFIX _name)
  if(NOT ${PREFIX}_${_name}_LIBRARY)
    if(${PREFIX}_LIBRARY)
      # Search the user provided libraries for _name
      foreach(_lib ${${PREFIX}_LIBRARY})
        get_filename_component(_lib_name ${_lib} NAME)
        string(FIND ${_lib_name} ${_name} _lib_found)
    
        if(NOT _lib_found EQUAL -1)
          # Set the component library list
          set(${PREFIX}_${_name}_LIBRARY ${_lib})
          break()
        endif()
      endforeach()
      
    else()
      set(${PREFIX}_LIB_SERACH_DIRS)
      if(${PREFIX}_ROOT_DIR)
        set(${PREFIX}_LIB_SERACH_DIRS "${${PREFIX}_ROOT_DIR}/lib")
      endif()

      # Search for the library
      find_library(${PREFIX}_${_name}_LIBRARY ${_name}
          HINTS ${${PREFIX}_PKGCONF_LIBRARY_DIRS} ${${PREFIX}_LIB_SERACH_DIRS}
          NO_CMAKE_SYSTEM_PATH)

    endif()
  endif()

  # Check that it exists and set the found variable
  if(${PREFIX}_${_name}_LIBRARY AND EXISTS ${${PREFIX}_${_name}_LIBRARY})
    set(${PREFIX}_${_name}_FOUND TRUE)
  else()
    set(${PREFIX}_${_name}_FOUND FALSE)
  endif()
  
  mark_as_advanced(${PREFIX}_${_name}_LIBRARY)

endmacro()

