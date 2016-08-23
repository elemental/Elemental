# Try to find the MPC library
# See http://www.multiprecision.org/index.php?prog=mpc&page=home
#
# This module supports requiring a minimum version, e.g. you can do
#   find_package(MPC 1.0.3)
# to require version 1.0.3 to newer of MPC.
#
# Once done this will define
#
#  MPC_FOUND - system has MPC lib with correct version
#  MPC_INCLUDES - the MPC include directory
#  MPC_LIBRARIES - the MPC library
#  MPC_VERSION - MPC version

# Set MPC_INCLUDES

find_path(MPC_INCLUDES
  NAMES
  mpc.h
  PATHS
  $ENV{GMPDIR}
  $ENV{MPFRDIR}
  $ENV{MPCDIR}
  ${INCLUDE_INSTALL_DIR}
)

# Set MPC_FIND_VERSION to 1.0.0 if no minimum version is specified

if(NOT MPC_FIND_VERSION)
  if(NOT MPC_FIND_VERSION_MAJOR)
    set(MPC_FIND_VERSION_MAJOR 1)
  endif(NOT MPC_FIND_VERSION_MAJOR)
  if(NOT MPC_FIND_VERSION_MINOR)
    set(MPC_FIND_VERSION_MINOR 0)
  endif(NOT MPC_FIND_VERSION_MINOR)
  if(NOT MPC_FIND_VERSION_PATCH)
    set(MPC_FIND_VERSION_PATCH 0)
  endif(NOT MPC_FIND_VERSION_PATCH)

  set(MPC_FIND_VERSION "${MPC_FIND_VERSION_MAJOR}.${MPC_FIND_VERSION_MINOR}.${MPC_FIND_VERSION_PATCH}")
endif(NOT MPC_FIND_VERSION)


if(MPC_INCLUDES)

  # Set MPC_VERSION

  file(READ "${MPC_INCLUDES}/mpc.h" _mpc_version_header)

  string(REGEX MATCH "define[ \t]+MPC_VERSION_MAJOR[ \t]+([0-9]+)" _mpc_major_version_match "${_mpc_version_header}")
  set(MPC_MAJOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPC_VERSION_MINOR[ \t]+([0-9]+)" _mpc_minor_version_match "${_mpc_version_header}")
  set(MPC_MINOR_VERSION "${CMAKE_MATCH_1}")
  string(REGEX MATCH "define[ \t]+MPC_VERSION_PATCHLEVEL[ \t]+([0-9]+)" _mpc_patchlevel_version_match "${_mpc_version_header}")
  set(MPC_PATCHLEVEL_VERSION "${CMAKE_MATCH_1}")

  set(MPC_VERSION ${MPC_MAJOR_VERSION}.${MPC_MINOR_VERSION}.${MPC_PATCHLEVEL_VERSION})

  # Check whether found version exceeds minimum version

  if(${MPC_VERSION} VERSION_LESS ${MPC_FIND_VERSION})
    set(MPC_VERSION_OK FALSE)
    message(STATUS "MPC version ${MPC_VERSION} found in ${MPC_INCLUDES}, "
                   "but at least version ${MPC_FIND_VERSION} is required")
  else(${MPC_VERSION} VERSION_LESS ${MPC_FIND_VERSION})
    set(MPC_VERSION_OK TRUE)
  endif(${MPC_VERSION} VERSION_LESS ${MPC_FIND_VERSION})

endif(MPC_INCLUDES)

# Set MPC_LIBRARIES

find_library(MPC_LIBRARIES mpc
  PATHS $ENV{GMPDIR} $ENV{MPFRDIR} $ENV{MPCDIR} ${LIB_INSTALL_DIR})

# Epilogue

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(MPC DEFAULT_MSG
                                  MPC_INCLUDES MPC_LIBRARIES MPC_VERSION_OK)
mark_as_advanced(MPC_INCLUDES MPC_LIBRARIES)
