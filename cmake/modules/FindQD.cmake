# Try to find the QD library
# See http://crd-legacy.lbl.gov/~dhbailey/mpdist/index.html
#
# Once done this will define
#
#  QD_FOUND - system has QD lib with correct version
#  QD_INCLUDES - the QD include directory
#  QD_LIBRARIES - the QD library

# Copyright (c) 2016 Jack Poulson, <jack.poulson@gmail.com>
# Redistribution and use is allowed according to the terms of the BSD license.

# Set QD_INCLUDES
find_path(QD_INCLUDES
  NAMES
  qd/dd_real.h
  qd/qd_real.h
  PATHS
  $ENV{QDDIR}
  ${INCLUDE_INSTALL_DIR}
)

# Set QD_LIBRARIES
find_library(QD_LIBRARIES qd PATHS $ENV{QDDIR} ${LIB_INSTALL_DIR})

# Epilogue
include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(QD DEFAULT_MSG QD_INCLUDES QD_LIBRARIES)
mark_as_advanced(QD_INCLUDES QD_LIBRARIES)
