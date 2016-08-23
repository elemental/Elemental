#
#  Copyright 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
set(EL_HAVE_QT5 FALSE)
if(Qt5_DIR)
  set(Qt5_LIBDIR "${Qt5_DIR}/lib")
endif()
if(Qt5_LIBDIR)
  set(Qt5Widgets_DIR "${Qt5_LIBDIR}/cmake/Qt5Widgets")
endif()
if(EL_USE_QT5 OR Qt5Widgets_DIR)
  # Search for Qt5
  find_package(Qt5Widgets)
  if(Qt5Widgets_FOUND)
    set(EL_HAVE_QT5 TRUE)
    message(STATUS "Found Qt5")
  else()
    message(STATUS "Did NOT find Qt5")
  endif()
endif()
