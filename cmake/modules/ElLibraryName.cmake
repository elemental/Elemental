#
#  Copyright 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
macro(EL_LIBRARY_NAME FULLNAME BASENAME)
  if(BUILD_SHARED_LIBS)
    set(${FULLNAME}
      ${CMAKE_SHARED_LIBRARY_PREFIX}${BASENAME}${CMAKE_SHARED_LIBRARY_SUFFIX})
  else()
    set(${FULLNAME}
      ${CMAKE_STATIC_LIBRARY_PREFIX}${BASENAME}${CMAKE_STATIC_LIBRARY_SUFFIX})
  endif()
endmacro()
