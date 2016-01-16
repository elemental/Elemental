#
#  Copyright 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
macro(EL_SET_PARENT_SCOPE)
  foreach(VAR ${ARGN})
    set(VAR ${VAR} PARENT_SCOPE)
  endforeach()
endmacro()
