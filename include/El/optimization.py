#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

from ctypes import CFUNCTYPE

# Optimization
# ************

# Utilities
# =========

# Coherence
# ---------
def Coherence(A):
  value = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCoherence_s(A.obj,pointer(value))
    elif A.tag == dTag: lib.ElCoherence_d(A.obj,pointer(value))
    elif A.tag == cTag: lib.ElCoherence_c(A.obj,pointer(value))
    elif A.tag == zTag: lib.ElCoherence_z(A.obj,pointer(value))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElCoherenceDist_s(A.obj,pointer(value))
    elif A.tag == dTag: lib.ElCoherenceDist_d(A.obj,pointer(value))
    elif A.tag == cTag: lib.ElCoherenceDist_c(A.obj,pointer(value))
    elif A.tag == zTag: lib.ElCoherenceDist_z(A.obj,pointer(value))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
  return value
