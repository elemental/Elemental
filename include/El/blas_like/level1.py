#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
import ctypes, numpy

# TODO: Switch to a different boolean type if appropriate
from ctypes import c_int    as bType
# TODO: Switch from c_int if Elemental was configured for 64-bit integers
from ctypes import c_int    as iType
from ctypes import c_float  as sType
from ctypes import c_double as dType
from ctypes import pointer
from ctypes import POINTER

# BLAS 1
# ======
# TODO: Move into a separate submodule

def Copy(A,B):
  CheckTag(A.tag)
  if A.tag != B.tag:
    raise Exception('Copying between datatypes is not yet supported in Python')
  if type(A) is Matrix and type(B) is Matrix:
    if   B.tag == iTag: lib.ElMatrixCopy_i(A.obj,B.obj)
    elif B.tag == sTag: lib.ElMatrixCopy_s(A.obj,B.obj)
    elif B.tag == dTag: lib.ElMatrixCopy_d(A.obj,B.obj)
    elif B.tag == cTag: lib.ElMatrixCopy_c(A.obj,B.obj)
    elif B.tag == zTag: lib.ElMatrixCopy_z(A.obj,B.obj)
  else: raise Exception('Unsupported matrix types')
