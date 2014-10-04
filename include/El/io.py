#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

# Input/Output
# ************
def Print(A,s):
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElPrint_i(A.obj,ctypes.c_char_p(s))
    elif A.tag == sTag: lib.ElPrint_s(A.obj,ctypes.c_char_p(s))
    elif A.tag == dTag: lib.ElPrint_d(A.obj,ctypes.c_char_p(s))
    elif A.tag == cTag: lib.ElPrint_c(A.obj,ctypes.c_char_p(s))
    elif A.tag == zTag: lib.ElPrint_z(A.obj,ctypes.c_char_p(s))
  else: print 'Unsupported matrix type'
