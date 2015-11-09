#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

# Lattice
# *******

# LLL
# ===
lib.ElLLL_s.argtypes = [c_void_p,sType,sType,sType,sType]
lib.ElLLL_d.argtypes = [c_void_p,dType,dType,dType,dType]

def LLL(B,delta,eta,theta,innerTol):
  args = [B.obj,delta,eta,theta,innerTol]
  if type(B) is Matrix:
    if   B.tag == sTag: lib.ElLLL_s(*args)
    elif B.tag == dTag: lib.ElLLL_d(*args)
    else: DataExcept()
  else: TypeExcept()
