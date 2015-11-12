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

# Gram-Schmidt
# ============
lib.ElLatticeGramSchmidt_s.argtypes = \
lib.ElLatticeGramSchmidt_d.argtypes = \
lib.ElLatticeGramSchmidt_c.argtypes = \
lib.ElLatticeGramSchmidt_z.argtypes = \
  [c_void_p,c_void_p]

def LatticeGramSchmidt(B):
  if type(B) is Matrix:
    G = Matrix(B.tag)
    M = Matrix(B.tag)
    args = [B.obj,G.obj,M.obj]
    if   B.tag == sTag: lib.ElLatticeGramSchmidt_s(*args)
    elif B.tag == dTag: lib.ElLatticeGramSchmidt_d(*args)
    elif B.tag == cTag: lib.ElLatticeGramSchmidt_c(*args)
    elif B.tag == zTag: lib.ElLatticeGramSchmidt_z(*args)
    else: DataExcept()
    return G, M
  else: TypeExcept()

# LLL
# ===
lib.ElLLL_s.argtypes = \
lib.ElLLL_c.argtypes = \
  [c_void_p,c_void_p,sType,sType,bType,bType,bType,POINTER(iType)]
lib.ElLLL_d.argtypes = \
lib.ElLLL_z.argtypes = \
  [c_void_p,c_void_p,dType,dType,bType,bType,bType,POINTER(iType)]

def LLL(B,delta,innerTol=0,presort=False,smallestFirst=True,progress=False):
  numBacktrack = iType()
  QRMat = Matrix(B.tag)
  args = [B.obj,QRMat.obj,delta,innerTol,
          presort,smallestFirst,progress,pointer(numBacktrack)]
  if type(B) is Matrix:
    if   B.tag == sTag: lib.ElLLL_s(*args)
    elif B.tag == dTag: lib.ElLLL_d(*args)
    elif B.tag == cTag: lib.ElLLL_c(*args)
    elif B.tag == zTag: lib.ElLLL_z(*args)
    else: DataExcept()
    return QRMat, numBacktrack.value
  else: TypeExcept()

lib.ElLLLDelta_s.argtypes = \
lib.ElLLLDelta_c.argtypes = \
  [c_void_p,POINTER(sType)]
lib.ElLLLDelta_d.argtypes = \
lib.ElLLLDelta_z.argtypes = \
  [c_void_p,POINTER(dType)]

def LLLDelta(QRMat):
  delta = TagToType(Base(QRMat.tag))()
  args = [QRMat.obj,pointer(delta)]
  if type(QRMat) is Matrix:
    if   QRMat.tag == sTag: lib.ElLLLDelta_s(*args)
    elif QRMat.tag == dTag: lib.ElLLLDelta_d(*args)
    elif QRMat.tag == cTag: lib.ElLLLDelta_c(*args)
    elif QRMat.tag == zTag: lib.ElLLLDelta_z(*args)
    else: DataExcept()
    return delta.value
  else: TypeExcept()
