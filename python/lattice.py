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

class LLLInfo(ctypes.Structure):
  _fields_ = [("nullity",iType),("numSwaps",iType)]

lib.ElLLLCtrlDefault_s.argtypes = \
lib.ElLLLCtrlDefault_d.argtypes = \
  [c_void_p]
class LLLCtrl_s(ctypes.Structure):
  _fields_ = [("delta",sType),
              ("weak",bType),
              ("presort",bType),
              ("smallestFirst",bType),
              ("reorthogTol",sType),
              ("zeroTol",sType),
              ("progress",bType),
              ("time",bType)]
  def __init__(self):
    lib.ElLLLCtrlDefault_s(pointer(self))
class LLLCtrl_d(ctypes.Structure):
  _fields_ = [("delta",dType),
              ("weak",bType),
              ("presort",bType),
              ("smallestFirst",bType),
              ("reorthogTol",dType),
              ("zeroTol",dType),
              ("progress",bType),
              ("time",bType)]
  def __init__(self):
    lib.ElLLLCtrlDefault_d(pointer(self))

lib.ElLLL_s.argtypes = \
lib.ElLLL_c.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_s,POINTER(LLLInfo)]
lib.ElLLL_d.argtypes = \
lib.ElLLL_z.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_d,POINTER(LLLInfo)]

lib.ElLLLFull_s.argtypes = \
lib.ElLLLFull_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,
   LLLCtrl_s,POINTER(LLLInfo)]
lib.ElLLLFull_d.argtypes = \
lib.ElLLLFull_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,
   LLLCtrl_d,POINTER(LLLInfo)]

def LLL(B,full=True,ctrl=None):
  info = LLLInfo()
  if ctrl==None:
    if   B.tag == sTag: ctrl = LLLCtrl_s()
    elif B.tag == dTag: ctrl = LLLCtrl_d()
    else: DataExcept()

  U = Matrix(B.tag)
  UInv = Matrix(B.tag)
  QRMat = Matrix(B.tag)

  args = [B.obj,QRMat.obj,ctrl,pointer(info)]
  argsFull = [B.obj,U.obj,UInv.obj,QRMat.obj,ctrl,pointer(info)]
  if type(B) is Matrix:
    if   B.tag == sTag:
      if full: lib.ElLLLFull_s(*argsFull)
      else:    lib.ElLLL_s(*args)
    elif B.tag == dTag:
      if full: lib.ElLLLFull_d(*argsFull)
      else:    lib.ElLLL_d(*args)
    elif B.tag == cTag: 
      if full: lib.ElLLLFull_c(*argsFull)
      else:    lib.ElLLL_c(*args)
    elif B.tag == zTag:
      if full: lib.ElLLLFull_z(*argsFull)
      else:    lib.ElLLL_z(*args)
    else: DataExcept()
    if full:
      return U, UInv, QRMat, info
    else:
      return QRMat, info
  else: TypeExcept()

lib.ElLLLDelta_s.argtypes = \
lib.ElLLLDelta_c.argtypes = \
  [c_void_p,LLLCtrl_s,POINTER(sType)]
lib.ElLLLDelta_d.argtypes = \
lib.ElLLLDelta_z.argtypes = \
  [c_void_p,LLLCtrl_d,POINTER(dType)]

def LLLDelta(QRMat,ctrl=None):
  if ctrl==None:
    if   B.tag == sTag: ctrl = LLLCtrl_s()
    elif B.tag == dTag: ctrl = LLLCtrl_d()
    else: DataExcept()

  delta = TagToType(Base(QRMat.tag))()
  args = [QRMat.obj,ctrl,pointer(delta)]
  if type(QRMat) is Matrix:
    if   QRMat.tag == sTag: lib.ElLLLDelta_s(*args)
    elif QRMat.tag == dTag: lib.ElLLLDelta_d(*args)
    elif QRMat.tag == cTag: lib.ElLLLDelta_c(*args)
    elif QRMat.tag == zTag: lib.ElLLLDelta_z(*args)
    else: DataExcept()
    return delta.value
  else: TypeExcept()
