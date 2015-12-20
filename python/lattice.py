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

class LLLInfo_s(ctypes.Structure):
  _fields_ = [("delta",sType),
              ("eta",sType),
              ("nullity",iType),
              ("numSwaps",iType)]
class LLLInfo_d(ctypes.Structure):
  _fields_ = [("delta",dType),
              ("eta",dType),
              ("nullity",iType),
              ("numSwaps",iType)]

lib.ElLLLCtrlDefault_s.argtypes = \
lib.ElLLLCtrlDefault_d.argtypes = \
  [c_void_p]
class LLLCtrl_s(ctypes.Structure):
  _fields_ = [("delta",sType),
              ("eta",sType),
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
              ("eta",dType),
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
  [c_void_p,LLLCtrl_s,POINTER(LLLInfo_s)]
lib.ElLLL_d.argtypes = \
lib.ElLLL_z.argtypes = \
  [c_void_p,LLLCtrl_d,POINTER(LLLInfo_d)]

lib.ElLLLFormR_s.argtypes = \
lib.ElLLLFormR_c.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_s,POINTER(LLLInfo_s)]
lib.ElLLLFormR_d.argtypes = \
lib.ElLLLFormR_z.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_d,POINTER(LLLInfo_d)]

lib.ElLLLFull_s.argtypes = \
lib.ElLLLFull_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,
   LLLCtrl_s,POINTER(LLLInfo_s)]
lib.ElLLLFull_d.argtypes = \
lib.ElLLLFull_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,
   LLLCtrl_d,POINTER(LLLInfo_d)]

(LLL_LATTICE_ONLY,LLL_FORM_R,LLL_FULL)=(0,1,2)

def LLL(B,mode=LLL_LATTICE_ONLY,ctrl=None):
  if   B.tag == sTag or B.tag == cTag: info = LLLInfo_s()
  elif B.tag == dTag or B.tag == zTag: info = LLLInfo_d()
  else: DataExcept()

  if ctrl==None:
    if   B.tag == sTag or B.tag == cTag: ctrl = LLLCtrl_s()
    elif B.tag == dTag or B.tag == zTag: ctrl = LLLCtrl_d()
    else: DataExcept()

  if type(B) is Matrix:
    U = Matrix(B.tag)
    UInv = Matrix(B.tag)
    R = Matrix(B.tag)

    args = [B.obj,ctrl,pointer(info)]
    argsFormR = [B.obj,R.obj,ctrl,pointer(info)]
    argsFull = [B.obj,U.obj,UInv.obj,R.obj,ctrl,pointer(info)]

    if   B.tag == sTag:
      if   mode==LLL_FULL:   lib.ElLLLFull_s(*argsFull)
      elif mode==LLL_FORM_R: lib.ElLLLFormR_s(*argsFormR)
      else:                  lib.ElLLL_s(*args)
    elif B.tag == dTag:
      if   mode==LLL_FULL:   lib.ElLLLFull_d(*argsFull)
      elif mode==LLL_FORM_R: lib.ElLLLFormR_d(*argsFormR)
      else:                  lib.ElLLL_d(*args)
    elif B.tag == cTag: 
      if   mode==LLL_FULL:   lib.ElLLLFull_c(*argsFull)
      elif mode==LLL_FORM_R: lib.ElLLLFormR_c(*argsFormR)
      else:                  lib.ElLLL_c(*args)
    elif B.tag == zTag:
      if   mode==LLL_FULL:   lib.ElLLLFull_z(*argsFull)
      elif mode==LLL_FORM_R: lib.ElLLLFormR_z(*argsFormR)
      else:                  lib.ElLLL_z(*args)
    else: DataExcept()

    if mode==LLL_FULL:
      return U, UInv, R, info
    elif mode==LLL_FORM_R:
      return R, info
    else:
      return info

  else: TypeExcept()

# Lattice image/kernel decomposition
# ==================================

lib.ElLatticeImageAndKernel_s.argtypes = \
lib.ElLatticeImageAndKernel_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,LLLCtrl_s]
lib.ElLatticeImageAndKernel_d.argtypes = \
lib.ElLatticeImageAndKernel_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,LLLCtrl_d]

def LatticeImageAndKernel(B,ctrl=None):
  if ctrl==None:
    if   B.tag == sTag or B.tag == cTag: ctrl = LLLCtrl_s()
    elif B.tag == dTag or B.tag == zTag: ctrl = LLLCtrl_d()

  if type(B) is Matrix:
    M = Matrix(B.tag)
    K = Matrix(B.tag)
    args = [B.obj,M.obj,K.obj,ctrl]
    if   B.tag == sTag: lib.ElLatticeImageAndKernel_s(*args)
    elif B.tag == dTag: lib.ElLatticeImageAndKernel_d(*args)
    elif B.tag == cTag: lib.ElLatticeImageAndKernel_c(*args)
    elif B.tag == zTag: lib.ElLatticeImageAndKernel_z(*args)
    else: DataExcept()
    return M, K
  else: TypeExcept()

lib.ElLatticeKernel_s.argtypes = \
lib.ElLatticeKernel_c.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_s]
lib.ElLatticeKernel_d.argtypes = \
lib.ElLatticeKernel_z.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_d]

def LatticeKernel(B,ctrl=None):
  if ctrl==None:
    if   B.tag == sTag or B.tag == cTag: ctrl = LLLCtrl_s()
    elif B.tag == dTag or B.tag == zTag: ctrl = LLLCtrl_d()

  if type(B) is Matrix:
    K = Matrix(B.tag)
    args = [B.obj,K.obj,ctrl]
    if   B.tag == sTag: lib.ElLatticeKernel_s(*args)
    elif B.tag == dTag: lib.ElLatticeKernel_d(*args)
    elif B.tag == cTag: lib.ElLatticeKernel_c(*args)
    elif B.tag == zTag: lib.ElLatticeKernel_z(*args)
    else: DataExcept()
    return K
  else: TypeExcept()

# Z-dependence search
# ===================

lib.ElZDependenceSearch_s.argtypes = \
lib.ElZDependenceSearch_c.argtypes = \
  [c_void_p,sType,c_void_p,c_void_p,LLLCtrl_s,POINTER(iType)]
lib.ElZDependenceSearch_d.argtypes = \
lib.ElZDependenceSearch_z.argtypes = \
  [c_void_p,dType,c_void_p,c_void_p,LLLCtrl_d,POINTER(iType)]

def ZDependenceSearch(z,NSqrt,ctrl=None):
  if ctrl==None:
    if   z.tag == sTag or z.tag == cTag: ctrl = LLLCtrl_s()
    elif z.tag == dTag or z.tag == zTag: ctrl = LLLCtrl_d()

  if type(z) is Matrix:
    numFound = iType()
    B = Matrix(z.tag)
    U = Matrix(z.tag)
    args = [z.obj,NSqrt,B.obj,U.obj,ctrl,pointer(numFound)]
    if   z.tag == sTag: lib.ElZDependenceSearch_s(*args)
    elif z.tag == dTag: lib.ElZDependenceSearch_d(*args)
    elif z.tag == cTag: lib.ElZDependenceSearch_c(*args)
    elif z.tag == zTag: lib.ElZDependenceSearch_z(*args)
    else: DataExcept()
    return numFound.value, B, U
  else: TypeExcept()
