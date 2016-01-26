#
#  Copyright (c) 2009-2016, Jack Poulson
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
              ("rank",iType),
              ("nullity",iType),
              ("numSwaps",iType),
              ("logVol",sType)]
class LLLInfo_d(ctypes.Structure):
  _fields_ = [("delta",dType),
              ("eta",dType),
              ("rank",iType),
              ("nullity",iType),
              ("numSwaps",iType),
              ("logVol",dType)]

(LLL_WEAK,LLL_NORMAL,LLL_DEEP,LLL_DEEP_REDUCE)=(0,1,2,3)

lib.ElLLLCtrlDefault_s.argtypes = \
lib.ElLLLCtrlDefault_d.argtypes = \
  [c_void_p]
class LLLCtrl_s(ctypes.Structure):
  _fields_ = [("delta",sType),
              ("eta",sType),
              ("variant",c_uint),
              ("recursive",bType),
              ("cutoff",iType),
              ("presort",bType),
              ("smallestFirst",bType),
              ("reorthogTol",sType),
              ("numOrthog",iType),
              ("zeroTol",sType),
              ("blockingThresh",sType),
              ("progress",bType),
              ("time",bType),
              ("jumpstart",bType),
              ("startCol",iType)]
  def __init__(self):
    lib.ElLLLCtrlDefault_s(pointer(self))
class LLLCtrl_d(ctypes.Structure):
  _fields_ = [("delta",dType),
              ("eta",dType),
              ("variant",c_uint),
              ("recursive",bType),
              ("cutoff",iType),
              ("presort",bType),
              ("smallestFirst",bType),
              ("reorthogTol",dType),
              ("numOrthog",iType),
              ("zeroTol",dType),
              ("blockingThresh",sType),
              ("progress",bType),
              ("time",bType),
              ("jumpstart",bType),
              ("startCol",iType)]
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
  [c_void_p,c_void_p,c_void_p,
   LLLCtrl_s,POINTER(LLLInfo_s)]
lib.ElLLLFull_d.argtypes = \
lib.ElLLLFull_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,
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
    R = Matrix(B.tag)

    args = [B.obj,ctrl,pointer(info)]
    argsFormR = [B.obj,R.obj,ctrl,pointer(info)]
    argsFull = [B.obj,U.obj,R.obj,ctrl,pointer(info)]

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
      return U, R, info
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

lib.ElLatticeImage_s.argtypes = \
lib.ElLatticeImage_c.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_s]
lib.ElLatticeImage_d.argtypes = \
lib.ElLatticeImage_z.argtypes = \
  [c_void_p,c_void_p,LLLCtrl_d]

def LatticeImage(B,ctrl=None):
  if ctrl==None:
    if   B.tag == sTag or B.tag == cTag: ctrl = LLLCtrl_s()
    elif B.tag == dTag or B.tag == zTag: ctrl = LLLCtrl_d()

  if type(B) is Matrix:
    M = Matrix(B.tag)
    args = [B.obj,M.obj,ctrl]
    if   B.tag == sTag: lib.ElLatticeImage_s(*args)
    elif B.tag == dTag: lib.ElLatticeImage_d(*args)
    elif B.tag == cTag: lib.ElLatticeImage_c(*args)
    elif B.tag == zTag: lib.ElLatticeImage_z(*args)
    else: DataExcept()
    return M
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
    numExact = iType()
    B = Matrix(z.tag)
    U = Matrix(z.tag)
    args = [z.obj,NSqrt,B.obj,U.obj,ctrl,pointer(numExact)]
    if   z.tag == sTag: lib.ElZDependenceSearch_s(*args)
    elif z.tag == dTag: lib.ElZDependenceSearch_d(*args)
    elif z.tag == cTag: lib.ElZDependenceSearch_c(*args)
    elif z.tag == zTag: lib.ElZDependenceSearch_z(*args)
    else: DataExcept()
    return numExact.value, B, U
  else: TypeExcept()

# Algebraic relation search
# =========================

lib.ElAlgebraicRelationSearch_s.argtypes = \
  [sType,iType,sType,c_void_p,c_void_p,LLLCtrl_s,POINTER(iType)]
lib.ElAlgebraicRelationSearch_d.argtypes = \
  [dType,iType,dType,c_void_p,c_void_p,LLLCtrl_d,POINTER(iType)]
lib.ElAlgebraicRelationSearch_c.argtypes = \
  [cType,iType,sType,c_void_p,c_void_p,LLLCtrl_s,POINTER(iType)]
lib.ElAlgebraicRelationSearch_z.argtypes = \
  [zType,iType,dType,c_void_p,c_void_p,LLLCtrl_d,POINTER(iType)]

def AlgebraicRelationSearch(alpha,n,NSqrt,ctrl=None):
  if   type(alpha) == sType:
    tag = sTag
  elif type(alpha) == dType:
    tag = dTag
  elif type(alpha) == cType:
    tag = cTag
  elif type(alpha) == zType:
    tag = zTag

  if ctrl==None:
    if   tag == sTag or tag == cTag: ctrl = LLLCtrl_s()
    elif tag == dTag or tag == zTag: ctrl = LLLCtrl_d()

  # TODO: Provide a way to pick a different return type from Matrix
  #if type(z) is Matrix:
  if True:
    numExact = iType()
    B = Matrix(tag)
    U = Matrix(tag)
    args = [alpha,n,NSqrt,B.obj,U.obj,ctrl,pointer(numExact)]
    if   type(alpha) == sType: lib.ElAlgebraicRelationSearch_s(*args)
    elif type(alpha) == dType: lib.ElAlgebraicRelationSearch_d(*args)
    elif type(alpha) == cType: lib.ElAlgebraicRelationSearch_c(*args)
    elif type(alpha) == zType: lib.ElAlgebraicRelationSearch_z(*args)
    else: DataExcept()
    return numExact.value, B, U
  else: TypeExcept()
