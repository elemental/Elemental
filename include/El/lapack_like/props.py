#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
import ctypes

# Emulate an enum for the sign scaling
(SIGN_SCALE_NONE,SIGN_SCALE_DET,SIGN_SCALE_FROB)=(0,1,2)

lib.ElSignCtrlDefault_s.argtypes = [c_void_p]
lib.ElSignCtrlDefault_s.restype = c_uint
class SignCtrl_s(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",sType),
              ("power",sType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_s(pointer(self))

lib.ElSignCtrlDefault_d.argtypes = [c_void_p]
lib.ElSignCtrlDefault_d.restype = c_uint
class SignCtrl_d(ctypes.Structure):
  _fields_ = [("maxIts",iType),
              ("tol",dType),
              ("power",dType),
              ("scaling",c_uint)]
  def __init__(self):
    lib.ElSignCtrlDefault_d(pointer(self))

lib.ElHessQRCtrlDefault.argtypes = [c_void_p]
lib.ElHessQRCtrlDefault.restype = c_uint
class HessQRCtrl(ctypes.Structure):
  _fields_ = [("distAED",bType),
              ("blockHeight",iType),("blockWidth",iType)]
  def __init__(self):
    lib.ElHessQRCtrlDefault(pointer(self))

lib.ElSDCCtrlDefault_s.argtypes = [c_void_p]
lib.ElSDCCtrlDefault_s.restype = c_uint
class SDCCtrl_s(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",sType),
              ("spreadFactor",sType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_s)]
  def __init__(self):
    lib.ElSDCCtrlDefault_s(pointer(self))

lib.ElSDCCtrlDefault_d.argtypes = [c_void_p]
lib.ElSDCCtrlDefault_d.restype = c_uint
class SDCCtrl_d(ctypes.Structure):
  _fields_ = [("cutoff",iType),
              ("maxInnerIts",iType),("maxOuterIts",iType),
              ("tol",dType),
              ("spreadFactor",dType),
              ("random",bType),
              ("progress",bType),
              ("signCtrl",SignCtrl_d)]
  def __init__(self):
    lib.ElSDCCtrlDefault_d(pointer(self))

lib.ElSchurCtrlDefault_s.argtypes = [c_void_p]
lib.ElSchurCtrlDefault_s.restype = c_uint
class SchurCtrl_s(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_s)]
  def __init__(self):
    lib.ElSchurCtrlDefault_s(pointer(self))

lib.ElSchurCtrlDefault_d.argtypes = [c_void_p]
lib.ElSchurCtrlDefault_d.restype = c_uint
class SchurCtrl_d(ctypes.Structure):
  _fields_ = [("useSDC",bType),
              ("qrCtrl",HessQRCtrl),
              ("sdcCtrl",SDCCtrl_d)]
  def __init__(self):
    lib.ElSchurCtrlDefault_d(pointer(self))

# Condition number
# ================
lib.ElCondition_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElCondition_s.restype = c_uint
lib.ElCondition_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElCondition_d.restype = c_uint
lib.ElCondition_c.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElCondition_c.restype = c_uint
lib.ElCondition_z.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElCondition_z.restype = c_uint
lib.ElConditionDist_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElConditionDist_s.restype = c_uint
lib.ElConditionDist_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElConditionDist_d.restype = c_uint
lib.ElConditionDist_c.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElConditionDist_c.restype = c_uint
lib.ElConditionDist_z.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElConditionDist_z.restype = c_uint
def Condition(A,normType=FROBENIUS_NORM):
  cond = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCondition_s(A.obj,normType,pointer(cond))
    elif A.tag == dTag: lib.ElCondition_d(A.obj,normType,pointer(cond))
    elif A.tag == cTag: lib.ElCondition_c(A.obj,normType,pointer(cond))
    elif A.tag == zTag: lib.ElCondition_z(A.obj,normType,pointer(cond))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElConditionDist_s(A.obj,normType,pointer(cond))
    elif A.tag == dTag: lib.ElConditionDist_d(A.obj,normType,pointer(cond))
    elif A.tag == cTag: lib.ElConditionDist_c(A.obj,normType,pointer(cond))
    elif A.tag == zTag: lib.ElConditionDist_z(A.obj,normType,pointer(cond))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
  return cond

lib.ElFrobeniusCondition_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusCondition_s.restype = c_uint
lib.ElFrobeniusCondition_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusCondition_d.restype = c_uint
lib.ElFrobeniusCondition_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusCondition_c.restype = c_uint
lib.ElFrobeniusCondition_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusCondition_z.restype = c_uint
lib.ElFrobeniusConditionDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusConditionDist_s.restype = c_uint
lib.ElFrobeniusConditionDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusConditionDist_d.restype = c_uint
lib.ElFrobeniusConditionDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusConditionDist_c.restype = c_uint
lib.ElFrobeniusConditionDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusConditionDist_z.restype = c_uint
def FrobeniusCondition(A):
  cond = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElFrobeniusCondition_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElFrobeniusCondition_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElFrobeniusCondition_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElFrobeniusCondition_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElFrobeniusConditionDist_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElFrobeniusConditionDist_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElFrobeniusConditionDist_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElFrobeniusConditionDist_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElInfinityCondition_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityCondition_s.restype = c_uint
lib.ElInfinityCondition_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityCondition_d.restype = c_uint
lib.ElInfinityCondition_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityCondition_c.restype = c_uint
lib.ElInfinityCondition_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityCondition_z.restype = c_uint
lib.ElInfinityConditionDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityConditionDist_s.restype = c_uint
lib.ElInfinityConditionDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityConditionDist_d.restype = c_uint
lib.ElInfinityConditionDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityConditionDist_c.restype = c_uint
lib.ElInfinityConditionDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityConditionDist_z.restype = c_uint
def InfinityCondition(A):
  cond = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElInfinityCondition_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElInfinityCondition_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElInfinityCondition_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElInfinityCondition_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElInfinityConditionDist_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElInfinityConditionDist_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElInfinityConditionDist_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElInfinityConditionDist_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElMaxCondition_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxCondition_s.restype = c_uint
lib.ElMaxCondition_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxCondition_d.restype = c_uint
lib.ElMaxCondition_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxCondition_c.restype = c_uint
lib.ElMaxCondition_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxCondition_z.restype = c_uint
lib.ElMaxConditionDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxConditionDist_s.restype = c_uint
lib.ElMaxConditionDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxConditionDist_d.restype = c_uint
lib.ElMaxConditionDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxConditionDist_c.restype = c_uint
lib.ElMaxConditionDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxConditionDist_z.restype = c_uint
def MaxCondition(A):
  cond = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElMaxCondition_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElMaxCondition_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElMaxCondition_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElMaxCondition_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElMaxConditionDist_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElMaxConditionDist_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElMaxConditionDist_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElMaxConditionDist_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElOneCondition_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneCondition_s.restype = c_uint
lib.ElOneCondition_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneCondition_d.restype = c_uint
lib.ElOneCondition_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneCondition_c.restype = c_uint
lib.ElOneCondition_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneCondition_z.restype = c_uint
lib.ElOneConditionDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneConditionDist_s.restype = c_uint
lib.ElOneConditionDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneConditionDist_d.restype = c_uint
lib.ElOneConditionDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneConditionDist_c.restype = c_uint
lib.ElOneConditionDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneConditionDist_z.restype = c_uint
def OneCondition(A):
  cond = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElOneCondition_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElOneCondition_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElOneCondition_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElOneCondition_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElOneConditionDist_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElOneConditionDist_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElOneConditionDist_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElOneConditionDist_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElTwoCondition_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoCondition_s.restype = c_uint
lib.ElTwoCondition_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoCondition_d.restype = c_uint
lib.ElTwoCondition_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoCondition_c.restype = c_uint
lib.ElTwoCondition_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoCondition_z.restype = c_uint
lib.ElTwoConditionDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoConditionDist_s.restype = c_uint
lib.ElTwoConditionDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoConditionDist_d.restype = c_uint
lib.ElTwoConditionDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoConditionDist_c.restype = c_uint
lib.ElTwoConditionDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoConditionDist_z.restype = c_uint
def TwoCondition(A):
  cond = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoCondition_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElTwoCondition_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElTwoCondition_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElTwoCondition_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoConditionDist_s(A.obj,pointer(cond))
    elif A.tag == dTag: lib.ElTwoConditionDist_d(A.obj,pointer(cond))
    elif A.tag == cTag: lib.ElTwoConditionDist_c(A.obj,pointer(cond))
    elif A.tag == zTag: lib.ElTwoConditionDist_z(A.obj,pointer(cond))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Determinant
# ===========

# Return the result in a safer, expanded format
# ---------------------------------------------
lib.ElSafeDeterminant_s.argtypes = [c_void_p,POINTER(SafeProduct_s)]
lib.ElSafeDeterminant_s.restype = c_uint
lib.ElSafeDeterminant_d.argtypes = [c_void_p,POINTER(SafeProduct_d)]
lib.ElSafeDeterminant_d.restype = c_uint
lib.ElSafeDeterminant_c.argtypes = [c_void_p,POINTER(SafeProduct_c)]
lib.ElSafeDeterminant_c.restype = c_uint
lib.ElSafeDeterminant_z.argtypes = [c_void_p,POINTER(SafeProduct_z)]
lib.ElSafeDeterminant_z.restype = c_uint
lib.ElSafeDeterminantDist_s.argtypes = [c_void_p,POINTER(SafeProduct_s)]
lib.ElSafeDeterminantDist_s.restype = c_uint
lib.ElSafeDeterminantDist_d.argtypes = [c_void_p,POINTER(SafeProduct_d)]
lib.ElSafeDeterminantDist_d.restype = c_uint
lib.ElSafeDeterminantDist_c.argtypes = [c_void_p,POINTER(SafeProduct_c)]
lib.ElSafeDeterminantDist_c.restype = c_uint
lib.ElSafeDeterminantDist_z.argtypes = [c_void_p,POINTER(SafeProduct_z)]
lib.ElSafeDeterminantDist_z.restype = c_uint
def SafeDeterminant(A):
  if type(A) is Matrix:
    if   A.tag == sTag: 
      safeProd = SafeProduct_s()
      lib.ElSafeDeterminant_s(A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == dTag: 
      safeProd = SafeProduct_d()
      lib.ElSafeDeterminant_d(A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == cTag: 
      safeProd = SafeProduct_c()
      lib.ElSafeDeterminant_c(A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == zTag: 
      safeProd = SafeProduct_z()
      lib.ElSafeDeterminant_z(A.obj,pointer(safeProd))
      return safeProd
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      safeProd = SafeProduct_s()
      lib.ElSafeDeterminantDist_s(A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == dTag: 
      safeProd = SafeProduct_d()
      lib.ElSafeDeterminantDist_d(A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == cTag: 
      safeProd = SafeProduct_c()
      lib.ElSafeDeterminantDist_c(A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == zTag: 
      safeProd = SafeProduct_z()
      lib.ElSafeDeterminantDist_z(A.obj,pointer(safeProd))
      return safeProd
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElSafeHPDDeterminant_s.argtypes = [c_uint,c_void_p,POINTER(SafeProduct_s)]
lib.ElSafeHPDDeterminant_s.restype = c_uint
lib.ElSafeHPDDeterminant_d.argtypes = [c_uint,c_void_p,POINTER(SafeProduct_d)]
lib.ElSafeHPDDeterminant_d.restype = c_uint
lib.ElSafeHPDDeterminant_c.argtypes = [c_uint,c_void_p,POINTER(SafeProduct_s)]
lib.ElSafeHPDDeterminant_c.restype = c_uint
lib.ElSafeHPDDeterminant_z.argtypes = [c_uint,c_void_p,POINTER(SafeProduct_d)]
lib.ElSafeHPDDeterminant_z.restype = c_uint
lib.ElSafeHPDDeterminantDist_s.argtypes = \
  [c_uint,c_void_p,POINTER(SafeProduct_s)]
lib.ElSafeHPDDeterminantDist_s.restype = c_uint
lib.ElSafeHPDDeterminantDist_d.argtypes = \
  [c_uint,c_void_p,POINTER(SafeProduct_d)]
lib.ElSafeHPDDeterminantDist_d.restype = c_uint
lib.ElSafeHPDDeterminantDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(SafeProduct_s)]
lib.ElSafeHPDDeterminantDist_c.restype = c_uint
lib.ElSafeHPDDeterminantDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(SafeProduct_d)]
lib.ElSafeHPDDeterminantDist_z.restype = c_uint
def SafeHPDDeterminant(uplo,A):
  if type(A) is Matrix:
    if   A.tag == sTag: 
      safeProd = SafeProduct_s()
      lib.ElSafeHPDDeterminant_s(uplo,A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == dTag: 
      safeProd = SafeProduct_d()
      lib.ElSafeHPDDeterminant_d(uplo,A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == cTag: 
      safeProd = SafeProduct_s()
      lib.ElSafeHPDDeterminant_c(uplo,A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == zTag: 
      safeProd = SafeProduct_d()
      lib.ElSafeHPDDeterminant_z(uplo,A.obj,pointer(safeProd))
      return safeProd
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      safeProd = SafeProduct_s()
      lib.ElSafeHPDDeterminantDist_s(uplo,A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == dTag: 
      safeProd = SafeProduct_d()
      lib.ElSafeHPDDeterminantDist_d(uplo,A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == cTag: 
      safeProd = SafeProduct_s()
      lib.ElSafeHPDDeterminantDist_c(uplo,A.obj,pointer(safeProd))
      return safeProd
    elif A.tag == zTag: 
      safeProd = SafeProduct_d()
      lib.ElSafeHPDDeterminantDist_z(uplo,A.obj,pointer(safeProd))
      return safeProd
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Directly return the result (warning: their may be under-/over-flow)
# -------------------------------------------------------------------
lib.ElDeterminant_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElDeterminant_s.restype = c_uint
lib.ElDeterminant_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElDeterminant_d.restype = c_uint
lib.ElDeterminant_c.argtypes = [c_void_p,POINTER(cType)]
lib.ElDeterminant_c.restype = c_uint
lib.ElDeterminant_z.argtypes = [c_void_p,POINTER(zType)]
lib.ElDeterminant_z.restype = c_uint
lib.ElDeterminantDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElDeterminantDist_s.restype = c_uint
lib.ElDeterminantDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElDeterminantDist_d.restype = c_uint
lib.ElDeterminantDist_c.argtypes = [c_void_p,POINTER(cType)]
lib.ElDeterminantDist_c.restype = c_uint
lib.ElDeterminantDist_z.argtypes = [c_void_p,POINTER(zType)]
lib.ElDeterminantDist_z.restype = c_uint
def Determinant(A):
  prod = TagToType(A.tag)()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElDeterminant_s(A.obj,pointer(prod))
    elif A.tag == dTag: lib.ElDeterminant_d(A.obj,pointer(prod))
    elif A.tag == cTag: lib.ElDeterminant_c(A.obj,pointer(prod))
    elif A.tag == zTag: lib.ElDeterminant_z(A.obj,pointer(prod))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElDeterminantDist_s(A.obj,pointer(prod))
    elif A.tag == dTag: lib.ElDeterminantDist_d(A.obj,pointer(prod))
    elif A.tag == cTag: lib.ElDeterminantDist_c(A.obj,pointer(prod))
    elif A.tag == zTag: lib.ElDeterminantDist_z(A.obj,pointer(prod))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
  return prod

lib.ElHPDDeterminant_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHPDDeterminant_s.restype = c_uint
lib.ElHPDDeterminant_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHPDDeterminant_d.restype = c_uint
lib.ElHPDDeterminant_c.argtypes = [c_uint,c_void_p,POINTER(cType)]
lib.ElHPDDeterminant_c.restype = c_uint
lib.ElHPDDeterminant_z.argtypes = [c_uint,c_void_p,POINTER(zType)]
lib.ElHPDDeterminant_z.restype = c_uint
lib.ElHPDDeterminantDist_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHPDDeterminantDist_s.restype = c_uint
lib.ElHPDDeterminantDist_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHPDDeterminantDist_d.restype = c_uint
lib.ElHPDDeterminantDist_c.argtypes = [c_uint,c_void_p,POINTER(cType)]
lib.ElHPDDeterminantDist_c.restype = c_uint
lib.ElHPDDeterminantDist_z.argtypes = [c_uint,c_void_p,POINTER(zType)]
lib.ElHPDDeterminantDist_z.restype = c_uint
def HPDDeterminant(uplo,A):
  prod = TagToType(Base(A.tag))()
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPDDeterminant_s(uplo,A.obj,pointer(prod))
    elif A.tag == dTag: lib.ElHPDDeterminant_d(uplo,A.obj,pointer(prod))
    elif A.tag == cTag: lib.ElHPDDeterminant_c(uplo,A.obj,pointer(prod))
    elif A.tag == zTag: lib.ElHPDDeterminant_z(uplo,A.obj,pointer(prod))
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPDDeterminantDist_s(uplo,A.obj,pointer(prod))
    elif A.tag == dTag: lib.ElHPDDeterminantDist_d(uplo,A.obj,pointer(prod))
    elif A.tag == cTag: lib.ElHPDDeterminantDist_c(uplo,A.obj,pointer(prod))
    elif A.tag == zTag: lib.ElHPDDeterminantDist_z(uplo,A.obj,pointer(prod))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
  return prod

# Inertia
# =======
# TODO

# Norm
# ====
# TODO

# Pseudospectra
# =============
lib.ElSnapshotCtrlDefault.argtypes = [c_void_p]
lib.ElSnapshotCtrlDefault.restype = c_uint
lib.ElSnapshotCtrlDestroy.argtypes = [c_void_p]
lib.ElSnapshotCtrlDestroy.restype = c_uint
class SnapshotCtrl(ctypes.Structure):
  _fields_ = [("realSize",iType),("imagSize",iType),
              ("imgSaveFreq",iType),("numSaveFreq",iType),
              ("imgDispFreq",iType),
              ("imgSaveCount",iType),("numSaveCount",iType),
              ("imgDispCount",iType),
              ("imgBase",c_char_p),("numBase",c_char_p),
              ("imgFormat",c_uint),("numFormat",c_uint),
              ("itCounts",bType)]
  def __init__(self):
    lib.ElSnaphsotCtrlDefault(pointer(self)) 
  def Destroy(self):
    lib.ElSnapshotCtrlDestroy(pointer(self))

# Emulate an enum for the pseudospectral norm
(PS_TWO_NORM,PS_ONE_NORM)=(0,1)

lib.ElPseudospecCtrlDefault_s.argtypes = [c_void_p]
lib.ElPseudospecCtrlDefault_s.restype = c_uint
lib.ElPseudospecCtrlDestroy_s.argtypes = [c_void_p]
lib.ElPseudospecCtrlDestroy_s.restype = c_uint
class PseudospecCtrl_s(ctypes.Structure):
  _fields_ = [("norm",c_uint),
              ("blockWidth",iType),
              ("schur",bType),
              ("forceComplexSchur",bType),
              ("forceComplexPs",bType),
              ("schurCtrl",SchurCtrl_s),
              ("maxIts",iType),
              ("tol",sType),
              ("deflate",bType),
              ("arnoldi",bType),
              ("basisSize",iType),
              ("reorthog",bType),
              ("progress",bType),
              ("snapCtrl",SnapshotCtrl)]
  def __init__(self):
    lib.ElPseudospecCtrlDefault_s(pointer(self))
  def Destroy(self):
    lib.ElPseudospecCtrlDestroy_s(pointer(self))

lib.ElPseudospecCtrlDefault_d.argtypes = [c_void_p]
lib.ElPseudospecCtrlDefault_d.restype = c_uint
lib.ElPseudospecCtrlDestroy_d.argtypes = [c_void_p]
lib.ElPseudospecCtrlDestroy_d.restype = c_uint
class PseudospecCtrl_d(ctypes.Structure):
  _fields_ = [("norm",c_uint),
              ("blockWidth",iType),
              ("schur",bType),
              ("forceComplexSchur",bType),
              ("forceComplexPs",bType),
              ("schurCtrl",SchurCtrl_d),
              ("maxIts",iType),
              ("tol",dType),
              ("deflate",bType),
              ("arnoldi",bType),
              ("basisSize",iType),
              ("reorthog",bType),
              ("progress",bType),
              ("snapCtrl",SnapshotCtrl)]
  def __init__(self):
    lib.ElPseudospecCtrlDefault_d(pointer(self))
  def Destroy(self):
    lib.ElPseudospecCtrlDestroy_d(pointer(self))

# Automatically determine the window of the complex plane
# -------------------------------------------------------
# The choice is based upon a few different norms of the Schur factor, as simply
# using the spectral radius would be insufficient for highly non-normal 
# matrices, e.g., a Jordan block with eigenvalue zero
lib.ElPseudospectralAutoWindow_s.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindow_s.restype = c_uint
lib.ElPseudospectralAutoWindow_d.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindow_d.restype = c_uint
lib.ElPseudospectralAutoWindow_c.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindow_c.restype = c_uint
lib.ElPseudospectralAutoWindow_z.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindow_z.restype = c_uint
lib.ElPseudospectralAutoWindowDist_s.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindowDist_s.restype = c_uint
lib.ElPseudospectralAutoWindowDist_d.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindowDist_d.restype = c_uint
lib.ElPseudospectralAutoWindowDist_c.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindowDist_c.restype = c_uint
lib.ElPseudospectralAutoWindowDist_z.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElPseudospectralAutoWindowDist_z.restype = c_uint
lib.ElPseudospectralAutoWindowX_s.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowX_s.restype = c_uint
lib.ElPseudospectralAutoWindowX_d.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowX_d.restype = c_uint
lib.ElPseudospectralAutoWindowX_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowX_c.restype = c_uint
lib.ElPseudospectralAutoWindowX_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowX_z.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_s.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowXDist_s.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_d.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowXDist_d.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralAutoWindowXDist_c.restype = c_uint
lib.ElPseudospectralAutoWindowXDist_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralAutoWindowXDist_z.restype = c_uint
def PseudospectralAutoWindow(A,realSize=200,imagSize=200,ctrl=None):
  if type(A) is Matrix:
    invNormMap = Matrix(Base(A.tag))
    if   A.tag == sTag: 
      if ctrl == None:
        lib.ElPseudospectralAutoWindow_s \
        (A.obj,invNormMap.obj,realSize,imagSize)
      else:
        lib.ElPseudospectralAutoWindowX_s \
        (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      lib.ElPseudospectralAutoWindowX_d \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      lib.ElPseudospectralAutoWindowX_c \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      lib.ElPseudospectralAutoWindowX_z \
      (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
    return invNormMap
  elif type(A) is DistMatrix:
    invNormMap = DistMatrix(Base(A.tag),MC,MR,A.Grid())
    if   A.tag == sTag: 
      if ctrl == None:
        lib.ElPseudospectralAutoWindowDist_s \
        (A.obj,invNormMap.obj,realSize,imagSize)
      else:
        lib.ElPseudospectralAutoWindowXDist_s \
        (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      if ctrl == None:
        lib.ElPseudospectralAutoWindowDist_d \
        (A.obj,invNormMap.obj,realSize,imagSize)
      else:
        lib.ElPseudospectralAutoWindowXDist_d \
        (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      if ctrl == None:
        lib.ElPseudospectralAutoWindowDist_c \
        (A.obj,invNormMap.obj,realSize,imagSize)
      else:
        lib.ElPseudospectralAutoWindowXDist_c \
        (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      if ctrl == None:
        lib.ElPseudospectralAutoWindowDist_z \
        (A.obj,invNormMap.obj,realSize,imagSize)
      else:
        lib.ElPseudospectralAutoWindowXDist_z \
        (A.obj,invNormMap.obj,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
    return invNormMap
  else: raise Exception('Unsupported matrix type')

# Investigate a particular window of the complex plane
# ----------------------------------------------------
lib.ElPseudospectralWindow_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType]
lib.ElPseudospectralWindow_s.restype = c_uint
lib.ElPseudospectralWindow_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType]
lib.ElPseudospectralWindow_d.restype = c_uint
lib.ElPseudospectralWindow_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType]
lib.ElPseudospectralWindow_c.restype = c_uint
lib.ElPseudospectralWindow_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType]
lib.ElPseudospectralWindow_z.restype = c_uint
lib.ElPseudospectralWindowDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType]
lib.ElPseudospectralWindowDist_s.restype = c_uint
lib.ElPseudospectralWindowDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType]
lib.ElPseudospectralWindowDist_d.restype = c_uint
lib.ElPseudospectralWindowDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType]
lib.ElPseudospectralWindowDist_c.restype = c_uint
lib.ElPseudospectralWindowDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType]
lib.ElPseudospectralWindowDist_z.restype = c_uint
lib.ElPseudospectralWindowX_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowX_s.restype = c_uint
lib.ElPseudospectralWindowX_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowX_d.restype = c_uint
lib.ElPseudospectralWindowX_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowX_c.restype = c_uint
lib.ElPseudospectralWindowX_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowX_z.restype = c_uint
lib.ElPseudospectralWindowXDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowXDist_s.restype = c_uint
lib.ElPseudospectralWindowXDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowXDist_d.restype = c_uint
lib.ElPseudospectralWindowXDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElPseudospectralWindowXDist_c.restype = c_uint
lib.ElPseudospectralWindowXDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElPseudospectralWindowXDist_z.restype = c_uint
def PseudospectralWindow \
    (A,centerPre,realWidth,imagWidth,realSize=200,imagSize=200,ctrl=None):
  center = TagToType(A.tag)(centerPre)
  if type(A) is Matrix:
    invNormMap = Matrix(Base(A.tag))
    if   A.tag == sTag: 
      if ctrl == None:
        lib.ElPseudospectralWindow_s \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowX_s \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      if ctrl == None:
        lib.ElPseudospectralWindow_d \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowX_d \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      if ctrl == None:
        lib.ElPseudospectralWindow_c \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowX_c \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      if ctrl == None:
        lib.ElPseudospectralWindow_z \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowX_z \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
    return invNormMap
  elif type(A) is DistMatrix:
    invNormMap = DistMatrix(Base(A.tag),MC,MR,A.Grid())
    if   A.tag == sTag: 
      if ctrl == None:
        lib.ElPseudospectralWindowDist_s \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowXDist_s \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == dTag:
      if ctrl == None:
        lib.ElPseudospectralWindowDist_d \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowXDist_d \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == cTag:
      if ctrl == None:
        lib.ElPseudospectralWindowDist_c \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowXDist_c \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    elif A.tag == zTag:
      if ctrl == None:
        lib.ElPseudospectralWindowDist_z \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize)
      else:
        lib.ElPseudospectralWindowXDist_z \
        (A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize,ctrl)
    else: raise Exception('Unsupported datatype')
    return invNormMap
  else: raise Exception('Unsupported matrix type')

# Investigate a specific set of points in the complex plane
# ---------------------------------------------------------
lib.ElPseudospectralCloud_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloud_s.restype = c_uint
lib.ElPseudospectralCloud_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloud_d.restype = c_uint
lib.ElPseudospectralCloud_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloud_c.restype = c_uint
lib.ElPseudospectralCloud_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloud_z.restype = c_uint
lib.ElPseudospectralCloudDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloudDist_s.restype = c_uint
lib.ElPseudospectralCloudDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloudDist_d.restype = c_uint
lib.ElPseudospectralCloudDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloudDist_c.restype = c_uint
lib.ElPseudospectralCloudDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPseudospectralCloudDist_z.restype = c_uint
lib.ElPseudospectralCloudX_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudX_s.restype = c_uint
lib.ElPseudospectralCloudX_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudX_d.restype = c_uint
lib.ElPseudospectralCloudX_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudX_c.restype = c_uint
lib.ElPseudospectralCloudX_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudX_z.restype = c_uint
lib.ElPseudospectralCloudXDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudXDist_s.restype = c_uint
lib.ElPseudospectralCloudXDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudXDist_d.restype = c_uint
lib.ElPseudospectralCloudXDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElPseudospectralCloudXDist_c.restype = c_uint
lib.ElPseudospectralCloudXDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElPseudospectralCloudXDist_z.restype = c_uint
def PseudospectralCloud(A,shifts,ctrl=None):
  if type(A) is Matrix:
    invNorms = Matrix(Base(A.tag)) 
    if   A.tag == sTag: 
      if ctrl == None:
        lib.ElPseudospectralCloud_s(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudX_s(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == dTag:
      if ctrl == None:
        lib.ElPseudospectralCloud_d(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudX_d(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == cTag:
      if ctrl == None:
        lib.ElPseudospectralCloud_c(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudX_c(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == zTag:
      if ctrl == None:
        lib.ElPseudospectralCloud_z(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudX_z(A.obj,shifts.obj,invNorms.obj,ctrl)
    else: raise Exception('Unsupported datatype')
    return invNorms
  elif type(A) is DistMatrix:
    invNorms = DistMatrix(Base(A.tag),VR,STAR,A.Grid())
    if   A.tag == sTag: 
      if ctrl == None:
        lib.ElPseudospectralCloudDist_s(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudXDist_s(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == dTag:
      if ctrl == None:
        lib.ElPseudospectralCloudDist_d(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudXDist_d(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == cTag:
      if ctrl == None:
        lib.ElPseudospectralCloudDist_c(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudXDist_c(A.obj,shifts.obj,invNorms.obj,ctrl)
    elif A.tag == zTag:
      if ctrl == None:
        lib.ElPseudospectralCloudDist_z(A.obj,shifts.obj,invNorms.obj)
      else:
        lib.ElPseudospectralCloudXDist_z(A.obj,shifts.obj,invNorms.obj,ctrl)
    else: raise Exception('Unsupported datatype')
    return invNorms
  else: raise Exception('Unsupported matrix type')

# Trace
# =====
# TODO
