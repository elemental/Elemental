#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
from factor import *
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
  args = [A.obj,normType,pointer(cond)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCondition_s(*args)
    elif A.tag == dTag: lib.ElCondition_d(*args)
    elif A.tag == cTag: lib.ElCondition_c(*args)
    elif A.tag == zTag: lib.ElCondition_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElConditionDist_s(*args)
    elif A.tag == dTag: lib.ElConditionDist_d(*args)
    elif A.tag == cTag: lib.ElConditionDist_c(*args)
    elif A.tag == zTag: lib.ElConditionDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
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
  args = [A.obj,pointer(cond)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElFrobeniusCondition_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusCondition_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusCondition_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusCondition_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElFrobeniusConditionDist_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusConditionDist_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusConditionDist_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusConditionDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
  args = [A.obj,pointer(cond)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElInfinityCondition_s(*args)
    elif A.tag == dTag: lib.ElInfinityCondition_d(*args)
    elif A.tag == cTag: lib.ElInfinityCondition_c(*args)
    elif A.tag == zTag: lib.ElInfinityCondition_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElInfinityConditionDist_s(*args)
    elif A.tag == dTag: lib.ElInfinityConditionDist_d(*args)
    elif A.tag == cTag: lib.ElInfinityConditionDist_c(*args)
    elif A.tag == zTag: lib.ElInfinityConditionDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
  args = [A.obj,pointer(cond)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElMaxCondition_s(*args)
    elif A.tag == dTag: lib.ElMaxCondition_d(*args)
    elif A.tag == cTag: lib.ElMaxCondition_c(*args)
    elif A.tag == zTag: lib.ElMaxCondition_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElMaxConditionDist_s(*args)
    elif A.tag == dTag: lib.ElMaxConditionDist_d(*args)
    elif A.tag == cTag: lib.ElMaxConditionDist_c(*args)
    elif A.tag == zTag: lib.ElMaxConditionDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
  args = [A.obj,pointer(cond)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElOneCondition_s(*args)
    elif A.tag == dTag: lib.ElOneCondition_d(*args)
    elif A.tag == cTag: lib.ElOneCondition_c(*args)
    elif A.tag == zTag: lib.ElOneCondition_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElOneConditionDist_s(*args)
    elif A.tag == dTag: lib.ElOneConditionDist_d(*args)
    elif A.tag == cTag: lib.ElOneConditionDist_c(*args)
    elif A.tag == zTag: lib.ElOneConditionDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
  args = [A.obj,pointer(cond)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoCondition_s(*args)
    elif A.tag == dTag: lib.ElTwoCondition_d(*args)
    elif A.tag == cTag: lib.ElTwoCondition_c(*args)
    elif A.tag == zTag: lib.ElTwoCondition_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoConditionDist_s(*args)
    elif A.tag == dTag: lib.ElTwoConditionDist_d(*args)
    elif A.tag == cTag: lib.ElTwoConditionDist_c(*args)
    elif A.tag == zTag: lib.ElTwoConditionDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
  safeProd = TagToSafeProduct(A.tag)
  args = [A.obj,pointer(safeProd)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSafeDeterminant_s(*args)
    elif A.tag == dTag: lib.ElSafeDeterminant_d(*args)
    elif A.tag == cTag: lib.ElSafeDeterminant_c(*args)
    elif A.tag == zTag: lib.ElSafeDeterminant_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSafeDeterminantDist_s(*args)
    elif A.tag == dTag: lib.ElSafeDeterminantDist_d(*args)
    elif A.tag == cTag: lib.ElSafeDeterminantDist_c(*args)
    elif A.tag == zTag: lib.ElSafeDeterminantDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return safeProd

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
  safeProd = TagToSafeProd(Base(A.tag))
  args = [uplo,A.obj,pointer(safeProd)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSafeHPDDeterminant_s(*args)
    elif A.tag == dTag: lib.ElSafeHPDDeterminant_d(*args)
    elif A.tag == cTag: lib.ElSafeHPDDeterminant_c(*args)
    elif A.tag == zTag: lib.ElSafeHPDDeterminant_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSafeHPDDeterminantDist_s(*args)
    elif A.tag == dTag: lib.ElSafeHPDDeterminantDist_d(*args)
    elif A.tag == cTag: lib.ElSafeHPDDeterminantDist_c(*args)
    elif A.tag == zTag: lib.ElSafeHPDDeterminantDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return safeProd

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
  args = [A.obj,pointer(prod)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElDeterminant_s(*args)
    elif A.tag == dTag: lib.ElDeterminant_d(*args)
    elif A.tag == cTag: lib.ElDeterminant_c(*args)
    elif A.tag == zTag: lib.ElDeterminant_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElDeterminantDist_s(*args)
    elif A.tag == dTag: lib.ElDeterminantDist_d(*args)
    elif A.tag == cTag: lib.ElDeterminantDist_c(*args)
    elif A.tag == zTag: lib.ElDeterminantDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
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
  args = [uplo,A.obj,pointer(prod)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPDDeterminant_s(*args)
    elif A.tag == dTag: lib.ElHPDDeterminant_d(*args)
    elif A.tag == cTag: lib.ElHPDDeterminant_c(*args)
    elif A.tag == zTag: lib.ElHPDDeterminant_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPDDeterminantDist_s(*args)
    elif A.tag == dTag: lib.ElHPDDeterminantDist_d(*args)
    elif A.tag == cTag: lib.ElHPDDeterminantDist_c(*args)
    elif A.tag == zTag: lib.ElHPDDeterminantDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return prod

# Inertia
# =======
lib.ElInertia_s.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertia_s.restype = c_uint
lib.ElInertia_d.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertia_d.restype = c_uint
lib.ElInertia_c.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertia_c.restype = c_uint
lib.ElInertia_z.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertia_z.restype = c_uint
lib.ElInertiaDist_s.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertiaDist_s.restype = c_uint
lib.ElInertiaDist_d.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertiaDist_d.restype = c_uint
lib.ElInertiaDist_c.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertiaDist_c.restype = c_uint
lib.ElInertiaDist_z.argtypes = [c_uint,c_void_p,c_uint,POINTER(InertiaType)]
lib.ElInertiaDist_z.restype = c_uint
def Inertia(uplo,A,pivType=BUNCH_KAUFMAN_A):
  inertia = InertiaType()
  args = [uplo,A.obj,pivType,pointer(inertia)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElInertia_s(*args)
    elif A.tag == dTag: lib.ElInertia_d(*args)
    elif A.tag == cTag: lib.ElInertia_c(*args)
    elif A.tag == zTag: lib.ElInertia_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElInertiaDist_s(*args)
    elif A.tag == dTag: lib.ElInertiaDist_d(*args)
    elif A.tag == cTag: lib.ElInertiaDist_c(*args)
    elif A.tag == zTag: lib.ElInertiaDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return inertia

# Norm
# ====
lib.ElNorm_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElNorm_s.restype = c_uint
lib.ElNorm_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElNorm_d.restype = c_uint
lib.ElNorm_c.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElNorm_c.restype = c_uint
lib.ElNorm_z.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElNorm_z.restype = c_uint
lib.ElNormDist_s.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElNormDist_s.restype = c_uint
lib.ElNormDist_d.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElNormDist_d.restype = c_uint
lib.ElNormDist_c.argtypes = [c_void_p,c_uint,POINTER(sType)]
lib.ElNormDist_c.restype = c_uint
lib.ElNormDist_z.argtypes = [c_void_p,c_uint,POINTER(dType)]
lib.ElNormDist_z.restype = c_uint
def Norm(A,normType=FROBENIUS_NORM):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,normType,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElNorm_s(*args)
    elif A.tag == dTag: lib.ElNorm_d(*args)
    elif A.tag == cTag: lib.ElNorm_c(*args)
    elif A.tag == zTag: lib.ElNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElNormDist_s(*args)
    elif A.tag == dTag: lib.ElNormDist_d(*args)
    elif A.tag == cTag: lib.ElNormDist_c(*args)
    elif A.tag == zTag: lib.ElNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricNorm_s.argtypes = [c_uint,c_void_p,c_uint,POINTER(sType)]
lib.ElSymmetricNorm_s.restype = c_uint
lib.ElSymmetricNorm_d.argtypes = [c_uint,c_void_p,c_uint,POINTER(dType)]
lib.ElSymmetricNorm_d.restype = c_uint
lib.ElSymmetricNorm_c.argtypes = [c_uint,c_void_p,c_uint,POINTER(sType)]
lib.ElSymmetricNorm_c.restype = c_uint
lib.ElSymmetricNorm_z.argtypes = [c_uint,c_void_p,c_uint,POINTER(dType)]
lib.ElSymmetricNorm_z.restype = c_uint
lib.ElSymmetricNormDist_s.argtypes = [c_uint,c_void_p,c_uint,POINTER(sType)]
lib.ElSymmetricNormDist_s.restype = c_uint
lib.ElSymmetricNormDist_d.argtypes = [c_uint,c_void_p,c_uint,POINTER(dType)]
lib.ElSymmetricNormDist_d.restype = c_uint
lib.ElSymmetricNormDist_c.argtypes = [c_uint,c_void_p,c_uint,POINTER(sType)]
lib.ElSymmetricNormDist_c.restype = c_uint
lib.ElSymmetricNormDist_z.argtypes = [c_uint,c_void_p,c_uint,POINTER(dType)]
lib.ElSymmetricNormDist_z.restype = c_uint
lib.ElHermitianNorm_c.argtypes = [c_uint,c_void_p,c_uint,POINTER(sType)]
lib.ElHermitianNorm_c.restype = c_uint
lib.ElHermitianNorm_z.argtypes = [c_uint,c_void_p,c_uint,POINTER(dType)]
lib.ElHermitianNorm_z.restype = c_uint
lib.ElHermitianNormDist_c.argtypes = [c_uint,c_void_p,c_uint,POINTER(sType)]
lib.ElHermitianNormDist_c.restype = c_uint
lib.ElHermitianNormDist_z.argtypes = [c_uint,c_void_p,c_uint,POINTER(dType)]
lib.ElHermitianNormDist_z.restype = c_uint
def SymmetricNorm(uplo,A,normType=FROBENIUS_NORM,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,normType,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianNorm_c(*args)
      else:         lib.ElSymmetricNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianNorm_z(*args)
      else:         lib.ElSymmetricNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianNormDist_c(*args)
      else:         lib.ElSymmetricNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianNormDist_z(*args)
      else:         lib.ElSymmetricNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianNorm(uplo,A,normType=FROBENIUS_NORM):
  return SymmetricNorm(uplo,A,normType,True)

lib.ElEntrywiseNorm_s.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElEntrywiseNorm_s.restype = c_uint
lib.ElEntrywiseNorm_d.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElEntrywiseNorm_d.restype = c_uint
lib.ElEntrywiseNorm_c.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElEntrywiseNorm_c.restype = c_uint
lib.ElEntrywiseNorm_z.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElEntrywiseNorm_z.restype = c_uint
lib.ElEntrywiseNormDist_s.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElEntrywiseNormDist_s.restype = c_uint
lib.ElEntrywiseNormDist_d.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElEntrywiseNormDist_d.restype = c_uint
lib.ElEntrywiseNormDist_c.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElEntrywiseNormDist_c.restype = c_uint
lib.ElEntrywiseNormDist_z.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElEntrywiseNormDist_z.restype = c_uint
def EntrywiseNorm(A,p):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,p,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElEntrywiseNorm_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseNorm_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseNorm_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElEntrywiseNormDist_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseNormDist_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseNormDist_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricEntrywiseNorm_s.argtypes = [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricEntrywiseNorm_s.restype = c_uint
lib.ElSymmetricEntrywiseNorm_d.argtypes = [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricEntrywiseNorm_d.restype = c_uint
lib.ElSymmetricEntrywiseNorm_c.argtypes = [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricEntrywiseNorm_c.restype = c_uint
lib.ElSymmetricEntrywiseNorm_z.argtypes = [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricEntrywiseNorm_z.restype = c_uint
lib.ElSymmetricEntrywiseNormDist_s.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricEntrywiseNormDist_s.restype = c_uint
lib.ElSymmetricEntrywiseNormDist_d.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricEntrywiseNormDist_d.restype = c_uint
lib.ElSymmetricEntrywiseNormDist_c.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricEntrywiseNormDist_c.restype = c_uint
lib.ElSymmetricEntrywiseNormDist_z.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricEntrywiseNormDist_z.restype = c_uint
lib.ElHermitianEntrywiseNorm_c.argtypes = [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElHermitianEntrywiseNorm_c.restype = c_uint
lib.ElHermitianEntrywiseNorm_z.argtypes = [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElHermitianEntrywiseNorm_z.restype = c_uint
lib.ElHermitianEntrywiseNormDist_c.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElHermitianEntrywiseNormDist_c.restype = c_uint
lib.ElHermitianEntrywiseNormDist_z.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElHermitianEntrywiseNormDist_z.restype = c_uint
def SymmetricEntrywiseNorm(uplo,A,p,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,p,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianEntrywiseNorm_c(*args)
      else:         lib.ElSymmetricEntrywiseNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianEntrywiseNorm_z(*args)
      else:         lib.ElSymmetricEntrywiseNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianEntrywiseNormDist_c(*args)
      else:         lib.ElSymmetricEntrywiseNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianEntrywiseNormDist_z(*args)
      else:         lib.ElSymmetricEntrywiseNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianEntrywiseNorm(uplo,A,p):
  return SymmetricEntrywiseNorm(uplo,A,p,True)

lib.ElEntrywiseOneNorm_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElEntrywiseOneNorm_s.restype = c_uint
lib.ElEntrywiseOneNorm_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElEntrywiseOneNorm_d.restype = c_uint
lib.ElEntrywiseOneNorm_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElEntrywiseOneNorm_c.restype = c_uint
lib.ElEntrywiseOneNorm_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElEntrywiseOneNorm_z.restype = c_uint
lib.ElEntrywiseOneNormDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElEntrywiseOneNormDist_s.restype = c_uint
lib.ElEntrywiseOneNormDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElEntrywiseOneNormDist_d.restype = c_uint
lib.ElEntrywiseOneNormDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElEntrywiseOneNormDist_c.restype = c_uint
lib.ElEntrywiseOneNormDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElEntrywiseOneNormDist_z.restype = c_uint
def EntrywiseOneNorm(A):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElEntrywiseOneNorm_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseOneNorm_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseOneNorm_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseOneNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElEntrywiseOneNormDist_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseOneNormDist_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseOneNormDist_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseOneNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricEntrywiseOneNorm_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricEntrywiseOneNorm_s.restype = c_uint
lib.ElSymmetricEntrywiseOneNorm_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricEntrywiseOneNorm_d.restype = c_uint
lib.ElSymmetricEntrywiseOneNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricEntrywiseOneNorm_c.restype = c_uint
lib.ElSymmetricEntrywiseOneNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricEntrywiseOneNorm_z.restype = c_uint
lib.ElSymmetricEntrywiseOneNormDist_s.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricEntrywiseOneNormDist_s.restype = c_uint
lib.ElSymmetricEntrywiseOneNormDist_d.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricEntrywiseOneNormDist_d.restype = c_uint
lib.ElSymmetricEntrywiseOneNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricEntrywiseOneNormDist_c.restype = c_uint
lib.ElSymmetricEntrywiseOneNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricEntrywiseOneNormDist_z.restype = c_uint
lib.ElHermitianEntrywiseOneNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianEntrywiseOneNorm_c.restype = c_uint
lib.ElHermitianEntrywiseOneNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianEntrywiseOneNorm_z.restype = c_uint
lib.ElHermitianEntrywiseOneNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianEntrywiseOneNormDist_c.restype = c_uint
lib.ElHermitianEntrywiseOneNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianEntrywiseOneNormDist_z.restype = c_uint
def SymmetricEntrywiseOneNorm(uplo,A,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseOneNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseOneNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianEntrywiseOneNorm_c(*args)
      else:         lib.ElSymmetricEntrywiseOneNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianEntrywiseOneNorm_z(*args)
      else:         lib.ElSymmetricEntrywiseOneNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseOneNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseOneNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianEntrywiseOneNormDist_c(*args)
      else:         lib.ElSymmetricEntrywiseOneNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianEntrywiseOneNormDist_z(*args)
      else:         lib.ElSymmetricEntrywiseOneNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianEntrywiseOneNorm(uplo,A):
  return SymmetricEntrywiseOneNorm(uplo,A,True)

lib.ElFrobeniusNorm_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusNorm_s.restype = c_uint
lib.ElFrobeniusNorm_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusNorm_d.restype = c_uint
lib.ElFrobeniusNorm_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusNorm_c.restype = c_uint
lib.ElFrobeniusNorm_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusNorm_z.restype = c_uint
lib.ElFrobeniusNormDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusNormDist_s.restype = c_uint
lib.ElFrobeniusNormDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusNormDist_d.restype = c_uint
lib.ElFrobeniusNormDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElFrobeniusNormDist_c.restype = c_uint
lib.ElFrobeniusNormDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElFrobeniusNormDist_z.restype = c_uint
def FrobeniusNorm(A):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElFrobeniusNorm_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusNorm_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusNorm_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElFrobeniusNormDist_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusNormDist_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusNormDist_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricFrobeniusNorm_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricFrobeniusNorm_s.restype = c_uint
lib.ElSymmetricFrobeniusNorm_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricFrobeniusNorm_d.restype = c_uint
lib.ElSymmetricFrobeniusNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricFrobeniusNorm_c.restype = c_uint
lib.ElSymmetricFrobeniusNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricFrobeniusNorm_z.restype = c_uint
lib.ElSymmetricFrobeniusNormDist_s.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricFrobeniusNormDist_s.restype = c_uint
lib.ElSymmetricFrobeniusNormDist_d.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricFrobeniusNormDist_d.restype = c_uint
lib.ElSymmetricFrobeniusNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricFrobeniusNormDist_c.restype = c_uint
lib.ElSymmetricFrobeniusNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricFrobeniusNormDist_z.restype = c_uint
lib.ElHermitianFrobeniusNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianFrobeniusNorm_c.restype = c_uint
lib.ElHermitianFrobeniusNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianFrobeniusNorm_z.restype = c_uint
lib.ElHermitianFrobeniusNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianFrobeniusNormDist_c.restype = c_uint
lib.ElHermitianFrobeniusNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianFrobeniusNormDist_z.restype = c_uint
def SymmetricFrobeniusNorm(uplo,A,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricFrobeniusNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricFrobeniusNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianFrobeniusNorm_c(*args)
      else:         lib.ElSymmetricFrobeniusNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianFrobeniusNorm_z(*args)
      else:         lib.ElSymmetricFrobeniusNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricFrobeniusNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricFrobeniusNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianFrobeniusNormDist_c(*args)
      else:         lib.ElSymmetricFrobeniusNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianFrobeniusNormDist_z(*args)
      else:         lib.ElSymmetricFrobeniusNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianFrobeniusNorm(uplo,A):
  return SymmetricFrobeniusNorm(uplo,A,True)

lib.ElInfinityNorm_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityNorm_s.restype = c_uint
lib.ElInfinityNorm_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityNorm_d.restype = c_uint
lib.ElInfinityNorm_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityNorm_c.restype = c_uint
lib.ElInfinityNorm_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityNorm_z.restype = c_uint
lib.ElInfinityNormDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityNormDist_s.restype = c_uint
lib.ElInfinityNormDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityNormDist_d.restype = c_uint
lib.ElInfinityNormDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElInfinityNormDist_c.restype = c_uint
lib.ElInfinityNormDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElInfinityNormDist_z.restype = c_uint
def InfinityNorm(A):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElInfinityNorm_s(*args)
    elif A.tag == dTag: lib.ElInfinityNorm_d(*args)
    elif A.tag == cTag: lib.ElInfinityNorm_c(*args)
    elif A.tag == zTag: lib.ElInfinityNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElInfinityNormDist_s(*args)
    elif A.tag == dTag: lib.ElInfinityNormDist_d(*args)
    elif A.tag == cTag: lib.ElInfinityNormDist_c(*args)
    elif A.tag == zTag: lib.ElInfinityNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricInfinityNorm_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricInfinityNorm_s.restype = c_uint
lib.ElSymmetricInfinityNorm_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricInfinityNorm_d.restype = c_uint
lib.ElSymmetricInfinityNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricInfinityNorm_c.restype = c_uint
lib.ElSymmetricInfinityNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricInfinityNorm_z.restype = c_uint
lib.ElSymmetricInfinityNormDist_s.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricInfinityNormDist_s.restype = c_uint
lib.ElSymmetricInfinityNormDist_d.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricInfinityNormDist_d.restype = c_uint
lib.ElSymmetricInfinityNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricInfinityNormDist_c.restype = c_uint
lib.ElSymmetricInfinityNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricInfinityNormDist_z.restype = c_uint
lib.ElHermitianInfinityNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianInfinityNorm_c.restype = c_uint
lib.ElHermitianInfinityNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianInfinityNorm_z.restype = c_uint
lib.ElHermitianInfinityNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianInfinityNormDist_c.restype = c_uint
lib.ElHermitianInfinityNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianInfinityNormDist_z.restype = c_uint
def SymmetricInfinityNorm(uplo,A,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricInfinityNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricInfinityNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianInfinityNorm_c(*args)
      else:         lib.ElSymmetricInfinityNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianInfinityNorm_z(*args)
      else:         lib.ElSymmetricInfinityNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricInfinityNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricInfinityNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianInfinityNormDist_c(*args)
      else:         lib.ElSymmetricInfinityNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianInfinityNormDist_z(*args)
      else:         lib.ElSymmetricInfinityNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianInfinityNorm(uplo,A):
  return SymmetricInfinityNorm(uplo,A,True)

lib.ElKyFanNorm_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElKyFanNorm_s.restype = c_uint
lib.ElKyFanNorm_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElKyFanNorm_d.restype = c_uint
lib.ElKyFanNorm_c.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElKyFanNorm_c.restype = c_uint
lib.ElKyFanNorm_z.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElKyFanNorm_z.restype = c_uint
lib.ElKyFanNormDist_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElKyFanNormDist_s.restype = c_uint
lib.ElKyFanNormDist_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElKyFanNormDist_d.restype = c_uint
lib.ElKyFanNormDist_c.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElKyFanNormDist_c.restype = c_uint
lib.ElKyFanNormDist_z.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElKyFanNormDist_z.restype = c_uint
def KyFanNorm(A,k):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,k,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElKyFanNorm_s(*args)
    elif A.tag == dTag: lib.ElKyFanNorm_d(*args)
    elif A.tag == cTag: lib.ElKyFanNorm_c(*args)
    elif A.tag == zTag: lib.ElKyFanNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElKyFanNormDist_s(*args)
    elif A.tag == dTag: lib.ElKyFanNormDist_d(*args)
    elif A.tag == cTag: lib.ElKyFanNormDist_c(*args)
    elif A.tag == zTag: lib.ElKyFanNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricKyFanNorm_s.argtypes = [c_uint,c_void_p,iType,POINTER(sType)]
lib.ElSymmetricKyFanNorm_s.restype = c_uint
lib.ElSymmetricKyFanNorm_d.argtypes = [c_uint,c_void_p,iType,POINTER(dType)]
lib.ElSymmetricKyFanNorm_d.restype = c_uint
lib.ElSymmetricKyFanNorm_c.argtypes = [c_uint,c_void_p,iType,POINTER(sType)]
lib.ElSymmetricKyFanNorm_c.restype = c_uint
lib.ElSymmetricKyFanNorm_z.argtypes = [c_uint,c_void_p,iType,POINTER(dType)]
lib.ElSymmetricKyFanNorm_z.restype = c_uint
lib.ElSymmetricKyFanNormDist_s.argtypes = \
  [c_uint,c_void_p,iType,POINTER(sType)]
lib.ElSymmetricKyFanNormDist_s.restype = c_uint
lib.ElSymmetricKyFanNormDist_d.argtypes = \
  [c_uint,c_void_p,iType,POINTER(dType)]
lib.ElSymmetricKyFanNormDist_d.restype = c_uint
lib.ElSymmetricKyFanNormDist_c.argtypes = \
  [c_uint,c_void_p,iType,POINTER(sType)]
lib.ElSymmetricKyFanNormDist_c.restype = c_uint
lib.ElSymmetricKyFanNormDist_z.argtypes = \
  [c_uint,c_void_p,iType,POINTER(dType)]
lib.ElSymmetricKyFanNormDist_z.restype = c_uint
lib.ElHermitianKyFanNorm_c.argtypes = [c_uint,c_void_p,iType,POINTER(sType)]
lib.ElHermitianKyFanNorm_c.restype = c_uint
lib.ElHermitianKyFanNorm_z.argtypes = [c_uint,c_void_p,iType,POINTER(dType)]
lib.ElHermitianKyFanNorm_z.restype = c_uint
lib.ElHermitianKyFanNormDist_c.argtypes = \
  [c_uint,c_void_p,iType,POINTER(sType)]
lib.ElHermitianKyFanNormDist_c.restype = c_uint
lib.ElHermitianKyFanNormDist_z.argtypes = \
  [c_uint,c_void_p,iType,POINTER(dType)]
lib.ElHermitianKyFanNormDist_z.restype = c_uint
def SymmetricKyFanNorm(uplo,A,k,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,k,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricKyFanNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricKyFanNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianKyFanNorm_c(*args)
      else:         lib.ElSymmetricKyFanNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianKyFanNorm_z(*args)
      else:         lib.ElSymmetricKyFanNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricKyFanNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricKyFanNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianKyFanNormDist_c(*args)
      else:         lib.ElSymmetricKyFanNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianKyFanNormDist_z(*args)
      else:         lib.ElSymmetricKyFanNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianKyFanNorm(uplo,A,k):
  return SymmetricKyFanNorm(uplo,A,k,True)

lib.ElKyFanSchattenNorm_s.argtypes = [c_void_p,iType,sType,POINTER(sType)]
lib.ElKyFanSchattenNorm_s.restype = c_uint
lib.ElKyFanSchattenNorm_d.argtypes = [c_void_p,iType,dType,POINTER(dType)]
lib.ElKyFanSchattenNorm_d.restype = c_uint
lib.ElKyFanSchattenNorm_c.argtypes = [c_void_p,iType,sType,POINTER(sType)]
lib.ElKyFanSchattenNorm_c.restype = c_uint
lib.ElKyFanSchattenNorm_z.argtypes = [c_void_p,iType,dType,POINTER(dType)]
lib.ElKyFanSchattenNorm_z.restype = c_uint
lib.ElKyFanSchattenNormDist_s.argtypes = [c_void_p,iType,sType,POINTER(sType)]
lib.ElKyFanSchattenNormDist_s.restype = c_uint
lib.ElKyFanSchattenNormDist_d.argtypes = [c_void_p,iType,dType,POINTER(dType)]
lib.ElKyFanSchattenNormDist_d.restype = c_uint
lib.ElKyFanSchattenNormDist_c.argtypes = [c_void_p,iType,sType,POINTER(sType)]
lib.ElKyFanSchattenNormDist_c.restype = c_uint
lib.ElKyFanSchattenNormDist_z.argtypes = [c_void_p,iType,dType,POINTER(dType)]
lib.ElKyFanSchattenNormDist_z.restype = c_uint
def KyFanSchattenNorm(A,k,p):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,k,p,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElKyFanSchattenNorm_s(*args)
    elif A.tag == dTag: lib.ElKyFanSchattenNorm_d(*args)
    elif A.tag == cTag: lib.ElKyFanSchattenNorm_c(*args)
    elif A.tag == zTag: lib.ElKyFanSchattenNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElKyFanSchattenNormDist_s(*args)
    elif A.tag == dTag: lib.ElKyFanSchattenNormDist_d(*args)
    elif A.tag == cTag: lib.ElKyFanSchattenNormDist_c(*args)
    elif A.tag == zTag: lib.ElKyFanSchattenNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricKyFanSchattenNorm_s.argtypes = \
  [c_uint,c_void_p,iType,sType,POINTER(sType)]
lib.ElSymmetricKyFanSchattenNorm_s.restype = c_uint
lib.ElSymmetricKyFanSchattenNorm_d.argtypes = \
  [c_uint,c_void_p,iType,dType,POINTER(dType)]
lib.ElSymmetricKyFanSchattenNorm_d.restype = c_uint
lib.ElSymmetricKyFanSchattenNorm_c.argtypes = \
  [c_uint,c_void_p,iType,sType,POINTER(sType)]
lib.ElSymmetricKyFanSchattenNorm_c.restype = c_uint
lib.ElSymmetricKyFanSchattenNorm_z.argtypes = \
  [c_uint,c_void_p,iType,dType,POINTER(dType)]
lib.ElSymmetricKyFanSchattenNorm_z.restype = c_uint
lib.ElSymmetricKyFanSchattenNormDist_s.argtypes = \
  [c_uint,c_void_p,iType,sType,POINTER(sType)]
lib.ElSymmetricKyFanSchattenNormDist_s.restype = c_uint
lib.ElSymmetricKyFanSchattenNormDist_d.argtypes = \
  [c_uint,c_void_p,iType,dType,POINTER(dType)]
lib.ElSymmetricKyFanSchattenNormDist_d.restype = c_uint
lib.ElSymmetricKyFanSchattenNormDist_c.argtypes = \
  [c_uint,c_void_p,iType,sType,POINTER(sType)]
lib.ElSymmetricKyFanSchattenNormDist_c.restype = c_uint
lib.ElSymmetricKyFanSchattenNormDist_z.argtypes = \
  [c_uint,c_void_p,iType,dType,POINTER(dType)]
lib.ElSymmetricKyFanSchattenNormDist_z.restype = c_uint
lib.ElHermitianKyFanSchattenNorm_c.argtypes = \
  [c_uint,c_void_p,iType,sType,POINTER(sType)]
lib.ElHermitianKyFanSchattenNorm_c.restype = c_uint
lib.ElHermitianKyFanSchattenNorm_z.argtypes = \
  [c_uint,c_void_p,iType,dType,POINTER(dType)]
lib.ElHermitianKyFanSchattenNorm_z.restype = c_uint
lib.ElHermitianKyFanSchattenNormDist_c.argtypes = \
  [c_uint,c_void_p,iType,sType,POINTER(sType)]
lib.ElHermitianKyFanSchattenNormDist_c.restype = c_uint
lib.ElHermitianKyFanSchattenNormDist_z.argtypes = \
  [c_uint,c_void_p,iType,dType,POINTER(dType)]
lib.ElHermitianKyFanSchattenNormDist_z.restype = c_uint
def SymmetricKyFanSchattenNorm(uplo,A,k,p,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,k,p,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricKyFanSchattenNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricKyFanSchattenNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianKyFanSchattenNorm_c(*args)
      else:         lib.ElSymmetricKyFanSchattenNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianKyFanSchattenNorm_z(*args)
      else:         lib.ElSymmetricKyFanSchattenNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricKyFanSchattenNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricKyFanSchattenNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianKyFanSchattenNormDist_c(*args)
      else:         lib.ElSymmetricKyFanSchattenNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianKyFanSchattenNormDist_z(*args)
      else:         lib.ElSymmetricKyFanSchattenNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianKyFanSchattenNorm(uplo,A,k,p):
  return SymmetricKyFanSchattenNorm(uplo,A,k,p,True)

lib.ElMaxNorm_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMaxNorm_i.restype = c_uint
lib.ElMaxNorm_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxNorm_s.restype = c_uint
lib.ElMaxNorm_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxNorm_d.restype = c_uint
lib.ElMaxNorm_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxNorm_c.restype = c_uint
lib.ElMaxNorm_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxNorm_z.restype = c_uint
lib.ElMaxNormDist_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMaxNormDist_i.restype = c_uint
lib.ElMaxNormDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxNormDist_s.restype = c_uint
lib.ElMaxNormDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxNormDist_d.restype = c_uint
lib.ElMaxNormDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElMaxNormDist_c.restype = c_uint
lib.ElMaxNormDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElMaxNormDist_z.restype = c_uint
def MaxNorm(A):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMaxNorm_i(*args)
    elif A.tag == sTag: lib.ElMaxNorm_s(*args)
    elif A.tag == dTag: lib.ElMaxNorm_d(*args)
    elif A.tag == cTag: lib.ElMaxNorm_c(*args)
    elif A.tag == zTag: lib.ElMaxNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMaxNormDist_i(*args)
    elif A.tag == sTag: lib.ElMaxNormDist_s(*args)
    elif A.tag == dTag: lib.ElMaxNormDist_d(*args)
    elif A.tag == cTag: lib.ElMaxNormDist_c(*args)
    elif A.tag == zTag: lib.ElMaxNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricMaxNorm_i.argtypes = [c_uint,c_void_p,POINTER(iType)]
lib.ElSymmetricMaxNorm_i.restype = c_uint
lib.ElSymmetricMaxNorm_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricMaxNorm_s.restype = c_uint
lib.ElSymmetricMaxNorm_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricMaxNorm_d.restype = c_uint
lib.ElSymmetricMaxNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricMaxNorm_c.restype = c_uint
lib.ElSymmetricMaxNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricMaxNorm_z.restype = c_uint
lib.ElSymmetricMaxNormDist_i.argtypes = [c_uint,c_void_p,POINTER(iType)]
lib.ElSymmetricMaxNormDist_i.restype = c_uint
lib.ElSymmetricMaxNormDist_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricMaxNormDist_s.restype = c_uint
lib.ElSymmetricMaxNormDist_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricMaxNormDist_d.restype = c_uint
lib.ElSymmetricMaxNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricMaxNormDist_c.restype = c_uint
lib.ElSymmetricMaxNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricMaxNormDist_z.restype = c_uint
lib.ElHermitianMaxNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianMaxNorm_c.restype = c_uint
lib.ElHermitianMaxNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianMaxNorm_z.restype = c_uint
lib.ElHermitianMaxNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianMaxNormDist_c.restype = c_uint
lib.ElHermitianMaxNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianMaxNormDist_z.restype = c_uint
def SymmetricMaxNorm(uplo,A,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSymmetricMaxNorm_i(*args)
    elif A.tag == sTag: lib.ElSymmetricMaxNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricMaxNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianMaxNorm_c(*args)
      else:         lib.ElSymmetricMaxNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianMaxNorm_z(*args)
      else:         lib.ElSymmetricMaxNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSymmetricMaxNormDist_i(*args)
    elif A.tag == sTag: lib.ElSymmetricMaxNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricMaxNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianMaxNormDist_c(*args)
      else:         lib.ElSymmetricMaxNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianMaxNormDist_z(*args)
      else:         lib.ElSymmetricMaxNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianMaxNorm(uplo,A):
  return SymmetricMaxNorm(uplo,A,True)

lib.ElNuclearNorm_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElNuclearNorm_s.restype = c_uint
lib.ElNuclearNorm_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElNuclearNorm_d.restype = c_uint
lib.ElNuclearNorm_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElNuclearNorm_c.restype = c_uint
lib.ElNuclearNorm_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElNuclearNorm_z.restype = c_uint
lib.ElNuclearNormDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElNuclearNormDist_s.restype = c_uint
lib.ElNuclearNormDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElNuclearNormDist_d.restype = c_uint
lib.ElNuclearNormDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElNuclearNormDist_c.restype = c_uint
lib.ElNuclearNormDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElNuclearNormDist_z.restype = c_uint
def NuclearNorm(A):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElNuclearNorm_s(*args)
    elif A.tag == dTag: lib.ElNuclearNorm_d(*args)
    elif A.tag == cTag: lib.ElNuclearNorm_c(*args)
    elif A.tag == zTag: lib.ElNuclearNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElNuclearNormDist_s(*args)
    elif A.tag == dTag: lib.ElNuclearNormDist_d(*args)
    elif A.tag == cTag: lib.ElNuclearNormDist_c(*args)
    elif A.tag == zTag: lib.ElNuclearNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricNuclearNorm_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricNuclearNorm_s.restype = c_uint
lib.ElSymmetricNuclearNorm_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricNuclearNorm_d.restype = c_uint
lib.ElSymmetricNuclearNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricNuclearNorm_c.restype = c_uint
lib.ElSymmetricNuclearNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricNuclearNorm_z.restype = c_uint
lib.ElSymmetricNuclearNormDist_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricNuclearNormDist_s.restype = c_uint
lib.ElSymmetricNuclearNormDist_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricNuclearNormDist_d.restype = c_uint
lib.ElSymmetricNuclearNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricNuclearNormDist_c.restype = c_uint
lib.ElSymmetricNuclearNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricNuclearNormDist_z.restype = c_uint
lib.ElHermitianNuclearNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianNuclearNorm_c.restype = c_uint
lib.ElHermitianNuclearNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianNuclearNorm_z.restype = c_uint
lib.ElHermitianNuclearNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianNuclearNormDist_c.restype = c_uint
lib.ElHermitianNuclearNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianNuclearNormDist_z.restype = c_uint
def SymmetricNuclearNorm(uplo,A,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricNuclearNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricNuclearNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianNuclearNorm_c(*args)
      else:         lib.ElSymmetricNuclearNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianNuclearNorm_z(*args)
      else:         lib.ElSymmetricNuclearNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricNuclearNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricNuclearNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianNuclearNormDist_c(*args)
      else:         lib.ElSymmetricNuclearNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianNuclearNormDist_z(*args)
      else:         lib.ElSymmetricNuclearNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianNuclearNorm(uplo,A):
  return SymmetricNuclearNorm(uplo,A,True)

lib.ElOneNorm_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneNorm_s.restype = c_uint
lib.ElOneNorm_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneNorm_d.restype = c_uint
lib.ElOneNorm_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneNorm_c.restype = c_uint
lib.ElOneNorm_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneNorm_z.restype = c_uint
lib.ElOneNormDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneNormDist_s.restype = c_uint
lib.ElOneNormDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneNormDist_d.restype = c_uint
lib.ElOneNormDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElOneNormDist_c.restype = c_uint
lib.ElOneNormDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElOneNormDist_z.restype = c_uint
def OneNorm(A):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElOneNorm_s(*args)
    elif A.tag == dTag: lib.ElOneNorm_d(*args)
    elif A.tag == cTag: lib.ElOneNorm_c(*args)
    elif A.tag == zTag: lib.ElOneNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElOneNormDist_s(*args)
    elif A.tag == dTag: lib.ElOneNormDist_d(*args)
    elif A.tag == cTag: lib.ElOneNormDist_c(*args)
    elif A.tag == zTag: lib.ElOneNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricOneNorm_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricOneNorm_s.restype = c_uint
lib.ElSymmetricOneNorm_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricOneNorm_d.restype = c_uint
lib.ElSymmetricOneNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricOneNorm_c.restype = c_uint
lib.ElSymmetricOneNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricOneNorm_z.restype = c_uint
lib.ElSymmetricOneNormDist_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricOneNormDist_s.restype = c_uint
lib.ElSymmetricOneNormDist_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricOneNormDist_d.restype = c_uint
lib.ElSymmetricOneNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricOneNormDist_c.restype = c_uint
lib.ElSymmetricOneNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricOneNormDist_z.restype = c_uint
lib.ElHermitianOneNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianOneNorm_c.restype = c_uint
lib.ElHermitianOneNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianOneNorm_z.restype = c_uint
lib.ElHermitianOneNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianOneNormDist_c.restype = c_uint
lib.ElHermitianOneNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianOneNormDist_z.restype = c_uint
def SymmetricOneNorm(uplo,A,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricOneNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricOneNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianOneNorm_c(*args)
      else:         lib.ElSymmetricOneNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianOneNorm_z(*args)
      else:         lib.ElSymmetricOneNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricOneNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricOneNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianOneNormDist_c(*args)
      else:         lib.ElSymmetricOneNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianOneNormDist_z(*args)
      else:         lib.ElSymmetricOneNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianOneNorm(uplo,A):
  return SymmetricOneNorm(uplo,A,True)

lib.ElSchattenNorm_s.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElSchattenNorm_s.restype = c_uint
lib.ElSchattenNorm_d.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElSchattenNorm_d.restype = c_uint
lib.ElSchattenNorm_c.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElSchattenNorm_c.restype = c_uint
lib.ElSchattenNorm_z.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElSchattenNorm_z.restype = c_uint
lib.ElSchattenNormDist_s.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElSchattenNormDist_s.restype = c_uint
lib.ElSchattenNormDist_d.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElSchattenNormDist_d.restype = c_uint
lib.ElSchattenNormDist_c.argtypes = [c_void_p,sType,POINTER(sType)]
lib.ElSchattenNormDist_c.restype = c_uint
lib.ElSchattenNormDist_z.argtypes = [c_void_p,dType,POINTER(dType)]
lib.ElSchattenNormDist_z.restype = c_uint
def SchattenNorm(A,p):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,p,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSchattenNorm_s(*args)
    elif A.tag == dTag: lib.ElSchattenNorm_d(*args)
    elif A.tag == cTag: lib.ElSchattenNorm_c(*args)
    elif A.tag == zTag: lib.ElSchattenNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSchattenNormDist_s(*args)
    elif A.tag == dTag: lib.ElSchattenNormDist_d(*args)
    elif A.tag == cTag: lib.ElSchattenNormDist_c(*args)
    elif A.tag == zTag: lib.ElSchattenNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricSchattenNorm_s.argtypes = [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricSchattenNorm_s.restype = c_uint
lib.ElSymmetricSchattenNorm_d.argtypes = [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricSchattenNorm_d.restype = c_uint
lib.ElSymmetricSchattenNorm_c.argtypes = [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricSchattenNorm_c.restype = c_uint
lib.ElSymmetricSchattenNorm_z.argtypes = [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricSchattenNorm_z.restype = c_uint
lib.ElSymmetricSchattenNormDist_s.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricSchattenNormDist_s.restype = c_uint
lib.ElSymmetricSchattenNormDist_d.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricSchattenNormDist_d.restype = c_uint
lib.ElSymmetricSchattenNormDist_c.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElSymmetricSchattenNormDist_c.restype = c_uint
lib.ElSymmetricSchattenNormDist_z.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElSymmetricSchattenNormDist_z.restype = c_uint
lib.ElHermitianSchattenNorm_c.argtypes = [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElHermitianSchattenNorm_c.restype = c_uint
lib.ElHermitianSchattenNorm_z.argtypes = [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElHermitianSchattenNorm_z.restype = c_uint
lib.ElHermitianSchattenNormDist_c.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]
lib.ElHermitianSchattenNormDist_c.restype = c_uint
lib.ElHermitianSchattenNormDist_z.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]
lib.ElHermitianSchattenNormDist_z.restype = c_uint
def SymmetricSchattenNorm(uplo,A,p,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,p,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricSchattenNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSchattenNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianSchattenNorm_c(*args)
      else:         lib.ElSymmetricSchattenNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianSchattenNorm_z(*args)
      else:         lib.ElSymmetricSchattenNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricSchattenNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSchattenNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianSchattenNormDist_c(*args)
      else:         lib.ElSymmetricSchattenNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianSchattenNormDist_z(*args)
      else:         lib.ElSymmetricSchattenNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianSchattenNorm(uplo,A,p):
  return SymmetricSchattenNorm(uplo,A,p,True)

lib.ElTwoNorm_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoNorm_s.restype = c_uint
lib.ElTwoNorm_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoNorm_d.restype = c_uint
lib.ElTwoNorm_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoNorm_c.restype = c_uint
lib.ElTwoNorm_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoNorm_z.restype = c_uint
lib.ElTwoNormDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoNormDist_s.restype = c_uint
lib.ElTwoNormDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoNormDist_d.restype = c_uint
lib.ElTwoNormDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElTwoNormDist_c.restype = c_uint
lib.ElTwoNormDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElTwoNormDist_z.restype = c_uint
def TwoNorm(A):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoNorm_s(*args)
    elif A.tag == dTag: lib.ElTwoNorm_d(*args)
    elif A.tag == cTag: lib.ElTwoNorm_c(*args)
    elif A.tag == zTag: lib.ElTwoNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoNormDist_s(*args)
    elif A.tag == dTag: lib.ElTwoNormDist_d(*args)
    elif A.tag == cTag: lib.ElTwoNormDist_c(*args)
    elif A.tag == zTag: lib.ElTwoNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricTwoNorm_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricTwoNorm_s.restype = c_uint
lib.ElSymmetricTwoNorm_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricTwoNorm_d.restype = c_uint
lib.ElSymmetricTwoNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricTwoNorm_c.restype = c_uint
lib.ElSymmetricTwoNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricTwoNorm_z.restype = c_uint
lib.ElSymmetricTwoNormDist_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricTwoNormDist_s.restype = c_uint
lib.ElSymmetricTwoNormDist_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricTwoNormDist_d.restype = c_uint
lib.ElSymmetricTwoNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricTwoNormDist_c.restype = c_uint
lib.ElSymmetricTwoNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricTwoNormDist_z.restype = c_uint
lib.ElHermitianTwoNorm_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianTwoNorm_c.restype = c_uint
lib.ElHermitianTwoNorm_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianTwoNorm_z.restype = c_uint
lib.ElHermitianTwoNormDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElHermitianTwoNormDist_c.restype = c_uint
lib.ElHermitianTwoNormDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElHermitianTwoNormDist_z.restype = c_uint
def SymmetricTwoNorm(uplo,A,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricTwoNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricTwoNorm_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianTwoNorm_c(*args)
      else:         lib.ElSymmetricTwoNorm_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianTwoNorm_z(*args)
      else:         lib.ElSymmetricTwoNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricTwoNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricTwoNormDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianTwoNormDist_c(*args)
      else:         lib.ElSymmetricTwoNormDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianTwoNormDist_z(*args)
      else:         lib.ElSymmetricTwoNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianTwoNorm(uplo,A):
  return SymmetricTwoNorm(uplo,A,True)

lib.ElZeroNorm_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElZeroNorm_i.restype = c_uint
lib.ElZeroNorm_s.argtypes = [c_void_p,sType,POINTER(iType)]
lib.ElZeroNorm_s.restype = c_uint
lib.ElZeroNorm_d.argtypes = [c_void_p,dType,POINTER(iType)]
lib.ElZeroNorm_d.restype = c_uint
lib.ElZeroNorm_c.argtypes = [c_void_p,sType,POINTER(iType)]
lib.ElZeroNorm_c.restype = c_uint
lib.ElZeroNorm_z.argtypes = [c_void_p,dType,POINTER(iType)]
lib.ElZeroNorm_z.restype = c_uint
lib.ElZeroNormDist_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElZeroNormDist_i.restype = c_uint
lib.ElZeroNormDist_s.argtypes = [c_void_p,sType,POINTER(iType)]
lib.ElZeroNormDist_s.restype = c_uint
lib.ElZeroNormDist_d.argtypes = [c_void_p,dType,POINTER(iType)]
lib.ElZeroNormDist_d.restype = c_uint
lib.ElZeroNormDist_c.argtypes = [c_void_p,sType,POINTER(iType)]
lib.ElZeroNormDist_c.restype = c_uint
lib.ElZeroNormDist_z.argtypes = [c_void_p,cType,POINTER(iType)]
lib.ElZeroNormDist_z.restype = c_uint
def ZeroNorm(A,tol=0):
  norm = iType()
  args = [A.obj,tol,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElZeroNorm_i(*args) 
    elif A.tag == sTag: lib.ElZeroNorm_s(*args)
    elif A.tag == dTag: lib.ElZeroNorm_d(*args)
    elif A.tag == cTag: lib.ElZeroNorm_c(*args)
    elif A.tag == zTag: lib.ElZeroNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElZeroNormDist_i(*args)
    elif A.tag == sTag: lib.ElZeroNormDist_s(*args)
    elif A.tag == dTag: lib.ElZeroNormDist_d(*args)
    elif A.tag == cTag: lib.ElZeroNormDist_c(*args)
    elif A.tag == zTag: lib.ElZeroNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElTwoNormEstimate_s.argtypes = [c_void_p,sType,iType,POINTER(sType)]
lib.ElTwoNormEstimate_s.restype = c_uint
lib.ElTwoNormEstimate_d.argtypes = [c_void_p,dType,iType,POINTER(dType)]
lib.ElTwoNormEstimate_d.restype = c_uint
lib.ElTwoNormEstimate_c.argtypes = [c_void_p,sType,iType,POINTER(sType)]
lib.ElTwoNormEstimate_c.restype = c_uint
lib.ElTwoNormEstimate_z.argtypes = [c_void_p,dType,iType,POINTER(dType)]
lib.ElTwoNormEstimate_z.restype = c_uint
lib.ElTwoNormEstimateDist_s.argtypes = [c_void_p,sType,iType,POINTER(sType)]
lib.ElTwoNormEstimateDist_s.restype = c_uint
lib.ElTwoNormEstimateDist_d.argtypes = [c_void_p,dType,iType,POINTER(dType)]
lib.ElTwoNormEstimateDist_d.restype = c_uint
lib.ElTwoNormEstimateDist_c.argtypes = [c_void_p,sType,iType,POINTER(sType)]
lib.ElTwoNormEstimateDist_c.restype = c_uint
lib.ElTwoNormEstimateDist_z.argtypes = [c_void_p,dType,iType,POINTER(dType)]
lib.ElTwoNormEstimateDist_z.restype = c_uint
def TwoNormEstimate(A,tol=1e-6,maxIts=100):
  norm = TagToType(Base(A.tag))()
  args = [A.obj,tol,maxIts,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoNormEstimate_s(*args)
    elif A.tag == dTag: lib.ElTwoNormEstimate_d(*args)
    elif A.tag == cTag: lib.ElTwoNormEstimate_c(*args)
    elif A.tag == zTag: lib.ElTwoNormEstimate_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoNormEstimateDist_s(*args)
    elif A.tag == dTag: lib.ElTwoNormEstimateDist_d(*args)
    elif A.tag == cTag: lib.ElTwoNormEstimateDist_c(*args)
    elif A.tag == zTag: lib.ElTwoNormEstimateDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm

lib.ElSymmetricTwoNormEstimate_s.argtypes = \
  [c_uint,c_void_p,sType,iType,POINTER(sType)]
lib.ElSymmetricTwoNormEstimate_s.restype = c_uint
lib.ElSymmetricTwoNormEstimate_d.argtypes = \
  [c_uint,c_void_p,dType,iType,POINTER(dType)]
lib.ElSymmetricTwoNormEstimate_d.restype = c_uint
lib.ElSymmetricTwoNormEstimate_c.argtypes = \
  [c_uint,c_void_p,sType,iType,POINTER(sType)]
lib.ElSymmetricTwoNormEstimate_c.restype = c_uint
lib.ElSymmetricTwoNormEstimate_z.argtypes = \
  [c_uint,c_void_p,dType,iType,POINTER(dType)]
lib.ElSymmetricTwoNormEstimate_z.restype = c_uint
lib.ElSymmetricTwoNormEstimateDist_s.argtypes = \
  [c_uint,c_void_p,sType,iType,POINTER(sType)]
lib.ElSymmetricTwoNormEstimateDist_s.restype = c_uint
lib.ElSymmetricTwoNormEstimateDist_d.argtypes = \
  [c_uint,c_void_p,dType,iType,POINTER(dType)]
lib.ElSymmetricTwoNormEstimateDist_d.restype = c_uint
lib.ElHermitianTwoNormEstimateDist_c.argtypes = \
  [c_uint,c_void_p,sType,iType,POINTER(sType)]
lib.ElHermitianTwoNormEstimateDist_c.restype = c_uint
lib.ElHermitianTwoNormEstimateDist_z.argtypes = \
  [c_uint,c_void_p,dType,iType,POINTER(dType)]
lib.ElHermitianTwoNormEstimateDist_z.restype = c_uint
lib.ElHermitianTwoNormEstimate_c.argtypes = \
  [c_uint,c_void_p,sType,iType,POINTER(sType)]
lib.ElHermitianTwoNormEstimate_c.restype = c_uint
lib.ElHermitianTwoNormEstimate_z.argtypes = \
  [c_uint,c_void_p,dType,iType,POINTER(dType)]
lib.ElHermitianTwoNormEstimate_z.restype = c_uint
lib.ElHermitianTwoNormEstimateDist_c.argtypes = \
  [c_uint,c_void_p,sType,iType,POINTER(sType)]
lib.ElHermitianTwoNormEstimateDist_c.restype = c_uint
lib.ElHermitianTwoNormEstimateDist_z.argtypes = \
  [c_uint,c_void_p,dType,iType,POINTER(dType)]
lib.ElHermitianTwoNormEstimateDist_z.restype = c_uint
def SymmetricTwoNormEstimate(uplo,A,tol=1e-6,maxIts=100,conjugate=False):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,tol,maxIts,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricTwoNormEstimate_s(*args)
    elif A.tag == dTag: lib.ElSymmetricTwoNormEstimate_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianTwoNormEstimate_c(*args)
      else:         lib.ElSymmetricTwoNormEstimate_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianTwoNormEstimate_z(*args)
      else:         lib.ElSymmetricTwoNormEstimate_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricTwoNormEstimateDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricTwoNormEstimateDist_d(*args)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianTwoNormEstimateDist_c(*args)
      else:         lib.ElSymmetricTwoNormEstimateDist_c(*args)
    elif A.tag == zTag: 
      if conjugate: lib.ElHermitianTwoNormEstimateDist_z(*args)
      else:         lib.ElSymmetricTwoNormEstimateDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm
def HermitianTwoNormEstimate(uplo,A,tol=1e-6,maxIts=100):
  return HermitianTwoNormEstimate(uplo,A,tol,maxits,True)

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
lib.ElSpectralPortrait_s.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortrait_s.restype = c_uint
lib.ElSpectralPortrait_d.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortrait_d.restype = c_uint
lib.ElSpectralPortrait_c.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortrait_c.restype = c_uint
lib.ElSpectralPortrait_z.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortrait_z.restype = c_uint
lib.ElSpectralPortraitDist_s.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortraitDist_s.restype = c_uint
lib.ElSpectralPortraitDist_d.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortraitDist_d.restype = c_uint
lib.ElSpectralPortraitDist_c.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortraitDist_c.restype = c_uint
lib.ElSpectralPortraitDist_z.argtypes = [c_void_p,c_void_p,iType,iType]
lib.ElSpectralPortraitDist_z.restype = c_uint
lib.ElSpectralPortraitX_s.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElSpectralPortraitX_s.restype = c_uint
lib.ElSpectralPortraitX_d.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElSpectralPortraitX_d.restype = c_uint
lib.ElSpectralPortraitX_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElSpectralPortraitX_c.restype = c_uint
lib.ElSpectralPortraitX_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElSpectralPortraitX_z.restype = c_uint
lib.ElSpectralPortraitXDist_s.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElSpectralPortraitXDist_s.restype = c_uint
lib.ElSpectralPortraitXDist_d.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElSpectralPortraitXDist_d.restype = c_uint
lib.ElSpectralPortraitXDist_c.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_s]
lib.ElSpectralPortraitXDist_c.restype = c_uint
lib.ElSpectralPortraitXDist_z.argtypes = \
  [c_void_p,c_void_p,iType,iType,PseudospecCtrl_d]
lib.ElSpectralPortraitXDist_z.restype = c_uint
def SpectralPortrait(A,realSize=200,imagSize=200,ctrl=None):
  if type(A) is Matrix:
    invNormMap = Matrix(Base(A.tag))
    args = [A.obj,invNormMap.obj,realSize,imagSize]
    argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSpectralPortrait_s(*args)
      else:            lib.ElSpectralPortraitX_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralPortrait_d(*args)
      else:            lib.ElSpectralPortraitX_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralPortrait_c(*args)
      else:            lib.ElSpectralPortraitX_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralPortrait_z(*args)
      else:            lib.ElSpectralPortraitX_z(*argsCtrl)
    else: DataExcept()
    return invNormMap
  elif type(A) is DistMatrix:
    invNormMap = DistMatrix(Base(A.tag),MC,MR,A.Grid())
    args = [A.obj,invNormMap.obj,realSize,imagSize]
    argsCtrl = [A.obj,invNormMap.obj,realSize,imagSize,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSpectralPortraitDist_s(*args)
      else:            lib.ElSpectralPortraitXDist_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralPortraitDist_d(*args)
      else:            lib.ElSpectralPortraitXDist_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralPortraitDist_c(*args)
      else:            lib.ElSpectralPortraitXDist_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralPortraitDist_z(*args)
      else:            lib.ElSpectralPortraitXDist_z(*argsCtrl)
    else: DataExcept()
    return invNormMap
  else: TypeExcept()

# Investigate a particular window of the complex plane
# ----------------------------------------------------
lib.ElSpectralWindow_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType]
lib.ElSpectralWindow_s.restype = c_uint
lib.ElSpectralWindow_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType]
lib.ElSpectralWindow_d.restype = c_uint
lib.ElSpectralWindow_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType]
lib.ElSpectralWindow_c.restype = c_uint
lib.ElSpectralWindow_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType]
lib.ElSpectralWindow_z.restype = c_uint
lib.ElSpectralWindowDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType]
lib.ElSpectralWindowDist_s.restype = c_uint
lib.ElSpectralWindowDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType]
lib.ElSpectralWindowDist_d.restype = c_uint
lib.ElSpectralWindowDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType]
lib.ElSpectralWindowDist_c.restype = c_uint
lib.ElSpectralWindowDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType]
lib.ElSpectralWindowDist_z.restype = c_uint
lib.ElSpectralWindowX_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElSpectralWindowX_s.restype = c_uint
lib.ElSpectralWindowX_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElSpectralWindowX_d.restype = c_uint
lib.ElSpectralWindowX_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElSpectralWindowX_c.restype = c_uint
lib.ElSpectralWindowX_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElSpectralWindowX_z.restype = c_uint
lib.ElSpectralWindowXDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElSpectralWindowXDist_s.restype = c_uint
lib.ElSpectralWindowXDist_d.argtypes = \
  [c_void_p,c_void_p,dType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElSpectralWindowXDist_d.restype = c_uint
lib.ElSpectralWindowXDist_c.argtypes = \
  [c_void_p,c_void_p,cType,sType,sType,iType,iType,PseudospecCtrl_s]
lib.ElSpectralWindowXDist_c.restype = c_uint
lib.ElSpectralWindowXDist_z.argtypes = \
  [c_void_p,c_void_p,zType,dType,dType,iType,iType,PseudospecCtrl_d]
lib.ElSpectralWindowXDist_z.restype = c_uint
def SpectralWindow \
    (A,centerPre,realWidth,imagWidth,realSize=200,imagSize=200,ctrl=None):
  center = TagToType(A.tag)(centerPre)
  if type(A) is Matrix:
    invNormMap = Matrix(Base(A.tag))
    args = [A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize]
    argsCtrl = [A.obj,invNormMap.obj,center,realWidth,imagWidth,
                realSize,imagSize,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSpectralWindow_s(*args)
      else:            lib.ElSpectralWindowX_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralWindow_d(*args)
      else:            lib.ElSpectralWindowX_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralWindow_c(*args)
      else:            lib.ElSpectralWindowX_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralWindow_z(*args)
      else:            lib.ElSpectralWindowX_z(*argsCtrl)
    else: DataExcept()
    return invNormMap
  elif type(A) is DistMatrix:
    invNormMap = DistMatrix(Base(A.tag),MC,MR,A.Grid())
    args = [A.obj,invNormMap.obj,center,realWidth,imagWidth,realSize,imagSize]
    argsCtrl = [A.obj,invNormMap.obj,center,realWidth,imagWidth,
                realSize,imagSize,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSpectralWindowDist_s(*args)
      else:            lib.ElSpectralWindowXDist_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralWindowDist_d(*args)
      else:            lib.ElSpectralWindowXDist_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralWindowDist_c(*args)
      else:            lib.ElSpectralWindowXDist_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralWindowDist_z(*args)
      else:            lib.ElSpectralWindowXDist_z(*argsCtrl)
    else: DataExcept()
    return invNormMap
  else: TypeExcept()

# Investigate a specific set of points in the complex plane
# ---------------------------------------------------------
lib.ElSpectralCloud_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloud_s.restype = c_uint
lib.ElSpectralCloud_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloud_d.restype = c_uint
lib.ElSpectralCloud_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloud_c.restype = c_uint
lib.ElSpectralCloud_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloud_z.restype = c_uint
lib.ElSpectralCloudDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloudDist_s.restype = c_uint
lib.ElSpectralCloudDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloudDist_d.restype = c_uint
lib.ElSpectralCloudDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloudDist_c.restype = c_uint
lib.ElSpectralCloudDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElSpectralCloudDist_z.restype = c_uint
lib.ElSpectralCloudX_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElSpectralCloudX_s.restype = c_uint
lib.ElSpectralCloudX_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElSpectralCloudX_d.restype = c_uint
lib.ElSpectralCloudX_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElSpectralCloudX_c.restype = c_uint
lib.ElSpectralCloudX_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElSpectralCloudX_z.restype = c_uint
lib.ElSpectralCloudXDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElSpectralCloudXDist_s.restype = c_uint
lib.ElSpectralCloudXDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElSpectralCloudXDist_d.restype = c_uint
lib.ElSpectralCloudXDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_s]
lib.ElSpectralCloudXDist_c.restype = c_uint
lib.ElSpectralCloudXDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,PseudospecCtrl_d]
lib.ElSpectralCloudXDist_z.restype = c_uint
def SpectralCloud(A,shifts,ctrl=None):
  if type(A) is Matrix:
    invNorms = Matrix(Base(A.tag)) 
    args = [A.obj,shifts.obj,invNorms.obj]
    argsCtrl = [A.obj,shifts.obj,invNorms.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSpectralCloud_s(*args)
      else:            lib.ElSpectralCloudX_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralCloud_d(*args)
      else:            lib.ElSpectralCloudX_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralCloud_c(*args)
      else:            lib.ElSpectralCloudX_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralCloud_z(*args)
      else:            lib.ElSpectralCloudX_z(*argsCtrl)
    else: DataExcept()
    return invNorms
  elif type(A) is DistMatrix:
    invNorms = DistMatrix(Base(A.tag),VR,STAR,A.Grid())
    args = [A.obj,shifts.obj,invNorms.obj]
    argsCtrl = [A.obj,shifts.obj,invNorms.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl == None: lib.ElSpectralCloudDist_s(*args)
      else:            lib.ElSpectralCloudXDist_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl == None: lib.ElSpectralCloudDist_d(*args)
      else:            lib.ElSpectralCloudXDist_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl == None: lib.ElSpectralCloudDist_c(*args)
      else:            lib.ElSpectralCloudXDist_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl == None: lib.ElSpectralCloudDist_z(*args)
      else:            lib.ElSpectralCloudXDist_z(*argsCtrl)
    else: DataExcept()
    return invNorms
  else: TypeExcept()

# Trace
# =====
lib.ElTrace_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElTrace_s.restype = c_uint
lib.ElTrace_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElTrace_d.restype = c_uint
lib.ElTrace_c.argtypes = [c_void_p,POINTER(cType)]
lib.ElTrace_c.restype = c_uint
lib.ElTrace_z.argtypes = [c_void_p,POINTER(zType)]
lib.ElTrace_z.restype = c_uint
lib.ElTraceDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElTraceDist_s.restype = c_uint
lib.ElTraceDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElTraceDist_d.restype = c_uint
lib.ElTraceDist_c.argtypes = [c_void_p,POINTER(cType)]
lib.ElTraceDist_c.restype = c_uint
lib.ElTraceDist_z.argtypes = [c_void_p,POINTER(zType)]
lib.ElTraceDist_z.restype = c_uint
def Trace(A):
  trace = TagToType(A.tag)()
  args = [A.obj,pointer(trace)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrace_s(*args)
    elif A.tag == dTag: lib.ElTrace_d(*args)
    elif A.tag == cTag: lib.ElTrace_c(*args)
    elif A.tag == zTag: lib.ElTrace_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTraceDist_s(*args)
    elif A.tag == dTag: lib.ElTraceDist_d(*args)
    elif A.tag == cTag: lib.ElTraceDist_c(*args)
    elif A.tag == zTag: lib.ElTraceDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return trace
