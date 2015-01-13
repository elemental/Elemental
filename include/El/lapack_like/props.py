#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
from factor import *
import ctypes

# Condition number
# ================
lib.ElCondition_s.argtypes = \
lib.ElCondition_c.argtypes = \
lib.ElConditionDist_s.argtypes = \
lib.ElConditionDist_c.argtypes = \
  [c_void_p,c_uint,POINTER(sType)]

lib.ElCondition_d.argtypes = \
lib.ElCondition_z.argtypes = \
lib.ElConditionDist_d.argtypes = \
lib.ElConditionDist_z.argtypes = \
  [c_void_p,c_uint,POINTER(dType)]

lib.ElCondition_s.restype = \
lib.ElCondition_d.restype = \
lib.ElCondition_c.restype = \
lib.ElCondition_z.restype = \
lib.ElConditionDist_s.restype = \
lib.ElConditionDist_d.restype = \
lib.ElConditionDist_c.restype = \
lib.ElConditionDist_z.restype = \
  c_uint

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
  return cond.value

lib.ElFrobeniusCondition_s.argtypes = \
lib.ElFrobeniusCondition_c.argtypes = \
lib.ElFrobeniusConditionDist_s.argtypes = \
lib.ElFrobeniusConditionDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElFrobeniusCondition_d.argtypes = \
lib.ElFrobeniusCondition_z.argtypes = \
lib.ElFrobeniusConditionDist_d.argtypes = \
lib.ElFrobeniusConditionDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElFrobeniusCondition_s.restype = \
lib.ElFrobeniusCondition_d.restype = \
lib.ElFrobeniusCondition_c.restype = \
lib.ElFrobeniusCondition_z.restype = \
lib.ElFrobeniusConditionDist_s.restype = \
lib.ElFrobeniusConditionDist_d.restype = \
lib.ElFrobeniusConditionDist_c.restype = \
lib.ElFrobeniusConditionDist_z.restype = \
  c_uint

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

lib.ElInfinityCondition_s.argtypes = \
lib.ElInfinityCondition_c.argtypes = \
lib.ElInfinityConditionDist_s.argtypes = \
lib.ElInfinityConditionDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElInfinityCondition_d.argtypes = \
lib.ElInfinityCondition_z.argtypes = \
lib.ElInfinityConditionDist_d.argtypes = \
lib.ElInfinityConditionDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElInfinityCondition_s.restype = \
lib.ElInfinityCondition_d.restype = \
lib.ElInfinityCondition_c.restype = \
lib.ElInfinityCondition_z.restype = \
lib.ElInfinityConditionDist_s.restype = \
lib.ElInfinityConditionDist_d.restype = \
lib.ElInfinityConditionDist_c.restype = \
lib.ElInfinityConditionDist_z.restype = \
  c_uint

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

lib.ElMaxCondition_s.argtypes = \
lib.ElMaxCondition_c.argtypes = \
lib.ElMaxConditionDist_s.argtypes = \
lib.ElMaxConditionDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElMaxCondition_d.argtypes = \
lib.ElMaxCondition_z.argtypes = \
lib.ElMaxConditionDist_d.argtypes = \
lib.ElMaxConditionDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElMaxCondition_s.restype = \
lib.ElMaxCondition_d.restype = \
lib.ElMaxCondition_c.restype = \
lib.ElMaxCondition_z.restype = \
lib.ElMaxConditionDist_s.restype = \
lib.ElMaxConditionDist_d.restype = \
lib.ElMaxConditionDist_c.restype = \
lib.ElMaxConditionDist_z.restype = \
  c_uint

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

lib.ElOneCondition_s.argtypes = \
lib.ElOneCondition_c.argtypes = \
lib.ElOneConditionDist_s.argtypes = \
lib.ElOneConditionDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElOneCondition_d.argtypes = \
lib.ElOneCondition_z.argtypes = \
lib.ElOneConditionDist_d.argtypes = \
lib.ElOneConditionDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElOneCondition_s.restype = \
lib.ElOneCondition_d.restype = \
lib.ElOneCondition_c.restype = \
lib.ElOneCondition_z.restype = \
lib.ElOneConditionDist_s.restype = \
lib.ElOneConditionDist_d.restype = \
lib.ElOneConditionDist_c.restype = \
lib.ElOneConditionDist_z.restype = \
  c_uint

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

lib.ElTwoCondition_s.argtypes = \
lib.ElTwoCondition_c.argtypes = \
lib.ElTwoConditionDist_s.argtypes = \
lib.ElTwoConditionDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElTwoCondition_d.argtypes = \
lib.ElTwoCondition_z.argtypes = \
lib.ElTwoConditionDist_d.argtypes = \
lib.ElTwoConditionDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElTwoCondition_s.restype = \
lib.ElTwoCondition_d.restype = \
lib.ElTwoCondition_c.restype = \
lib.ElTwoCondition_z.restype = \
lib.ElTwoConditionDist_s.restype = \
lib.ElTwoConditionDist_d.restype = \
lib.ElTwoConditionDist_c.restype = \
lib.ElTwoConditionDist_z.restype = \
  c_uint

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
lib.ElSafeDeterminant_s.argtypes = \
lib.ElSafeDeterminantDist_s.argtypes = \
  [c_void_p,POINTER(SafeProduct_s)]

lib.ElSafeDeterminant_c.argtypes = \
lib.ElSafeDeterminantDist_c.argtypes = \
  [c_void_p,POINTER(SafeProduct_c)]

lib.ElSafeDeterminant_d.argtypes = \
lib.ElSafeDeterminantDist_d.argtypes = \
  [c_void_p,POINTER(SafeProduct_d)]

lib.ElSafeDeterminant_z.argtypes = \
lib.ElSafeDeterminantDist_z.argtypes = \
  [c_void_p,POINTER(SafeProduct_z)]

lib.ElSafeDeterminant_s.restype = \
lib.ElSafeDeterminant_d.restype = \
lib.ElSafeDeterminant_c.restype = \
lib.ElSafeDeterminant_z.restype = \
lib.ElSafeDeterminantDist_s.restype = \
lib.ElSafeDeterminantDist_d.restype = \
lib.ElSafeDeterminantDist_c.restype = \
lib.ElSafeDeterminantDist_z.restype = \
  c_uint

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

lib.ElSafeHPDDeterminant_s.argtypes = \
lib.ElSafeHPDDeterminant_c.argtypes = \
lib.ElSafeHPDDeterminantDist_s.argtypes = \
lib.ElSafeHPDDeterminantDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(SafeProduct_s)]

lib.ElSafeHPDDeterminant_d.argtypes = \
lib.ElSafeHPDDeterminant_z.argtypes = \
lib.ElSafeHPDDeterminantDist_d.argtypes = \
lib.ElSafeHPDDeterminantDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(SafeProduct_d)]

lib.ElSafeHPDDeterminant_s.restype = \
lib.ElSafeHPDDeterminant_d.restype = \
lib.ElSafeHPDDeterminant_c.restype = \
lib.ElSafeHPDDeterminant_z.restype = \
lib.ElSafeHPDDeterminantDist_s.restype = \
lib.ElSafeHPDDeterminantDist_d.restype = \
lib.ElSafeHPDDeterminantDist_c.restype = \
lib.ElSafeHPDDeterminantDist_z.restype = \
  c_uint

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
lib.ElDeterminant_s.argtypes = \
lib.ElDeterminantDist_s.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElDeterminant_d.argtypes = \
lib.ElDeterminantDist_d.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElDeterminant_c.argtypes = \
lib.ElDeterminantDist_c.argtypes = \
  [c_void_p,POINTER(cType)]

lib.ElDeterminant_z.argtypes = \
lib.ElDeterminantDist_z.argtypes = \
  [c_void_p,POINTER(zType)]

lib.ElDeterminant_s.restype = \
lib.ElDeterminant_d.restype = \
lib.ElDeterminant_c.restype = \
lib.ElDeterminant_z.restype = \
lib.ElDeterminantDist_s.restype = \
lib.ElDeterminantDist_d.restype = \
lib.ElDeterminantDist_c.restype = \
lib.ElDeterminantDist_z.restype = \
  c_uint

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
  return prod.value

lib.ElHPDDeterminant_s.argtypes = \
lib.ElHPDDeterminantDist_s.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]

lib.ElHPDDeterminant_d.argtypes = \
lib.ElHPDDeterminantDist_d.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]

lib.ElHPDDeterminant_c.argtypes = \
lib.ElHPDDeterminantDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(cType)]

lib.ElHPDDeterminant_z.argtypes = \
lib.ElHPDDeterminantDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(zType)]

lib.ElHPDDeterminant_s.restype = \
lib.ElHPDDeterminant_d.restype = \
lib.ElHPDDeterminant_c.restype = \
lib.ElHPDDeterminant_z.restype = \
lib.ElHPDDeterminantDist_s.restype = \
lib.ElHPDDeterminantDist_d.restype = \
lib.ElHPDDeterminantDist_c.restype = \
lib.ElHPDDeterminantDist_z.restype = \
  c_uint

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
  return prod.value

# Inertia
# =======
lib.ElInertia_s.argtypes = \
lib.ElInertia_d.argtypes = \
lib.ElInertia_c.argtypes = \
lib.ElInertia_z.argtypes = \
lib.ElInertiaDist_s.argtypes = \
lib.ElInertiaDist_d.argtypes = \
lib.ElInertiaDist_c.argtypes = \
lib.ElInertiaDist_z.argtypes = \
  [c_uint,c_void_p,c_uint,POINTER(InertiaType)]

lib.ElInertia_s.restype = \
lib.ElInertia_d.restype = \
lib.ElInertia_c.restype = \
lib.ElInertia_z.restype = \
lib.ElInertiaDist_s.restype = \
lib.ElInertiaDist_d.restype = \
lib.ElInertiaDist_c.restype = \
lib.ElInertiaDist_z.restype = \
  c_uint

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
lib.ElNorm_s.argtypes = \
lib.ElNorm_c.argtypes = \
lib.ElNormDist_s.argtypes = \
lib.ElNormDist_c.argtypes = \
  [c_void_p,c_uint,POINTER(sType)]

lib.ElNorm_d.argtypes = \
lib.ElNorm_z.argtypes = \
lib.ElNormDist_d.argtypes = \
lib.ElNormDist_z.argtypes = \
  [c_void_p,c_uint,POINTER(dType)]

lib.ElNorm_s.restype = \
lib.ElNorm_d.restype = \
lib.ElNorm_c.restype = \
lib.ElNorm_z.restype = \
lib.ElNormDist_s.restype = \
lib.ElNormDist_d.restype = \
lib.ElNormDist_c.restype = \
lib.ElNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricNorm_s.argtypes = \
lib.ElSymmetricNorm_c.argtypes = \
lib.ElSymmetricNormDist_s.argtypes = \
lib.ElSymmetricNormDist_c.argtypes = \
lib.ElHermitianNorm_c.argtypes = \
lib.ElHermitianNormDist_c.argtypes = \
  [c_uint,c_void_p,c_uint,POINTER(sType)]

lib.ElSymmetricNorm_d.argtypes = \
lib.ElSymmetricNorm_z.argtypes = \
lib.ElSymmetricNormDist_d.argtypes = \
lib.ElSymmetricNormDist_z.argtypes = \
lib.ElHermitianNorm_z.argtypes = \
lib.ElHermitianNormDist_z.argtypes = \
  [c_uint,c_void_p,c_uint,POINTER(dType)]

lib.ElSymmetricNorm_s.restype = \
lib.ElSymmetricNorm_d.restype = \
lib.ElSymmetricNorm_c.restype = \
lib.ElSymmetricNorm_z.restype = \
lib.ElSymmetricNormDist_s.restype = \
lib.ElSymmetricNormDist_d.restype = \
lib.ElSymmetricNormDist_c.restype = \
lib.ElSymmetricNormDist_z.restype = \
lib.ElHermitianNorm_c.restype = \
lib.ElHermitianNorm_z.restype = \
lib.ElHermitianNormDist_c.restype = \
lib.ElHermitianNormDist_z.restype = \
  c_uint

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
  return norm.value
def HermitianNorm(uplo,A,normType=FROBENIUS_NORM):
  return SymmetricNorm(uplo,A,normType,True)

lib.ElEntrywiseNorm_s.argtypes = \
lib.ElEntrywiseNorm_c.argtypes = \
lib.ElEntrywiseNormDist_s.argtypes = \
lib.ElEntrywiseNormDist_c.argtypes = \
lib.ElEntrywiseNormSparse_s.argtypes = \
lib.ElEntrywiseNormSparse_c.argtypes = \
lib.ElEntrywiseNormDistSparse_s.argtypes = \
lib.ElEntrywiseNormDistSparse_c.argtypes = \
lib.ElEntrywiseNormDistMultiVec_s.argtypes = \
lib.ElEntrywiseNormDistMultiVec_c.argtypes = \
  [c_void_p,sType,POINTER(sType)]

lib.ElEntrywiseNorm_d.argtypes = \
lib.ElEntrywiseNorm_z.argtypes = \
lib.ElEntrywiseNormDist_d.argtypes = \
lib.ElEntrywiseNormDist_z.argtypes = \
lib.ElEntrywiseNormSparse_d.argtypes = \
lib.ElEntrywiseNormSparse_z.argtypes = \
lib.ElEntrywiseNormDistSparse_d.argtypes = \
lib.ElEntrywiseNormDistSparse_z.argtypes = \
lib.ElEntrywiseNormDistMultiVec_d.argtypes = \
lib.ElEntrywiseNormDistMultiVec_z.argtypes = \
  [c_void_p,dType,POINTER(dType)]

lib.ElEntrywiseNorm_s.restype = \
lib.ElEntrywiseNorm_d.restype = \
lib.ElEntrywiseNorm_c.restype = \
lib.ElEntrywiseNorm_z.restype = \
lib.ElEntrywiseNormDist_s.restype = \
lib.ElEntrywiseNormDist_d.restype = \
lib.ElEntrywiseNormDist_c.restype = \
lib.ElEntrywiseNormDist_z.restype = \
lib.ElEntrywiseNormSparse_s.restype = \
lib.ElEntrywiseNormSparse_d.restype = \
lib.ElEntrywiseNormSparse_c.restype = \
lib.ElEntrywiseNormSparse_z.restype = \
lib.ElEntrywiseNormDistSparse_s.restype = \
lib.ElEntrywiseNormDistSparse_d.restype = \
lib.ElEntrywiseNormDistSparse_c.restype = \
lib.ElEntrywiseNormDistSparse_z.restype = \
lib.ElEntrywiseNormDistMultiVec_s.restype = \
lib.ElEntrywiseNormDistMultiVec_d.restype = \
lib.ElEntrywiseNormDistMultiVec_c.restype = \
lib.ElEntrywiseNormDistMultiVec_z.restype = \
  c_uint

def EntrywiseNorm(A,p=1):
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
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElEntrywiseNormSparse_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseNormSparse_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseNormSparse_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseNormSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElEntrywiseNormDistSparse_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseNormDistSparse_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseNormDistSparse_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseNormDistSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == sTag: lib.ElEntrywiseNormDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseNormDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseNormDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseNormDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value

lib.ElSymmetricEntrywiseNorm_s.argtypes = \
lib.ElSymmetricEntrywiseNorm_c.argtypes = \
lib.ElSymmetricEntrywiseNormDist_s.argtypes = \
lib.ElSymmetricEntrywiseNormDist_c.argtypes = \
lib.ElSymmetricEntrywiseNormSparse_s.argtypes = \
lib.ElSymmetricEntrywiseNormSparse_c.argtypes = \
lib.ElSymmetricEntrywiseNormDistSparse_s.argtypes = \
lib.ElSymmetricEntrywiseNormDistSparse_c.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]

lib.ElSymmetricEntrywiseNorm_d.argtypes = \
lib.ElSymmetricEntrywiseNorm_z.argtypes = \
lib.ElSymmetricEntrywiseNormDist_d.argtypes = \
lib.ElSymmetricEntrywiseNormDist_z.argtypes = \
lib.ElSymmetricEntrywiseNormSparse_d.argtypes = \
lib.ElSymmetricEntrywiseNormSparse_z.argtypes = \
lib.ElSymmetricEntrywiseNormDistSparse_d.argtypes = \
lib.ElSymmetricEntrywiseNormDistSparse_z.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]

lib.ElSymmetricEntrywiseNorm_s.restype = \
lib.ElSymmetricEntrywiseNorm_d.restype = \
lib.ElSymmetricEntrywiseNorm_c.restype = \
lib.ElSymmetricEntrywiseNorm_z.restype = \
lib.ElSymmetricEntrywiseNormDist_s.restype = \
lib.ElSymmetricEntrywiseNormDist_d.restype = \
lib.ElSymmetricEntrywiseNormDist_c.restype = \
lib.ElSymmetricEntrywiseNormDist_z.restype = \
lib.ElSymmetricEntrywiseNormSparse_s.restype = \
lib.ElSymmetricEntrywiseNormSparse_d.restype = \
lib.ElSymmetricEntrywiseNormSparse_c.restype = \
lib.ElSymmetricEntrywiseNormSparse_z.restype = \
lib.ElSymmetricEntrywiseNormDistSparse_s.restype = \
lib.ElSymmetricEntrywiseNormDistSparse_d.restype = \
lib.ElSymmetricEntrywiseNormDistSparse_c.restype = \
lib.ElSymmetricEntrywiseNormDistSparse_z.restype = \
  c_uint

def SymmetricEntrywiseNorm(uplo,A,p=1):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,p,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseNorm_d(*args)
    elif A.tag == cTag: lib.ElSymmetricEntrywiseNorm_c(*args)
    elif A.tag == zTag: lib.ElSymmetricEntrywiseNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseNormDist_d(*args)
    elif A.tag == cTag: lib.ElSymmetricEntrywiseNormDist_c(*args)
    elif A.tag == zTag: lib.ElSymmetricEntrywiseNormDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseNormSparse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseNormSparse_d(*args)
    elif A.tag == cTag: lib.ElSymmetricEntrywiseNormSparse_c(*args)
    elif A.tag == zTag: lib.ElSymmetricEntrywiseNormSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElSymmetricEntrywiseNormDistSparse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricEntrywiseNormDistSparse_d(*args)
    elif A.tag == cTag: lib.ElSymmetricEntrywiseNormDistSparse_c(*args)
    elif A.tag == zTag: lib.ElSymmetricEntrywiseNormDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value
def HermitianEntrywiseNorm(uplo,A,p=1):
  return SymmetricEntrywiseNorm(uplo,A,p)

lib.ElFrobeniusNorm_s.argtypes = \
lib.ElFrobeniusNorm_c.argtypes = \
lib.ElFrobeniusNormDist_s.argtypes = \
lib.ElFrobeniusNormDist_c.argtypes = \
lib.ElFrobeniusNormSparse_s.argtypes = \
lib.ElFrobeniusNormSparse_c.argtypes = \
lib.ElFrobeniusNormDistSparse_s.argtypes = \
lib.ElFrobeniusNormDistSparse_c.argtypes = \
lib.ElFrobeniusNormDistMultiVec_s.argtypes = \
lib.ElFrobeniusNormDistMultiVec_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElFrobeniusNorm_d.argtypes = \
lib.ElFrobeniusNorm_z.argtypes = \
lib.ElFrobeniusNormDist_d.argtypes = \
lib.ElFrobeniusNormDist_z.argtypes = \
lib.ElFrobeniusNormSparse_d.argtypes = \
lib.ElFrobeniusNormSparse_z.argtypes = \
lib.ElFrobeniusNormDistSparse_d.argtypes = \
lib.ElFrobeniusNormDistSparse_z.argtypes = \
lib.ElFrobeniusNormDistMultiVec_d.argtypes = \
lib.ElFrobeniusNormDistMultiVec_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElFrobeniusNorm_s.restype = \
lib.ElFrobeniusNorm_d.restype = \
lib.ElFrobeniusNorm_c.restype = \
lib.ElFrobeniusNorm_z.restype = \
lib.ElFrobeniusNormDist_s.restype = \
lib.ElFrobeniusNormDist_d.restype = \
lib.ElFrobeniusNormDist_c.restype = \
lib.ElFrobeniusNormDist_z.restype = \
lib.ElFrobeniusNormSparse_s.restype = \
lib.ElFrobeniusNormSparse_d.restype = \
lib.ElFrobeniusNormSparse_c.restype = \
lib.ElFrobeniusNormSparse_z.restype = \
lib.ElFrobeniusNormDistSparse_s.restype = \
lib.ElFrobeniusNormDistSparse_d.restype = \
lib.ElFrobeniusNormDistSparse_c.restype = \
lib.ElFrobeniusNormDistSparse_z.restype = \
lib.ElFrobeniusNormDistMultiVec_s.restype = \
lib.ElFrobeniusNormDistMultiVec_d.restype = \
lib.ElFrobeniusNormDistMultiVec_c.restype = \
lib.ElFrobeniusNormDistMultiVec_z.restype = \
  c_uint

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
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElFrobeniusNormSparse_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusNormSparse_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusNormSparse_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusNormSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElFrobeniusNormDistSparse_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusNormDistSparse_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusNormDistSparse_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusNormDistSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == sTag: lib.ElFrobeniusNormDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusNormDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusNormDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusNormDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value

lib.ElSymmetricFrobeniusNorm_s.argtypes = \
lib.ElSymmetricFrobeniusNorm_c.argtypes = \
lib.ElSymmetricFrobeniusNormDist_s.argtypes = \
lib.ElSymmetricFrobeniusNormDist_c.argtypes = \
lib.ElSymmetricFrobeniusNormSparse_s.argtypes = \
lib.ElSymmetricFrobeniusNormSparse_c.argtypes = \
lib.ElSymmetricFrobeniusNormDistSparse_s.argtypes = \
lib.ElSymmetricFrobeniusNormDistSparse_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]

lib.ElSymmetricFrobeniusNorm_d.argtypes = \
lib.ElSymmetricFrobeniusNorm_z.argtypes = \
lib.ElSymmetricFrobeniusNormDist_d.argtypes = \
lib.ElSymmetricFrobeniusNormDist_z.argtypes = \
lib.ElSymmetricFrobeniusNormSparse_d.argtypes = \
lib.ElSymmetricFrobeniusNormSparse_z.argtypes = \
lib.ElSymmetricFrobeniusNormDistSparse_d.argtypes = \
lib.ElSymmetricFrobeniusNormDistSparse_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]

lib.ElSymmetricFrobeniusNorm_s.restype = \
lib.ElSymmetricFrobeniusNorm_d.restype = \
lib.ElSymmetricFrobeniusNorm_c.restype = \
lib.ElSymmetricFrobeniusNorm_z.restype = \
lib.ElSymmetricFrobeniusNormDist_s.restype = \
lib.ElSymmetricFrobeniusNormDist_d.restype = \
lib.ElSymmetricFrobeniusNormDist_c.restype = \
lib.ElSymmetricFrobeniusNormDist_z.restype = \
lib.ElSymmetricFrobeniusNormSparse_s.restype = \
lib.ElSymmetricFrobeniusNormSparse_d.restype = \
lib.ElSymmetricFrobeniusNormSparse_c.restype = \
lib.ElSymmetricFrobeniusNormSparse_z.restype = \
lib.ElSymmetricFrobeniusNormDistSparse_s.restype = \
lib.ElSymmetricFrobeniusNormDistSparse_d.restype = \
lib.ElSymmetricFrobeniusNormDistSparse_c.restype = \
lib.ElSymmetricFrobeniusNormDistSparse_z.restype = \
  c_uint

def SymmetricFrobeniusNorm(uplo,A):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricFrobeniusNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricFrobeniusNorm_d(*args)
    elif A.tag == cTag: lib.ElSymmetricFrobeniusNorm_c(*args)
    elif A.tag == zTag: lib.ElSymmetricFrobeniusNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricFrobeniusNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricFrobeniusNormDist_d(*args)
    elif A.tag == cTag: lib.ElSymmetricFrobeniusNormDist_c(*args)
    elif A.tag == zTag: lib.ElSymmetricFrobeniusNormDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElSymmetricFrobeniusNormSparse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricFrobeniusNormSparse_d(*args)
    elif A.tag == cTag: lib.ElSymmetricFrobeniusNormSparse_c(*args)
    elif A.tag == zTag: lib.ElSymmetricFrobeniusNormSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElSymmetricFrobeniusNormDistSparse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricFrobeniusNormDistSparse_d(*args)
    elif A.tag == cTag: lib.ElSymmetricFrobeniusNormDistSparse_c(*args)
    elif A.tag == zTag: lib.ElSymmetricFrobeniusNormDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value
def HermitianFrobeniusNorm(uplo,A):
  return SymmetricFrobeniusNorm(uplo,A)

lib.ElInfinityNorm_s.argtypes = \
lib.ElInfinityNorm_c.argtypes = \
lib.ElInfinityNormDist_s.argtypes = \
lib.ElInfinityNormDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElInfinityNorm_d.argtypes = \
lib.ElInfinityNorm_z.argtypes = \
lib.ElInfinityNormDist_d.argtypes = \
lib.ElInfinityNormDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElInfinityNorm_s.restype = \
lib.ElInfinityNorm_d.restype = \
lib.ElInfinityNorm_c.restype = \
lib.ElInfinityNorm_z.restype = \
lib.ElInfinityNormDist_s.restype = \
lib.ElInfinityNormDist_d.restype = \
lib.ElInfinityNormDist_c.restype = \
lib.ElInfinityNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricInfinityNorm_s.argtypes = \
lib.ElSymmetricInfinityNorm_c.argtypes = \
lib.ElSymmetricInfinityNormDist_s.argtypes = \
lib.ElSymmetricInfinityNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]

lib.ElSymmetricInfinityNorm_d.argtypes = \
lib.ElSymmetricInfinityNorm_z.argtypes = \
lib.ElSymmetricInfinityNormDist_d.argtypes = \
lib.ElSymmetricInfinityNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]

lib.ElSymmetricInfinityNorm_s.restype = \
lib.ElSymmetricInfinityNorm_d.restype = \
lib.ElSymmetricInfinityNorm_c.restype = \
lib.ElSymmetricInfinityNorm_z.restype = \
lib.ElSymmetricInfinityNormDist_s.restype = \
lib.ElSymmetricInfinityNormDist_d.restype = \
lib.ElSymmetricInfinityNormDist_c.restype = \
lib.ElSymmetricInfinityNormDist_z.restype = \
  c_uint

def SymmetricInfinityNorm(uplo,A):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricInfinityNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricInfinityNorm_d(*args)
    elif A.tag == cTag: lib.ElSymmetricInfinityNorm_c(*args)
    elif A.tag == zTag: lib.ElSymmetricInfinityNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricInfinityNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricInfinityNormDist_d(*args)
    elif A.tag == cTag: lib.ElSymmetricInfinityNormDist_c(*args)
    elif A.tag == zTag: lib.ElSymmetricInfinityNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value
def HermitianInfinityNorm(uplo,A):
  return SymmetricInfinityNorm(uplo,A)

lib.ElKyFanNorm_s.argtypes = \
lib.ElKyFanNorm_c.argtypes = \
lib.ElKyFanNormDist_s.argtypes = \
lib.ElKyFanNormDist_c.argtypes = \
  [c_void_p,iType,POINTER(sType)]

lib.ElKyFanNorm_d.argtypes = \
lib.ElKyFanNorm_z.argtypes = \
lib.ElKyFanNormDist_d.argtypes = \
lib.ElKyFanNormDist_z.argtypes = \
  [c_void_p,iType,POINTER(dType)]

lib.ElKyFanNorm_s.restype = \
lib.ElKyFanNorm_d.restype = \
lib.ElKyFanNorm_c.restype = \
lib.ElKyFanNorm_z.restype = \
lib.ElKyFanNormDist_s.restype = \
lib.ElKyFanNormDist_d.restype = \
lib.ElKyFanNormDist_c.restype = \
lib.ElKyFanNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricKyFanNorm_s.argtypes = \
lib.ElSymmetricKyFanNorm_c.argtypes = \
lib.ElSymmetricKyFanNormDist_s.argtypes = \
lib.ElSymmetricKyFanNormDist_c.argtypes = \
lib.ElHermitianKyFanNorm_c.argtypes = \
lib.ElHermitianKyFanNormDist_c.argtypes = \
  [c_uint,c_void_p,iType,POINTER(sType)]

lib.ElSymmetricKyFanNorm_d.argtypes = \
lib.ElSymmetricKyFanNorm_z.argtypes = \
lib.ElSymmetricKyFanNormDist_d.argtypes = \
lib.ElSymmetricKyFanNormDist_z.argtypes = \
lib.ElHermitianKyFanNorm_z.argtypes = \
lib.ElHermitianKyFanNormDist_z.argtypes = \
  [c_uint,c_void_p,iType,POINTER(dType)]

lib.ElSymmetricKyFanNorm_s.restype = \
lib.ElSymmetricKyFanNorm_d.restype = \
lib.ElSymmetricKyFanNorm_c.restype = \
lib.ElSymmetricKyFanNorm_z.restype = \
lib.ElSymmetricKyFanNormDist_s.restype = \
lib.ElSymmetricKyFanNormDist_d.restype = \
lib.ElSymmetricKyFanNormDist_c.restype = \
lib.ElSymmetricKyFanNormDist_z.restype = \
lib.ElHermitianKyFanNorm_c.restype = \
lib.ElHermitianKyFanNorm_z.restype = \
lib.ElHermitianKyFanNormDist_c.restype = \
lib.ElHermitianKyFanNormDist_z.restype = \
  c_uint

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
  return norm.value
def HermitianKyFanNorm(uplo,A,k):
  return SymmetricKyFanNorm(uplo,A,k,True)

lib.ElKyFanSchattenNorm_s.argtypes = \
lib.ElKyFanSchattenNorm_c.argtypes = \
lib.ElKyFanSchattenNormDist_s.argtypes = \
lib.ElKyFanSchattenNormDist_c.argtypes = \
  [c_void_p,iType,sType,POINTER(sType)]

lib.ElKyFanSchattenNorm_d.argtypes = \
lib.ElKyFanSchattenNorm_z.argtypes = \
lib.ElKyFanSchattenNormDist_d.argtypes = \
lib.ElKyFanSchattenNormDist_z.argtypes = \
  [c_void_p,iType,dType,POINTER(dType)]

lib.ElKyFanSchattenNorm_s.restype = \
lib.ElKyFanSchattenNorm_d.restype = \
lib.ElKyFanSchattenNorm_c.restype = \
lib.ElKyFanSchattenNorm_z.restype = \
lib.ElKyFanSchattenNormDist_s.restype = \
lib.ElKyFanSchattenNormDist_d.restype = \
lib.ElKyFanSchattenNormDist_c.restype = \
lib.ElKyFanSchattenNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricKyFanSchattenNorm_s.argtypes = \
lib.ElSymmetricKyFanSchattenNorm_c.argtypes = \
lib.ElSymmetricKyFanSchattenNormDist_s.argtypes = \
lib.ElSymmetricKyFanSchattenNormDist_c.argtypes = \
lib.ElHermitianKyFanSchattenNorm_c.argtypes = \
lib.ElHermitianKyFanSchattenNormDist_c.argtypes = \
  [c_uint,c_void_p,iType,sType,POINTER(sType)]

lib.ElSymmetricKyFanSchattenNorm_d.argtypes = \
lib.ElSymmetricKyFanSchattenNorm_z.argtypes = \
lib.ElSymmetricKyFanSchattenNormDist_d.argtypes = \
lib.ElSymmetricKyFanSchattenNormDist_z.argtypes = \
lib.ElHermitianKyFanSchattenNorm_z.argtypes = \
lib.ElHermitianKyFanSchattenNormDist_z.argtypes = \
  [c_uint,c_void_p,iType,dType,POINTER(dType)]

lib.ElSymmetricKyFanSchattenNorm_s.restype = \
lib.ElSymmetricKyFanSchattenNorm_d.restype = \
lib.ElSymmetricKyFanSchattenNorm_c.restype = \
lib.ElSymmetricKyFanSchattenNorm_z.restype = \
lib.ElSymmetricKyFanSchattenNormDist_s.restype = \
lib.ElSymmetricKyFanSchattenNormDist_d.restype = \
lib.ElSymmetricKyFanSchattenNormDist_c.restype = \
lib.ElSymmetricKyFanSchattenNormDist_z.restype = \
lib.ElHermitianKyFanSchattenNorm_c.restype = \
lib.ElHermitianKyFanSchattenNorm_z.restype = \
lib.ElHermitianKyFanSchattenNormDist_c.restype = \
lib.ElHermitianKyFanSchattenNormDist_z.restype = \
  c_uint

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
  return norm.value
def HermitianKyFanSchattenNorm(uplo,A,k,p):
  return SymmetricKyFanSchattenNorm(uplo,A,k,p,True)

lib.ElMaxNorm_i.argtypes = \
lib.ElMaxNormDist_i.argtypes = \
lib.ElMaxNormSparse_i.argtypes = \
lib.ElMaxNormDistSparse_i.argtypes = \
lib.ElMaxNormDistMultiVec_i.argtypes = \
  [c_void_p,POINTER(iType)]

lib.ElMaxNorm_s.argtypes = \
lib.ElMaxNorm_c.argtypes = \
lib.ElMaxNormDist_s.argtypes = \
lib.ElMaxNormDist_c.argtypes = \
lib.ElMaxNormSparse_s.argtypes = \
lib.ElMaxNormSparse_c.argtypes = \
lib.ElMaxNormDistSparse_s.argtypes = \
lib.ElMaxNormDistSparse_c.argtypes = \
lib.ElMaxNormDistMultiVec_s.argtypes = \
lib.ElMaxNormDistMultiVec_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElMaxNorm_d.argtypes = \
lib.ElMaxNorm_z.argtypes = \
lib.ElMaxNormDist_d.argtypes = \
lib.ElMaxNormDist_z.argtypes = \
lib.ElMaxNormSparse_d.argtypes = \
lib.ElMaxNormSparse_z.argtypes = \
lib.ElMaxNormDistSparse_d.argtypes = \
lib.ElMaxNormDistSparse_z.argtypes = \
lib.ElMaxNormDistMultiVec_d.argtypes = \
lib.ElMaxNormDistMultiVec_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElMaxNorm_i.restype = \
lib.ElMaxNorm_s.restype = \
lib.ElMaxNorm_d.restype = \
lib.ElMaxNorm_c.restype = \
lib.ElMaxNorm_z.restype = \
lib.ElMaxNormDist_i.restype = \
lib.ElMaxNormDist_s.restype = \
lib.ElMaxNormDist_d.restype = \
lib.ElMaxNormDist_c.restype = \
lib.ElMaxNormDist_z.restype = \
lib.ElMaxNormSparse_i.restype = \
lib.ElMaxNormSparse_s.restype = \
lib.ElMaxNormSparse_d.restype = \
lib.ElMaxNormSparse_c.restype = \
lib.ElMaxNormSparse_z.restype = \
lib.ElMaxNormDistSparse_i.restype = \
lib.ElMaxNormDistSparse_s.restype = \
lib.ElMaxNormDistSparse_d.restype = \
lib.ElMaxNormDistSparse_c.restype = \
lib.ElMaxNormDistSparse_z.restype = \
lib.ElMaxNormDistMultiVec_i.restype = \
lib.ElMaxNormDistMultiVec_s.restype = \
lib.ElMaxNormDistMultiVec_d.restype = \
lib.ElMaxNormDistMultiVec_c.restype = \
lib.ElMaxNormDistMultiVec_z.restype = \
  c_uint

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
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElMaxNormSparse_i(*args)
    elif A.tag == sTag: lib.ElMaxNormSparse_s(*args)
    elif A.tag == dTag: lib.ElMaxNormSparse_d(*args)
    elif A.tag == cTag: lib.ElMaxNormSparse_c(*args)
    elif A.tag == zTag: lib.ElMaxNormSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElMaxNormDistSparse_i(*args)
    elif A.tag == sTag: lib.ElMaxNormDistSparse_s(*args)
    elif A.tag == dTag: lib.ElMaxNormDistSparse_d(*args)
    elif A.tag == cTag: lib.ElMaxNormDistSparse_c(*args)
    elif A.tag == zTag: lib.ElMaxNormDistSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElMaxNormDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElMaxNormDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElMaxNormDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElMaxNormDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElMaxNormDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value

lib.ElSymmetricMaxNorm_i.argtypes = \
lib.ElSymmetricMaxNormDist_i.argtypes = \
lib.ElSymmetricMaxNormSparse_i.argtypes = \
lib.ElSymmetricMaxNormDistSparse_i.argtypes = \
  [c_uint,c_void_p,POINTER(iType)]
lib.ElSymmetricMaxNorm_s.argtypes = \
lib.ElSymmetricMaxNorm_c.argtypes = \
lib.ElSymmetricMaxNormDist_s.argtypes = \
lib.ElSymmetricMaxNormDist_c.argtypes = \
lib.ElSymmetricMaxNormSparse_s.argtypes = \
lib.ElSymmetricMaxNormSparse_c.argtypes = \
lib.ElSymmetricMaxNormDistSparse_s.argtypes = \
lib.ElSymmetricMaxNormDistSparse_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricMaxNorm_d.argtypes = \
lib.ElSymmetricMaxNorm_z.argtypes = \
lib.ElSymmetricMaxNormDist_d.argtypes = \
lib.ElSymmetricMaxNormDist_z.argtypes = \
lib.ElSymmetricMaxNormSparse_d.argtypes = \
lib.ElSymmetricMaxNormSparse_z.argtypes = \
lib.ElSymmetricMaxNormDistSparse_d.argtypes = \
lib.ElSymmetricMaxNormDistSparse_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]

lib.ElSymmetricMaxNorm_i.restype = \
lib.ElSymmetricMaxNorm_s.restype = \
lib.ElSymmetricMaxNorm_d.restype = \
lib.ElSymmetricMaxNorm_c.restype = \
lib.ElSymmetricMaxNorm_z.restype = \
lib.ElSymmetricMaxNormDist_i.restype = \
lib.ElSymmetricMaxNormDist_s.restype = \
lib.ElSymmetricMaxNormDist_d.restype = \
lib.ElSymmetricMaxNormDist_c.restype = \
lib.ElSymmetricMaxNormDist_z.restype = \
lib.ElSymmetricMaxNormSparse_i.restype = \
lib.ElSymmetricMaxNormSparse_s.restype = \
lib.ElSymmetricMaxNormSparse_d.restype = \
lib.ElSymmetricMaxNormSparse_c.restype = \
lib.ElSymmetricMaxNormSparse_z.restype = \
lib.ElSymmetricMaxNormDistSparse_i.restype = \
lib.ElSymmetricMaxNormDistSparse_s.restype = \
lib.ElSymmetricMaxNormDistSparse_d.restype = \
lib.ElSymmetricMaxNormDistSparse_c.restype = \
lib.ElSymmetricMaxNormDistSparse_z.restype = \
  c_uint

def SymmetricMaxNorm(uplo,A):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSymmetricMaxNorm_i(*args)
    elif A.tag == sTag: lib.ElSymmetricMaxNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricMaxNorm_d(*args)
    elif A.tag == cTag: lib.ElSymmetricMaxNorm_c(*args)
    elif A.tag == zTag: lib.ElSymmetricMaxNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSymmetricMaxNormDist_i(*args)
    elif A.tag == sTag: lib.ElSymmetricMaxNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricMaxNormDist_d(*args)
    elif A.tag == cTag: lib.ElSymmetricMaxNormDist_c(*args)
    elif A.tag == zTag: lib.ElSymmetricMaxNormDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElSymmetricMaxNormSparse_i(*args)
    elif A.tag == sTag: lib.ElSymmetricMaxNormSparse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricMaxNormSparse_d(*args)
    elif A.tag == cTag: lib.ElSymmetricMaxNormSparse_c(*args)
    elif A.tag == zTag: lib.ElSymmetricMaxNormSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElSymmetricMaxNormDistSparse_i(*args)
    elif A.tag == sTag: lib.ElSymmetricMaxNormDistSparse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricMaxNormDistSparse_d(*args)
    elif A.tag == cTag: lib.ElSymmetricMaxNormDistSparse_c(*args)
    elif A.tag == zTag: lib.ElSymmetricMaxNormDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value
def HermitianMaxNorm(uplo,A):
  return SymmetricMaxNorm(uplo,A)

lib.ElNuclearNorm_s.argtypes = \
lib.ElNuclearNorm_c.argtypes = \
lib.ElNuclearNormDist_s.argtypes = \
lib.ElNuclearNormDist_c.argtypes = \
  [c_void_p,POINTER(sType)]
lib.ElNuclearNorm_d.argtypes = \
lib.ElNuclearNorm_z.argtypes = \
lib.ElNuclearNormDist_d.argtypes = \
lib.ElNuclearNormDist_z.argtypes = \
  [c_void_p,POINTER(dType)]
lib.ElNuclearNorm_s.restype = \
lib.ElNuclearNorm_d.restype = \
lib.ElNuclearNorm_c.restype = \
lib.ElNuclearNorm_z.restype = \
lib.ElNuclearNormDist_s.restype = \
lib.ElNuclearNormDist_d.restype = \
lib.ElNuclearNormDist_c.restype = \
lib.ElNuclearNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricNuclearNorm_s.argtypes = \
lib.ElSymmetricNuclearNorm_c.argtypes = \
lib.ElSymmetricNuclearNormDist_s.argtypes = \
lib.ElSymmetricNuclearNormDist_c.argtypes = \
lib.ElHermitianNuclearNorm_c.argtypes = \
lib.ElHermitianNuclearNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]
lib.ElSymmetricNuclearNorm_d.argtypes = \
lib.ElSymmetricNuclearNorm_z.argtypes = \
lib.ElSymmetricNuclearNormDist_d.argtypes = \
lib.ElSymmetricNuclearNormDist_z.argtypes = \
lib.ElHermitianNuclearNorm_z.argtypes = \
lib.ElHermitianNuclearNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]
lib.ElSymmetricNuclearNorm_s.restype = \
lib.ElSymmetricNuclearNorm_d.restype = \
lib.ElSymmetricNuclearNorm_c.restype = \
lib.ElSymmetricNuclearNorm_z.restype = \
lib.ElSymmetricNuclearNormDist_s.restype = \
lib.ElSymmetricNuclearNormDist_d.restype = \
lib.ElSymmetricNuclearNormDist_c.restype = \
lib.ElSymmetricNuclearNormDist_z.restype = \
lib.ElHermitianNuclearNorm_c.restype = \
lib.ElHermitianNuclearNorm_z.restype = \
lib.ElHermitianNuclearNormDist_c.restype = \
lib.ElHermitianNuclearNormDist_z.restype = \
  c_uint

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
  return norm.value
def HermitianNuclearNorm(uplo,A):
  return SymmetricNuclearNorm(uplo,A,True)

lib.ElOneNorm_s.argtypes = \
lib.ElOneNorm_c.argtypes = \
lib.ElOneNormDist_s.argtypes = \
lib.ElOneNormDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElOneNorm_d.argtypes = \
lib.ElOneNorm_z.argtypes = \
lib.ElOneNormDist_d.argtypes = \
lib.ElOneNormDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElOneNorm_s.restype = \
lib.ElOneNorm_d.restype = \
lib.ElOneNorm_c.restype = \
lib.ElOneNorm_z.restype = \
lib.ElOneNormDist_s.restype = \
lib.ElOneNormDist_d.restype = \
lib.ElOneNormDist_c.restype = \
lib.ElOneNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricOneNorm_s.argtypes = \
lib.ElSymmetricOneNorm_c.argtypes = \
lib.ElSymmetricOneNormDist_s.argtypes = \
lib.ElSymmetricOneNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]

lib.ElSymmetricOneNorm_d.argtypes = \
lib.ElSymmetricOneNorm_z.argtypes = \
lib.ElSymmetricOneNormDist_d.argtypes = \
lib.ElSymmetricOneNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]

lib.ElSymmetricOneNorm_s.restype = \
lib.ElSymmetricOneNorm_d.restype = \
lib.ElSymmetricOneNorm_c.restype = \
lib.ElSymmetricOneNorm_z.restype = \
lib.ElSymmetricOneNormDist_s.restype = \
lib.ElSymmetricOneNormDist_d.restype = \
lib.ElSymmetricOneNormDist_c.restype = \
lib.ElSymmetricOneNormDist_z.restype = \
  c_uint

def SymmetricOneNorm(uplo,A):
  norm = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(norm)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricOneNorm_s(*args)
    elif A.tag == dTag: lib.ElSymmetricOneNorm_d(*args)
    elif A.tag == cTag: lib.ElSymmetricOneNorm_c(*args)
    elif A.tag == zTag: lib.ElSymmetricOneNorm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricOneNormDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricOneNormDist_d(*args)
    elif A.tag == cTag: lib.ElSymmetricOneNormDist_c(*args)
    elif A.tag == zTag: lib.ElSymmetricOneNormDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value
def HermitianOneNorm(uplo,A):
  return SymmetricOneNorm(uplo,A)

lib.ElSchattenNorm_s.argtypes = \
lib.ElSchattenNorm_c.argtypes = \
lib.ElSchattenNormDist_s.argtypes = \
lib.ElSchattenNormDist_c.argtypes = \
  [c_void_p,sType,POINTER(sType)]

lib.ElSchattenNorm_d.argtypes = \
lib.ElSchattenNorm_z.argtypes = \
lib.ElSchattenNormDist_d.argtypes = \
lib.ElSchattenNormDist_z.argtypes = \
  [c_void_p,dType,POINTER(dType)]

lib.ElSchattenNorm_s.restype = \
lib.ElSchattenNorm_d.restype = \
lib.ElSchattenNorm_c.restype = \
lib.ElSchattenNorm_z.restype = \
lib.ElSchattenNormDist_s.restype = \
lib.ElSchattenNormDist_d.restype = \
lib.ElSchattenNormDist_c.restype = \
lib.ElSchattenNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricSchattenNorm_s.argtypes = \
lib.ElSymmetricSchattenNorm_c.argtypes = \
lib.ElSymmetricSchattenNormDist_s.argtypes = \
lib.ElSymmetricSchattenNormDist_c.argtypes = \
lib.ElHermitianSchattenNorm_c.argtypes = \
lib.ElHermitianSchattenNormDist_c.argtypes = \
  [c_uint,c_void_p,sType,POINTER(sType)]

lib.ElSymmetricSchattenNorm_d.argtypes = \
lib.ElSymmetricSchattenNorm_z.argtypes = \
lib.ElSymmetricSchattenNormDist_d.argtypes = \
lib.ElSymmetricSchattenNormDist_z.argtypes = \
lib.ElHermitianSchattenNorm_z.argtypes = \
lib.ElHermitianSchattenNormDist_z.argtypes = \
  [c_uint,c_void_p,dType,POINTER(dType)]

lib.ElSymmetricSchattenNorm_s.restype = \
lib.ElSymmetricSchattenNorm_d.restype = \
lib.ElSymmetricSchattenNorm_c.restype = \
lib.ElSymmetricSchattenNorm_z.restype = \
lib.ElSymmetricSchattenNormDist_s.restype = \
lib.ElSymmetricSchattenNormDist_d.restype = \
lib.ElSymmetricSchattenNormDist_c.restype = \
lib.ElSymmetricSchattenNormDist_z.restype = \
lib.ElHermitianSchattenNorm_c.restype = \
lib.ElHermitianSchattenNorm_z.restype = \
lib.ElHermitianSchattenNormDist_c.restype = \
lib.ElHermitianSchattenNormDist_z.restype = \
  c_uint

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
  return norm.value
def HermitianSchattenNorm(uplo,A,p):
  return SymmetricSchattenNorm(uplo,A,p,True)

lib.ElTwoNorm_s.argtypes = \
lib.ElTwoNorm_c.argtypes = \
lib.ElTwoNormDist_s.argtypes = \
lib.ElTwoNormDist_c.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElTwoNorm_d.argtypes = \
lib.ElTwoNorm_z.argtypes = \
lib.ElTwoNormDist_d.argtypes = \
lib.ElTwoNormDist_z.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElTwoNorm_s.restype = \
lib.ElTwoNorm_d.restype = \
lib.ElTwoNorm_c.restype = \
lib.ElTwoNorm_z.restype = \
lib.ElTwoNormDist_s.restype = \
lib.ElTwoNormDist_d.restype = \
lib.ElTwoNormDist_c.restype = \
lib.ElTwoNormDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricTwoNorm_s.argtypes = \
lib.ElSymmetricTwoNorm_c.argtypes = \
lib.ElSymmetricTwoNormDist_s.argtypes = \
lib.ElSymmetricTwoNormDist_c.argtypes = \
lib.ElHermitianTwoNorm_c.argtypes = \
lib.ElHermitianTwoNormDist_c.argtypes = \
  [c_uint,c_void_p,POINTER(sType)]

lib.ElSymmetricTwoNorm_d.argtypes = \
lib.ElSymmetricTwoNorm_z.argtypes = \
lib.ElSymmetricTwoNormDist_d.argtypes = \
lib.ElSymmetricTwoNormDist_z.argtypes = \
lib.ElHermitianTwoNorm_z.argtypes = \
lib.ElHermitianTwoNormDist_z.argtypes = \
  [c_uint,c_void_p,POINTER(dType)]

lib.ElSymmetricTwoNorm_s.restype = \
lib.ElSymmetricTwoNorm_d.restype = \
lib.ElSymmetricTwoNorm_c.restype = \
lib.ElSymmetricTwoNorm_z.restype = \
lib.ElSymmetricTwoNormDist_s.restype = \
lib.ElSymmetricTwoNormDist_d.restype = \
lib.ElSymmetricTwoNormDist_c.restype = \
lib.ElSymmetricTwoNormDist_z.restype = \
lib.ElHermitianTwoNorm_c.restype = \
lib.ElHermitianTwoNorm_z.restype = \
lib.ElHermitianTwoNormDist_c.restype = \
lib.ElHermitianTwoNormDist_z.restype = \
  c_uint

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
  return norm.value
def HermitianTwoNorm(uplo,A):
  return SymmetricTwoNorm(uplo,A,True)

lib.ElZeroNorm_i.argtypes = \
lib.ElZeroNormDist_i.argtypes = \
lib.ElZeroNormSparse_i.argtypes = \
lib.ElZeroNormDistSparse_i.argtypes = \
  [c_void_p,iType,POINTER(iType)]

lib.ElZeroNorm_s.argtypes = \
lib.ElZeroNormDist_s.argtypes = \
lib.ElZeroNormSparse_s.argtypes = \
lib.ElZeroNormDistSparse_s.argtypes = \
  [c_void_p,sType,POINTER(iType)]

lib.ElZeroNorm_d.argtypes = \
lib.ElZeroNormDist_d.argtypes = \
lib.ElZeroNormSparse_d.argtypes = \
lib.ElZeroNormDistSparse_d.argtypes = \
  [c_void_p,dType,POINTER(iType)]

lib.ElZeroNorm_c.argtypes = \
lib.ElZeroNormDist_c.argtypes = \
lib.ElZeroNormSparse_c.argtypes = \
lib.ElZeroNormDistSparse_c.argtypes = \
  [c_void_p,sType,POINTER(iType)]

lib.ElZeroNorm_z.argtypes = \
lib.ElZeroNormDist_z.argtypes = \
lib.ElZeroNormSparse_z.argtypes = \
lib.ElZeroNormDistSparse_z.argtypes = \
  [c_void_p,cType,POINTER(iType)]

lib.ElZeroNorm_i.restype = \
lib.ElZeroNorm_s.restype = \
lib.ElZeroNorm_d.restype = \
lib.ElZeroNorm_c.restype = \
lib.ElZeroNorm_z.restype = \
lib.ElZeroNormDist_i.restype = \
lib.ElZeroNormDist_s.restype = \
lib.ElZeroNormDist_d.restype = \
lib.ElZeroNormDist_c.restype = \
lib.ElZeroNormDist_z.restype = \
lib.ElZeroNormSparse_i.restype = \
lib.ElZeroNormSparse_s.restype = \
lib.ElZeroNormSparse_d.restype = \
lib.ElZeroNormSparse_c.restype = \
lib.ElZeroNormSparse_z.restype = \
lib.ElZeroNormDistSparse_i.restype = \
lib.ElZeroNormDistSparse_s.restype = \
lib.ElZeroNormDistSparse_d.restype = \
lib.ElZeroNormDistSparse_c.restype = \
lib.ElZeroNormDistSparse_z.restype = \
  c_uint

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
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElZeroNormSparse_i(*args)
    elif A.tag == sTag: lib.ElZeroNormSparse_s(*args)
    elif A.tag == dTag: lib.ElZeroNormSparse_d(*args)
    elif A.tag == cTag: lib.ElZeroNormSparse_c(*args)
    elif A.tag == zTag: lib.ElZeroNormSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElZeroNormDistSparse_i(*args)
    elif A.tag == sTag: lib.ElZeroNormDistSparse_s(*args)
    elif A.tag == dTag: lib.ElZeroNormDistSparse_d(*args)
    elif A.tag == cTag: lib.ElZeroNormDistSparse_c(*args)
    elif A.tag == zTag: lib.ElZeroNormDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return norm.value

lib.ElTwoNormEstimate_s.argtypes = \
lib.ElTwoNormEstimate_c.argtypes = \
lib.ElTwoNormEstimateDist_s.argtypes = \
lib.ElTwoNormEstimateDist_c.argtypes = \
  [c_void_p,sType,iType,POINTER(sType)]

lib.ElTwoNormEstimate_d.argtypes = \
lib.ElTwoNormEstimate_z.argtypes = \
lib.ElTwoNormEstimateDist_d.argtypes = \
lib.ElTwoNormEstimateDist_z.argtypes = \
  [c_void_p,dType,iType,POINTER(dType)]

lib.ElTwoNormEstimate_s.restype = \
lib.ElTwoNormEstimate_d.restype = \
lib.ElTwoNormEstimate_c.restype = \
lib.ElTwoNormEstimate_z.restype = \
lib.ElTwoNormEstimateDist_s.restype = \
lib.ElTwoNormEstimateDist_d.restype = \
lib.ElTwoNormEstimateDist_c.restype = \
lib.ElTwoNormEstimateDist_z.restype = \
  c_uint

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
  return norm.value

lib.ElSymmetricTwoNormEstimate_s.argtypes = \
lib.ElSymmetricTwoNormEstimate_c.argtypes = \
lib.ElSymmetricTwoNormEstimateDist_s.argtypes = \
lib.ElHermitianTwoNormEstimateDist_c.argtypes = \
lib.ElHermitianTwoNormEstimate_c.argtypes = \
lib.ElHermitianTwoNormEstimateDist_c.argtypes = \
  [c_uint,c_void_p,sType,iType,POINTER(sType)]

lib.ElSymmetricTwoNormEstimate_d.argtypes = \
lib.ElSymmetricTwoNormEstimate_z.argtypes = \
lib.ElSymmetricTwoNormEstimateDist_d.argtypes = \
lib.ElHermitianTwoNormEstimateDist_z.argtypes = \
lib.ElHermitianTwoNormEstimate_z.argtypes = \
lib.ElHermitianTwoNormEstimateDist_z.argtypes = \
  [c_uint,c_void_p,dType,iType,POINTER(dType)]

lib.ElSymmetricTwoNormEstimate_s.restype = \
lib.ElSymmetricTwoNormEstimate_d.restype = \
lib.ElSymmetricTwoNormEstimate_c.restype = \
lib.ElSymmetricTwoNormEstimate_z.restype = \
lib.ElSymmetricTwoNormEstimateDist_s.restype = \
lib.ElSymmetricTwoNormEstimateDist_d.restype = \
lib.ElHermitianTwoNormEstimateDist_c.restype = \
lib.ElHermitianTwoNormEstimateDist_z.restype = \
lib.ElHermitianTwoNormEstimate_c.restype = \
lib.ElHermitianTwoNormEstimate_z.restype = \
lib.ElHermitianTwoNormEstimateDist_c.restype = \
lib.ElHermitianTwoNormEstimateDist_z.restype = \
  c_uint

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
  return norm.value
def HermitianTwoNormEstimate(uplo,A,tol=1e-6,maxIts=100):
  return HermitianTwoNormEstimate(uplo,A,tol,maxits,True)

# Trace
# =====
lib.ElTrace_s.argtypes = \
lib.ElTraceDist_s.argtypes = \
  [c_void_p,POINTER(sType)]

lib.ElTrace_d.argtypes = \
lib.ElTraceDist_d.argtypes = \
  [c_void_p,POINTER(dType)]

lib.ElTrace_c.argtypes = \
lib.ElTraceDist_c.argtypes = \
  [c_void_p,POINTER(cType)]

lib.ElTrace_z.argtypes = \
lib.ElTraceDist_z.argtypes = \
  [c_void_p,POINTER(zType)]

lib.ElTrace_s.restype = \
lib.ElTrace_d.restype = \
lib.ElTrace_c.restype = \
lib.ElTrace_z.restype = \
lib.ElTraceDist_s.restype = \
lib.ElTraceDist_d.restype = \
lib.ElTraceDist_c.restype = \
lib.ElTraceDist_z.restype = \
  c_uint

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
  return trace.value
