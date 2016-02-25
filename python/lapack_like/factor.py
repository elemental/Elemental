#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
from perm import *
import ctypes

# Cholesky factorization
# ======================
lib.ElCholesky_s.argtypes = \
lib.ElCholesky_d.argtypes = \
lib.ElCholesky_c.argtypes = \
lib.ElCholesky_z.argtypes = \
lib.ElCholeskyDist_s.argtypes = \
lib.ElCholeskyDist_d.argtypes = \
lib.ElCholeskyDist_c.argtypes = \
lib.ElCholeskyDist_z.argtypes = \
  [c_uint,c_void_p]

lib.ElCholeskyPiv_s.argtypes = \
lib.ElCholeskyPiv_d.argtypes = \
lib.ElCholeskyPiv_c.argtypes = \
lib.ElCholeskyPiv_z.argtypes = \
lib.ElCholeskyPivDist_s.argtypes = \
lib.ElCholeskyPivDist_d.argtypes = \
lib.ElCholeskyPivDist_c.argtypes = \
lib.ElCholeskyPivDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]

def Cholesky(uplo,A,piv=False):
  if type(A) is Matrix:
    if piv:
      p = Permutation()
      args = [uplo,A.obj,p.obj]
      if   A.tag == sTag: lib.ElCholeskyPiv_s(*args)
      elif A.tag == dTag: lib.ElCholeskyPiv_d(*args)
      elif A.tag == cTag: lib.ElCholeskyPiv_c(*args)
      elif A.tag == zTag: lib.ElCholeskyPiv_z(*args)
      else: DataExcept()
      return p
    else:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElCholesky_s(*args)
      elif A.tag == dTag: lib.ElCholesky_d(*args)
      elif A.tag == cTag: lib.ElCholesky_c(*args)
      elif A.tag == zTag: lib.ElCholesky_z(*args)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if piv:
      p = DistPermutation(A.Grid())
      args = [uplo,A.obj,p.obj]
      if   A.tag == sTag: lib.ElCholeskyPivDist_s(*args)
      elif A.tag == dTag: lib.ElCholeskyPivDist_d(*args)
      elif A.tag == cTag: lib.ElCholeskyPivDist_c(*args)
      elif A.tag == zTag: lib.ElCholeskyPivDist_z(*args)
      else: DataExcept()
      return p
    else:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElCholeskyDist_s(*args)
      elif A.tag == dTag: lib.ElCholeskyDist_d(*args)
      elif A.tag == cTag: lib.ElCholeskyDist_c(*args)
      elif A.tag == zTag: lib.ElCholeskyDist_z(*args)
      else: DataExcept()
  else: TypeExcept()

lib.ElCholeskyMod_s.argtypes = \
lib.ElCholeskyMod_c.argtypes = \
lib.ElCholeskyModDist_s.argtypes = \
lib.ElCholeskyModDist_c.argtypes = \
  [c_uint,c_void_p,sType,c_void_p]

lib.ElCholeskyMod_d.argtypes = \
lib.ElCholeskyMod_z.argtypes = \
lib.ElCholeskyModDist_d.argtypes = \
lib.ElCholeskyModDist_z.argtypes = \
  [c_uint,c_void_p,dType,c_void_p]

def CholeskyMod(uplo,T,alpha,V):
  if type(T) is not type(V):
    raise Exception('Types of T and V must match')
  if T.tag != V.tag:
    raise Exception('Datatypes of T and V must match')
  args = [uplo,T.obj,alpha,V.obj]
  if type(T) is Matrix:
    if   T.tag == sTag: lib.ElCholeskyMod_s(*args)
    elif T.tag == dTag: lib.ElCholeskyMod_d(*args)
    elif T.tag == cTag: lib.ElCholeskyMod_c(*args)
    elif T.tag == zTag: lib.ElCholeskyMod_z(*args)
    else: DataExcept()
  elif type(T) is DistMatrix:
    if   T.tag == sTag: lib.ElCholeskyModDist_s(*args)
    elif T.tag == dTag: lib.ElCholeskyModDist_d(*args)
    elif T.tag == cTag: lib.ElCholeskyModDist_c(*args)
    elif T.tag == zTag: lib.ElCholeskyModDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElHPSDCholesky_s.argtypes = \
lib.ElHPSDCholesky_d.argtypes = \
lib.ElHPSDCholesky_c.argtypes = \
lib.ElHPSDCholesky_z.argtypes = \
lib.ElHPSDCholeskyDist_s.argtypes = \
lib.ElHPSDCholeskyDist_d.argtypes = \
lib.ElHPSDCholeskyDist_c.argtypes = \
lib.ElHPSDCholeskyDist_z.argtypes = \
  [c_uint,c_void_p]

def HPSDCholesky(uplo,A):
  args = [uplo,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPSDCholesky_s(*args)
    elif A.tag == dTag: lib.ElHPSDCholesky_d(*args)
    elif A.tag == cTag: lib.ElHPSDCholesky_c(*args)
    elif A.tag == zTag: lib.ElHPSDCholesky_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPSDCholeskyDist_s(*args)
    elif A.tag == dTag: lib.ElHPSDCholeskyDist_d(*args)
    elif A.tag == cTag: lib.ElHPSDCholeskyDist_c(*args)
    elif A.tag == zTag: lib.ElHPSDCholeskyDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElSolveAfterCholesky_s.argtypes = \
lib.ElSolveAfterCholesky_d.argtypes = \
lib.ElSolveAfterCholesky_c.argtypes = \
lib.ElSolveAfterCholesky_z.argtypes = \
lib.ElSolveAfterCholeskyDist_s.argtypes = \
lib.ElSolveAfterCholeskyDist_d.argtypes = \
lib.ElSolveAfterCholeskyDist_c.argtypes = \
lib.ElSolveAfterCholeskyDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p]

lib.ElSolveAfterCholeskyPiv_s.argtypes = \
lib.ElSolveAfterCholeskyPiv_d.argtypes = \
lib.ElSolveAfterCholeskyPiv_c.argtypes = \
lib.ElSolveAfterCholeskyPiv_z.argtypes = \
lib.ElSolveAfterCholeskyPivDist_s.argtypes = \
lib.ElSolveAfterCholeskyPivDist_d.argtypes = \
lib.ElSolveAfterCholeskyPivDist_c.argtypes = \
lib.ElSolveAfterCholeskyPivDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]

def SolveAfterCholesky(uplo,orient,A,B):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [uplo,orient,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterCholesky_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterCholesky_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterCholesky_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterCholesky_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterCholeskyDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterCholeskyDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterCholeskyDist_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterCholeskyDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

def SolveAfterCholeskyPiv(uplo,orient,A,p,B):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [uplo,orient,A.obj,p.obj,B.obj]
  if type(A) is Matrix:
    if type(p) is not Permutation:
      raise Exception('p must be a Permutation')
    if   A.tag == sTag: lib.ElSolveAfterCholeskyPiv_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterCholeskyPiv_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterCholeskyPiv_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterCholeskyPiv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(p) is not DistPermutation:
      raise Exception('p must be a DistPermutation')
    if   A.tag == sTag: lib.ElSolveAfterCholeskyPivDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterCholeskyPivDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterCholeskyPivDist_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterCholeskyPivDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# LDL
# ===

# Emulate an enum for LDL pivot types
(BUNCH_KAUFMAN_A,BUNCH_KAUFMAN_C,BUNCH_KAUFMAN_D,BUNCH_KAUFMAN_BOUNDED,
 BUNCH_PARLETT,LDL_WITHOUT_PIVOTING)=(0,1,2,3,4,5)

class LDLPivot(ctypes.Structure):
  _fields_ = [("nb",iType),("from",(iType*2))]

lib.ElLDLPivotConstant_s.argtypes = [c_uint,POINTER(sType)]
lib.ElLDLPivotConstant_d.argtypes = [c_uint,POINTER(dType)]
def LDLPivotConstant_s(pivotType):
  gamma = sType()
  lib.ElLDLPivotConstant_s(pivotType,pointer(gamma))
  return gamma
def LDLPivotConstant_d(pivotType):
  gamma = dType()
  lib.ElLDLPivotConstant_d(pivotType,pointer(gamma))
  return gamma

class LDLPivotCtrl_s(ctypes.Structure):
  _fields_ = [("pivotType",c_uint),("gamma",sType)]
  def __init__(self,pivType=BUNCH_KAUFMAN_A):
    pivotType = pivType
    gamma = LDLPivotConstant_s(pivType)
class LDLPivotCtrl_d(ctypes.Structure):
  _fields_ = [("pivotType",c_uint),("gamma",dType)]
  def __init__(self,pivType=BUNCH_KAUFMAN_A):
    pivotType = pivType
    gamma = LDLPivotConstant_d(pivType)

def TagToPivotCtrl(tag,pivType=BUNCH_KAUFMAN_A):
  if   tag == sTag: return LDLPivotCtrl_s(pivType)
  elif tag == dTag: return LDLPivotCtrl_d(pivType)
  elif tag == cTag: return LDLPivotCtrl_s(pivType)
  elif tag == zTag: return LDLPivotCtrl_d(pivType)
  else: DataExcept()

lib.ElLDL_s.argtypes = \
lib.ElLDL_d.argtypes = \
lib.ElLDLDist_s.argtypes = \
lib.ElLDLDist_d.argtypes = \
  [c_void_p]

lib.ElLDL_c.argtypes = \
lib.ElLDL_z.argtypes = \
lib.ElLDLDist_c.argtypes = \
lib.ElLDLDist_z.argtypes = \
  [c_void_p,bType]

lib.ElLDLPiv_s.argtypes = \
lib.ElLDLPivDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,LDLPivotCtrl_s]

lib.ElLDLPiv_d.argtypes = \
lib.ElLDLPivDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,LDLPivotCtrl_d]

lib.ElLDLPiv_c.argtypes = \
lib.ElLDLPivDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,bType,LDLPivotCtrl_s]

lib.ElLDLPiv_z.argtypes = \
lib.ElLDLPivDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,bType,LDLPivotCtrl_d]

def LDL(A,conjugate=True,pivType=BUNCH_KAUFMAN_A):
  if type(A) is Matrix:
    if pivType == LDL_WITHOUT_PIVOTING:
      args = [A.obj]
      argsCpx = [A.obj,conjugate]
      if   A.tag == sTag: lib.ElLDL_s(*args)
      elif A.tag == dTag: lib.ElLDL_d(*args)
      elif A.tag == cTag: lib.ElLDL_c(*argsCpx)
      elif A.tag == zTag: lib.ElLDL_z(*argsCpx)
      else: DataExcept()
    else:
      dSub = Matrix(A.tag)
      p = Permutation()
      ctrl = TagToPivotCtrl(A.tag,pivType)
      args = [A.obj,dSub.obj,p.obj,ctrl]
      argsCpx = [A.obj,dSub.obj,p.obj,conjugate,ctrl]
      if   A.tag == sTag: lib.ElLDLPiv_s(*args)
      elif A.tag == dTag: lib.ElLDLPiv_d(*args)
      elif A.tag == cTag: lib.ElLDLPiv_c(*argsCpx)
      elif A.tag == zTag: lib.ElLDLPiv_z(*argsCpx)
      else: DataExcept()
      return dSub, p
  elif type(A) is DistMatrix:
    if pivType == LDL_WITHOUT_PIVOTING:
      args = [A.obj]
      argsCpx = [A.obj,conjugate]
      if   A.tag == sTag: lib.ElLDLDist_s(*args)
      elif A.tag == dTag: lib.ElLDLDist_d(*args)
      elif A.tag == cTag: lib.ElLDLDist_c(*argsCpx)
      elif A.tag == zTag: lib.ElLDLDist_z(*argsCpx)
      else: DataExcept()
    else:
      dSub = DistMatrix(A.tag,VC,STAR,A.Grid())
      p = DistPermutation(A.Grid())
      ctrl = TagToPivotCtrl(A.tag,pivType)
      args = [A.obj,dSub.obj,p.obj,ctrl]
      argsCpx = [A.obj,dSub.obj,p.obj,conjugate,ctrl]
      if   A.tag == sTag: lib.ElLDLPivDist_s(*args)
      elif A.tag == dTag: lib.ElLDLPivDist_d(*args)
      elif A.tag == cTag: lib.ElLDLPivDist_c(*argsCpx)
      elif A.tag == zTag: lib.ElLDLPivDist_z(*argsCpx)
      else: DataExcept()
      return dSub, p
  else: TypeExcept()

lib.ElInertiaAfterLDL_s.argtypes = \
lib.ElInertiaAfterLDL_d.argtypes = \
lib.ElInertiaAfterLDL_c.argtypes = \
lib.ElInertiaAfterLDL_z.argtypes = \
lib.ElInertiaAfterLDLDist_s.argtypes = \
lib.ElInertiaAfterLDLDist_d.argtypes = \
lib.ElInertiaAfterLDLDist_c.argtypes = \
lib.ElInertiaAfterLDLDist_z.argtypes = \
  [c_void_p,c_void_p,POINTER(InertiaType)]

def InertiaAfterLDL(d,dSub):
  inertia = InertiaType()
  if type(d) is not type(dSub):
    raise Exception('Types of d and dSub must match')
  if d.tag != Base(dSub.tag):
    raise Exception('Datatype of d must be the base of that of dSub')
  args = [d.obj,dSub.obj,pointer(inertia)]
  if type(d) is Matrix:
    if   dSub.tag == sTag: lib.ElInertiaAfterLDL_s(*args)
    elif dSub.tag == dTag: lib.ElInertiaAfterLDL_d(*args)
    elif dSub.tag == cTag: lib.ElInertiaAfterLDL_c(*args)
    elif dSub.tag == zTag: lib.ElInertiaAfterLDL_z(*args)
    else: DataExcept()
  elif type(d) is DistMatrix:
    if   dSub.tag == sTag: lib.ElInertiaAfterLDLDist_s(*args)
    elif dSub.tag == dTag: lib.ElInertiaAfterLDLDist_d(*args)
    elif dSub.tag == cTag: lib.ElInertiaAfterLDLDist_c(*args)
    elif dSub.tag == zTag: lib.ElInertiaAfterLDLDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return inertia

lib.ElSolveAfterLDL_s.argtypes = \
lib.ElSolveAfterLDL_d.argtypes = \
lib.ElSolveAfterLDLDist_s.argtypes = \
lib.ElSolveAfterLDLDist_d.argtypes = \
  [c_void_p,c_void_p]

lib.ElSolveAfterLDL_c.argtypes = \
lib.ElSolveAfterLDL_z.argtypes = \
lib.ElSolveAfterLDLDist_c.argtypes = \
lib.ElSolveAfterLDLDist_z.argtypes = \
  [c_void_p,c_void_p,bType]

def SolveAfterLDL(A,B,conjugate=True):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [A.obj,B.obj]
  argsCpx = [A.obj,B.obj,conjugate]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterLDL_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLDL_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLDL_c(*argsCpx)
    elif A.tag == zTag: lib.ElSolveAfterLDL_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterLDLDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLDLDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLDLDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElSolveAfterLDLDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

def SolveAfterLDLPiv(A,dSub,p,B,conjugate=True):
  if type(A) is not type(dSub) or type(dSub) is not type(B):
    raise Exception('Types of {A,dSub,B} must match')
  if A.tag != dSub.tag or dSub.tag != B.tag:
    raise Exception('Datatypes of {A,dSub,B} must match')
  args = [A.obj,dSub.obj,p.obj,B.obj]
  argsCpx = [A.obj,dSub.obj,p.obj,B.obj,conjugate]
  if type(A) is Matrix:
    if type(p) is not Permutation:
      raise Exception('Expected p to be a Permutation')
    if   A.tag == sTag: lib.ElSolveAfterLDLPiv_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLDLPiv_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLDLPiv_c(*argsCpx)
    elif A.tag == zTag: lib.ElSolveAfterLDLPiv_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(p) is not DistPermutation:
      raise Exception('Expected p to be a DistPermutation')
    if   A.tag == sTag: lib.ElSolveAfterLDLPivDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLDLPivDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLDLPivDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElSolveAfterLDLPivDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

lib.ElMultiplyAfterLDL_s.argtypes = \
lib.ElMultiplyAfterLDL_d.argtypes = \
lib.ElMultiplyAfterLDLDist_s.argtypes = \
lib.ElMultiplyAfterLDLDist_d.argtypes = \
  [c_void_p,c_void_p]

lib.ElMultiplyAfterLDL_c.argtypes = \
lib.ElMultiplyAfterLDL_z.argtypes = \
lib.ElMultiplyAfterLDLDist_c.argtypes = \
lib.ElMultiplyAfterLDLDist_z.argtypes = \
  [c_void_p,c_void_p,bType]

def MultiplyAfterLDL(A,B,conjugate=True):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [A.obj,B.obj]
  argsCpx = [A.obj,B.obj,conjugate]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElMultiplyAfterLDL_s(*args)
    elif A.tag == dTag: lib.ElMultiplyAfterLDL_d(*args)
    elif A.tag == cTag: lib.ElMultiplyAfterLDL_c(*argsCpx)
    elif A.tag == zTag: lib.ElMultiplyAfterLDL_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElMultiplyAfterLDLDist_s(*args)
    elif A.tag == dTag: lib.ElMultiplyAfterLDLDist_d(*args)
    elif A.tag == cTag: lib.ElMultiplyAfterLDLDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElMultiplyAfterLDLDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

lib.ElMultiplyAfterLDLPiv_s.argtypes = \
lib.ElMultiplyAfterLDLPiv_d.argtypes = \
lib.ElMultiplyAfterLDLPivDist_s.argtypes = \
lib.ElMultiplyAfterLDLPivDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElMultiplyAfterLDLPiv_c.argtypes = \
lib.ElMultiplyAfterLDLPiv_z.argtypes = \
lib.ElMultiplyAfterLDLPivDist_c.argtypes = \
lib.ElMultiplyAfterLDLPivDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,bType]

def MultiplyAfterLDLPiv(A,dSub,p,B,conjugate=True):
  if type(A) is not type(dSub) or type(dSub) is not type(B):
    raise Exception('Types of {A,dSub,B} must match')
  if A.tag != dSub.tag or dSub.tag != B.tag:
    raise Exception('Datatypes of {A,dSub,B} must match')
  args = [A.obj,dSub.obj,p.obj,B.obj]
  argsCpx = [A.obj,dSub.obj,p.obj,B.obj,conjugate]
  if type(A) is Matrix:
    if type(p) is not Permutation:
      raise Exception('Expected p to be a Permutation')
    if   A.tag == sTag: lib.ElMultiplyAfterLDLPiv_s(*args)
    elif A.tag == dTag: lib.ElMultiplyAfterLDLPiv_d(*args)
    elif A.tag == cTag: lib.ElMultiplyAfterLDLPiv_c(*argsCpx)
    elif A.tag == zTag: lib.ElMultiplyAfterLDLPiv_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(p) is not DistPermutation:
      raise Exception('Expected p to be a DistPermutation')
    if   A.tag == sTag: lib.ElMultiplyAfterLDLPivDist_s(*args)
    elif A.tag == dTag: lib.ElMultiplyAfterLDLPivDist_d(*args)
    elif A.tag == cTag: lib.ElMultiplyAfterLDLPivDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElMultiplyAfterLDLPivDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# Regularized Quasi-Semidefinite LDL^H factorization
# ==================================================

# Emulate an enum for the Regularized LDL solve/refinement
(REG_SOLVE_FGMRES,REG_SOLVE_LGMRES)=(0,1)

lib.ElRegSolveCtrlDefault_s.argtypes = \
lib.ElRegSolveCtrlDefault_d.argtypes = \
  [c_void_p]
class RegSolveCtrl_s(ctypes.Structure):
  _fields_ = [("alg",c_uint),("relTol",sType),("relTolRefine",sType),
              ("maxIts",iType),("maxRefineIts",iType),("restart",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElRegSolveCtrlDefault_s(pointer(self))
class RegSolveCtrl_d(ctypes.Structure):
  _fields_ = [("alg",c_uint),("relTol",dType),("relTolRefine",dType),
              ("maxIts",iType),("maxRefineIts",iType),("restart",iType),
              ("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElRegSolveCtrlDefault_d(pointer(self))

# TODO: Wrappers for the factorization and solve

# LU factorization
# ================

# Emulate an enum for the pivot type for LU factorization
(LU_PARTIAL,LU_FULL,LU_ROOK,LU_WITHOUT_PIVOTING)=(0,1,2,3)

lib.ElLU_s.argtypes = \
lib.ElLU_d.argtypes = \
lib.ElLU_c.argtypes = \
lib.ElLU_z.argtypes = \
lib.ElLUDist_s.argtypes = \
lib.ElLUDist_d.argtypes = \
lib.ElLUDist_c.argtypes = \
lib.ElLUDist_z.argtypes = \
  [c_void_p]

lib.ElLUPartialPiv_s.argtypes = \
lib.ElLUPartialPiv_d.argtypes = \
lib.ElLUPartialPiv_c.argtypes = \
lib.ElLUPartialPiv_z.argtypes = \
lib.ElLUPartialPivDist_s.argtypes = \
lib.ElLUPartialPivDist_d.argtypes = \
lib.ElLUPartialPivDist_c.argtypes = \
lib.ElLUPartialPivDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElLUFullPiv_s.argtypes = \
lib.ElLUFullPiv_d.argtypes = \
lib.ElLUFullPiv_c.argtypes = \
lib.ElLUFullPiv_z.argtypes = \
lib.ElLUFullPivDist_s.argtypes = \
lib.ElLUFullPivDist_d.argtypes = \
lib.ElLUFullPivDist_c.argtypes = \
lib.ElLUFullPivDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

def LU(A,pivType=LU_PARTIAL):
  if type(A) is Matrix: 
    if pivType == LU_WITHOUT_PIVOTING:
      args = [A.obj]
      if   A.tag == sTag: lib.ElLU_s(*args)
      elif A.tag == dTag: lib.ElLU_d(*args)
      elif A.tag == cTag: lib.ElLU_c(*args)
      elif A.tag == zTag: lib.ElLU_z(*args)
      else: DataExcept()
    elif pivType == LU_PARTIAL:
      P = Permutation()
      args = [A.obj,P.obj]
      if   A.tag == sTag: lib.ElLUPartialPiv_s(*args)
      elif A.tag == dTag: lib.ElLUPartialPiv_d(*args)
      elif A.tag == cTag: lib.ElLUPartialPiv_c(*args)
      elif A.tag == zTag: lib.ElLUPartialPiv_z(*args)
      else: DataExcept()
      return P
    elif pivType == LU_FULL:
      P = Permutation()
      Q = Permutation()
      args = [A.obj,P.obj,Q.obj]
      if   A.tag == sTag: lib.ElLUFullPiv_s(*args)
      elif A.tag == dTag: lib.ElLUFullPiv_d(*args)
      elif A.tag == cTag: lib.ElLUFullPiv_c(*args)
      elif A.tag == zTag: lib.ElLUFullPiv_z(*args)
      else: DataExcept()
      return P, Q
    else: raise Exception('Unsupported pivot type')
  elif type(A) is DistMatrix:
    if pivType == LU_WITHOUT_PIVOTING:
      args = [A.obj]
      if   A.tag == sTag: lib.ElLUDist_s(*args)
      elif A.tag == dTag: lib.ElLUDist_d(*args)
      elif A.tag == cTag: lib.ElLUDist_c(*args)
      elif A.tag == zTag: lib.ElLUDist_z(*args)
      else: DataExcept()
    elif pivType == LU_PARTIAL:
      P = DistPermutation(A.Grid())
      args = [A.obj,P.obj]
      if   A.tag == sTag: lib.ElLUPartialPivDist_s(*args)
      elif A.tag == dTag: lib.ElLUPartialPivDist_d(*args)
      elif A.tag == cTag: lib.ElLUPartialPivDist_c(*args)
      elif A.tag == zTag: lib.ElLUPartialPivDist_z(*args)
      else: DataExcept()
      return P
    elif pivType == LU_FULL:
      P = DistPermutation(A.Grid())
      Q = DistPermutation(A.Grid())
      args = [A.obj,P.obj,Q.obj]
      if   A.tag == sTag: lib.ElLUFullPivDist_s(*args)
      elif A.tag == dTag: lib.ElLUFullPivDist_d(*args)
      elif A.tag == cTag: lib.ElLUFullPivDist_c(*args)
      elif A.tag == zTag: lib.ElLUFullPivDist_z(*args)
      else: DataExcept()
      return P, Q
    else: raise Exception('Unsupported pivot type')
  else: TypeExcept()

# LEFT OFF HERE

lib.ElLUMod_s.argtypes = \
lib.ElLUMod_c.argtypes = \
lib.ElLUModDist_s.argtypes = \
lib.ElLUModDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,sType]

lib.ElLUMod_d.argtypes = \
lib.ElLUMod_z.argtypes = \
lib.ElLUModDist_d.argtypes = \
lib.ElLUModDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,dType]

def LUMod(A,P,u,v,conjugate=True,tau=0.1):
  if type(A) is not type(u) or type(u) is not type(v):
    raise Exception('Types of {A,u,v} must be equal')
  if A.tag != u.tag or u.tag != v.tag:
    raise Exception('Datatypes of {A,u,v} must be equal')
  args = [A.obj,P.obj,u.obj,v.obj,tau]
  argsCpx = [A.obj,P.obj,u.obj,v.obj,conjugate,tau]
  if type(A) is Matrix:
    if type(P) is not Permutation:
      raise Exception('Expected P to be a Permutation')
    if   A.tag == sTag: lib.ElLUMod_s(*args)
    elif A.tag == dTag: lib.ElLUMod_d(*args)
    elif A.tag == cTag: lib.ElLUMod_c(*argsCpx)
    elif A.tag == zTag: lib.ElLUMod_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(P) is not DistPermutation:
      raise Exception('Expected P to be a DistPermutation')
    if   A.tag == sTag: lib.ElLUModDist_s(*args)
    elif A.tag == dTag: lib.ElLUModDist_d(*args)
    elif A.tag == cTag: lib.ElLUModDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElLUModDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

lib.ElSolveAfterLU_s.argtypes = \
lib.ElSolveAfterLU_d.argtypes = \
lib.ElSolveAfterLU_c.argtypes = \
lib.ElSolveAfterLU_z.argtypes = \
lib.ElSolveAfterLUDist_s.argtypes = \
lib.ElSolveAfterLUDist_d.argtypes = \
lib.ElSolveAfterLUDist_c.argtypes = \
lib.ElSolveAfterLUDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]

def SolveAfterLU(orient,A,B):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [orient,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterLU_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLU_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLU_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLU_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterLUDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLUDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLUDist_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLUDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElSolveAfterLUPartialPiv_s.argtypes = \
lib.ElSolveAfterLUPartialPiv_d.argtypes = \
lib.ElSolveAfterLUPartialPiv_c.argtypes = \
lib.ElSolveAfterLUPartialPiv_z.argtypes = \
lib.ElSolveAfterLUPartialPivDist_d.argtypes = \
lib.ElSolveAfterLUPartialPivDist_c.argtypes = \
lib.ElSolveAfterLUPartialPivDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p]

def SolveAfterLUPartialPiv(orient,A,P,B):
  if type(A) is not type(B):
    raise Exception('Types of {A,B} must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [orient,A.obj,P.obj,B.obj]
  if type(A) is Matrix:
    if type(P) is not Permutation:
      raise Exception('Expected P to be a Permutation')
    if   A.tag == sTag: lib.ElSolveAfterLUPartialPiv_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLUPartialPiv_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLUPartialPiv_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLUPartialPiv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(P) is not DistPermutation:
      raise Exception('Expected P to be a DistPermutation')
    if   A.tag == sTag: lib.ElSolveAfterLUPartialPivDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLUPartialPivDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLUPartialPivDist_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLUPartialPivDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElSolveAfterLUFullPiv_s.argtypes = \
lib.ElSolveAfterLUFullPiv_d.argtypes = \
lib.ElSolveAfterLUFullPiv_c.argtypes = \
lib.ElSolveAfterLUFullPiv_z.argtypes = \
lib.ElSolveAfterLUFullPivDist_d.argtypes = \
lib.ElSolveAfterLUFullPivDist_c.argtypes = \
lib.ElSolveAfterLUFullPivDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

def SolveAfterLUFullPiv(orient,A,P,Q,B):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [orient,A.obj,P.obj,Q.obj,B.obj]
  if type(A) is Matrix:
    if type(P) is not Permutation or type(Q) is not Permutation:
      raise Exception('Expected P and Q to be Permutations')
    if   A.tag == sTag: lib.ElSolveAfterLUFullPiv_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLUFullPiv_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLUFullPiv_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLUFullPiv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(P) is not DistPermutation or type(Q) is not DistPermutation:
      raise Exception('Expected P and Q to be DistPermutations')
    if   A.tag == sTag: lib.ElSolveAfterLUFullPivDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLUFullPivDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLUFullPivDist_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLUFullPivDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# LQ factorization
# ================
(LQ_IMPLICIT,LQ_EXPLICIT,LQ_EXPLICIT_TRIANG,LQ_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElLQ_s.argtypes = \
lib.ElLQ_d.argtypes = \
lib.ElLQ_c.argtypes = \
lib.ElLQ_z.argtypes = \
lib.ElLQDist_s.argtypes = \
lib.ElLQDist_d.argtypes = \
lib.ElLQDist_c.argtypes = \
lib.ElLQDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElLQExplicit_s.argtypes = \
lib.ElLQExplicit_d.argtypes = \
lib.ElLQExplicit_c.argtypes = \
lib.ElLQExplicit_z.argtypes = \
lib.ElLQExplicitDist_s.argtypes = \
lib.ElLQExplicitDist_d.argtypes = \
lib.ElLQExplicitDist_c.argtypes = \
lib.ElLQExplicitDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElLQExplicitUnitary_s.argtypes = \
lib.ElLQExplicitUnitary_d.argtypes = \
lib.ElLQExplicitUnitary_c.argtypes = \
lib.ElLQExplicitUnitary_z.argtypes = \
lib.ElLQExplicitUnitaryDist_s.argtypes = \
lib.ElLQExplicitUnitaryDist_d.argtypes = \
lib.ElLQExplicitUnitaryDist_c.argtypes = \
lib.ElLQExplicitUnitaryDist_z.argtypes = \
lib.ElLQExplicitTriang_s.argtypes = \
lib.ElLQExplicitTriang_d.argtypes = \
lib.ElLQExplicitTriang_c.argtypes = \
lib.ElLQExplicitTriang_z.argtypes = \
lib.ElLQExplicitTriangDist_s.argtypes = \
lib.ElLQExplicitTriangDist_d.argtypes = \
lib.ElLQExplicitTriangDist_c.argtypes = \
lib.ElLQExplicitTriangDist_z.argtypes = \
  [c_void_p]

def LQ(A,factType=LQ_IMPLICIT):
  if type(A) is Matrix:
    if factType == LQ_IMPLICIT:
      t = Matrix(A.tag)
      d = Matrix(Base(A.tag))
      args = [A.obj,t.obj,d.obj]
      if   A.tag == sTag: lib.ElLQ_s(*args)
      elif A.tag == dTag: lib.ElLQ_d(*args)
      elif A.tag == cTag: lib.ElLQ_c(*args)
      elif A.tag == zTag: lib.ElLQ_z(*args)
      else: DataExcept()
      return t, d
    elif factType == LQ_EXPLICIT:
      L = Matrix(A.tag)
      args = [L.obj,A.obj]
      if   A.tag == sTag: lib.ElLQExplicit_s(*args)
      elif A.tag == dTag: lib.ElLQExplicit_d(*args)
      elif A.tag == cTag: lib.ElLQExplicit_c(*args)
      elif A.tag == zTag: lib.ElLQExplicit_z(*args)
      else: DataExcept()
      return L
    elif factType == LQ_EXPLICIT_UNITARY:
      args = [A.obj]
      if   A.tag == sTag: lib.ElLQExplicitUnitary_s(*args)
      elif A.tag == dTag: lib.ElLQExplicitUnitary_d(*args)
      elif A.tag == cTag: lib.ElLQExplicitUnitary_c(*args)
      elif A.tag == zTag: lib.ElLQExplicitUnitary_z(*args)
      else: DataExcept()
    elif factType == LQ_EXPLICIT_TRIANG:
      args = [A.obj]
      if   A.tag == sTag: lib.ElLQExplicitTriang_s(*args)
      elif A.tag == dTag: lib.ElLQExplicitTriang_d(*args)
      elif A.tag == cTag: lib.ElLQExplicitTriang_c(*args)
      elif A.tag == zTag: lib.ElLQExplicitTriang_z(*args)
      else: DataExcept()
    else: raise Exception('Unsupported LQ factorization type')
  elif type(A) is DistMatrix:
    if factType == LQ_IMPLICIT:
      t = DistMatrix(A.tag,MC,STAR,A.Grid())
      d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
      args = [A.obj,t.obj,d.obj]
      if   A.tag == sTag: lib.ElLQDist_s(*args)
      elif A.tag == dTag: lib.ElLQDist_d(*args)
      elif A.tag == cTag: lib.ElLQDist_c(*args)
      elif A.tag == zTag: lib.ElLQDist_z(*args)
      else: DataExcept()
      return t, d
    elif factType == LQ_EXPLICIT:
      L = DistMatrix(A.tag,MC,MR,A.Grid())
      args = [L.obj,A.obj]
      if   A.tag == sTag: lib.ElLQExplicitDist_s(*args)
      elif A.tag == dTag: lib.ElLQExplicitDist_d(*args)
      elif A.tag == cTag: lib.ElLQExplicitDist_c(*args)
      elif A.tag == zTag: lib.ElLQExplicitDist_z(*args)
      else: DataExcept()
      return L
    elif factType == LQ_EXPLICIT_UNITARY:
      args = [A.obj]
      if   A.tag == sTag: lib.ElLQExplicitUnitaryDist_s(*args)
      elif A.tag == dTag: lib.ElLQExplicitUnitaryDist_d(*args)
      elif A.tag == cTag: lib.ElLQExplicitUnitaryDist_c(*args)
      elif A.tag == zTag: lib.ElLQExplicitUnitaryDist_z(*args)
      else: DataExcept()
    elif factType == LQ_EXPLICIT_TRIANG:
      args = [A.obj]
      if   A.tag == sTag: lib.ElLQExplicitTriangDist_s(*args)
      elif A.tag == dTag: lib.ElLQExplicitTriangDist_d(*args)
      elif A.tag == cTag: lib.ElLQExplicitTriangDist_c(*args)
      elif A.tag == zTag: lib.ElLQExplicitTriangDist_z(*args)
      else: DataExcept()
    else: raise Exception('Unsupported LQ factorization type')
  else: TypeExcept()

lib.ElApplyQAfterLQ_s.argtypes = \
lib.ElApplyQAfterLQ_d.argtypes = \
lib.ElApplyQAfterLQ_c.argtypes = \
lib.ElApplyQAfterLQ_z.argtypes = \
lib.ElApplyQAfterLQDist_s.argtypes = \
lib.ElApplyQAfterLQDist_d.argtypes = \
lib.ElApplyQAfterLQDist_c.argtypes = \
lib.ElApplyQAfterLQDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

def ApplyQAfterLQ(side,orient,A,t,d,B):
  if type(A) is not type(t) or type(t) is not type(d) or type(d) is not type(B):
    raise Exception('Matrix types of {A,t,d,B} must match')
  if A.tag != t.tag or t.tag != B.tag or d.tag != Base(A.tag):
    raise Exception('Datatypes of {A,t,B} must match and d must have base type')
  args = [side,orient,A.obj,t.obj,d.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElApplyQAfterLQ_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterLQ_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterLQ_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterLQ_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElApplyQAfterLQDist_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterLQDist_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterLQDist_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterLQDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElSolveAfterLQ_s.argtypes = \
lib.ElSolveAfterLQ_d.argtypes = \
lib.ElSolveAfterLQ_c.argtypes = \
lib.ElSolveAfterLQ_z.argtypes = \
lib.ElSolveAfterLQDist_s.argtypes = \
lib.ElSolveAfterLQDist_d.argtypes = \
lib.ElSolveAfterLQDist_c.argtypes = \
lib.ElSolveAfterLQDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]

def SolveAfterLQ(orient,A,t,d,B):
  if type(A) is Matrix:
    X = Matrix(A.tag)
    args = [orient,A.obj,t.obj,d.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElSolveAfterLQ_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLQ_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLQ_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLQ_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [orient,A.obj,t.obj,d.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElSolveAfterLQDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterLQDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterLQDist_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterLQDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

# QR factorization
# ================
lib.ElQRCtrlDefault_s.argtypes = \
lib.ElQRCtrlDefault_d.argtypes = \
  [c_void_p]
class QRCtrl_s(ctypes.Structure):
  _fields_ = [("colPiv",bType),("boundRank",bType),("maxRank",iType),
              ("adaptive",bType),("tol",sType),("alwaysRecomputeNorms",bType),
              ("smallestFirst",bType)]
  def __init__(self):
    lib.ElQRCtrlDefault_s(pointer(self))
class QRCtrl_d(ctypes.Structure):
  _fields_ = [("colPiv",bType),("boundRank",bType),("maxRank",iType),
              ("adaptive",bType),("tol",dType),("alwaysRecomputeNorms",bType),
              ("smallestFirst",bType)]
  def __init__(self):
    lib.ElQRCtrlDefault_d(pointer(self))

(QR_IMPLICIT,QR_EXPLICIT,QR_EXPLICIT_TRIANG,QR_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElQR_s.argtypes = \
lib.ElQR_d.argtypes = \
lib.ElQR_c.argtypes = \
lib.ElQR_z.argtypes = \
lib.ElQRDist_s.argtypes = \
lib.ElQRDist_d.argtypes = \
lib.ElQRDist_c.argtypes = \
lib.ElQRDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

lib.ElQRExplicit_s.argtypes = \
lib.ElQRExplicit_d.argtypes = \
lib.ElQRExplicit_c.argtypes = \
lib.ElQRExplicit_z.argtypes = \
lib.ElQRExplicitDist_s.argtypes = \
lib.ElQRExplicitDist_d.argtypes = \
lib.ElQRExplicitDist_c.argtypes = \
lib.ElQRExplicitDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElQRExplicitTriang_s.argtypes = \
lib.ElQRExplicitTriang_d.argtypes = \
lib.ElQRExplicitTriang_c.argtypes = \
lib.ElQRExplicitTriang_z.argtypes = \
lib.ElQRExplicitTriangDist_s.argtypes = \
lib.ElQRExplicitTriangDist_d.argtypes = \
lib.ElQRExplicitTriangDist_c.argtypes = \
lib.ElQRExplicitTriangDist_z.argtypes = \
  [c_void_p]

lib.ElQRExplicitUnitary_s.argtypes = \
lib.ElQRExplicitUnitary_d.argtypes = \
lib.ElQRExplicitUnitary_c.argtypes = \
lib.ElQRExplicitUnitary_z.argtypes = \
lib.ElQRExplicitUnitaryDist_s.argtypes = \
lib.ElQRExplicitUnitaryDist_d.argtypes = \
lib.ElQRExplicitUnitaryDist_c.argtypes = \
lib.ElQRExplicitUnitaryDist_z.argtypes = \
  [c_void_p]

lib.ElQRColPiv_s.argtypes = \
lib.ElQRColPiv_d.argtypes = \
lib.ElQRColPiv_c.argtypes = \
lib.ElQRColPiv_z.argtypes = \
lib.ElQRColPivDist_s.argtypes = \
lib.ElQRColPivDist_d.argtypes = \
lib.ElQRColPivDist_c.argtypes = \
lib.ElQRColPivDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElQRColPivX_s.argtypes = \
lib.ElQRColPivX_c.argtypes = \
lib.ElQRColPivXDist_s.argtypes = \
lib.ElQRColPivXDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]

lib.ElQRColPivX_d.argtypes = \
lib.ElQRColPivX_z.argtypes = \
lib.ElQRColPivXDist_d.argtypes = \
lib.ElQRColPivXDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]

lib.ElQRColPivExplicit_s.argtypes = \
lib.ElQRColPivExplicit_d.argtypes = \
lib.ElQRColPivExplicit_c.argtypes = \
lib.ElQRColPivExplicit_z.argtypes = \
lib.ElQRColPivExplicitDist_s.argtypes = \
lib.ElQRColPivExplicitDist_d.argtypes = \
lib.ElQRColPivExplicitDist_c.argtypes = \
lib.ElQRColPivExplicitDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

def QR(A,piv=False,factType=QR_IMPLICIT,ctrl=None):
  if type(A) is Matrix:
    if piv:
      if factType == QR_IMPLICIT:  
        t = Matrix(A.tag)
        d = Matrix(Base(A.tag))
        Omega = Permutation()
        args = [A.obj,t.obj,d.obj,Omega.obj]
        argsCtrl = [A.obj,t.obj,d.obj,Omega.obj,ctrl]
        if   A.tag == sTag:
          if ctrl == None: lib.ElQRColPivDist_s(*args)
          else:            lib.ElQRColPivXDist_s(*argsCtrl)
        elif A.tag == dTag:
          if ctrl == None: lib.ElQRColPivDist_d(*args)
          else:            lib.ElQRColPivXDist_d(*argsCtrl)
        elif A.tag == cTag:
          if ctrl == None: lib.ElQRColPivDist_c(*args)
          else:            lib.ElQRColPivXDist_c(*argsCtrl)
        elif A.tag == zTag:
          if ctrl == None: lib.ElQRColPivDist_z(*args)
          else:            lib.ElQRColPivXDist_z(*argsCtrl)
        else: DataExcept()
        return t, d, Omega
      elif factType == QR_EXPLICIT:
        if ctrl != None: 
          raise Exception('\'ctrl\' not yet supported for explicit piv fact\'s')
        R = Matrix(A.tag)
        Omega = Matrix(iTag)
        args = [A.obj,R.obj,Omega.obj]
        if   A.tag == sTag: lib.ElQRColPivExplicit_s(*args)
        elif A.tag == dTag: lib.ElQRColPivExplicit_d(*args)
        elif A.tag == cTag: lib.ElQRColPivExplicit_c(*args)
        elif A.tag == zTag: lib.ElQRColPivExplicit_z(*args)
        else: DataExcept()
        return R, Omega
      else: 
        raise Exception('Partial pivoted explicit fact\'s not yet supported')
    else:
      if factType == QR_IMPLICIT:
        t = Matrix(A.tag)
        d = Matrix(Base(A.tag))
        args = [A.obj,t.obj,d.obj]
        if   A.tag == sTag: lib.ElQR_s(*args)
        elif A.tag == dTag: lib.ElQR_d(*args)
        elif A.tag == cTag: lib.ElQR_c(*args)
        elif A.tag == zTag: lib.ElQR_z(*args)
        else: DataExcept()
        return t, d
      elif factType == QR_EXPLICIT:
        R = Matrix(A.tag)
        args = [A.obj,R.obj]
        if   A.tag == sTag: lib.ElQRExplicit_s(*args)
        elif A.tag == dTag: lib.ElQRExplicit_d(*args)
        elif A.tag == cTag: lib.ElQRExplicit_c(*args)
        elif A.tag == zTag: lib.ElQRExplicit_z(*args)
        else: DataExcept()
        return R
      elif factType == QR_EXPLICIT_TRIANG:
        args = [A.obj]
        if   A.tag == sTag: lib.ElQRExplicitTriang_s(*args)
        elif A.tag == dTag: lib.ElQRExplicitTriang_d(*args)
        elif A.tag == cTag: lib.ElQRExplicitTriang_c(*args)
        elif A.tag == zTag: lib.ElQRExplicitTriang_z(*args)
        else: DataExcept()
      elif factType == QR_EXPLICIT_UNITARY:
        args = [A.obj]
        if   A.tag == sTag: lib.ElQRExplicitUnitary_s(*args)
        elif A.tag == dTag: lib.ElQRExplicitUnitary_d(*args)
        elif A.tag == cTag: lib.ElQRExplicitUnitary_c(*args)
        elif A.tag == zTag: lib.ElQRExplicitUnitary_z(*args)
        else: DataExcept()
      else: raise Exception('Unsupported QR factorization type')
  elif type(A) is DistMatrix:
    if piv:
      if factType == QR_IMPLICIT:  
        t = DistMatrix(A.tag,MC,STAR,A.Grid())
        d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
        Omega = DistPermutation(A.Grid())
        args = [A.obj,t.obj,d.obj,Omega.obj]
        argsCtrl = [A.obj,t.obj,d.obj,Omega.obj,ctrl]
        if   A.tag == sTag:
          if ctrl == None: lib.ElQRColPivDist_s(*args)
          else:            lib.ElQRColPivXDist_s(*argsCtrl)
        elif A.tag == dTag:
          if ctrl == None: lib.ElQRColPivDist_d(*args)
          else:            lib.ElQRColPivXDist_d(*argsCtrl)
        elif A.tag == cTag:
          if ctrl == None: lib.ElQRColPivDist_c(*args)
          else:            lib.ElQRColPivXDist_c(*argsCtrl)
        elif A.tag == zTag:
          if ctrl == None: lib.ElQRColPivDist_z(*args)
          else:            lib.ElQRColPivXDist_z(*argsCtrl)
        else: DataExcept()
        return t, d, Omega
      elif factType == QR_EXPLICIT:
        if ctrl != None: 
          raise Exception('\'ctrl\' not yet supported for explicit piv fact\'s')
        R = DistMatrix(A.tag,MC,MR,A.Grid())
        Omega = DistMatrix(iTag,MC,MR,A.Grid())
        args = [A.obj,R.obj,Omega.obj]
        if   A.tag == sTag: lib.ElQRColPivExplicitDist_s(*args)
        elif A.tag == dTag: lib.ElQRColPivExplicitDist_d(*args)
        elif A.tag == cTag: lib.ElQRColPivExplicitDist_c(*args)
        elif A.tag == zTag: lib.ElQRColPivExplicitDist_z(*args)
        else: DataExcept()
        return R, Omega
      else: 
        raise Exception('Partial pivoted explicit fact\'s not yet supported')
    else:
      if factType == QR_IMPLICIT:
        t = DistMatrix(A.tag,MC,STAR,A.Grid())
        d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
        args = [A.obj,t.obj,d.obj]
        if   A.tag == sTag: lib.ElQRDist_s(*args)
        elif A.tag == dTag: lib.ElQRDist_d(*args)
        elif A.tag == cTag: lib.ElQRDist_c(*args)
        elif A.tag == zTag: lib.ElQRDist_z(*args)
        else: DataExcept()
        return t, d
      elif factType == QR_EXPLICIT:
        R = Matrix(A.tag)
        args = [A.obj,R.obj]
        if   A.tag == sTag: lib.ElQRExplicitDist_s(*args)
        elif A.tag == dTag: lib.ElQRExplicitDist_d(*args)
        elif A.tag == cTag: lib.ElQRExplicitDist_c(*args)
        elif A.tag == zTag: lib.ElQRExplicitDist_z(*args)
        else: DataExcept()
        return R
      elif factType == QR_EXPLICIT_TRIANG:
        args = [A.obj]
        if   A.tag == sTag: lib.ElQRExplicitTriangDist_s(*args)
        elif A.tag == dTag: lib.ElQRExplicitTriangDist_d(*args)
        elif A.tag == cTag: lib.ElQRExplicitTriangDist_c(*args)
        elif A.tag == zTag: lib.ElQRExplicitTriangDist_z(*args)
        else: DataExcept()
      elif factType == QR_EXPLICIT_UNITARY:
        args = [A.obj]
        if   A.tag == sTag: lib.ElQRExplicitUnitaryDist_s(*args)
        elif A.tag == dTag: lib.ElQRExplicitUnitaryDist_d(*args)
        elif A.tag == cTag: lib.ElQRExplicitUnitaryDist_c(*args)
        elif A.tag == zTag: lib.ElQRExplicitUnitaryDist_z(*args)
        else: DataExcept()
      else: raise Exception('Unsupported QR factorization type')
  else: TypeExcept()

lib.ElCholeskyQR_s.argtypes = \
lib.ElCholeskyQR_d.argtypes = \
lib.ElCholeskyQR_c.argtypes = \
lib.ElCholeskyQR_z.argtypes = \
lib.ElCholeskyQRDist_s.argtypes = \
lib.ElCholeskyQRDist_d.argtypes = \
lib.ElCholeskyQRDist_c.argtypes = \
lib.ElCholeskyQRDist_z.argtypes = \
  [c_void_p,c_void_p]

def CholeskyQR(A):
  if type(A) is Matrix:
    R = Matrix(A.tag)
    args = [A.obj,R.obj]
    if   A.tag == sTag: lib.ElCholeskyQR_s(*args)
    elif A.tag == dTag: lib.ElCholeskyQR_d(*args)
    elif A.tag == cTag: lib.ElCholeskyQR_c(*args)
    elif A.tag == zTag: lib.ElCholeskyQR_z(*args)
    else: DataExcept()
    return R
  elif type(A) is DistMatrix:
    R = DistMatrix(A.tag,STAR,STAR,A.Grid())
    args = [A.obj,R.obj]
    if   A.tag == sTag: lib.ElCholeskyQRDist_s(*args)
    elif A.tag == dTag: lib.ElCholeskyQRDist_d(*args)
    elif A.tag == cTag: lib.ElCholeskyQRDist_c(*args)
    elif A.tag == zTag: lib.ElCholeskyQRDist_z(*args)
    else: DataExcept()
    return R
  else: TypeExcept()

lib.ElApplyQAfterQR_s.argtypes = \
lib.ElApplyQAfterQR_d.argtypes = \
lib.ElApplyQAfterQR_c.argtypes = \
lib.ElApplyQAfterQR_z.argtypes = \
lib.ElApplyQAfterQRDist_s.argtypes = \
lib.ElApplyQAfterQRDist_d.argtypes = \
lib.ElApplyQAfterQRDist_c.argtypes = \
lib.ElApplyQAfterQRDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

def ApplyQAfterQR(side,orient,A,t,d,B):
  if type(A) is not type(t) or type(t) is not type(d) or type(d) is not type(B):
    raise Exception('Matrix types of {A,t,d,B} must match')
  if A.tag != t.tag or t.tag != B.tag:
    raise Exception('Datatypes of {A,t,B} must match')
  if d.tag != Base(A.tag):
    raise Exception('Base type of A must match that of d')
  args = [side,orient,A.obj,t.obj,d.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElApplyQAfterQR_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterQR_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterQR_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterQR_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElApplyQAfterQRDist_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterQRDist_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterQRDist_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterQRDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# TODO: TSQR
# TODO: ExplicitTSQR

# RQ factorization
# ================
(RQ_IMPLICIT,RQ_EXPLICIT,RQ_EXPLICIT_TRIANG,RQ_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElRQ_s.argtypes = \
lib.ElRQ_d.argtypes = \
lib.ElRQ_c.argtypes = \
lib.ElRQ_z.argtypes = \
lib.ElRQDist_s.argtypes = \
lib.ElRQDist_d.argtypes = \
lib.ElRQDist_c.argtypes = \
lib.ElRQDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

#lib.ElRQExplicit_s.argtypes = \
#lib.ElRQExplicit_d.argtypes = \
#lib.ElRQExplicit_c.argtypes = \
#lib.ElRQExplicit_z.argtypes = \
#lib.ElRQExplicitDist_s.argtypes = \
#lib.ElRQExplicitDist_d.argtypes = \
#lib.ElRQExplicitDist_c.argtypes = \
#lib.ElRQExplicitDist_z.argtypes = \
#  [c_void_p,c_void_p]

#lib.ElRQExplicitUnitary_s.argtypes = \
#lib.ElRQExplicitUnitary_d.argtypes = \
#lib.ElRQExplicitUnitary_c.argtypes = \
#lib.ElRQExplicitUnitary_z.argtypes = \
#lib.ElRQExplicitUnitaryDist_s.argtypes = \
#lib.ElRQExplicitUnitaryDist_d.argtypes = \
#lib.ElRQExplicitUnitaryDist_c.argtypes = \
#lib.ElRQExplicitUnitaryDist_z.argtypes = \
#  [c_void_p]

lib.ElRQExplicitTriang_s.argtypes = \
lib.ElRQExplicitTriang_d.argtypes = \
lib.ElRQExplicitTriang_c.argtypes = \
lib.ElRQExplicitTriang_z.argtypes = \
lib.ElRQExplicitTriangDist_s.argtypes = \
lib.ElRQExplicitTriangDist_d.argtypes = \
lib.ElRQExplicitTriangDist_c.argtypes = \
lib.ElRQExplicitTriangDist_z.argtypes = \
  [c_void_p]

def RQ(A,factType=RQ_IMPLICIT):
  if type(A) is Matrix:
    if factType == RQ_IMPLICIT:
      t = Matrix(A.tag)
      d = Matrix(Base(A.tag))
      args = [A.obj,t.obj,d.obj]
      if   A.tag == sTag: lib.ElRQ_s(*args)
      elif A.tag == dTag: lib.ElRQ_d(*args)
      elif A.tag == cTag: lib.ElRQ_c(*args)
      elif A.tag == zTag: lib.ElRQ_z(*args)
      else: DataExcept()
      return t, d
    elif factType == RQ_EXPLICIT:
      raise Exception('Explicit RQ factorization not yet supported')
      #R = Matrix(A.tag)
      #if   A.tag == sTag: lib.ElRQExplicit_s(R.obj,A.obj)
      #elif A.tag == dTag: lib.ElRQExplicit_d(R.obj,A.obj)
      #elif A.tag == cTag: lib.ElRQExplicit_c(R.obj,A.obj)
      #elif A.tag == zTag: lib.ElRQExplicit_z(R.obj,A.obj)
      #else: raise Exception('Unsupported datatype')
      #return R
    elif factType == RQ_EXPLICIT_UNITARY:
      raise Exception('Explicit unitary factor from RQ not yet supported')
      #if   A.tag == sTag: lib.ElRQExplicitUnitary_s(A.obj)
      #elif A.tag == dTag: lib.ElRQExplicitUnitary_d(A.obj)
      #elif A.tag == cTag: lib.ElRQExplicitUnitary_c(A.obj)
      #elif A.tag == zTag: lib.ElRQExplicitUnitary_z(A.obj)
      #else: raise Exception('Unsupported datatype')
    elif factType == RQ_EXPLICIT_TRIANG:
      args = [A.obj]
      if   A.tag == sTag: lib.ElRQExplicitTriang_s(*args)
      elif A.tag == dTag: lib.ElRQExplicitTriang_d(*args)
      elif A.tag == cTag: lib.ElRQExplicitTriang_c(*args)
      elif A.tag == zTag: lib.ElRQExplicitTriang_z(*args)
      else: DataExcept()
    else: raise Exception('Unsupported RQ factorization type')
  elif type(A) is DistMatrix:
    if factType == RQ_IMPLICIT:
      t = DistMatrix(A.tag,MC,STAR,A.Grid())
      d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
      args = [A.obj,t.obj,d.obj]
      if   A.tag == sTag: lib.ElRQDist_s(*args)
      elif A.tag == dTag: lib.ElRQDist_d(*args)
      elif A.tag == cTag: lib.ElRQDist_c(*args)
      elif A.tag == zTag: lib.ElRQDist_z(*args)
      else: DataExcept()
      return t, d
    elif factType == RQ_EXPLICIT:
      raise Exception('Explicit RQ factorization not yet supported')
      #R = DistMatrix(A.tag,MC,MR,A.Grid())
      #if   A.tag == sTag: lib.ElRQExplicitDist_s(R.obj,A.obj)
      #elif A.tag == dTag: lib.ElRQExplicitDist_d(R.obj,A.obj)
      #elif A.tag == cTag: lib.ElRQExplicitDist_c(R.obj,A.obj)
      #elif A.tag == zTag: lib.ElRQExplicitDist_z(R.obj,A.obj)
      #else: raise Exception('Unsupported datatype')
      #return R
    elif factType == RQ_EXPLICIT_UNITARY:
      raise Exception('Explicit unitary factor from RQ not yet supported')
      #if   A.tag == sTag: lib.ElRQExplicitUnitaryDist_s(A.obj)
      #elif A.tag == dTag: lib.ElRQExplicitUnitaryDist_d(A.obj)
      #elif A.tag == cTag: lib.ElRQExplicitUnitaryDist_c(A.obj)
      #elif A.tag == zTag: lib.ElRQExplicitUnitaryDist_z(A.obj)
      #else: raise Exception('Unsupported datatype')
    elif factType == RQ_EXPLICIT_TRIANG:
      args = [A.obj]
      if   A.tag == sTag: lib.ElRQExplicitTriangDist_s(*args)
      elif A.tag == dTag: lib.ElRQExplicitTriangDist_d(*args)
      elif A.tag == cTag: lib.ElRQExplicitTriangDist_c(*args)
      elif A.tag == zTag: lib.ElRQExplicitTriangDist_z(*args)
      else: DataExcept()
    else: raise Exception('Unsupported RQ factorization type')
  else: TypeExcept()

lib.ElApplyQAfterRQ_s.argtypes = \
lib.ElApplyQAfterRQ_d.argtypes = \
lib.ElApplyQAfterRQ_c.argtypes = \
lib.ElApplyQAfterRQ_z.argtypes = \
lib.ElApplyQAfterRQDist_s.argtypes = \
lib.ElApplyQAfterRQDist_d.argtypes = \
lib.ElApplyQAfterRQDist_c.argtypes = \
lib.ElApplyQAfterRQDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

def ApplyQAfterRQ(side,orient,A,t,d,B):
  if type(A) is not type(t) or type(t) is not type(d) or type(d) is not type(B):
    raise Exception('Matrix types of {A,t,d,B} must match')
  if A.tag != t.tag or t.tag != B.tag or d.tag != Base(A.tag):
    raise Exception('Datatypes of {A,t,B} must match and d must have base type')
  args = [side,orient,A.obj,t.obj,d.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElApplyQAfterRQ_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterRQ_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterRQ_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterRQ_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElApplyQAfterRQDist_s(*args)
    elif A.tag == dTag: lib.ElApplyQAfterRQDist_d(*args)
    elif A.tag == cTag: lib.ElApplyQAfterRQDist_c(*args)
    elif A.tag == zTag: lib.ElApplyQAfterRQDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElSolveAfterRQ_s.argtypes = \
lib.ElSolveAfterRQ_d.argtypes = \
lib.ElSolveAfterRQ_c.argtypes = \
lib.ElSolveAfterRQ_z.argtypes = \
lib.ElSolveAfterRQDist_s.argtypes = \
lib.ElSolveAfterRQDist_d.argtypes = \
lib.ElSolveAfterRQDist_c.argtypes = \
lib.ElSolveAfterRQDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]

def SolveAfterRQ(orient,A,t,d,B):
  if type(A) is Matrix:
    X = Matrix(A.tag)
    args = [orient,A.obj,t.obj,d.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElSolveAfterRQ_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterRQ_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterRQ_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterRQ_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [orient,A.obj,t.obj,d.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElSolveAfterRQDist_s(*args)
    elif A.tag == dTag: lib.ElSolveAfterRQDist_d(*args)
    elif A.tag == cTag: lib.ElSolveAfterRQDist_c(*args)
    elif A.tag == zTag: lib.ElSolveAfterRQDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

# Generalized QR factorization
# ============================
(GQR_IMPLICIT,GQR_EXPLICIT,GQR_EXPLICIT_TRIANG,GQR_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElGQR_s.argtypes = \
lib.ElGQR_d.argtypes = \
lib.ElGQR_c.argtypes = \
lib.ElGQR_z.argtypes = \
lib.ElGQRDist_s.argtypes = \
lib.ElGQRDist_d.argtypes = \
lib.ElGQRDist_c.argtypes = \
lib.ElGQRDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElGQRExplicitTriang_s.argtypes = \
lib.ElGQRExplicitTriang_d.argtypes = \
lib.ElGQRExplicitTriang_c.argtypes = \
lib.ElGQRExplicitTriang_z.argtypes = \
lib.ElGQRExplicitTriangDist_s.argtypes = \
lib.ElGQRExplicitTriangDist_d.argtypes = \
lib.ElGQRExplicitTriangDist_c.argtypes = \
lib.ElGQRExplicitTriangDist_z.argtypes = \
  [c_void_p,c_void_p]

def GQR(A,B,factType=GQR_IMPLICIT):
  if type(A) is Matrix:
    if factType == GQR_IMPLICIT:
      tA = Matrix(A.tag)
      dA = Matrix(Base(A.tag))
      tB = Matrix(A.tag)
      dB = Matrix(Base(A.tag))
      args = [A.obj,tA.obj,dA.obj,tB.obj,dB.obj]
      if   A.tag == sTag: lib.ElGQR_s(*args)
      elif A.tag == dTag: lib.ElGQR_d(*args)
      elif A.tag == cTag: lib.ElGQR_c(*args)
      elif A.tag == zTag: lib.ElGQR_z(*args)
      else: DataExcept()
      return tA, dA, tB, dB
    elif factType == GQR_EXPLICIT:
      raise Exception('GQR_EXPLICIT not yet supported')
    elif factType == GQR_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      args = [A.obj,B.obj]
      if   A.tag == sTag: lib.ElGQRExplicitTriang_s(*args)
      elif A.tag == dTag: lib.ElGQRExplicitTriang_d(*args)
      elif A.tag == cTag: lib.ElGQRExplicitTriang_c(*args)
      elif A.tag == zTag: lib.ElGQRExplicitTriang_z(*args)
      else: DataExcept()
      return B
    elif factType == GQR_EXPLICIT_UNITARY:
      raise Exception('GQR_EXPLICIT_UNITARY not yet supported') 
    else: raise Exception('Unsupported GQR factorization type')
  elif type(A) is DistMatrix:
    if factType == GQR_IMPLICIT:
      tA = Matrix(A.tag)
      dA = Matrix(Base(A.tag))
      tB = Matrix(A.tag)
      dB = Matrix(Base(A.tag))
      args = [A.obj,tA.obj,dA.obj,tB.obj,dB.obj]
      if   A.tag == sTag: lib.ElGQRDist_s(*args)
      elif A.tag == dTag: lib.ElGQRDist_d(*args)
      elif A.tag == cTag: lib.ElGQRDist_c(*args)
      elif A.tag == zTag: lib.ElGQRDist_z(*args)
      else: DataExcept()
      return tA, dA, tB, dB
    elif factType == GQR_EXPLICIT:
      raise Exception('GQR_EXPLICIT not yet supported')
    elif factType == GQR_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      args = [A.obj,B.obj]
      if   A.tag == sTag: lib.ElGQRExplicitTriangDist_s(*args)
      elif A.tag == dTag: lib.ElGQRExplicitTriangDist_d(*args)
      elif A.tag == cTag: lib.ElGQRExplicitTriangDist_c(*args)
      elif A.tag == zTag: lib.ElGQRExplicitTriangDist_z(*args)
      else: DataExcept()
      return B
    elif factType == GQR_EXPLICIT_UNITARY:
      raise Exception('GQR_EXPLICIT_UNITARY not yet supported') 
    else: raise Exception('Unsupported GQR factorization type')
  else: TypeExcept()

# Generalized RQ factorization
# ============================
(GRQ_IMPLICIT,GRQ_EXPLICIT,GRQ_EXPLICIT_TRIANG,GRQ_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElGRQ_s.argtypes = \
lib.ElGRQ_d.argtypes = \
lib.ElGRQ_c.argtypes = \
lib.ElGRQ_z.argtypes = \
lib.ElGRQDist_s.argtypes = \
lib.ElGRQDist_d.argtypes = \
lib.ElGRQDist_c.argtypes = \
lib.ElGRQDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElGRQExplicitTriang_s.argtypes = \
lib.ElGRQExplicitTriang_d.argtypes = \
lib.ElGRQExplicitTriang_c.argtypes = \
lib.ElGRQExplicitTriang_z.argtypes = \
lib.ElGRQExplicitTriangDist_s.argtypes = \
lib.ElGRQExplicitTriangDist_d.argtypes = \
lib.ElGRQExplicitTriangDist_c.argtypes = \
lib.ElGRQExplicitTriangDist_z.argtypes = \
  [c_void_p,c_void_p]

def GRQ(A,B,factType=GRQ_IMPLICIT):
  if type(A) is Matrix:
    if factType == GRQ_IMPLICIT:
      tA = Matrix(A.tag)
      dA = Matrix(Base(A.tag))
      tB = Matrix(A.tag)
      dB = Matrix(Base(A.tag))
      args = [A.obj,tA.obj,dA.obj,tB.obj,dB.obj]
      if   A.tag == sTag: lib.ElGRQ_s(*args)
      elif A.tag == dTag: lib.ElGRQ_d(*args)
      elif A.tag == cTag: lib.ElGRQ_c(*args)
      elif A.tag == zTag: lib.ElGRQ_z(*args)
      else: DataExcept()
      return tA, dA, tB, dB
    elif factType == GRQ_EXPLICIT:
      raise Exception('GRQ_EXPLICIT not yet supported')
    elif factType == GRQ_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      args = [A.obj,B.obj]
      if   A.tag == sTag: lib.ElGRQExplicitTriang_s(*args)
      elif A.tag == dTag: lib.ElGRQExplicitTriang_d(*args)
      elif A.tag == cTag: lib.ElGRQExplicitTriang_c(*args)
      elif A.tag == zTag: lib.ElGRQExplicitTriang_z(*args)
      else: DataExcept()
      return B
    elif factType == GRQ_EXPLICIT_UNITARY:
      raise Exception('GRQ_EXPLICIT_UNITARY not yet supported') 
    else: raise Exception('Unsupported GRQ factorization type')
  elif type(A) is DistMatrix:
    if factType == GRQ_IMPLICIT:
      tA = Matrix(A.tag)
      dA = Matrix(Base(A.tag))
      tB = Matrix(A.tag)
      dB = Matrix(Base(A.tag))
      args = [A.obj,tA.obj,dA.obj,tB.obj,dB.obj]
      if   A.tag == sTag: lib.ElGRQDist_s(*args)
      elif A.tag == dTag: lib.ElGRQDist_d(*args)
      elif A.tag == cTag: lib.ElGRQDist_c(*args)
      elif A.tag == zTag: lib.ElGRQDist_z(*args)
      else: DataExcept()
      return tA, dA, tB, dB
    elif factType == GRQ_EXPLICIT:
      raise Exception('GRQ_EXPLICIT not yet supported')
    elif factType == GRQ_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      args = [A.obj,B.obj]
      if   A.tag == sTag: lib.ElGRQExplicitTriangDist_s(*args)
      elif A.tag == dTag: lib.ElGRQExplicitTriangDist_d(*args)
      elif A.tag == cTag: lib.ElGRQExplicitTriangDist_c(*args)
      elif A.tag == zTag: lib.ElGRQExplicitTriangDist_z(*args)
      else: DataExcept()
      return B
    elif factType == GRQ_EXPLICIT_UNITARY:
      raise Exception('GRQ_EXPLICIT_UNITARY not yet supported') 
    else: raise Exception('Unsupported GRQ factorization type')
  else: TypeExcept()

# Interpolative decomposition
# ===========================
lib.ElID_s.argtypes = \
lib.ElID_c.argtypes = \
lib.ElIDDist_s.argtypes = \
lib.ElIDDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,QRCtrl_s,bType]

lib.ElID_d.argtypes = \
lib.ElID_z.argtypes = \
lib.ElIDDist_d.argtypes = \
lib.ElIDDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,QRCtrl_d,bType]

def ID(A,ctrl,canOverwrite=False):
  if type(A) is Matrix: 
    Omega = Permutation()
    Z = Matrix(A.tag)
    args = [A.obj,Omega.obj,Z.obj,ctrl,canOverwrite]
    if   A.tag == sTag: lib.ElID_s(*args)
    elif A.tag == dTag: lib.ElID_d(*args)
    elif A.tag == cTag: lib.ElID_c(*args)
    elif A.tag == zTag: lib.ElID_z(*args)
    else: DataExcept()
    return Omega, Z
  elif type(A) is DistMatrix:
    Omega = DistPermutation(A.Grid())
    Z = DistMatrix(A.tag,STAR,VR,A.Grid())
    args = [A.obj,Omega.obj,Z.obj,ctrl,canOverwrite]
    if   A.tag == sTag: lib.ElIDDist_s(*args)
    elif A.tag == dTag: lib.ElIDDist_d(*args)
    elif A.tag == cTag: lib.ElIDDist_c(*args)
    elif A.tag == zTag: lib.ElIDDist_z(*args)
    else: DataExcept()
    return Omega, Z
  else: TypeExcept()

# Skeleton decomposition
# ======================
lib.ElSkeleton_s.argtypes = \
lib.ElSkeleton_c.argtypes = \
lib.ElSkeletonDist_s.argtypes = \
lib.ElSkeletonDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]

lib.ElSkeleton_d.argtypes = \
lib.ElSkeleton_z.argtypes = \
lib.ElSkeletonDist_d.argtypes = \
lib.ElSkeletonDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]

def Skeleton(A,ctrl):
  if type(A) is Matrix:
    PR = Permutation()
    PC = Permutation()
    Z = Matrix(A.tag)
    args = [A.obj,PR.obj,PC.obj,Z.obj,ctrl]
    if   A.tag == sTag: lib.ElSkeleton_s(*args)
    elif A.tag == dTag: lib.ElSkeleton_d(*args)
    elif A.tag == cTag: lib.ElSkeleton_c(*args)
    elif A.tag == zTag: lib.ElSkeleton_z(*args)
    else: DataExcept()
    return PR, PC, Z
  elif type(A) is DistMatrix:
    PR = DistPermutation(A.Grid())
    PC = DistPermutation(A.Grid())
    Z = DistMatrix(A.tag,STAR,VR,A.Grid())
    args = [A.obj,PR.obj,PC.obj,Z.obj,ctrl]
    if   A.tag == sTag: lib.ElSkeletonDist_s(*args)
    elif A.tag == dTag: lib.ElSkeletonDist_d(*args)
    elif A.tag == cTag: lib.ElSkeletonDist_c(*args)
    elif A.tag == zTag: lib.ElSkeletonDist_z(*args)
    else: DataExcept()
    return PR, PC, Z
  else: TypeExcept()
