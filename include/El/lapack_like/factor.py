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

# Cholesky factorization
# ======================
lib.ElCholesky_s.argtypes = [c_uint,c_void_p]
lib.ElCholesky_s.restype = c_uint
lib.ElCholesky_d.argtypes = [c_uint,c_void_p]
lib.ElCholesky_d.restype = c_uint
lib.ElCholesky_c.argtypes = [c_uint,c_void_p]
lib.ElCholesky_c.restype = c_uint
lib.ElCholesky_z.argtypes = [c_uint,c_void_p]
lib.ElCholesky_z.restype = c_uint
lib.ElCholeskyDist_s.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_s.restype = c_uint
lib.ElCholeskyDist_d.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_d.restype = c_uint
lib.ElCholeskyDist_c.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_c.restype = c_uint
lib.ElCholeskyDist_z.argtypes = [c_uint,c_void_p]
lib.ElCholeskyDist_z.restype = c_uint
lib.ElCholeskyPiv_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_s.restype = c_uint
lib.ElCholeskyPiv_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_d.restype = c_uint
lib.ElCholeskyPiv_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_c.restype = c_uint
lib.ElCholeskyPiv_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPiv_z.restype = c_uint
lib.ElCholeskyPivDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_s.restype = c_uint
lib.ElCholeskyPivDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_d.restype = c_uint
lib.ElCholeskyPivDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_c.restype = c_uint
lib.ElCholeskyPivDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElCholeskyPivDist_z.restype = c_uint

def Cholesky(uplo,A,piv=False):
  if type(A) is Matrix:
    if piv:
      p = Matrix(iTag)
      if   A.tag == sTag: lib.ElCholeskyPiv_s(uplo,A.obj,p.obj)
      elif A.tag == dTag: lib.ElCholeskyPiv_d(uplo,A.obj,p.obj)
      elif A.tag == cTag: lib.ElCholeskyPiv_c(uplo,A.obj,p.obj)
      elif A.tag == zTag: lib.ElCholeskyPiv_z(uplo,A.obj,p.obj)
      else: raise Exception('Unsupported datatype')
      return p
    else:
      if   A.tag == sTag: lib.ElCholesky_s(uplo,A.obj)
      elif A.tag == dTag: lib.ElCholesky_d(uplo,A.obj)
      elif A.tag == cTag: lib.ElCholesky_c(uplo,A.obj)
      elif A.tag == zTag: lib.ElCholesky_z(uplo,A.obj)
      else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if piv:
      p = DistMatrix(iTag,VC,STAR,A.Grid())
      if   A.tag == sTag: lib.ElCholeskyPivDist_s(uplo,A.obj,p.obj)
      elif A.tag == dTag: lib.ElCholeskyPivDist_d(uplo,A.obj,p.obj)
      elif A.tag == cTag: lib.ElCholeskyPivDist_c(uplo,A.obj,p.obj)
      elif A.tag == zTag: lib.ElCholeskyPivDist_z(uplo,A.obj,p.obj)
      else: raise Exception('Unsupported datatype')
      return p
    else:
      if   A.tag == sTag: lib.ElCholeskyDist_s(uplo,A.obj)
      elif A.tag == dTag: lib.ElCholeskyDist_d(uplo,A.obj)
      elif A.tag == cTag: lib.ElCholeskyDist_c(uplo,A.obj)
      elif A.tag == zTag: lib.ElCholeskyDist_z(uplo,A.obj)
      else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElCholeskyMod_s.argtypes = [c_uint,c_void_p,sType,c_void_p]
lib.ElCholeskyMod_s.restype = c_uint
lib.ElCholeskyMod_d.argtypes = [c_uint,c_void_p,dType,c_void_p]
lib.ElCholeskyMod_d.restype = c_uint
lib.ElCholeskyMod_c.argtypes = [c_uint,c_void_p,sType,c_void_p]
lib.ElCholeskyMod_c.restype = c_uint
lib.ElCholeskyMod_z.argtypes = [c_uint,c_void_p,dType,c_void_p]
lib.ElCholeskyMod_z.restype = c_uint
lib.ElCholeskyModDist_s.argtypes = [c_uint,c_void_p,sType,c_void_p]
lib.ElCholeskyModDist_s.restype = c_uint
lib.ElCholeskyModDist_d.argtypes = [c_uint,c_void_p,dType,c_void_p]
lib.ElCholeskyModDist_d.restype = c_uint
lib.ElCholeskyModDist_c.argtypes = [c_uint,c_void_p,sType,c_void_p]
lib.ElCholeskyModDist_c.restype = c_uint
lib.ElCholeskyModDist_z.argtypes = [c_uint,c_void_p,dType,c_void_p]
lib.ElCholeskyModDist_z.restype = c_uint
def CholeskyMod(uplo,T,alpha,V):
  if type(T) is not type(V):
    raise Exception('Types of T and V must match')
  if T.tag != V.tag:
    raise Exception('Datatypes of T and V must match')
  if type(T) is Matrix:
    if   T.tag == sTag: lib.ElCholeskyMod_s(uplo,T.obj,alpha,V.obj)
    elif T.tag == dTag: lib.ElCholeskyMod_d(uplo,T.obj,alpha,V.obj)
    elif T.tag == cTag: lib.ElCholeskyMod_c(uplo,T.obj,alpha,V.obj)
    elif T.tag == zTag: lib.ElCholeskyMod_z(uplo,T.obj,alpha,V.obj)
    else: raise Exception('Unsupported datatype')
  elif type(T) is DistMatrix:
    if   T.tag == sTag: lib.ElCholeskyModDist_s(uplo,T.obj,alpha,V.obj)
    elif T.tag == dTag: lib.ElCholeskyModDist_d(uplo,T.obj,alpha,V.obj)
    elif T.tag == cTag: lib.ElCholeskyModDist_c(uplo,T.obj,alpha,V.obj)
    elif T.tag == zTag: lib.ElCholeskyModDist_z(uplo,T.obj,alpha,V.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHPSDCholesky_s.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholesky_s.restype = c_uint
lib.ElHPSDCholesky_d.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholesky_d.restype = c_uint
lib.ElHPSDCholesky_c.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholesky_c.restype = c_uint
lib.ElHPSDCholesky_z.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholesky_z.restype = c_uint
lib.ElHPSDCholeskyDist_s.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholeskyDist_s.restype = c_uint
lib.ElHPSDCholeskyDist_d.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholeskyDist_d.restype = c_uint
lib.ElHPSDCholeskyDist_c.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholeskyDist_c.restype = c_uint
lib.ElHPSDCholeskyDist_z.argtypes = [c_uint,c_void_p]
lib.ElHPSDCholeskyDist_z.restype = c_uint
def HPSDCholesky(uplo,A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPSDCholesky_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHPSDCholesky_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHPSDCholesky_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHPSDCholesky_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPSDCholeskyDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHPSDCholeskyDist_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHPSDCholeskyDist_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHPSDCholeskyDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElSolveAfterCholesky_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_s.restype = c_uint
lib.ElSolveAfterCholesky_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_d.restype = c_uint
lib.ElSolveAfterCholesky_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_c.restype = c_uint
lib.ElSolveAfterCholesky_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholesky_z.restype = c_uint
lib.ElSolveAfterCholeskyDist_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_s.restype = c_uint
lib.ElSolveAfterCholeskyDist_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_d.restype = c_uint
lib.ElSolveAfterCholeskyDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_c.restype = c_uint
lib.ElSolveAfterCholeskyDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyDist_z.restype = c_uint
lib.ElSolveAfterCholeskyPiv_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPiv_s.restype = c_uint
lib.ElSolveAfterCholeskyPiv_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPiv_d.restype = c_uint
lib.ElSolveAfterCholeskyPiv_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPiv_c.restype = c_uint
lib.ElSolveAfterCholeskyPiv_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPiv_z.restype = c_uint
lib.ElSolveAfterCholeskyPivDist_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPivDist_s.restype = c_uint
lib.ElSolveAfterCholeskyPivDist_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPivDist_d.restype = c_uint
lib.ElSolveAfterCholeskyPivDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPivDist_c.restype = c_uint
lib.ElSolveAfterCholeskyPivDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterCholeskyPivDist_z.restype = c_uint

def SolveAfterCholesky(uplo,orient,A,B):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterCholesky_s(uplo,orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterCholesky_d(uplo,orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterCholesky_c(uplo,orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElSolveAfterCholesky_z(uplo,orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterCholeskyDist_s(uplo,orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterCholeskyDist_d(uplo,orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterCholeskyDist_c(uplo,orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElSolveAfterCholeskyDist_z(uplo,orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

def SolveAfterCholeskyPiv(uplo,orient,A,p,B):
  if type(A) is not type(p) or type(p) is not type(B):
    raise Exception('Types of {A,p,B} must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == sTag: 
      lib.ElSolveAfterCholeskyPiv_s(uplo,orient,A.obj,B.obj,p.obj)
    elif A.tag == dTag: 
      lib.ElSolveAfterCholeskyPiv_d(uplo,orient,A.obj,B.obj,p.obj)
    elif A.tag == cTag: 
      lib.ElSolveAfterCholeskyPiv_c(uplo,orient,A.obj,B.obj,p.obj)
    elif A.tag == zTag: 
      lib.ElSolveAfterCholeskyPiv_z(uplo,orient,A.obj,B.obj,p.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElSolveAfterCholeskyPivDist_s(uplo,orient,A.obj,B.obj,p.obj)
    elif A.tag == dTag: 
      lib.ElSolveAfterCholeskyPivDist_d(uplo,orient,A.obj,B.obj,p.obj)
    elif A.tag == cTag: 
      lib.ElSolveAfterCholeskyPivDist_c(uplo,orient,A.obj,B.obj,p.obj)
    elif A.tag == zTag: 
      lib.ElSolveAfterCholeskyPivDist_z(uplo,orient,A.obj,B.obj,p.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# LDL
# ===

# Emulate an enum for LDL pivot types
(BUNCH_KAUFMAN_A,BUNCH_KAUFMAN_C,BUNCH_KAUFMAN_D,BUNCH_KAUFMAN_BOUNDED,
 BUNCH_PARLETT,LDL_WITHOUT_PIVOTING)=(0,1,2,3,4,5)

class LDLPivot(ctypes.Structure):
  _fields_ = [("nb",iType),("from",(iType*2))]

lib.ElLDL_s.argtypes = [c_void_p]
lib.ElLDL_s.restype = c_uint
lib.ElLDL_d.argtypes = [c_void_p]
lib.ElLDL_d.restype = c_uint
lib.ElLDL_c.argtypes = [c_void_p,bType]
lib.ElLDL_c.restype = c_uint
lib.ElLDL_z.argtypes = [c_void_p,bType]
lib.ElLDL_z.restype = c_uint
lib.ElLDLDist_s.argtypes = [c_void_p]
lib.ElLDLDist_s.restype = c_uint
lib.ElLDLDist_d.argtypes = [c_void_p]
lib.ElLDLDist_d.restype = c_uint
lib.ElLDLDist_c.argtypes = [c_void_p,bType]
lib.ElLDLDist_c.restype = c_uint
lib.ElLDLDist_z.argtypes = [c_void_p,bType]
lib.ElLDLDist_z.restype = c_uint
lib.ElLDLPiv_s.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElLDLPiv_s.restype = c_uint
lib.ElLDLPiv_d.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElLDLPiv_d.restype = c_uint
lib.ElLDLPiv_c.argtypes = [c_void_p,c_void_p,c_void_p,bType,c_uint]
lib.ElLDLPiv_c.restype = c_uint
lib.ElLDLPiv_z.argtypes = [c_void_p,c_void_p,c_void_p,bType,c_uint]
lib.ElLDLPiv_z.restype = c_uint
lib.ElLDLPivDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElLDLPivDist_s.restype = c_uint
lib.ElLDLPivDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_uint]
lib.ElLDLPivDist_d.restype = c_uint
lib.ElLDLPivDist_c.argtypes = [c_void_p,c_void_p,c_void_p,bType,c_uint]
lib.ElLDLPivDist_c.restype = c_uint
lib.ElLDLPivDist_z.argtypes = [c_void_p,c_void_p,c_void_p,bType,c_uint]
lib.ElLDLPivDist_z.restype = c_uint

def LDL(A,conjugate=True,pivType=BUNCH_KAUFMAN_A):
  if type(A) is Matrix:
    if pivType == LDL_WITHOUT_PIVOTING:
      if   A.tag == sTag: lib.ElLDL_s(A.obj)
      elif A.tag == dTag: lib.ElLDL_d(A.obj)
      elif A.tag == cTag: lib.ElLDL_c(A.obj,conjugate)
      elif A.tag == zTag: lib.ElLDL_z(A.obj,conjugate)
      else: raise Exception('Unsupported datatype')
    else:
      dSub = Matrix(A.tag)
      p = Matrix(iTag)
      if   A.tag == sTag: lib.ElLDLPiv_s(A.obj,dSub.obj,p.obj,pivType)
      elif A.tag == dTag: lib.ElLDLPiv_d(A.obj,dSub.obj,p.obj,pivType)
      elif A.tag == cTag: lib.ElLDLPiv_c(A.obj,dSub.obj,p.obj,conjugate,pivType)
      elif A.tag == zTag: lib.ElLDLPiv_z(A.obj,dSub.obj,p.obj,conjugate,pivType)
      else: raise Exception('Unsupported datatype')
      return dSub, p
  elif type(A) is DistMatrix:
    if pivType == LDL_WITHOUT_PIVOTING:
      if   A.tag == sTag: lib.ElLDLDist_s(A.obj)
      elif A.tag == dTag: lib.ElLDLDist_d(A.obj)
      elif A.tag == cTag: lib.ElLDLDist_c(A.obj,conjugate)
      elif A.tag == zTag: lib.ElLDLDist_z(A.obj,conjugate)
      else: raise Exception('Unsupported datatype')
    else:
      dSub = DistMatrix(A.tag,VC,STAR,A.Grid())
      p = DistMatrix(iTag,VC,STAR,A.Grid())
      if   A.tag == sTag: 
        lib.ElLDLPivDist_s(A.obj,dSub.obj,p.obj,pivType)
      elif A.tag == dTag: 
        lib.ElLDLPivDist_d(A.obj,dSub.obj,p.obj,pivType)
      elif A.tag == cTag: 
        lib.ElLDLPivDist_c(A.obj,dSub.obj,p.obj,conjugate,pivType)
      elif A.tag == zTag: 
        lib.ElLDLPivDist_z(A.obj,dSub.obj,p.obj,conjugate,pivType)
      else: raise Exception('Unsupported datatype')
      return dSub, p
  else: raise Exception('Unsupported matrix type')

class InertiaType(ctypes.Structure):
  _fields_ = [("numPositive",iType),("numNegative",iType),("numZero",iType)]

lib.ElInertiaAfterLDL_s.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDL_s.restype = c_uint
lib.ElInertiaAfterLDL_d.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDL_d.restype = c_uint
lib.ElInertiaAfterLDL_c.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDL_c.restype = c_uint
lib.ElInertiaAfterLDL_z.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDL_z.restype = c_uint
lib.ElInertiaAfterLDLDist_s.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDLDist_s.restype = c_uint
lib.ElInertiaAfterLDLDist_d.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDLDist_d.restype = c_uint
lib.ElInertiaAfterLDLDist_c.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDLDist_c.restype = c_uint
lib.ElInertiaAfterLDLDist_z.argtypes = [c_void_p,c_void_p,POINTER(InertiaType)]
lib.ElInertiaAfterLDLDist_z.restype = c_uint
def InertiaAfterLDL(d,dSub):
  inertia = InertiaType()
  if type(d) is not type(dSub):
    raise Exception('Types of d and dSub must match')
  if d.tag != Base(dSub.tag):
    raise Exception('Datatype of d must be the base of that of dSub')
  if type(d) is Matrix:
    if   dSub.tag == sTag: 
      lib.ElInertiaAfterLDL_s(d.obj,dSub.obj,pointer(inertia))
    elif dSub.tag == dTag:
      lib.ElInertiaAfterLDL_d(d.obj,dSub.obj,pointer(inertia))
    elif dSub.tag == cTag:
      lib.ElInertiaAfterLDL_c(d.obj,dSub.obj,pointer(inertia))
    elif dSub.tag == zTag:
      lib.ElInertiaAfterLDL_z(d.obj,dSub.obj,pointer(inertia))
    else: raise Exception('Unsupported datatype')
  elif type(d) is DistMatrix:
    if   dSub.tag == sTag: 
      lib.ElInertiaAfterLDLDist_s(d.obj,dSub.obj,pointer(inertia))
    elif dSub.tag == dTag:
      lib.ElInertiaAfterLDLDist_d(d.obj,dSub.obj,pointer(inertia))
    elif dSub.tag == cTag:
      lib.ElInertiaAfterLDLDist_c(d.obj,dSub.obj,pointer(inertia))
    elif dSub.tag == zTag:
      lib.ElInertiaAfterLDLDist_z(d.obj,dSub.obj,pointer(inertia))
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
  return inertia

lib.ElSolveAfterLDL_s.argtypes = [c_void_p,c_void_p]
lib.ElSolveAfterLDL_s.restype = c_uint
lib.ElSolveAfterLDL_d.argtypes = [c_void_p,c_void_p]
lib.ElSolveAfterLDL_d.restype = c_uint
lib.ElSolveAfterLDL_c.argtypes = [c_void_p,c_void_p,bType]
lib.ElSolveAfterLDL_c.restype = c_uint
lib.ElSolveAfterLDL_z.argtypes = [c_void_p,c_void_p,bType]
lib.ElSolveAfterLDL_z.restype = c_uint
lib.ElSolveAfterLDLDist_s.argtypes = [c_void_p,c_void_p]
lib.ElSolveAfterLDLDist_s.restype = c_uint
lib.ElSolveAfterLDLDist_d.argtypes = [c_void_p,c_void_p]
lib.ElSolveAfterLDLDist_d.restype = c_uint
lib.ElSolveAfterLDLDist_c.argtypes = [c_void_p,c_void_p,bType]
lib.ElSolveAfterLDLDist_c.restype = c_uint
lib.ElSolveAfterLDLDist_z.argtypes = [c_void_p,c_void_p,bType]
lib.ElSolveAfterLDLDist_z.restype = c_uint
def SolveAfterLDL(A,B,conjugate=True):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterLDL_s(A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterLDL_d(A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterLDL_c(A.obj,B.obj,conjugate)
    elif A.tag == zTag: lib.ElSolveAfterLDL_z(A.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterLDLDist_s(A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterLDLDist_d(A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterLDLDist_c(A.obj,B.obj,conjugate)
    elif A.tag == zTag: lib.ElSolveAfterLDLDist_z(A.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

def SolveAfterLDLPiv(A,dSub,p,B,conjugate=True):
  if type(A) is not type(dSub) or type(dSub) is not type(p) or \
     type(p) is not type(B):
    raise Exception('Types of {A,dSub,p,B} must match')
  if A.tag != dSub.tag or dSub.tag != B.tag:
    raise Exception('Datatypes of {A,dSub,B} must match')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterLDLPiv_s(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterLDLPiv_d(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == cTag:
      lib.ElSolveAfterLDLPiv_c(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    elif A.tag == zTag:
      lib.ElSolveAfterLDLPiv_z(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterLDLPivDist_s(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterLDLPivDist_d(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == cTag:
      lib.ElSolveAfterLDLPivDist_c(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    elif A.tag == zTag:
      lib.ElSolveAfterLDLPivDist_z(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported datatype')

lib.ElMultiplyAfterLDL_s.argtypes = [c_void_p,c_void_p]
lib.ElMultiplyAfterLDL_s.restype = c_uint
lib.ElMultiplyAfterLDL_d.argtypes = [c_void_p,c_void_p]
lib.ElMultiplyAfterLDL_d.restype = c_uint
lib.ElMultiplyAfterLDL_c.argtypes = [c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDL_c.restype = c_uint
lib.ElMultiplyAfterLDL_z.argtypes = [c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDL_z.restype = c_uint
lib.ElMultiplyAfterLDLDist_s.argtypes = [c_void_p,c_void_p]
lib.ElMultiplyAfterLDLDist_s.restype = c_uint
lib.ElMultiplyAfterLDLDist_d.argtypes = [c_void_p,c_void_p]
lib.ElMultiplyAfterLDLDist_d.restype = c_uint
lib.ElMultiplyAfterLDLDist_c.argtypes = [c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDLDist_c.restype = c_uint
lib.ElMultiplyAfterLDLDist_z.argtypes = [c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDLDist_z.restype = c_uint
def MultiplyAfterLDL(A,B,conjugate=True):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElMultiplyAfterLDL_s(A.obj,B.obj)
    elif A.tag == dTag: lib.ElMultiplyAfterLDL_d(A.obj,B.obj)
    elif A.tag == cTag: lib.ElMultiplyAfterLDL_c(A.obj,B.obj,conjugate)
    elif A.tag == zTag: lib.ElMultiplyAfterLDL_z(A.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElMultiplyAfterLDLDist_s(A.obj,B.obj)
    elif A.tag == dTag: lib.ElMultiplyAfterLDLDist_d(A.obj,B.obj)
    elif A.tag == cTag: lib.ElMultiplyAfterLDLDist_c(A.obj,B.obj,conjugate)
    elif A.tag == zTag: lib.ElMultiplyAfterLDLDist_z(A.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElMultiplyAfterLDLPiv_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElMultiplyAfterLDLPiv_s.restype = c_uint
lib.ElMultiplyAfterLDLPiv_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElMultiplyAfterLDLPiv_d.restype = c_uint
lib.ElMultiplyAfterLDLPiv_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDLPiv_c.restype = c_uint
lib.ElMultiplyAfterLDLPiv_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDLPiv_z.restype = c_uint
lib.ElMultiplyAfterLDLPivDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElMultiplyAfterLDLPivDist_s.restype = c_uint
lib.ElMultiplyAfterLDLPivDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElMultiplyAfterLDLPivDist_d.restype = c_uint
lib.ElMultiplyAfterLDLPivDist_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDLPivDist_c.restype = c_uint
lib.ElMultiplyAfterLDLPivDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,bType]
lib.ElMultiplyAfterLDLPivDist_z.restype = c_uint
def MultiplyAfterLDLPiv(A,dSub,p,B,conjugate=True):
  if type(A) is not type(dSub) or type(dSub) is not type(p) or \
     type(p) is not type(B):
    raise Exception('Types of {A,dSub,p,B} must match')
  if A.tag != dSub.tag or dSub.tag != B.tag:
    raise Exception('Datatypes of {A,dSub,B} must match')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElMultiplyAfterLDLPiv_s(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == dTag: lib.ElMultiplyAfterLDLPiv_d(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == cTag:
      lib.ElMultiplyAfterLDLPiv_c(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    elif A.tag == zTag:
      lib.ElMultiplyAfterLDLPiv_z(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElMultiplyAfterLDLPivDist_s(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == dTag: 
      lib.ElMultiplyAfterLDLPivDist_d(A.obj,dSub.obj,p.obj,B.obj)
    elif A.tag == cTag:
      lib.ElMultiplyAfterLDLPivDist_c(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    elif A.tag == zTag:
      lib.ElMultiplyAfterLDLPivDist_z(A.obj,dSub.obj,p.obj,B.obj,conjugate)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# LU factorization
# ================

# Emulate an enum for the pivot type for LU factorization
(LU_PARTIAL,LU_FULL,LU_ROOK,LU_WITHOUT_PIVOTING)=(0,1,2,3)

lib.ElLU_s.argtypes = [c_void_p]
lib.ElLU_s.restype = c_uint
lib.ElLU_d.argtypes = [c_void_p]
lib.ElLU_d.restype = c_uint
lib.ElLU_c.argtypes = [c_void_p]
lib.ElLU_c.restype = c_uint
lib.ElLU_z.argtypes = [c_void_p]
lib.ElLU_z.restype = c_uint
lib.ElLUDist_s.argtypes = [c_void_p]
lib.ElLUDist_s.restype = c_uint
lib.ElLUDist_d.argtypes = [c_void_p]
lib.ElLUDist_d.restype = c_uint
lib.ElLUDist_c.argtypes = [c_void_p]
lib.ElLUDist_c.restype = c_uint
lib.ElLUDist_z.argtypes = [c_void_p]
lib.ElLUDist_z.restype = c_uint
lib.ElLUPartialPiv_s.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPiv_s.restype = c_uint
lib.ElLUPartialPiv_d.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPiv_d.restype = c_uint
lib.ElLUPartialPiv_c.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPiv_c.restype = c_uint
lib.ElLUPartialPiv_z.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPiv_z.restype = c_uint
lib.ElLUPartialPivDist_s.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPivDist_s.restype = c_uint
lib.ElLUPartialPivDist_d.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPivDist_d.restype = c_uint
lib.ElLUPartialPivDist_c.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPivDist_c.restype = c_uint
lib.ElLUPartialPivDist_z.argtypes = [c_void_p,c_void_p]
lib.ElLUPartialPivDist_z.restype = c_uint
lib.ElLUFullPiv_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPiv_s.restype = c_uint
lib.ElLUFullPiv_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPiv_d.restype = c_uint
lib.ElLUFullPiv_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPiv_c.restype = c_uint
lib.ElLUFullPiv_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPiv_z.restype = c_uint
lib.ElLUFullPivDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPivDist_s.restype = c_uint
lib.ElLUFullPivDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPivDist_d.restype = c_uint
lib.ElLUFullPivDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPivDist_c.restype = c_uint
lib.ElLUFullPivDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLUFullPivDist_z.restype = c_uint

def LU(A,pivType=LU_PARTIAL):
  if type(A) is Matrix: 
    if pivType == LU_WITHOUT_PIVOTING:
      if   A.tag == sTag: lib.ElLU_s(A.obj)
      elif A.tag == dTag: lib.ElLU_d(A.obj)
      elif A.tag == cTag: lib.ElLU_c(A.obj)
      elif A.tag == zTag: lib.ElLU_z(A.obj)
      else: raise Exception('Unsupported datatype')
    elif pivType == LU_PARTIAL:
      p = Matrix(iTag)
      if   A.tag == sTag: lib.ElLUPartialPiv_s(A.obj,p.obj)
      elif A.tag == dTag: lib.ElLUPartialPiv_d(A.obj,p.obj)
      elif A.tag == cTag: lib.ElLUPartialPiv_c(A.obj,p.obj)
      elif A.tag == zTag: lib.ElLUPartialPiv_z(A.obj,p.obj)
      else: raise Exception('Unsupported datatype')
      return p
    elif pivType == LU_FULL:
      p = Matrix(iTag)
      q = Matrix(iTag)
      if   A.tag == sTag: lib.ElLUFullPiv_s(A.obj,p.obj,q.obj)
      elif A.tag == dTag: lib.ElLUFullPiv_d(A.obj,p.obj,q.obj)
      elif A.tag == cTag: lib.ElLUFullPiv_c(A.obj,p.obj,q.obj)
      elif A.tag == zTag: lib.ElLUFullPiv_z(A.obj,p.obj,q.obj)
      else: raise Exception('Unsupported datatype')
      return p, q
    else: raise Exception('Unsupported pivot type')
  elif type(A) is DistMatrix:
    if pivType == LU_WITHOUT_PIVOTING:
      if   A.tag == sTag: lib.ElLUDist_s(A.obj)
      elif A.tag == dTag: lib.ElLUDist_d(A.obj)
      elif A.tag == cTag: lib.ElLUDist_c(A.obj)
      elif A.tag == zTag: lib.ElLUDist_z(A.obj)
      else: raise Exception('Unsupported datatype')
    elif pivType == LU_PARTIAL:
      p = DistMatrix(iTag,VC,STAR,A.Grid())
      if   A.tag == sTag: lib.ElLUPartialPivDist_s(A.obj,p.obj)
      elif A.tag == dTag: lib.ElLUPartialPivDist_d(A.obj,p.obj)
      elif A.tag == cTag: lib.ElLUPartialPivDist_c(A.obj,p.obj)
      elif A.tag == zTag: lib.ElLUPartialPivDist_z(A.obj,p.obj)
      else: raise Exception('Unsupported datatype')
      return p
    elif pivType == LU_FULL:
      p = DistMatrix(iTag,VC,STAR,A.Grid())
      q = DistMatrix(iTag,VC,STAR,A.Grid())
      if   A.tag == sTag: lib.ElLUFullPivDist_s(A.obj,p.obj,q.obj)
      elif A.tag == dTag: lib.ElLUFullPivDist_d(A.obj,p.obj,q.obj)
      elif A.tag == cTag: lib.ElLUFullPivDist_c(A.obj,p.obj,q.obj)
      elif A.tag == zTag: lib.ElLUFullPivDist_z(A.obj,p.obj,q.obj)
      else: raise Exception('Unsupported datatype')
      return p, q
    else: raise Exception('Unsupported pivot type')
  else: raise Exception('Unsupported matrix type')

lib.ElLUMod_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,sType]
lib.ElLUMod_s.restype = c_uint
lib.ElLUMod_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,dType]
lib.ElLUMod_d.restype = c_uint
lib.ElLUMod_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,sType]
lib.ElLUMod_c.restype = c_uint
lib.ElLUMod_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,dType]
lib.ElLUMod_z.restype = c_uint
lib.ElLUModDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,sType]
lib.ElLUModDist_s.restype = c_uint
lib.ElLUModDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,dType]
lib.ElLUModDist_d.restype = c_uint
lib.ElLUModDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,sType]
lib.ElLUModDist_c.restype = c_uint
lib.ElLUModDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,dType]
lib.ElLUModDist_z.restype = c_uint
def LUMod(A,p,u,v,conjugate=True,tau=0.1):
  if type(A) is not type(p) or type(p) is not type(u) or type(u) is not type(v):
    raise Exception('Types of {A,p,u,v} must be equal')
  if A.tag != u.tag or u.tag != v.tag:
    raise Exception('Datatypes of {A,u,v} must be equal')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLUMod_s(A.obj,p.obj,u.obj,v.obj,tau)
    elif A.tag == dTag: lib.ElLUMod_d(A.obj,p.obj,u.obj,v.obj,tau)
    elif A.tag == cTag: lib.ElLUMod_c(A.obj,p.obj,u.obj,v.obj,conjugate,tau)
    elif A.tag == zTag: lib.ElLUMod_z(A.obj,p.obj,u.obj,v.obj,conjugate,tau)
    else: raise Exception('Unsupported datatype') 
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLUModDist_s(A.obj,p.obj,u.obj,v.obj,tau)
    elif A.tag == dTag: lib.ElLUModDist_d(A.obj,p.obj,u.obj,v.obj,tau)
    elif A.tag == cTag: lib.ElLUModDist_c(A.obj,p.obj,u.obj,v.obj,conjugate,tau)
    elif A.tag == zTag: lib.ElLUModDist_z(A.obj,p.obj,u.obj,v.obj,conjugate,tau)
    else: raise Exception('Unsupported datatype') 
  else: raise Exception('Unsupported matrix type')

lib.ElSolveAfterLU_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLU_s.restype = c_uint
lib.ElSolveAfterLU_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLU_d.restype = c_uint
lib.ElSolveAfterLU_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLU_c.restype = c_uint
lib.ElSolveAfterLU_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLU_z.restype = c_uint
lib.ElSolveAfterLUDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLUDist_s.restype = c_uint
lib.ElSolveAfterLUDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLUDist_d.restype = c_uint
lib.ElSolveAfterLUDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLUDist_c.restype = c_uint
lib.ElSolveAfterLUDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSolveAfterLUDist_z.restype = c_uint
def SolveAfterLU(orient,A,B):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterLU_s(orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterLU_d(orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterLU_c(orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElSolveAfterLU_z(orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSolveAfterLUDist_s(orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterLUDist_d(orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterLUDist_c(orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElSolveAfterLUDist_z(orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElSolveAfterLUPartialPiv_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUPartialPiv_s.restype = c_uint
lib.ElSolveAfterLUPartialPiv_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUPartialPiv_d.restype = c_uint
lib.ElSolveAfterLUPartialPiv_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUPartialPiv_c.restype = c_uint
lib.ElSolveAfterLUPartialPiv_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUPartialPiv_z.restype = c_uint
lib.ElSolveAfterLUPartialPivDist_s.restype = c_uint
lib.ElSolveAfterLUPartialPivDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUPartialPivDist_d.restype = c_uint
lib.ElSolveAfterLUPartialPivDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUPartialPivDist_c.restype = c_uint
lib.ElSolveAfterLUPartialPivDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUPartialPivDist_z.restype = c_uint
def SolveAfterLUPartialPiv(orient,A,p,B):
  if type(A) is not type(p) or type(p) is not type(B):
    raise Exception('Types of {A,p,B} must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSolveAfterLUPartialPiv_s(orient,A.obj,p.obj,B.obj)
    elif A.tag == dTag: lib.ElSolveAfterLUPartialPiv_d(orient,A.obj,p.obj,B.obj)
    elif A.tag == cTag: lib.ElSolveAfterLUPartialPiv_c(orient,A.obj,p.obj,B.obj)
    elif A.tag == zTag: lib.ElSolveAfterLUPartialPiv_z(orient,A.obj,p.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElSolveAfterLUPartialPivDist_s(orient,A.obj,p.obj,B.obj)
    elif A.tag == dTag: 
      lib.ElSolveAfterLUPartialPivDist_d(orient,A.obj,p.obj,B.obj)
    elif A.tag == cTag: 
      lib.ElSolveAfterLUPartialPivDist_c(orient,A.obj,p.obj,B.obj)
    elif A.tag == zTag: 
      lib.ElSolveAfterLUPartialPivDist_z(orient,A.obj,p.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElSolveAfterLUFullPiv_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUFullPiv_s.restype = c_uint
lib.ElSolveAfterLUFullPiv_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUFullPiv_d.restype = c_uint
lib.ElSolveAfterLUFullPiv_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUFullPiv_c.restype = c_uint
lib.ElSolveAfterLUFullPiv_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUFullPiv_z.restype = c_uint
lib.ElSolveAfterLUFullPivDist_s.restype = c_uint
lib.ElSolveAfterLUFullPivDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUFullPivDist_d.restype = c_uint
lib.ElSolveAfterLUFullPivDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUFullPivDist_c.restype = c_uint
lib.ElSolveAfterLUFullPivDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLUFullPivDist_z.restype = c_uint
def SolveAfterLUFullPiv(orient,A,p,q,B):
  if type(A) is not type(p) or type(p) is not type(q) or type(q) is not type(B):
    raise Exception('Types of {A,p,q,B} must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if p.tag != iTag or q.tag != iTag:
    raise Exception('p and q must be integral')
  if type(A) is Matrix:
    if   A.tag == sTag:
      lib.ElSolveAfterLUFullPiv_s(orient,A.obj,p.obj,q.obj,B.obj)
    elif A.tag == dTag:
      lib.ElSolveAfterLUFullPiv_d(orient,A.obj,p.obj,q.obj,B.obj)
    elif A.tag == cTag:
      lib.ElSolveAfterLUFullPiv_c(orient,A.obj,p.obj,q.obj,B.obj)
    elif A.tag == zTag:
      lib.ElSolveAfterLUFullPiv_z(orient,A.obj,p.obj,q.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElSolveAfterLUFullPivDist_s(orient,A.obj,p.obj,q.obj,B.obj)
    elif A.tag == dTag: 
      lib.ElSolveAfterLUFullPivDist_d(orient,A.obj,p.obj,q.obj,B.obj)
    elif A.tag == cTag: 
      lib.ElSolveAfterLUFullPivDist_c(orient,A.obj,p.obj,q.obj,B.obj)
    elif A.tag == zTag: 
      lib.ElSolveAfterLUFullPivDist_z(orient,A.obj,p.obj,q.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# LQ factorization
# ================
# TODO

# QR factorization
# ================
lib.ElQRCtrlFillDefault_s.argtypes = [c_void_p]
lib.ElQRCtrlFillDefault_s.restype = c_uint
lib.ElQRCtrlFillDefault_d.argtypes = [c_void_p]
lib.ElQRCtrlFillDefault_d.restype = c_uint
class QRCtrl_s(ctypes.Structure):
  _fields_ = [("colPiv",bType),("boundRank",bType),("maxRank",iType),
              ("adaptive",bType),("tol",sType),("alwaysRecomputeNorms",bType)]
  def __init__(self):
    lib.ElQRCtrlFillDefault_s(pointer(self))
class QRCtrl_d(ctypes.Structure):
  _fields_ = [("colPiv",bType),("boundRank",bType),("maxRank",iType),
              ("adaptive",bType),("tol",dType),("alwaysRecomputeNorms",bType)]
  def __init__(self):
    lib.ElQRCtrlFillDefault_d(pointer(self))
# TODO

# RQ factorization
# ================
# TODO
