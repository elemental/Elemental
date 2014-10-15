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
(LQ_IMPLICIT,LQ_EXPLICIT,LQ_EXPLICIT_TRIANG,LQ_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElLQ_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQ_s.restype = c_uint
lib.ElLQ_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQ_d.restype = c_uint
lib.ElLQ_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQ_c.restype = c_uint
lib.ElLQ_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQ_z.restype = c_uint
lib.ElLQDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQDist_s.restype = c_uint
lib.ElLQDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQDist_d.restype = c_uint
lib.ElLQDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQDist_c.restype = c_uint
lib.ElLQDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLQDist_z.restype = c_uint
lib.ElLQExplicit_s.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicit_s.restype = c_uint
lib.ElLQExplicit_d.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicit_d.restype = c_uint
lib.ElLQExplicit_c.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicit_c.restype = c_uint
lib.ElLQExplicit_z.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicit_z.restype = c_uint
lib.ElLQExplicitDist_s.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicitDist_s.restype = c_uint
lib.ElLQExplicitDist_d.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicitDist_d.restype = c_uint
lib.ElLQExplicitDist_c.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicitDist_c.restype = c_uint
lib.ElLQExplicitDist_z.argtypes = [c_void_p,c_void_p]
lib.ElLQExplicitDist_z.restype = c_uint
lib.ElLQExplicitUnitary_s.argtypes = [c_void_p]
lib.ElLQExplicitUnitary_s.restype = c_uint
lib.ElLQExplicitUnitary_d.argtypes = [c_void_p]
lib.ElLQExplicitUnitary_d.restype = c_uint
lib.ElLQExplicitUnitary_c.argtypes = [c_void_p]
lib.ElLQExplicitUnitary_c.restype = c_uint
lib.ElLQExplicitUnitary_z.argtypes = [c_void_p]
lib.ElLQExplicitUnitary_z.restype = c_uint
lib.ElLQExplicitUnitaryDist_s.argtypes = [c_void_p]
lib.ElLQExplicitUnitaryDist_s.restype = c_uint
lib.ElLQExplicitUnitaryDist_d.argtypes = [c_void_p]
lib.ElLQExplicitUnitaryDist_d.restype = c_uint
lib.ElLQExplicitUnitaryDist_c.argtypes = [c_void_p]
lib.ElLQExplicitUnitaryDist_c.restype = c_uint
lib.ElLQExplicitUnitaryDist_z.argtypes = [c_void_p]
lib.ElLQExplicitUnitaryDist_z.restype = c_uint
lib.ElLQExplicitTriang_s.argtypes = [c_void_p]
lib.ElLQExplicitTriang_s.restype = c_uint
lib.ElLQExplicitTriang_d.argtypes = [c_void_p]
lib.ElLQExplicitTriang_d.restype = c_uint
lib.ElLQExplicitTriang_c.argtypes = [c_void_p]
lib.ElLQExplicitTriang_c.restype = c_uint
lib.ElLQExplicitTriang_z.argtypes = [c_void_p]
lib.ElLQExplicitTriang_z.restype = c_uint
lib.ElLQExplicitTriangDist_s.argtypes = [c_void_p]
lib.ElLQExplicitTriangDist_s.restype = c_uint
lib.ElLQExplicitTriangDist_d.argtypes = [c_void_p]
lib.ElLQExplicitTriangDist_d.restype = c_uint
lib.ElLQExplicitTriangDist_c.argtypes = [c_void_p]
lib.ElLQExplicitTriangDist_c.restype = c_uint
lib.ElLQExplicitTriangDist_z.argtypes = [c_void_p]
lib.ElLQExplicitTriangDist_z.restype = c_uint
def LQ(A,factType=LQ_IMPLICIT):
  if type(A) is Matrix:
    if factType == LQ_IMPLICIT:
      t = Matrix(A.tag)
      d = Matrix(Base(A.tag))
      if   A.tag == sTag: lib.ElLQ_s(A.obj,t.obj,d.obj)
      elif A.tag == dTag: lib.ElLQ_d(A.obj,t.obj,d.obj)
      elif A.tag == cTag: lib.ElLQ_c(A.obj,t.obj,d.obj)
      elif A.tag == zTag: lib.ElLQ_z(A.obj,t.obj,d.obj)
      else: raise Exception('Unsupported datatype')
      return t, d
    elif factType == LQ_EXPLICIT:
      L = Matrix(A.tag)
      if   A.tag == sTag: lib.ElLQExplicit_s(L.obj,A.obj)
      elif A.tag == dTag: lib.ElLQExplicit_d(L.obj,A.obj)
      elif A.tag == cTag: lib.ElLQExplicit_c(L.obj,A.obj)
      elif A.tag == zTag: lib.ElLQExplicit_z(L.obj,A.obj)
      else: raise Exception('Unsupported datatype')
      return L
    elif factType == LQ_EXPLICIT_UNITARY:
      if   A.tag == sTag: lib.ElLQExplicitUnitary_s(A.obj)
      elif A.tag == dTag: lib.ElLQExplicitUnitary_d(A.obj)
      elif A.tag == cTag: lib.ElLQExplicitUnitary_c(A.obj)
      elif A.tag == zTag: lib.ElLQExplicitUnitary_z(A.obj)
      else: raise Exception('Unsupported datatype')
    elif factType == LQ_EXPLICIT_TRIANG:
      if   A.tag == sTag: lib.ElLQExplicitTriang_s(A.obj)
      elif A.tag == dTag: lib.ElLQExplicitTriang_d(A.obj)
      elif A.tag == cTag: lib.ElLQExplicitTriang_c(A.obj)
      elif A.tag == zTag: lib.ElLQExplicitTriang_z(A.obj)
      else: raise Exception('Unsupported datatype')
    else: raise Exception('Unsupported LQ factorization type')
  elif type(A) is DistMatrix:
    if factType == LQ_IMPLICIT:
      t = DistMatrix(A.tag,MC,STAR,A.Grid())
      d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
      if   A.tag == sTag: lib.ElLQDist_s(A.obj,t.obj,d.obj)
      elif A.tag == dTag: lib.ElLQDist_d(A.obj,t.obj,d.obj)
      elif A.tag == cTag: lib.ElLQDist_c(A.obj,t.obj,d.obj)
      elif A.tag == zTag: lib.ElLQDist_z(A.obj,t.obj,d.obj)
      else: raise Exception('Unsupported datatype')
      return t, d
    elif factType == LQ_EXPLICIT:
      L = DistMatrix(A.tag,MC,MR,A.Grid())
      if   A.tag == sTag: lib.ElLQExplicitDist_s(L.obj,A.obj)
      elif A.tag == dTag: lib.ElLQExplicitDist_d(L.obj,A.obj)
      elif A.tag == cTag: lib.ElLQExplicitDist_c(L.obj,A.obj)
      elif A.tag == zTag: lib.ElLQExplicitDist_z(L.obj,A.obj)
      else: raise Exception('Unsupported datatype')
      return L
    elif factType == LQ_EXPLICIT_UNITARY:
      if   A.tag == sTag: lib.ElLQExplicitUnitaryDist_s(A.obj)
      elif A.tag == dTag: lib.ElLQExplicitUnitaryDist_d(A.obj)
      elif A.tag == cTag: lib.ElLQExplicitUnitaryDist_c(A.obj)
      elif A.tag == zTag: lib.ElLQExplicitUnitaryDist_z(A.obj)
      else: raise Exception('Unsupported datatype')
    elif factType == LQ_EXPLICIT_TRIANG:
      if   A.tag == sTag: lib.ElLQExplicitTriangDist_s(A.obj)
      elif A.tag == dTag: lib.ElLQExplicitTriangDist_d(A.obj)
      elif A.tag == cTag: lib.ElLQExplicitTriangDist_c(A.obj)
      elif A.tag == zTag: lib.ElLQExplicitTriangDist_z(A.obj)
      else: raise Exception('Unsupported datatype')
    else: raise Exception('Unsupported LQ factorization type')
  else: raise Exception('Unsupported datatype')

lib.ElApplyQAfterLQ_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQ_s.restype = c_uint
lib.ElApplyQAfterLQ_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQ_d.restype = c_uint
lib.ElApplyQAfterLQ_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQ_c.restype = c_uint
lib.ElApplyQAfterLQ_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQ_z.restype = c_uint
lib.ElApplyQAfterLQDist_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQDist_s.restype = c_uint
lib.ElApplyQAfterLQDist_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQDist_d.restype = c_uint
lib.ElApplyQAfterLQDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQDist_c.restype = c_uint
lib.ElApplyQAfterLQDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterLQDist_z.restype = c_uint
def ApplyQAfterLQ(side,orient,A,t,d,B):
  if type(A) is not type(t) or type(t) is not type(d) or type(d) is not type(B):
    raise Exception('Matrix types of {A,t,d,B} must match')
  if A.tag != t.tag or t.tag != B.tag or d.tag != Base(A.tag):
    raise Exception('Datatypes of {A,t,B} must match and d must have base type')
  if type(A) is Matrix:
    if   A.tag == sTag: 
      lib.ElApplyQAfterLQ_s(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == dTag:
      lib.ElApplyQAfterLQ_d(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == cTag:
      lib.ElApplyQAfterLQ_c(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == zTag:
      lib.ElApplyQAfterLQ_z(side,orient,A.obj,t.obj,d.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElApplyQAfterLQDist_s(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == dTag:
      lib.ElApplyQAfterLQDist_d(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == cTag:
      lib.ElApplyQAfterLQDist_c(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == zTag:
      lib.ElApplyQAfterLQDist_z(side,orient,A.obj,t.obj,d.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElSolveAfterLQ_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQ_s.restype = c_uint
lib.ElSolveAfterLQ_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQ_d.restype = c_uint
lib.ElSolveAfterLQ_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQ_c.restype = c_uint
lib.ElSolveAfterLQ_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQ_z.restype = c_uint
lib.ElSolveAfterLQDist_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQDist_s.restype = c_uint
lib.ElSolveAfterLQDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQDist_d.restype = c_uint
lib.ElSolveAfterLQDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQDist_c.restype = c_uint
lib.ElSolveAfterLQDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterLQDist_z.restype = c_uint
def SolveAfterLQ(orient,A,t,d,B):
  if type(A) is Matrix:
    X = Matrix(A.tag)
    if   A.tag == sTag: 
      lib.ElSolveAfterLQ_s(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == dTag:
      lib.ElSolveAfterLQ_d(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == cTag:
      lib.ElSolveAfterLQ_c(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == zTag:
      lib.ElSolveAfterLQ_z(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    if   A.tag == sTag: 
      lib.ElSolveAfterLQDist_s(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == dTag:
      lib.ElSolveAfterLQDist_d(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == cTag:
      lib.ElSolveAfterLQDist_c(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == zTag:
      lib.ElSolveAfterLQDist_z(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  else: raise Exception('Unsupported matrix type')

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

(QR_IMPLICIT,QR_EXPLICIT,QR_EXPLICIT_TRIANG,QR_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElQR_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQR_s.restype = c_uint
lib.ElQR_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQR_d.restype = c_uint
lib.ElQR_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQR_c.restype = c_uint
lib.ElQR_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQR_z.restype = c_uint
lib.ElQRDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRDist_s.restype = c_uint
lib.ElQRDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRDist_d.restype = c_uint
lib.ElQRDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRDist_c.restype = c_uint
lib.ElQRDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRDist_z.restype = c_uint

lib.ElQRExplicit_s.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicit_s.restype = c_uint
lib.ElQRExplicit_d.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicit_d.restype = c_uint
lib.ElQRExplicit_c.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicit_c.restype = c_uint
lib.ElQRExplicit_z.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicit_z.restype = c_uint
lib.ElQRExplicitDist_s.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicitDist_s.restype = c_uint
lib.ElQRExplicitDist_d.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicitDist_d.restype = c_uint
lib.ElQRExplicitDist_c.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicitDist_c.restype = c_uint
lib.ElQRExplicitDist_z.argtypes = [c_void_p,c_void_p]
lib.ElQRExplicitDist_z.restype = c_uint

lib.ElQRExplicitTriang_s.argtypes = [c_void_p]
lib.ElQRExplicitTriang_s.restype = c_uint
lib.ElQRExplicitTriang_d.argtypes = [c_void_p]
lib.ElQRExplicitTriang_d.restype = c_uint
lib.ElQRExplicitTriang_c.argtypes = [c_void_p]
lib.ElQRExplicitTriang_c.restype = c_uint
lib.ElQRExplicitTriang_z.argtypes = [c_void_p]
lib.ElQRExplicitTriang_z.restype = c_uint
lib.ElQRExplicitTriangDist_s.argtypes = [c_void_p]
lib.ElQRExplicitTriangDist_s.restype = c_uint
lib.ElQRExplicitTriangDist_d.argtypes = [c_void_p]
lib.ElQRExplicitTriangDist_d.restype = c_uint
lib.ElQRExplicitTriangDist_c.argtypes = [c_void_p]
lib.ElQRExplicitTriangDist_c.restype = c_uint
lib.ElQRExplicitTriangDist_z.argtypes = [c_void_p]
lib.ElQRExplicitTriangDist_z.restype = c_uint

lib.ElQRExplicitUnitary_s.argtypes = [c_void_p]
lib.ElQRExplicitUnitary_s.restype = c_uint
lib.ElQRExplicitUnitary_d.argtypes = [c_void_p]
lib.ElQRExplicitUnitary_d.restype = c_uint
lib.ElQRExplicitUnitary_c.argtypes = [c_void_p]
lib.ElQRExplicitUnitary_c.restype = c_uint
lib.ElQRExplicitUnitary_z.argtypes = [c_void_p]
lib.ElQRExplicitUnitary_z.restype = c_uint
lib.ElQRExplicitUnitaryDist_s.argtypes = [c_void_p]
lib.ElQRExplicitUnitaryDist_s.restype = c_uint
lib.ElQRExplicitUnitaryDist_d.argtypes = [c_void_p]
lib.ElQRExplicitUnitaryDist_d.restype = c_uint
lib.ElQRExplicitUnitaryDist_c.argtypes = [c_void_p]
lib.ElQRExplicitUnitaryDist_c.restype = c_uint
lib.ElQRExplicitUnitaryDist_z.argtypes = [c_void_p]
lib.ElQRExplicitUnitaryDist_z.restype = c_uint

lib.ElQRColPiv_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPiv_s.restype = c_uint
lib.ElQRColPiv_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPiv_d.restype = c_uint
lib.ElQRColPiv_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPiv_c.restype = c_uint
lib.ElQRColPiv_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPiv_z.restype = c_uint
lib.ElQRColPivDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPivDist_s.restype = c_uint
lib.ElQRColPivDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPivDist_d.restype = c_uint
lib.ElQRColPivDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPivDist_c.restype = c_uint
lib.ElQRColPivDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElQRColPivDist_z.restype = c_uint

lib.ElQRColPivX_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElQRColPivX_s.restype = c_uint
lib.ElQRColPivX_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElQRColPivX_d.restype = c_uint
lib.ElQRColPivX_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElQRColPivX_c.restype = c_uint
lib.ElQRColPivX_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElQRColPivX_z.restype = c_uint
lib.ElQRColPivXDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElQRColPivXDist_s.restype = c_uint
lib.ElQRColPivXDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElQRColPivXDist_d.restype = c_uint
lib.ElQRColPivXDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElQRColPivXDist_c.restype = c_uint
lib.ElQRColPivXDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElQRColPivXDist_z.restype = c_uint

lib.ElQRColPivExplicit_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicit_s.restype = c_uint
lib.ElQRColPivExplicit_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicit_d.restype = c_uint
lib.ElQRColPivExplicit_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicit_c.restype = c_uint
lib.ElQRColPivExplicit_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicit_z.restype = c_uint
lib.ElQRColPivExplicitDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicitDist_s.restype = c_uint
lib.ElQRColPivExplicitDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicitDist_d.restype = c_uint
lib.ElQRColPivExplicitDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicitDist_c.restype = c_uint
lib.ElQRColPivExplicitDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElQRColPivExplicitDist_z.restype = c_uint

def QR(A,piv=False,factType=QR_IMPLICIT,ctrl=None):
  if type(A) is Matrix:
    if piv:
      if factType == QR_IMPLICIT:  
        t = Matrix(A.tag)
        d = Matrix(Base(A.tag))
        p = Matrix(iTag)
        if   A.tag == sTag:
          if ctrl == None:
            lib.ElQRColPivDist_s(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_s(A.obj,t.obj,d.obj,p.obj,ctrl)
        elif A.tag == dTag:
          if ctrl == None:
            lib.ElQRColPivDist_d(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_d(A.obj,t.obj,d.obj,p.obj,ctrl)
        elif A.tag == cTag:
          if ctrl == None:
            lib.ElQRColPivDist_c(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_c(A.obj,t.obj,d.obj,p.obj,ctrl)
        elif A.tag == zTag:
          if ctrl == None:
            lib.ElQRColPivDist_z(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_z(A.obj,t.obj,d.obj,p.obj,ctrl)
        else: raise Exception('Unsupported datatype')
        return t, d, p
      elif factType == QR_EXPLICIT:
        if ctrl != None: 
          raise Exception('\'ctrl\' not yet supported for explicit piv fact\'s')
        R = Matrix(A.tag)
        P = Matrix(iTag)
        if   A.tag == sTag: lib.ElQRColPivExplicit_s(A.obj,R.obj,P.obj)
        elif A.tag == dTag: lib.ElQRColPivExplicit_d(A.obj,R.obj,P.obj)
        elif A.tag == cTag: lib.ElQRColPivExplicit_c(A.obj,R.obj,P.obj)
        elif A.tag == zTag: lib.ElQRColPivExplicit_z(A.obj,R.obj,P.obj)
        else: raise Exception('Unsupported datatype')
        return R, P
      else: 
        raise Exception('Partial pivoted explicit fact\'s not yet supported')
    else:
      if factType == QR_IMPLICIT:
        t = Matrix(A.tag)
        d = Matrix(Base(A.tag))
        if   A.tag == sTag: lib.ElQR_s(A.obj,t.obj,d.obj)
        elif A.tag == dTag: lib.ElQR_d(A.obj,t.obj,d.obj)
        elif A.tag == cTag: lib.ElQR_c(A.obj,t.obj,d.obj)
        elif A.tag == zTag: lib.ElQR_z(A.obj,t.obj,d.obj)
        else: raise Exception('Unsupported datatype')
        return t, d
      elif factType == QR_EXPLICIT:
        R = Matrix(A.tag)
        if   A.tag == sTag: lib.ElQRExplicit_s(A.obj,R.obj)
        elif A.tag == dTag: lib.ElQRExplicit_d(A.obj,R.obj)
        elif A.tag == cTag: lib.ElQRExplicit_c(A.obj,R.obj)
        elif A.tag == zTag: lib.ElQRExplicit_z(A.obj,R.obj)
        else: raise Exception('Unsupported datatype')
        return R
      elif factType == QR_EXPLICIT_TRIANG:
        if   A.tag == sTag: lib.ElQRExplicitTriang_s(A.obj)
        elif A.tag == dTag: lib.ElQRExplicitTriang_d(A.obj)
        elif A.tag == cTag: lib.ElQRExplicitTriang_c(A.obj)
        elif A.tag == zTag: lib.ElQRExplicitTriang_z(A.obj)
        else: raise Exception('Unsupported datatype')
      elif factType == QR_EXPLICIT_UNITARY:
        if   A.tag == sTag: lib.ElQRExplicitUnitary_s(A.obj)
        elif A.tag == dTag: lib.ElQRExplicitUnitary_d(A.obj)
        elif A.tag == cTag: lib.ElQRExplicitUnitary_c(A.obj)
        elif A.tag == zTag: lib.ElQRExplicitUnitary_z(A.obj)
        else: raise Exception('Unsupported datatype')
      else: raise Exception('Unsupported QR factorization type')
  elif type(A) is DistMatrix:
    if piv:
      if factType == QR_IMPLICIT:  
        t = DistMatrix(A.tag,MC,STAR,A.Grid())
        d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
        p = DistMatrix(iTag,MC,STAR,A.Grid())
        if   A.tag == sTag:
          if ctrl == None:
            lib.ElQRColPivDist_s(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_s(A.obj,t.obj,d.obj,p.obj,ctrl)
        elif A.tag == dTag:
          if ctrl == None:
            lib.ElQRColPivDist_d(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_d(A.obj,t.obj,d.obj,p.obj,ctrl)
        elif A.tag == cTag:
          if ctrl == None:
            lib.ElQRColPivDist_c(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_c(A.obj,t.obj,d.obj,p.obj,ctrl)
        elif A.tag == zTag:
          if ctrl == None:
            lib.ElQRColPivDist_z(A.obj,t.obj,d.obj,p.obj)
          else:
            lib.ElQRColPivXDist_z(A.obj,t.obj,d.obj,p.obj,ctrl)
        else: raise Exception('Unsupported datatype')
        return t, d, p
      elif factType == QR_EXPLICIT:
        if ctrl != None: 
          raise Exception('\'ctrl\' not yet supported for explicit piv fact\'s')
        R = DistMatrix(A.tag,MC,MR,A.Grid())
        P = DistMatrix(iTag,MC,MR,A.Grid())
        if   A.tag == sTag: lib.ElQRColPivExplicitDist_s(A.obj,R.obj,P.obj)
        elif A.tag == dTag: lib.ElQRColPivExplicitDist_d(A.obj,R.obj,P.obj)
        elif A.tag == cTag: lib.ElQRColPivExplicitDist_c(A.obj,R.obj,P.obj)
        elif A.tag == zTag: lib.ElQRColPivExplicitDist_z(A.obj,R.obj,P.obj)
        else: raise Exception('Unsupported datatype')
        return R, P
      else: 
        raise Exception('Partial pivoted explicit fact\'s not yet supported')
    else:
      if factType == QR_IMPLICIT:
        t = DistMatrix(A.tag,MC,STAR,A.Grid())
        d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
        if   A.tag == sTag: lib.ElQRDist_s(A.obj,t.obj,d.obj)
        elif A.tag == dTag: lib.ElQRDist_d(A.obj,t.obj,d.obj)
        elif A.tag == cTag: lib.ElQRDist_c(A.obj,t.obj,d.obj)
        elif A.tag == zTag: lib.ElQRDist_z(A.obj,t.obj,d.obj)
        else: raise Exception('Unsupported datatype')
        return t, d
      elif factType == QR_EXPLICIT:
        R = Matrix(A.tag)
        if   A.tag == sTag: lib.ElQRExplicitDist_s(A.obj,R.obj)
        elif A.tag == dTag: lib.ElQRExplicitDist_d(A.obj,R.obj)
        elif A.tag == cTag: lib.ElQRExplicitDist_c(A.obj,R.obj)
        elif A.tag == zTag: lib.ElQRExplicitDist_z(A.obj,R.obj)
        else: raise Exception('Unsupported datatype')
        return R
      elif factType == QR_EXPLICIT_TRIANG:
        if   A.tag == sTag: lib.ElQRExplicitTriangDist_s(A.obj)
        elif A.tag == dTag: lib.ElQRExplicitTriangDist_d(A.obj)
        elif A.tag == cTag: lib.ElQRExplicitTriangDist_c(A.obj)
        elif A.tag == zTag: lib.ElQRExplicitTriangDist_z(A.obj)
        else: raise Exception('Unsupported datatype')
      elif factType == QR_EXPLICIT_UNITARY:
        if   A.tag == sTag: lib.ElQRExplicitUnitaryDist_s(A.obj)
        elif A.tag == dTag: lib.ElQRExplicitUnitaryDist_d(A.obj)
        elif A.tag == cTag: lib.ElQRExplicitUnitaryDist_c(A.obj)
        elif A.tag == zTag: lib.ElQRExplicitUnitaryDist_z(A.obj)
        else: raise Exception('Unsupported datatype')
      else: raise Exception('Unsupported QR factorization type')
  else: raise Exception('Unsupported matrix type')

lib.ElCholeskyQR_s.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQR_s.restype = c_uint
lib.ElCholeskyQR_d.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQR_d.restype = c_uint
lib.ElCholeskyQR_c.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQR_c.restype = c_uint
lib.ElCholeskyQR_z.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQR_z.restype = c_uint
lib.ElCholeskyQRDist_s.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQRDist_s.restype = c_uint
lib.ElCholeskyQRDist_d.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQRDist_d.restype = c_uint
lib.ElCholeskyQRDist_c.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQRDist_c.restype = c_uint
lib.ElCholeskyQRDist_z.argtypes = [c_void_p,c_void_p]
lib.ElCholeskyQRDist_z.restype = c_uint
def CholeskyQR(A):
  if type(A) is Matrix:
    R = Matrix(A.tag)
    if   A.tag == sTag: lib.ElCholeskyQR_s(A.obj,R.obj)
    elif A.tag == dTag: lib.ElCholeskyQR_d(A.obj,R.obj)
    elif A.tag == cTag: lib.ElCholeskyQR_c(A.obj,R.obj)
    elif A.tag == zTag: lib.ElCholeskyQR_z(A.obj,R.obj)
    else: raise Exception('Unsupported datatype')
    return R
  elif type(A) is DistMatrix:
    R = DistMatrix(A.tag,STAR,STAR,A.Grid())
    if   A.tag == sTag: lib.ElCholeskyQRDist_s(A.obj,R.obj)
    elif A.tag == dTag: lib.ElCholeskyQRDist_d(A.obj,R.obj)
    elif A.tag == cTag: lib.ElCholeskyQRDist_c(A.obj,R.obj)
    elif A.tag == zTag: lib.ElCholeskyQRDist_z(A.obj,R.obj)
    else: raise Exception('Unsupported datatype')
    return R
  else: raise Exception('Unsupported matrix type')

lib.ElApplyQAfterQR_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQR_s.restype = c_uint
lib.ElApplyQAfterQR_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQR_d.restype = c_uint
lib.ElApplyQAfterQR_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQR_c.restype = c_uint
lib.ElApplyQAfterQR_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQR_z.restype = c_uint
lib.ElApplyQAfterQRDist_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQRDist_s.restype = c_uint
lib.ElApplyQAfterQRDist_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQRDist_d.restype = c_uint
lib.ElApplyQAfterQRDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQRDist_c.restype = c_uint
lib.ElApplyQAfterQRDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterQRDist_z.restype = c_uint
def ApplyQAfterQR(side,orient,A,t,d,B):
  if type(A) is not type(t) or type(t) is not type(d) or type(d) is not type(B):
    raise Exception('Matrix types of {A,t,d,B} must match')
  if A.tag != t.tag or t.tag != B.tag:
    raise Exception('Datatypes of {A,t,B} must match')
  if d.tag != Base(A.tag):
    raise Exception('Base type of A must match that of d')
  if type(A) is Matrix:
    if   A.tag == sTag:
      lib.ElApplyQAfterQR_s(side,orient,A.tag,t.tag,d.tag,B.tag)
    elif A.tag == dTag:
      lib.ElApplyQAfterQR_d(side,orient,A.tag,t.tag,d.tag,B.tag)
    elif A.tag == cTag:
      lib.ElApplyQAfterQR_c(side,orient,A.tag,t.tag,d.tag,B.tag)
    elif A.tag == zTag:
      lib.ElApplyQAfterQR_z(side,orient,A.tag,t.tag,d.tag,B.tag)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag:
      lib.ElApplyQAfterQRDist_s(side,orient,A.tag,t.tag,d.tag,B.tag)
    elif A.tag == dTag:
      lib.ElApplyQAfterQRDist_d(side,orient,A.tag,t.tag,d.tag,B.tag)
    elif A.tag == cTag:
      lib.ElApplyQAfterQRDist_c(side,orient,A.tag,t.tag,d.tag,B.tag)
    elif A.tag == zTag:
      lib.ElApplyQAfterQRDist_z(side,orient,A.tag,t.tag,d.tag,B.tag)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# TODO: TSQR
# TODO: ExplicitTSQR

# RQ factorization
# ================
(RQ_IMPLICIT,RQ_EXPLICIT,RQ_EXPLICIT_TRIANG,RQ_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElRQ_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQ_s.restype = c_uint
lib.ElRQ_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQ_d.restype = c_uint
lib.ElRQ_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQ_c.restype = c_uint
lib.ElRQ_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQ_z.restype = c_uint
lib.ElRQDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQDist_s.restype = c_uint
lib.ElRQDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQDist_d.restype = c_uint
lib.ElRQDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQDist_c.restype = c_uint
lib.ElRQDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRQDist_z.restype = c_uint
#lib.ElRQExplicit_s.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicit_s.restype = c_uint
#lib.ElRQExplicit_d.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicit_d.restype = c_uint
#lib.ElRQExplicit_c.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicit_c.restype = c_uint
#lib.ElRQExplicit_z.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicit_z.restype = c_uint
#lib.ElRQExplicitDist_s.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicitDist_s.restype = c_uint
#lib.ElRQExplicitDist_d.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicitDist_d.restype = c_uint
#lib.ElRQExplicitDist_c.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicitDist_c.restype = c_uint
#lib.ElRQExplicitDist_z.argtypes = [c_void_p,c_void_p]
#lib.ElRQExplicitDist_z.restype = c_uint
#lib.ElRQExplicitUnitary_s.argtypes = [c_void_p]
#lib.ElRQExplicitUnitary_s.restype = c_uint
#lib.ElRQExplicitUnitary_d.argtypes = [c_void_p]
#lib.ElRQExplicitUnitary_d.restype = c_uint
#lib.ElRQExplicitUnitary_c.argtypes = [c_void_p]
#lib.ElRQExplicitUnitary_c.restype = c_uint
#lib.ElRQExplicitUnitary_z.argtypes = [c_void_p]
#lib.ElRQExplicitUnitary_z.restype = c_uint
#lib.ElRQExplicitUnitaryDist_s.argtypes = [c_void_p]
#lib.ElRQExplicitUnitaryDist_s.restype = c_uint
#lib.ElRQExplicitUnitaryDist_d.argtypes = [c_void_p]
#lib.ElRQExplicitUnitaryDist_d.restype = c_uint
#lib.ElRQExplicitUnitaryDist_c.argtypes = [c_void_p]
#lib.ElRQExplicitUnitaryDist_c.restype = c_uint
#lib.ElRQExplicitUnitaryDist_z.argtypes = [c_void_p]
#lib.ElRQExplicitUnitaryDist_z.restype = c_uint
lib.ElRQExplicitTriang_s.argtypes = [c_void_p]
lib.ElRQExplicitTriang_s.restype = c_uint
lib.ElRQExplicitTriang_d.argtypes = [c_void_p]
lib.ElRQExplicitTriang_d.restype = c_uint
lib.ElRQExplicitTriang_c.argtypes = [c_void_p]
lib.ElRQExplicitTriang_c.restype = c_uint
lib.ElRQExplicitTriang_z.argtypes = [c_void_p]
lib.ElRQExplicitTriang_z.restype = c_uint
lib.ElRQExplicitTriangDist_s.argtypes = [c_void_p]
lib.ElRQExplicitTriangDist_s.restype = c_uint
lib.ElRQExplicitTriangDist_d.argtypes = [c_void_p]
lib.ElRQExplicitTriangDist_d.restype = c_uint
lib.ElRQExplicitTriangDist_c.argtypes = [c_void_p]
lib.ElRQExplicitTriangDist_c.restype = c_uint
lib.ElRQExplicitTriangDist_z.argtypes = [c_void_p]
lib.ElRQExplicitTriangDist_z.restype = c_uint
def RQ(A,factType=RQ_IMPLICIT):
  if type(A) is Matrix:
    if factType == RQ_IMPLICIT:
      t = Matrix(A.tag)
      d = Matrix(Base(A.tag))
      if   A.tag == sTag: lib.ElRQ_s(A.obj,t.obj,d.obj)
      elif A.tag == dTag: lib.ElRQ_d(A.obj,t.obj,d.obj)
      elif A.tag == cTag: lib.ElRQ_c(A.obj,t.obj,d.obj)
      elif A.tag == zTag: lib.ElRQ_z(A.obj,t.obj,d.obj)
      else: raise Exception('Unsupported datatype')
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
      if   A.tag == sTag: lib.ElRQExplicitTriang_s(A.obj)
      elif A.tag == dTag: lib.ElRQExplicitTriang_d(A.obj)
      elif A.tag == cTag: lib.ElRQExplicitTriang_c(A.obj)
      elif A.tag == zTag: lib.ElRQExplicitTriang_z(A.obj)
      else: raise Exception('Unsupported datatype')
    else: raise Exception('Unsupported RQ factorization type')
  elif type(A) is DistMatrix:
    if factType == RQ_IMPLICIT:
      t = DistMatrix(A.tag,MC,STAR,A.Grid())
      d = DistMatrix(Base(A.tag),MC,STAR,A.Grid())
      if   A.tag == sTag: lib.ElRQDist_s(A.obj,t.obj,d.obj)
      elif A.tag == dTag: lib.ElRQDist_d(A.obj,t.obj,d.obj)
      elif A.tag == cTag: lib.ElRQDist_c(A.obj,t.obj,d.obj)
      elif A.tag == zTag: lib.ElRQDist_z(A.obj,t.obj,d.obj)
      else: raise Exception('Unsupported datatype')
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
      if   A.tag == sTag: lib.ElRQExplicitTriangDist_s(A.obj)
      elif A.tag == dTag: lib.ElRQExplicitTriangDist_d(A.obj)
      elif A.tag == cTag: lib.ElRQExplicitTriangDist_c(A.obj)
      elif A.tag == zTag: lib.ElRQExplicitTriangDist_z(A.obj)
      else: raise Exception('Unsupported datatype')
    else: raise Exception('Unsupported RQ factorization type')
  else: raise Exception('Unsupported datatype')

lib.ElApplyQAfterRQ_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQ_s.restype = c_uint
lib.ElApplyQAfterRQ_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQ_d.restype = c_uint
lib.ElApplyQAfterRQ_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQ_c.restype = c_uint
lib.ElApplyQAfterRQ_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQ_z.restype = c_uint
lib.ElApplyQAfterRQDist_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQDist_s.restype = c_uint
lib.ElApplyQAfterRQDist_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQDist_d.restype = c_uint
lib.ElApplyQAfterRQDist_c.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQDist_c.restype = c_uint
lib.ElApplyQAfterRQDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElApplyQAfterRQDist_z.restype = c_uint
def ApplyQAfterRQ(side,orient,A,t,d,B):
  if type(A) is not type(t) or type(t) is not type(d) or type(d) is not type(B):
    raise Exception('Matrix types of {A,t,d,B} must match')
  if A.tag != t.tag or t.tag != B.tag or d.tag != Base(A.tag):
    raise Exception('Datatypes of {A,t,B} must match and d must have base type')
  if type(A) is Matrix:
    if   A.tag == sTag: 
      lib.ElApplyQAfterRQ_s(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == dTag:
      lib.ElApplyQAfterRQ_d(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == cTag:
      lib.ElApplyQAfterRQ_c(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == zTag:
      lib.ElApplyQAfterRQ_z(side,orient,A.obj,t.obj,d.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElApplyQAfterRQDist_s(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == dTag:
      lib.ElApplyQAfterRQDist_d(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == cTag:
      lib.ElApplyQAfterRQDist_c(side,orient,A.obj,t.obj,d.obj,B.obj)
    elif A.tag == zTag:
      lib.ElApplyQAfterRQDist_z(side,orient,A.obj,t.obj,d.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElSolveAfterRQ_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQ_s.restype = c_uint
lib.ElSolveAfterRQ_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQ_d.restype = c_uint
lib.ElSolveAfterRQ_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQ_c.restype = c_uint
lib.ElSolveAfterRQ_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQ_z.restype = c_uint
lib.ElSolveAfterRQDist_s.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQDist_s.restype = c_uint
lib.ElSolveAfterRQDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQDist_d.restype = c_uint
lib.ElSolveAfterRQDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQDist_c.restype = c_uint
lib.ElSolveAfterRQDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElSolveAfterRQDist_z.restype = c_uint
def SolveAfterRQ(orient,A,t,d,B):
  if type(A) is Matrix:
    X = Matrix(A.tag)
    if   A.tag == sTag: 
      lib.ElSolveAfterRQ_s(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == dTag:
      lib.ElSolveAfterRQ_d(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == cTag:
      lib.ElSolveAfterRQ_c(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == zTag:
      lib.ElSolveAfterRQ_z(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    if   A.tag == sTag: 
      lib.ElSolveAfterRQDist_s(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == dTag:
      lib.ElSolveAfterRQDist_d(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == cTag:
      lib.ElSolveAfterRQDist_c(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    elif A.tag == zTag:
      lib.ElSolveAfterRQDist_z(orient,A.obj,t.obj,d.obj,B.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  else: raise Exception('Unsupported matrix type')

# Generalized QR factorization
# ============================
(GQR_IMPLICIT,GQR_EXPLICIT,GQR_EXPLICIT_TRIANG,GQR_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElGQR_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQR_s.restype = c_uint
lib.ElGQR_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQR_d.restype = c_uint
lib.ElGQR_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQR_c.restype = c_uint
lib.ElGQR_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQR_z.restype = c_uint
lib.ElGQRDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQRDist_s.restype = c_uint
lib.ElGQRDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQRDist_d.restype = c_uint
lib.ElGQRDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQRDist_c.restype = c_uint
lib.ElGQRDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGQRDist_z.restype = c_uint

lib.ElGQRExplicitTriang_s.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriang_s.restype = c_uint
lib.ElGQRExplicitTriang_d.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriang_d.restype = c_uint
lib.ElGQRExplicitTriang_c.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriang_c.restype = c_uint
lib.ElGQRExplicitTriang_z.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriang_z.restype = c_uint
lib.ElGQRExplicitTriangDist_s.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriangDist_s.restype = c_uint
lib.ElGQRExplicitTriangDist_d.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriangDist_d.restype = c_uint
lib.ElGQRExplicitTriangDist_c.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriangDist_c.restype = c_uint
lib.ElGQRExplicitTriangDist_z.argtypes = [c_void_p,c_void_p]
lib.ElGQRExplicitTriangDist_z.restype = c_uint

def GQR(A,B,factType=GQR_IMPLICIT):
  if type(A) is Matrix:
    if factType == GQR_IMPLICIT:
      tA = Matrix(A.tag)
      dA = Matrix(Base(A.tag))
      tB = Matrix(A.tag)
      dB = Matrix(Base(A.tag))
      if   A.tag == sTag: lib.ElGQR_s(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == dTag: lib.ElGQR_d(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == cTag: lib.ElGQR_c(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == zTag: lib.ElGQR_z(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      else: raise Exception('Unsupported datatype')
      return tA, dA, tB, dB
    elif factType == GQR_EXPLICIT:
      raise Exception('GQR_EXPLICIT not yet supported')
    elif factType == GQR_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      if   A.tag == sTag: lib.ElGQRExplicitTriang_s(A.obj,B.obj)
      elif A.tag == dTag: lib.ElGQRExplicitTriang_d(A.obj,B.obj)
      elif A.tag == cTag: lib.ElGQRExplicitTriang_c(A.obj,B.obj)
      elif A.tag == zTag: lib.ElGQRExplicitTriang_z(A.obj,B.obj)
      else: raise Exception('Unsupported datatype')
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
      if   A.tag == sTag: lib.ElGQRDist_s(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == dTag: lib.ElGQRDist_d(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == cTag: lib.ElGQRDist_c(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == zTag: lib.ElGQRDist_z(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      else: raise Exception('Unsupported datatype')
      return tA, dA, tB, dB
    elif factType == GQR_EXPLICIT:
      raise Exception('GQR_EXPLICIT not yet supported')
    elif factType == GQR_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      if   A.tag == sTag: lib.ElGQRExplicitTriangDist_s(A.obj,B.obj)
      elif A.tag == dTag: lib.ElGQRExplicitTriangDist_d(A.obj,B.obj)
      elif A.tag == cTag: lib.ElGQRExplicitTriangDist_c(A.obj,B.obj)
      elif A.tag == zTag: lib.ElGQRExplicitTriangDist_z(A.obj,B.obj)
      else: raise Exception('Unsupported datatype')
      return B
    elif factType == GQR_EXPLICIT_UNITARY:
      raise Exception('GQR_EXPLICIT_UNITARY not yet supported') 
    else: raise Exception('Unsupported GQR factorization type')
  else: raise Exception('Unsupported matrix type')

# Generalized RQ factorization
# ============================
(GRQ_IMPLICIT,GRQ_EXPLICIT,GRQ_EXPLICIT_TRIANG,GRQ_EXPLICIT_UNITARY)=(0,1,2,3)

lib.ElGRQ_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQ_s.restype = c_uint
lib.ElGRQ_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQ_d.restype = c_uint
lib.ElGRQ_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQ_c.restype = c_uint
lib.ElGRQ_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQ_z.restype = c_uint
lib.ElGRQDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQDist_s.restype = c_uint
lib.ElGRQDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQDist_d.restype = c_uint
lib.ElGRQDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQDist_c.restype = c_uint
lib.ElGRQDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGRQDist_z.restype = c_uint

lib.ElGRQExplicitTriang_s.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriang_s.restype = c_uint
lib.ElGRQExplicitTriang_d.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriang_d.restype = c_uint
lib.ElGRQExplicitTriang_c.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriang_c.restype = c_uint
lib.ElGRQExplicitTriang_z.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriang_z.restype = c_uint
lib.ElGRQExplicitTriangDist_s.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriangDist_s.restype = c_uint
lib.ElGRQExplicitTriangDist_d.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriangDist_d.restype = c_uint
lib.ElGRQExplicitTriangDist_c.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriangDist_c.restype = c_uint
lib.ElGRQExplicitTriangDist_z.argtypes = [c_void_p,c_void_p]
lib.ElGRQExplicitTriangDist_z.restype = c_uint

def GRQ(A,B,factType=GRQ_IMPLICIT):
  if type(A) is Matrix:
    if factType == GRQ_IMPLICIT:
      tA = Matrix(A.tag)
      dA = Matrix(Base(A.tag))
      tB = Matrix(A.tag)
      dB = Matrix(Base(A.tag))
      if   A.tag == sTag: lib.ElGRQ_s(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == dTag: lib.ElGRQ_d(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == cTag: lib.ElGRQ_c(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == zTag: lib.ElGRQ_z(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      else: raise Exception('Unsupported datatype')
      return tA, dA, tB, dB
    elif factType == GRQ_EXPLICIT:
      raise Exception('GRQ_EXPLICIT not yet supported')
    elif factType == GRQ_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      if   A.tag == sTag: lib.ElGRQExplicitTriang_s(A.obj,B.obj)
      elif A.tag == dTag: lib.ElGRQExplicitTriang_d(A.obj,B.obj)
      elif A.tag == cTag: lib.ElGRQExplicitTriang_c(A.obj,B.obj)
      elif A.tag == zTag: lib.ElGRQExplicitTriang_z(A.obj,B.obj)
      else: raise Exception('Unsupported datatype')
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
      if   A.tag == sTag: lib.ElGRQDist_s(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == dTag: lib.ElGRQDist_d(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == cTag: lib.ElGRQDist_c(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      elif A.tag == zTag: lib.ElGRQDist_z(A.obj,tA.obj,dA.obj,tB.obj,dB.obj)
      else: raise Exception('Unsupported datatype')
      return tA, dA, tB, dB
    elif factType == GRQ_EXPLICIT:
      raise Exception('GRQ_EXPLICIT not yet supported')
    elif factType == GRQ_EXPLICIT_TRIANG:
      B = Matrix(A.tag)
      if   A.tag == sTag: lib.ElGRQExplicitTriangDist_s(A.obj,B.obj)
      elif A.tag == dTag: lib.ElGRQExplicitTriangDist_d(A.obj,B.obj)
      elif A.tag == cTag: lib.ElGRQExplicitTriangDist_c(A.obj,B.obj)
      elif A.tag == zTag: lib.ElGRQExplicitTriangDist_z(A.obj,B.obj)
      else: raise Exception('Unsupported datatype')
      return B
    elif factType == GRQ_EXPLICIT_UNITARY:
      raise Exception('GRQ_EXPLICIT_UNITARY not yet supported') 
    else: raise Exception('Unsupported GRQ factorization type')
  else: raise Exception('Unsupported matrix type')

# Interpolative decomposition
# ===========================
lib.ElID_s.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_s,bType]
lib.ElID_s.restype = c_uint
lib.ElID_d.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_d,bType]
lib.ElID_d.restype = c_uint
lib.ElID_c.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_s,bType]
lib.ElID_c.restype = c_uint
lib.ElID_z.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_d,bType]
lib.ElID_z.restype = c_uint
lib.ElIDDist_s.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_s,bType]
lib.ElIDDist_s.restype = c_uint
lib.ElIDDist_d.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_d,bType]
lib.ElIDDist_d.restype = c_uint
lib.ElIDDist_c.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_s,bType]
lib.ElIDDist_c.restype = c_uint
lib.ElIDDist_z.argtypes = [c_void_p,c_void_p,c_void_p,QRCtrl_d,bType]
lib.ElIDDist_z.restype = c_uint
def ID(A,ctrl,canOverwrite=False):
  if type(A) is Matrix: 
    p = Matrix(iTag)
    Z = Matrix(A.tag)
    if   A.tag == sTag: lib.ElID_s(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    elif A.tag == dTag: lib.ElID_d(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    elif A.tag == cTag: lib.ElID_c(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    elif A.tag == zTag: lib.ElID_z(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    else: raise Exception('Unsupported datatype')
    return p, Z
  elif type(A) is DistMatrix:
    p = DistMatrix(iTag,MC,STAR,A.Grid())
    Z = DistMatrix(A.tag,STAR,VR,A.Grid())
    if   A.tag == sTag: lib.ElIDDist_s(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    elif A.tag == dTag: lib.ElIDDist_d(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    elif A.tag == cTag: lib.ElIDDist_c(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    elif A.tag == zTag: lib.ElIDDist_z(A.obj,p.obj,Z.obj,ctrl,canOverwrite)
    else: raise Exception('Unsupported datatype')
    return p, Z
  else: raise Exception('Unsupported matrix type')

# Skeleton decomposition
# ======================
lib.ElSkeleton_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElSkeleton_s.restype = c_uint
lib.ElSkeleton_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElSkeleton_d.restype = c_uint
lib.ElSkeleton_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElSkeleton_c.restype = c_uint
lib.ElSkeleton_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElSkeleton_z.restype = c_uint
lib.ElSkeletonDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElSkeletonDist_s.restype = c_uint
lib.ElSkeletonDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElSkeletonDist_d.restype = c_uint
lib.ElSkeletonDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_s]
lib.ElSkeletonDist_c.restype = c_uint
lib.ElSkeletonDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,QRCtrl_d]
lib.ElSkeletonDist_z.restype = c_uint
def Skeleton(A,ctrl):
  if type(A) is Matrix:
    pR = Matrix(iTag)
    pC = Matrix(iTag)
    Z = Matrix(A.tag)
    if   A.tag == sTag: lib.ElSkeleton_s(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    elif A.tag == dTag: lib.ElSkeleton_d(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    elif A.tag == cTag: lib.ElSkeleton_c(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    elif A.tag == zTag: lib.ElSkeleton_z(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    else: raise Exception('Unsupported datatype')
    return pR, pC, Z
  elif type(A) is DistMatrix:
    pR = DistMatrix(iTag,MC,STAR,A.Grid())
    pC = DistMatrix(iTag,MC,STAR,A.Grid())
    Z = DistMatrix(A.tag,STAR,VR,A.Grid())
    if   A.tag == sTag: lib.ElSkeletonDist_s(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    elif A.tag == dTag: lib.ElSkeletonDist_d(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    elif A.tag == cTag: lib.ElSkeletonDist_c(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    elif A.tag == zTag: lib.ElSkeletonDist_z(A.obj,pR.obj,pC.obj,Z.obj,ctrl)
    else: raise Exception('Unsupported datatype')
    return pR, pC, Z
  else: raise Exception('Unsupported matrix type')
