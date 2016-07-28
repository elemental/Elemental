#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *

# BLAS 3
# ======

# Gemm
# ----

# Emulate an enum for the Gemm algorithm
(GEMM_DEFAULT,GEMM_SUMMA_A,GEMM_SUMMA_B,GEMM_SUMMA_C,GEMM_SUMMA_DOT,
 GEMM_CANNON)=(0,1,2,3,4,5)

lib.ElGemm_i.argtypes = [c_uint,c_uint,iType,c_void_p,c_void_p,iType,c_void_p]
lib.ElGemm_s.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElGemm_d.argtypes = [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElGemm_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElGemm_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElGemmXDist_i.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p,iType,c_void_p,c_uint]
lib.ElGemmXDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElGemmXDist_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p,c_uint]
lib.ElGemmXDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p,c_uint]
lib.ElGemmXDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p,c_uint]

def Gemm(orientA,orientB,alphaPre,A,B,betaPre,C,alg=GEMM_DEFAULT):
  if A.tag != B.tag or B.tag != C.tag:
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C):
    raise Exception('Matrix types of {A,B,C} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  args = [orientA,orientB,alpha,A.obj,B.obj,beta,C.obj]
  argsAlg = [orientA,orientB,alpha,A.obj,B.obj,beta,C.obj,alg]
  if type(A) is Matrix:
    if   B.tag == iTag: lib.ElGemm_i(*args)
    elif B.tag == sTag: lib.ElGemm_s(*args)
    elif B.tag == dTag: lib.ElGemm_d(*args)
    elif B.tag == cTag: lib.ElGemm_c(*args)
    elif B.tag == zTag: lib.ElGemm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   B.tag == iTag: lib.ElGemmXDist_i(*argsAlg)
    elif B.tag == sTag: lib.ElGemmXDist_s(*argsAlg)
    elif B.tag == dTag: lib.ElGemmXDist_d(*argsAlg)
    elif B.tag == cTag: lib.ElGemmXDist_c(*argsAlg)
    elif B.tag == zTag: lib.ElGemmXDist_z(*argsAlg)
    else: DataExcept()
  else: TypeExcept()

# Symm/Hemm
# ---------
lib.ElSymm_s.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSymm_d.argtypes = [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSymm_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSymm_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSymmDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSymmDist_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSymmDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSymmDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]

lib.ElHemm_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElHemm_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElHemmDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElHemmDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]

def Symm(side,uplo,alphaPre,A,B,betaPre,C,conj=False):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  args = [side,uplo,alpha,A.obj,B.obj,beta,C.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymm_s(*args)
    elif A.tag == dTag: lib.ElSymm_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElHemm_c(*args)
      else:    lib.ElSymm_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElHemm_z(*args)
      else:    lib.ElSymm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmDist_s(*args)
    elif A.tag == dTag: lib.ElSymmDist_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElHemmDist_c(*args)
      else:    lib.ElSymmDist_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElHemmDist_z(*args)
      else:    lib.ElSymmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

def Hemm(side,uplo,alpha,A,B,betaPre,C):
  Symm(side,uplo,alpha,A,B,betaPre,C,True) 

# Syrk/Herk
# ---------
lib.ElSyrk_s.argtypes = \
lib.ElSyrkDist_s.argtypes = \
lib.ElSyrkSparse_s.argtypes = \
lib.ElSyrkDistSparse_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,sType,c_void_p]

lib.ElSyrk_d.argtypes = \
lib.ElSyrkDist_d.argtypes = \
lib.ElSyrkSparse_d.argtypes = \
lib.ElSyrkDistSparse_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,dType,c_void_p]

lib.ElSyrk_c.argtypes = \
lib.ElSyrkDist_c.argtypes = \
lib.ElSyrkSparse_c.argtypes = \
lib.ElSyrkDistSparse_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,cType,c_void_p]

lib.ElSyrk_z.argtypes = \
lib.ElSyrkDist_z.argtypes = \
lib.ElSyrkSparse_z.argtypes = \
lib.ElSyrkDistSparse_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,zType,c_void_p]

lib.ElHerk_c.argtypes = \
lib.ElHerkDist_c.argtypes = \
lib.ElHerkSparse_c.argtypes = \
lib.ElHerkDistSparse_c.argtypes = \
  [c_uint,c_uint,sType,c_void_p,sType,c_void_p]

lib.ElHerk_z.argtypes = \
lib.ElHerkDist_z.argtypes = \
lib.ElHerkSparse_z.argtypes = \
lib.ElHerkDistSparse_z.argtypes = \
  [c_uint,c_uint,dType,c_void_p,dType,c_void_p]

def Syrk(uplo,orient,alphaPre,A,betaPre,C,conj=False):
  if A.tag != C.tag: raise Exception('Datatypes of A and C must match')
  if type(A) is not type(C): raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  args = [uplo,orient,alpha,A.obj,beta,C.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyrk_s(*args)
    elif A.tag == dTag: lib.ElSyrk_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElHerk_c(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrk_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElHerk_z(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrk_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSyrkDist_s(*args)
    elif A.tag == dTag: lib.ElSyrkDist_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElHerkDist_c(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkDist_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElHerkDist_z(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == sTag: lib.ElSyrkSparse_s(*args)
    elif A.tag == dTag: lib.ElSyrkSparse_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElHerkSparse_c(orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkSparse_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElHerkSparse_z(orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == sTag: lib.ElSyrkDistSparse_s(*args)
    elif A.tag == dTag: lib.ElSyrkDistSparse_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElHerkDistSparse_c(orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkDistSparse_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElHerkDistSparse_z(orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

def Herk(uplo,orient,alpha,A,beta,C):
  Syrk(uplo,orient,alpha,A,beta,C,True)

# Syr2k/Her2k
# -----------
lib.ElSyr2k_s.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSyr2k_d.argtypes = [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSyr2k_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSyr2k_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSyr2kDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSyr2kDist_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSyr2kDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSyr2kDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]

lib.ElHer2k_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHer2k_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,dType,c_void_p]
lib.ElHer2kDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHer2kDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,dType,c_void_p]

def Syr2k(uplo,orient,alphaPre,A,B,beta,C,conj=False):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [uplo,orient,alpha,A.obj,B.obj,beta,C.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyr2k_s(*args)
    elif A.tag == dTag: lib.ElSyr2k_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElHer2k_c(uplo,orient,alpha,A.obj,B.obj,beta.real,C.obj)
      else:    lib.ElSyr2k_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElHer2k_z(uplo,orient,alpha,A.obj,B.obj,beta.real,C.obj)
      else:    lib.ElSyr2k_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSyr2kDist_s(*args)
    elif A.tag == dTag: lib.ElSyr2kDist_d(*args)
    elif A.tag == cTag:
      if conj: 
        lib.ElHer2kDist_c(uplo,orient,alpha,A.obj,B.obj,beta.real,C.obj)
      else: lib.ElSyr2kDist_c(*args)
    elif A.tag == zTag:
      if conj:
        lib.ElHer2kDist_z(uplo,orient,alpha,A.obj,B.obj,beta.real,C.obj)
      else: lib.ElSyr2kDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

def Her2k(uplo,orient,alpha,A,B,beta,C):
  Syr2k(uplo,orient,alpha,A,B,beta,C,True)

# Multiply
# --------
lib.ElMultiply_i.argtypes = \
  [c_uint,iType,c_void_p,c_void_p,iType,c_void_p]
lib.ElMultiply_s.argtypes = \
  [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElMultiply_d.argtypes = \
  [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElMultiply_c.argtypes = \
  [c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElMultiply_z.argtypes = \
  [c_uint,zType,c_void_p,c_void_p,zType,c_void_p]

lib.ElMultiplyDist_i.argtypes = \
  [c_uint,iType,c_void_p,c_void_p,iType,c_void_p]
lib.ElMultiplyDist_s.argtypes = \
  [c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElMultiplyDist_d.argtypes = \
  [c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElMultiplyDist_c.argtypes = \
  [c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElMultiplyDist_z.argtypes = \
  [c_uint,zType,c_void_p,c_void_p,zType,c_void_p]

def Multiply(orient,alpha,A,X,beta,Y):
  if type(A) is SparseMatrix:
    if type(X) is not Matrix or type(Y) is not Matrix:
      raise Exception("Types of X and Y must match")
    if A.tag != X.tag or X.tag != Y.tag:
      raise Exception("Datatypes of {A,X,Y} must match")
    args = [orient,alpha,A.obj,X.obj,beta,Y.obj]
    if   A.tag == iTag: lib.ElMultiply_i(*args)
    elif A.tag == sTag: lib.ElMultiply_s(*args)
    elif A.tag == dTag: lib.ElMultiply_d(*args)
    elif A.tag == cTag: lib.ElMultiply_c(*args)
    elif A.tag == zTag: lib.ElMultiply_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(X) is not DistMultiVec or type(Y) is not DistMultiVec:
      raise Exception("Types of X and Y must match")
    if A.tag != X.tag or X.tag != Y.tag:
      raise Exception("Datatypes of {A,X,Y} must match")
    args = [orient,alpha,A.obj,X.obj,beta,Y.obj]
    if   A.tag == iTag: lib.ElMultiplyDist_i(*args)
    elif A.tag == sTag: lib.ElMultiplyDist_s(*args)
    elif A.tag == dTag: lib.ElMultiplyDist_d(*args)
    elif A.tag == cTag: lib.ElMultiplyDist_c(*args)
    elif A.tag == zTag: lib.ElMultiplyDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# MultiShiftQuasiTrsm
# -------------------
lib.ElMultiShiftQuasiTrsm_s.argtypes = \
lib.ElMultiShiftQuasiTrsmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]

lib.ElMultiShiftQuasiTrsm_d.argtypes = \
lib.ElMultiShiftQuasiTrsmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,dType,c_void_p,c_void_p,c_void_p]

lib.ElMultiShiftQuasiTrsm_c.argtypes = \
lib.ElMultiShiftQuasiTrsmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,cType,c_void_p,c_void_p,c_void_p]

lib.ElMultiShiftQuasiTrsm_z.argtypes = \
lib.ElMultiShiftQuasiTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,zType,c_void_p,c_void_p,c_void_p]

def MultiShiftQuasiTrsm(side,uplo,orient,alphaPre,A,shifts,B):
  if type(A) is not type(shifts) or type(shifts) is not type(B): 
    raise Exception('Types of A and B must match')
  if A.tag != shifts.tag or shifts.tag != B.tag: 
    raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [side,uplo,orient,alpha,A.obj,shifts.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElMultiShiftQuasiTrsm_s(*args)
    elif A.tag == dTag: lib.ElMultiShiftQuasiTrsm_d(*args)
    elif A.tag == cTag: lib.ElMultiShiftQuasiTrsm_c(*args)
    elif A.tag == zTag: lib.ElMultiShiftQuasiTrsm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElMultiShiftQuasiTrsmDist_s(*args)
    elif A.tag == dTag: lib.ElMultiShiftQuasiTrsmDist_d(*args)
    elif A.tag == cTag: lib.ElMultiShiftQuasiTrsmDist_c(*args)
    elif A.tag == zTag: lib.ElMultiShiftQuasiTrsmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# MultiShiftTrsm
# --------------
lib.ElMultiShiftTrsm_s.argtypes = \
lib.ElMultiShiftTrsmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]

lib.ElMultiShiftTrsm_d.argtypes = \
lib.ElMultiShiftTrsmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,dType,c_void_p,c_void_p,c_void_p]

lib.ElMultiShiftTrsm_c.argtypes = \
lib.ElMultiShiftTrsmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,cType,c_void_p,c_void_p,c_void_p]

lib.ElMultiShiftTrsm_z.argtypes = \
lib.ElMultiShiftTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,zType,c_void_p,c_void_p,c_void_p]

def MultiShiftTrsm(side,uplo,orient,alphaPre,A,shifts,B):
  if type(A) is not type(shifts) or type(shifts) is not type(B): 
    raise Exception('Types of A and B must match')
  if A.tag != shifts.tag or shifts.tag != B.tag: 
    raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [side,uplo,orient,alpha,A.obj,shifts.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElMultiShiftTrsm_s(*args)
    elif A.tag == dTag: lib.ElMultiShiftTrsm_d(*args)
    elif A.tag == cTag: lib.ElMultiShiftTrsm_c(*args)
    elif A.tag == zTag: lib.ElMultiShiftTrsm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElMultiShiftTrsmDist_s(*args)
    elif A.tag == dTag: lib.ElMultiShiftTrsmDist_d(*args)
    elif A.tag == cTag: lib.ElMultiShiftTrsmDist_c(*args)
    elif A.tag == zTag: lib.ElMultiShiftTrsmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# QuasiTrsm
# ---------
lib.ElQuasiTrsm_s.argtypes = \
lib.ElQuasiTrsmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]

lib.ElQuasiTrsm_d.argtypes = \
lib.ElQuasiTrsmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]

lib.ElQuasiTrsm_c.argtypes = \
lib.ElQuasiTrsmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]

lib.ElQuasiTrsm_z.argtypes = \
lib.ElQuasiTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]

def QuasiTrsm(side,uplo,orient,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [side,uplo,orient,alpha,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElQuasiTrsm_s(*args)
    elif A.tag == dTag: lib.ElQuasiTrsm_d(*args)
    elif A.tag == cTag: lib.ElQuasiTrsm_c(*args)
    elif A.tag == zTag: lib.ElQuasiTrsm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElQuasiTrsmDist_s(*args)
    elif A.tag == dTag: lib.ElQuasiTrsmDist_d(*args)
    elif A.tag == cTag: lib.ElQuasiTrsmDist_c(*args)
    elif A.tag == zTag: lib.ElQuasiTrsmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trdtrmm
# -------
lib.ElTrdtrmm_s.argtypes = \
lib.ElTrdtrmm_d.argtypes = \
lib.ElTrdtrmmDist_s.argtypes = \
lib.ElTrdtrmmDist_d.argtypes = \
  [c_uint,c_void_p]

lib.ElTrdtrmm_c.argtypes = \
lib.ElTrdtrmm_z.argtypes = \
lib.ElTrdtrmmDist_c.argtypes = \
lib.ElTrdtrmmDist_z.argtypes = \
  [c_uint,c_void_p,bType]

def Trdtrmm(uplo,A,conjugate=False):
  args = [uplo,A.obj]
  argsCpx = [uplo,A.obj,conjugate]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrdtrmm_s(*args)
    elif A.tag == dTag: lib.ElTrdtrmm_d(*args)
    elif A.tag == cTag: lib.ElTrdtrmm_c(*argsCpx)
    elif A.tag == zTag: lib.ElTrdtrmm_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrdtrmmDist_s(*args)
    elif A.tag == dTag: lib.ElTrdtrmmDist_d(*args)
    elif A.tag == cTag: lib.ElTrdtrmmDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElTrdtrmmDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# TrdtrmmQuasi
# ------------
lib.ElTrdtrmmQuasi_s.argtypes = \
lib.ElTrdtrmmQuasi_d.argtypes = \
lib.ElTrdtrmmQuasiDist_s.argtypes = \
lib.ElTrdtrmmQuasiDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p]

lib.ElTrdtrmmQuasi_c.argtypes = \
lib.ElTrdtrmmQuasi_z.argtypes = \
lib.ElTrdtrmmQuasiDist_c.argtypes = \
lib.ElTrdtrmmQuasiDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,bType]

def TrdtrmmQuasi(uplo,A,dOff,conjugate=False):
  if type(A) is not type(dOff): 
    raise Exception('Types of A and dOff must match')
  if A.tag != dOff.tag:
    raise Exception('Datatypes of A and dOff must match')
  args = [uplo,A.obj]
  argsCpx = [uplo,A.obj,conjugate]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrdtrmmQuasi_s(*args)
    elif A.tag == dTag: lib.ElTrdtrmmQuasi_d(*args)
    elif A.tag == cTag: lib.ElTrdtrmmQuasi_c(*argsCpx)
    elif A.tag == zTag: lib.ElTrdtrmmQuasi_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrdtrmmQuasiDist_s(*args)
    elif A.tag == dTag: lib.ElTrdtrmmQuasiDist_d(*args)
    elif A.tag == cTag: lib.ElTrdtrmmQuasiDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElTrdtrmmQuasiDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# Trmm
# ----
lib.ElTrmm_s.argtypes = \
lib.ElTrmmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]

lib.ElTrmm_d.argtypes = \
lib.ElTrmmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]

lib.ElTrmm_c.argtypes = \
lib.ElTrmmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]

lib.ElTrmm_z.argtypes = \
lib.ElTrmmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]

def Trmm(side,uplo,orient,diag,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [side,uplo,orient,diag,alpha,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrmm_s(*args)
    elif A.tag == dTag: lib.ElTrmm_d(*args)
    elif A.tag == cTag: lib.ElTrmm_c(*args)
    elif A.tag == zTag: lib.ElTrmm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrmmDist_s(*args)
    elif A.tag == dTag: lib.ElTrmmDist_d(*args)
    elif A.tag == cTag: lib.ElTrmmDist_c(*args)
    elif A.tag == zTag: lib.ElTrmmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trrk
# ----
lib.ElTrrk_s.argtypes = \
lib.ElTrrkDist_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]

lib.ElTrrk_d.argtypes = \
lib.ElTrrkDist_d.argtypes = \
  [c_uint,c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]

lib.ElTrrk_c.argtypes = \
lib.ElTrrkDist_c.argtypes = \
  [c_uint,c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]

lib.ElTrrk_z.argtypes = \
lib.ElTrrkDist_z.argtypes = \
  [c_uint,c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]

def Trrk(uplo,orientA,orientB,alphaPre,A,B,betaPre,C):
  if type(A) is not type(B) or type(B) is not type(C):
    raise Exception('Types of {A,B,C} must match')
  if A.tag != B.tag or B.tag != C.tag:
    raise Exception('Datatypes of {A,B,C} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  args = [uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrrk_s(*args)
    elif A.tag == dTag: lib.ElTrrk_d(*args)
    elif A.tag == cTag: lib.ElTrrk_c(*args)
    elif A.tag == zTag: lib.ElTrrk_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrrkDist_s(*args)
    elif A.tag == dTag: lib.ElTrrkDist_d(*args)
    elif A.tag == cTag: lib.ElTrrkDist_c(*args)
    elif A.tag == zTag: lib.ElTrrkDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trr2k
# -----
#lib.ElTrr2k_s.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   sType,c_void_p,c_void_p,sType,c_void_p,c_void_p,sType,c_void_p]
#lib.ElTrr2k_d.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   dType,c_void_p,c_void_p,dType,c_void_p,c_void_p,dType,c_void_p]
#lib.ElTrr2k_c.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   cType,c_void_p,c_void_p,cType,c_void_p,c_void_p,cType,c_void_p]
#lib.ElTrr2k_z.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   zType,c_void_p,c_void_p,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElTrr2kDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   sType,c_void_p,c_void_p,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElTrr2kDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   dType,c_void_p,c_void_p,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElTrr2kDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   cType,c_void_p,c_void_p,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElTrr2kDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   zType,c_void_p,c_void_p,zType,c_void_p,c_void_p,zType,c_void_p]
def Trr2k(uplo,orientA,orientB,orientC,orientD,
          alphaPre,A,B,betaPre,C,D,gammaPre,E):
  if type(A) is not type(B) or type(B) is not type(C) or \
     type(C) is not type(D) or type(D) is not type(E):
    raise Exception('Types of {A,B,C,D,E} must match')
  if A.tag != B.tag or B.tag != C.tag or C.tag != D.tag or D.tag != E.tag:
    raise Exception('Datatypes of {A,B,C,D,E} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  gamma = TagToType(A.tag)(gammaPre)
  args = [uplo,orientA,orientB,orientC,orientD,
          alpha,A.obj,B.obj,beta,C.obj,D.obj,gamma,E.obj]
  if type(A) is Matrix:
    raise Exception('Sequential implementation does not yet exist')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrr2kDist_s(*args)
    elif A.tag == dTag: lib.ElTrr2kDist_d(*args)
    elif A.tag == cTag: lib.ElTrr2kDist_c(*args)
    elif A.tag == zTag: lib.ElTrr2kDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trsm
# ----
lib.ElTrsm_s.argtypes = \
lib.ElTrsmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]

lib.ElTrsm_d.argtypes = \
lib.ElTrsmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]

lib.ElTrsm_c.argtypes = \
lib.ElTrsmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]

lib.ElTrsm_z.argtypes = \
lib.ElTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]

def Trsm(side,uplo,orient,diag,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [side,uplo,orient,diag,alpha,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrsm_s(*args)
    elif A.tag == dTag: lib.ElTrsm_d(*args)
    elif A.tag == cTag: lib.ElTrsm_c(*args)
    elif A.tag == zTag: lib.ElTrsm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrsmDist_s(*args)
    elif A.tag == dTag: lib.ElTrsmDist_d(*args)
    elif A.tag == cTag: lib.ElTrsmDist_c(*args)
    elif A.tag == zTag: lib.ElTrsmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trstrm
# ------
lib.ElTrstrm_s.argtypes = \
lib.ElTrstrmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]

lib.ElTrstrm_d.argtypes = \
lib.ElTrstrmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]

lib.ElTrstrm_c.argtypes = \
lib.ElTrstrmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]

lib.ElTrstrm_z.argtypes = \
lib.ElTrstrmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]

def Trstrm(side,uplo,orient,diag,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  args = [side,uplo,orient,diag,alpha,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrstrm_s(*args)
    elif A.tag == dTag: lib.ElTrstrm_d(*args)
    elif A.tag == cTag: lib.ElTrstrm_c(*args)
    elif A.tag == zTag: lib.ElTrstrm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrstrmDist_s(*args)
    elif A.tag == dTag: lib.ElTrstrmDist_d(*args)
    elif A.tag == cTag: lib.ElTrstrmDist_c(*args)
    elif A.tag == zTag: lib.ElTrstrmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Trtrmm
# ------
lib.ElTrtrmm_s.argtypes = \
lib.ElTrtrmm_d.argtypes = \
lib.ElTrtrmmDist_s.argtypes = \
lib.ElTrtrmmDist_d.argtypes = \
  [c_uint,c_void_p]

lib.ElTrtrmm_c.argtypes = \
lib.ElTrtrmm_z.argtypes = \
lib.ElTrtrmmDist_c.argtypes = \
lib.ElTrtrmmDist_z.argtypes = \
  [c_uint,c_void_p,bType]

def Trtrmm(uplo,A,conjugate=False):
  args = [uplo,A.obj]
  argsCpx = [uplo,A.obj,conjugate]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrtrmm_s(*args)
    elif A.tag == dTag: lib.ElTrtrmm_d(*args)
    elif A.tag == cTag: lib.ElTrtrmm_c(*argsCpx)
    elif A.tag == zTag: lib.ElTrtrmm_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrtrmmDist_s(*args)
    elif A.tag == dTag: lib.ElTrtrmmDist_d(*args)
    elif A.tag == cTag: lib.ElTrtrmmDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElTrtrmmDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# Two-sided Trmm
# --------------
lib.ElTwoSidedTrmm_s.argtypes = \
lib.ElTwoSidedTrmm_d.argtypes = \
lib.ElTwoSidedTrmm_c.argtypes = \
lib.ElTwoSidedTrmm_z.argtypes = \
lib.ElTwoSidedTrmmDist_s.argtypes = \
lib.ElTwoSidedTrmmDist_d.argtypes = \
lib.ElTwoSidedTrmmDist_c.argtypes = \
lib.ElTwoSidedTrmmDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p]

def TwoSidedTrmm(uplo,diag,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  args = [uplo,diag,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoSidedTrmm_s(*args)
    elif A.tag == dTag: lib.ElTwoSidedTrmm_d(*args)
    elif A.tag == cTag: lib.ElTwoSidedTrmm_c(*args)
    elif A.tag == zTag: lib.ElTwoSidedTrmm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoSidedTrmmDist_s(*args)
    elif A.tag == dTag: lib.ElTwoSidedTrmmDist_d(*args)
    elif A.tag == cTag: lib.ElTwoSidedTrmmDist_c(*args)
    elif A.tag == zTag: lib.ElTwoSidedTrmmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Two-sided Trsm
# --------------
lib.ElTwoSidedTrsm_s.argtypes = \
lib.ElTwoSidedTrsm_d.argtypes = \
lib.ElTwoSidedTrsm_c.argtypes = \
lib.ElTwoSidedTrsm_z.argtypes = \
lib.ElTwoSidedTrsmDist_s.argtypes = \
lib.ElTwoSidedTrsmDist_d.argtypes = \
lib.ElTwoSidedTrsmDist_c.argtypes = \
lib.ElTwoSidedTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p]

def TwoSidedTrsm(uplo,diag,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  args = [uplo,diag,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoSidedTrsm_s(*args)
    elif A.tag == dTag: lib.ElTwoSidedTrsm_d(*args)
    elif A.tag == cTag: lib.ElTwoSidedTrsm_c(*args)
    elif A.tag == zTag: lib.ElTwoSidedTrsm_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoSidedTrsmDist_s(*args)
    elif A.tag == dTag: lib.ElTwoSidedTrsmDist_d(*args)
    elif A.tag == cTag: lib.ElTwoSidedTrsmDist_c(*args)
    elif A.tag == zTag: lib.ElTwoSidedTrsmDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
