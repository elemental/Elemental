#
#  Copyright (c) 2009-2014, Jack Poulson
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
lib.ElGemm_i.restype = c_uint
lib.ElGemm_s.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElGemm_s.restype = c_uint
lib.ElGemm_d.argtypes = [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElGemm_d.restype = c_uint
lib.ElGemm_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElGemm_c.restype = c_uint
lib.ElGemm_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElGemm_z.restype = c_uint
lib.ElGemmXDist_i.argtypes = \
  [c_uint,c_uint,iType,c_void_p,c_void_p,iType,c_void_p,c_uint]
lib.ElGemmXDist_i.restype = c_uint
lib.ElGemmXDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElGemmXDist_s.restype = c_uint
lib.ElGemmXDist_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p,c_uint]
lib.ElGemmXDist_d.restype = c_uint
lib.ElGemmXDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p,c_uint]
lib.ElGemmXDist_c.restype = c_uint
lib.ElGemmXDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p,c_uint]
lib.ElGemmXDist_z.restype = c_uint
def Gemm(orientA,orientB,alphaPre,A,B,betaPre,C,alg=GEMM_DEFAULT):
  if A.tag != B.tag or B.tag != C.tag:
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C):
    raise Exception('Matrix types of {A,B,C} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   B.tag == iTag:
      lib.ElGemm_i(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif B.tag == sTag: 
      lib.ElGemm_s(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif B.tag == dTag:
      lib.ElGemm_d(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif B.tag == cTag:
      lib.ElGemm_c(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif B.tag == zTag:
      lib.ElGemm_z(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   B.tag == iTag: 
      lib.ElGemmXDist_i(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj,alg)
    elif B.tag == sTag:
      lib.ElGemmXDist_s(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj,alg)
    elif B.tag == dTag: 
      lib.ElGemmXDist_d(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj,alg)
    elif B.tag == cTag:
      lib.ElGemmXDist_c(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj,alg)
    elif B.tag == zTag:
      lib.ElGemmXDist_z(orientA,orientB,alpha,A.obj,B.obj,beta,C.obj,alg)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix types')

# Symm/Hemm
# ---------
lib.ElSymm_s.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSymm_s.restype = c_uint
lib.ElSymm_d.argtypes = [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSymm_d.restype = c_uint
lib.ElSymm_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSymm_c.restype = c_uint
lib.ElSymm_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSymm_z.restype = c_uint
lib.ElSymmDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSymmDist_s.restype = c_uint
lib.ElSymmDist_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSymmDist_d.restype = c_uint
lib.ElSymmDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSymmDist_c.restype = c_uint
lib.ElSymmDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSymmDist_z.restype = c_uint

lib.ElHemm_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElHemm_c.restype = c_uint
lib.ElHemm_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElHemm_z.restype = c_uint
lib.ElHemmDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElHemmDist_c.restype = c_uint
lib.ElHemmDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElHemmDist_z.restype = c_uint

def Symm(side,uplo,alphaPre,A,B,betaPre,C,conj=False):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymm_s(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSymm_d(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag:
      if conj: lib.ElHemm_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
      else:    lib.ElSymm_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag:
      if conj: lib.ElHemm_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
      else:    lib.ElSymm_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmDist_s(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSymmDist_d(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: 
      if conj: lib.ElHemmDist_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
      else:    lib.ElSymmDist_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag:
      if conj: lib.ElHemmDist_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
      else:    lib.ElSymmDist_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

def Hemm(side,uplo,alpha,A,B,betaPre,C):
  Symm(side,uplo,alpha,A,B,betaPre,C,True) 

# Syrk/Herk
# ---------
lib.ElSyrk_s.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElSyrk_s.restype = c_uint
lib.ElSyrk_d.argtypes = [c_uint,c_uint,dType,c_void_p,dType,c_void_p]
lib.ElSyrk_d.restype = c_uint
lib.ElSyrk_c.argtypes = [c_uint,c_uint,cType,c_void_p,cType,c_void_p]
lib.ElSyrk_c.restype = c_uint
lib.ElSyrk_z.argtypes = [c_uint,c_uint,zType,c_void_p,zType,c_void_p]
lib.ElSyrk_z.restype = c_uint
lib.ElSyrkDist_s.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElSyrkDist_s.restype = c_uint
lib.ElSyrkDist_d.argtypes = [c_uint,c_uint,dType,c_void_p,dType,c_void_p]
lib.ElSyrkDist_d.restype = c_uint
lib.ElSyrkDist_c.argtypes = [c_uint,c_uint,cType,c_void_p,cType,c_void_p]
lib.ElSyrkDist_c.restype = c_uint
lib.ElSyrkDist_z.argtypes = [c_uint,c_uint,zType,c_void_p,zType,c_void_p]
lib.ElSyrkDist_z.restype = c_uint

lib.ElHerk_c.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElHerk_c.restype = c_uint
lib.ElHerk_z.argtypes = [c_uint,c_uint,dType,c_void_p,dType,c_void_p]
lib.ElHerk_z.restype = c_uint
lib.ElHerkDist_c.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElHerkDist_c.restype = c_uint
lib.ElHerkDist_z.argtypes = [c_uint,c_uint,dType,c_void_p,dType,c_void_p]
lib.ElHerkDist_z.restype = c_uint

def Syrk(uplo,orient,alphaPre,A,betaPre,C,conj=False):
  if A.tag != C.tag: raise Exception('Datatypes of A and C must match')
  if type(A) is not type(C): raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyrk_s(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyrk_d(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == cTag: 
      if conj: lib.ElHerk_c(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrk_c(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == zTag:
      if conj: lib.ElHerk_z(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrk_z(uplo,orient,alpha,A.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSyrkDist_s(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyrkDist_d(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == cTag: 
      if conj: lib.ElHerkDist_c(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkDist_c(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == zTag:
      if conj: lib.ElHerkDist_z(uplo,orient,alpha.real,A.obj,beta.real,C.obj)
      else:    lib.ElSyrkDist_z(uplo,orient,alpha,A.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

def Herk(uplo,orient,alpha,A,beta,C):
  Syrk(uplo,orient,alpha,A,beta,C,True)

# Syr2k/Her2k
# -----------
lib.ElSyr2k_s.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSyr2k_s.restype = c_uint
lib.ElSyr2k_d.argtypes = [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSyr2k_d.restype = c_uint
lib.ElSyr2k_c.argtypes = [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSyr2k_c.restype = c_uint
lib.ElSyr2k_z.argtypes = [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSyr2k_z.restype = c_uint
lib.ElSyr2kDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSyr2kDist_s.restype = c_uint
lib.ElSyr2kDist_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSyr2kDist_d.restype = c_uint
lib.ElSyr2kDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSyr2kDist_c.restype = c_uint
lib.ElSyr2kDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSyr2kDist_z.restype = c_uint

lib.ElHer2k_c.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHer2k_c.restype = c_uint
lib.ElHer2k_z.argtypes = [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElHer2k_z.restype = c_uint
lib.ElHer2kDist_c.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHer2kDist_c.restype = c_uint
lib.ElHer2kDist_z.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElHer2kDist_z.restype = c_uint

def Syr2k(uplo,orient,alphaPre,A,B,betaPre,C,conj=False):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyr2k_s(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyr2k_d(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag:
      if conj: lib.ElHer2k_c(uplo,orient,alpha.real,A.obj,B.obj,beta.real,C.obj)
      else:    lib.ElSyr2k_c(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag:
      if conj: lib.ElHer2k_z(uplo,orient,alpha.real,A.obj,B.obj,beta.real,C.obj)
      else:    lib.ElSyr2k_z(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElSyr2kDist_s(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: 
      lib.ElSyr2kDist_d(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag:
      if conj: 
        lib.ElHer2kDist_c(uplo,orient,alpha.real,A.obj,B.obj,beta.real,C.obj)
      else:
        lib.ElSyr2kDist_c(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag:
      if conj:
        lib.ElHer2kDist_z(uplo,orient,alpha.real,A.obj,B.obj,beta.real,C.obj)
      else:
        lib.ElSyr2kDist_z(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

def Her2k(uplo,orient,alpha,A,B,beta,C):
  Syr2k(uplo,orient,alpha,A,B,beta,C,True)

# MultiShiftQuasiTrsm
# -------------------
lib.ElMultiShiftQuasiTrsm_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsm_s.restype = c_uint
lib.ElMultiShiftQuasiTrsm_d.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsm_d.restype = c_uint
lib.ElMultiShiftQuasiTrsm_c.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsm_c.restype = c_uint
lib.ElMultiShiftQuasiTrsm_z.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsm_z.restype = c_uint
lib.ElMultiShiftQuasiTrsmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsmDist_s.restype = c_uint
lib.ElMultiShiftQuasiTrsmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsmDist_d.restype = c_uint
lib.ElMultiShiftQuasiTrsmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsmDist_c.restype = c_uint
lib.ElMultiShiftQuasiTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftQuasiTrsmDist_z.restype = c_uint
def MultiShiftQuasiTrsm(side,uplo,orient,alphaPre,A,shifts,B):
  if type(A) is not type(shifts) or type(shifts) is not type(B): 
    raise Exception('Types of A and B must match')
  if A.tag != shifts.tag or shifts.tag != B.tag: 
    raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  if type(A) is Matrix:
    if A.tag == sTag: 
      lib.ElMultiShiftQuasiTrsm_s(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == dTag:
      lib.ElMultiShiftQuasiTrsm_d(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == cTag:
      lib.ElMultiShiftQuasiTrsm_c(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == zTag:
      lib.ElMultiShiftQuasiTrsm_z(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if A.tag == sTag: 
      lib.ElMultiShiftQuasiTrsmDist_s \
      (side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == dTag:
      lib.ElMultiShiftQuasiTrsmDist_d \
      (side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == cTag:
      lib.ElMultiShiftQuasiTrsmDist_c \
      (side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == zTag:
      lib.ElMultiShiftQuasiTrsmDist_z \
      (side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# MultiShiftTrsm
# --------------
lib.ElMultiShiftTrsm_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsm_s.restype = c_uint
lib.ElMultiShiftTrsm_d.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsm_d.restype = c_uint
lib.ElMultiShiftTrsm_c.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsm_c.restype = c_uint
lib.ElMultiShiftTrsm_z.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsm_z.restype = c_uint
lib.ElMultiShiftTrsmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsmDist_s.restype = c_uint
lib.ElMultiShiftTrsmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsmDist_d.restype = c_uint
lib.ElMultiShiftTrsmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsmDist_c.restype = c_uint
lib.ElMultiShiftTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftTrsmDist_z.restype = c_uint
def MultiShiftTrsm(side,uplo,orient,alphaPre,A,shifts,B):
  if type(A) is not type(shifts) or type(shifts) is not type(B): 
    raise Exception('Types of A and B must match')
  if A.tag != shifts.tag or shifts.tag != B.tag: 
    raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  if type(A) is Matrix:
    if A.tag == sTag: 
      lib.ElMultiShiftTrsm_s(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == dTag:
      lib.ElMultiShiftTrsm_d(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == cTag:
      lib.ElMultiShiftTrsm_c(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == zTag:
      lib.ElMultiShiftTrsm_z(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if A.tag == sTag: 
      lib.ElMultiShiftTrsmDist_s(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == dTag:
      lib.ElMultiShiftTrsmDist_d(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == cTag:
      lib.ElMultiShiftTrsmDist_c(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    elif A.tag == zTag:
      lib.ElMultiShiftTrsmDist_z(side,uplo,orient,alpha,A.obj,shifts.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# QuasiTrsm
# ---------
lib.ElQuasiTrsm_s.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsm_s.restype = c_uint
lib.ElQuasiTrsm_d.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsm_d.restype = c_uint
lib.ElQuasiTrsm_c.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsm_c.restype = c_uint
lib.ElQuasiTrsm_z.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsm_z.restype = c_uint
lib.ElQuasiTrsmDist_s.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsmDist_s.restype = c_uint
lib.ElQuasiTrsmDist_d.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsmDist_d.restype = c_uint
lib.ElQuasiTrsmDist_c.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsmDist_c.restype = c_uint
lib.ElQuasiTrsmDist_z.argtypes = [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElQuasiTrsmDist_z.restype = c_uint
def QuasiTrsm(side,uplo,orient,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  if type(A) is Matrix:
    if A.tag == sTag: 
      lib.ElQuasiTrsm_s(side,uplo,orient,alpha,A.obj,B.obj)
    elif A.tag == dTag:
      lib.ElQuasiTrsm_d(side,uplo,orient,alpha,A.obj,B.obj)
    elif A.tag == cTag:
      lib.ElQuasiTrsm_c(side,uplo,orient,alpha,A.obj,B.obj)
    elif A.tag == zTag:
      lib.ElQuasiTrsm_z(side,uplo,orient,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if A.tag == sTag: 
      lib.ElQuasiTrsmDist_s(side,uplo,orient,alpha,A.obj,B.obj)
    elif A.tag == dTag:
      lib.ElQuasiTrsmDist_d(side,uplo,orient,alpha,A.obj,B.obj)
    elif A.tag == cTag:
      lib.ElQuasiTrsmDist_c(side,uplo,orient,alpha,A.obj,B.obj)
    elif A.tag == zTag:
      lib.ElQuasiTrsmDist_z(side,uplo,orient,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trdtrmm
# -------
lib.ElTrdtrmm_s.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmm_s.restype = c_uint
lib.ElTrdtrmm_d.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmm_d.restype = c_uint
lib.ElTrdtrmm_c.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmm_c.restype = c_uint
lib.ElTrdtrmm_z.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmm_z.restype = c_uint
lib.ElTrdtrmmDist_s.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmmDist_s.restype = c_uint
lib.ElTrdtrmmDist_d.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmmDist_d.restype = c_uint
lib.ElTrdtrmmDist_c.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmmDist_c.restype = c_uint
lib.ElTrdtrmmDist_z.argtypes = [c_uint,c_void_p]
lib.ElTrdtrmmDist_z.restype = c_uint
def Trdtrmm(uplo,A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrdtrmm_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElTrdtrmm_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElTrdtrmm_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElTrdtrmm_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrdtrmmDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElTrdtrmmDist_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElTrdtrmmDist_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElTrdtrmmDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# TrdtrmmQuasi
# ------------
lib.ElTrdtrmmQuasi_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasi_s.restype = c_uint
lib.ElTrdtrmmQuasi_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasi_d.restype = c_uint
lib.ElTrdtrmmQuasi_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasi_c.restype = c_uint
lib.ElTrdtrmmQuasi_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasi_z.restype = c_uint
lib.ElTrdtrmmQuasiDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasiDist_s.restype = c_uint
lib.ElTrdtrmmQuasiDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasiDist_d.restype = c_uint
lib.ElTrdtrmmQuasiDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasiDist_c.restype = c_uint
lib.ElTrdtrmmQuasiDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElTrdtrmmQuasiDist_z.restype = c_uint
def TrdtrmmQuasi(uplo,A,dOff):
  if type(A) is not type(dOff): 
    raise Exception('Types of A and dOff must match')
  if A.tag != dOff.tag:
    raise Exception('Datatypes of A and dOff must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrdtrmm_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElTrdtrmm_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElTrdtrmm_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElTrdtrmm_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrdtrmmDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElTrdtrmmDist_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElTrdtrmmDist_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElTrdtrmmDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trmm
# ----
lib.ElTrmm_s.argtypes = [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElTrmm_s.restype = c_uint
lib.ElTrmm_d.argtypes = [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]
lib.ElTrmm_d.restype = c_uint
lib.ElTrmm_c.argtypes = [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]
lib.ElTrmm_c.restype = c_uint
lib.ElTrmm_z.argtypes = [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]
lib.ElTrmm_z.restype = c_uint
lib.ElTrmmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElTrmmDist_s.restype = c_uint
lib.ElTrmmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]
lib.ElTrmmDist_d.restype = c_uint
lib.ElTrmmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]
lib.ElTrmmDist_c.restype = c_uint
lib.ElTrmmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]
lib.ElTrmmDist_z.restype = c_uint
def Trmm(side,uplo,orient,diag,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrmm_s(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == dTag: lib.ElTrmm_d(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == cTag: lib.ElTrmm_c(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == zTag: lib.ElTrmm_z(side,uplo,orient,diag,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElTrmmDist_s(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == dTag: 
      lib.ElTrmmDist_d(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == cTag: 
      lib.ElTrmmDist_c(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == zTag: 
      lib.ElTrmmDist_z(side,uplo,orient,diag,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trrk
# ----
lib.ElTrrk_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElTrrk_s.restype = c_uint
lib.ElTrrk_d.argtypes = \
  [c_uint,c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElTrrk_d.restype = c_uint
lib.ElTrrk_c.argtypes = \
  [c_uint,c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElTrrk_c.restype = c_uint
lib.ElTrrk_z.argtypes = \
  [c_uint,c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElTrrk_z.restype = c_uint
lib.ElTrrkDist_s.argtypes = \
  [c_uint,c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElTrrkDist_s.restype = c_uint
lib.ElTrrkDist_d.argtypes = \
  [c_uint,c_uint,c_uint,dType,c_void_p,c_void_p,dType,c_void_p]
lib.ElTrrkDist_d.restype = c_uint
lib.ElTrrkDist_c.argtypes = \
  [c_uint,c_uint,c_uint,cType,c_void_p,c_void_p,cType,c_void_p]
lib.ElTrrkDist_c.restype = c_uint
lib.ElTrrkDist_z.argtypes = \
  [c_uint,c_uint,c_uint,zType,c_void_p,c_void_p,zType,c_void_p]
lib.ElTrrkDist_z.restype = c_uint
def Trrk(uplo,orientA,orientB,alphaPre,A,B,betaPre,C):
  if type(A) is not type(B) or type(B) is not type(C):
    raise Exception('Types of {A,B,C} must match')
  if A.tag != B.tag or B.tag != C.tag:
    raise Exception('Datatypes of {A,B,C} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag:
      lib.ElTrrk_s(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag:
      lib.ElTrrk_d(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag:
      lib.ElTrrk_c(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag:
      lib.ElTrrk_z(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElTrrkDist_s(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag:
      lib.ElTrrkDist_d(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag:
      lib.ElTrrkDist_c(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag:
      lib.ElTrrkDist_z(uplo,orientA,orientB,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trr2k
# -----
#lib.ElTrr2k_s.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   sType,c_void_p,c_void_p,c_void_p,c_void_p,sType,c_void_p]
#lib.ElTrr2k_s.restype = c_uint
#lib.ElTrr2k_d.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   dType,c_void_p,c_void_p,c_void_p,c_void_p,dType,c_void_p]
#lib.ElTrr2k_d.restype = c_uint
#lib.ElTrr2k_c.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   cType,c_void_p,c_void_p,c_void_p,c_void_p,cType,c_void_p]
#lib.ElTrr2k_c.restype = c_uint
#lib.ElTrr2k_z.argtypes = \
#  [c_uint,c_uint,c_uint,c_uint,c_uint,
#   zType,c_void_p,c_void_p,c_void_p,c_void_p,zType,c_void_p]
#lib.ElTrr2k_z.restype = c_uint
lib.ElTrr2kDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   sType,c_void_p,c_void_p,c_void_p,c_void_p,sType,c_void_p]
lib.ElTrr2kDist_s.restype = c_uint
lib.ElTrr2kDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   dType,c_void_p,c_void_p,c_void_p,c_void_p,dType,c_void_p]
lib.ElTrr2kDist_d.restype = c_uint
lib.ElTrr2kDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   cType,c_void_p,c_void_p,c_void_p,c_void_p,cType,c_void_p]
lib.ElTrr2kDist_c.restype = c_uint
lib.ElTrr2kDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,c_uint,
   zType,c_void_p,c_void_p,c_void_p,c_void_p,zType,c_void_p]
lib.ElTrr2kDist_z.restype = c_uint
def Trr2k(uplo,orientA,orientB,orientC,orientD,alphaPre,A,B,C,D,betaPre,E):
  if type(A) is not type(B) or type(B) is not type(C) or \
     type(C) is not type(D) or type(D) is not type(E):
    raise Exception('Types of {A,B,C,D,E} must match')
  if A.tag != B.tag or B.tag != C.tag or C.tag != D.tag or D.tag != E.tag:
    raise Exception('Datatypes of {A,B,C,D,E} must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    raise Exception('Sequential implementation does not yet exist')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElTrr2kDist_s \
      (uplo,orientA,orientB,orientC,orientD,
       alpha,A.obj,B.obj,C.obj,D.obj,beta,E.obj)
    elif A.tag == dTag:
      lib.ElTrr2kDist_d \
      (uplo,orientA,orientB,orientC,orientD,
       alpha,A.obj,B.obj,C.obj,D.obj,beta,E.obj)
    elif A.tag == cTag:
      lib.ElTrr2kDist_c \
      (uplo,orientA,orientB,orientC,orientD,
       alpha,A.obj,B.obj,C.obj,D.obj,beta,E.obj)
    elif A.tag == zTag:
      lib.ElTrr2kDist_z \
      (uplo,orientA,orientB,orientC,orientD,
       alpha,A.obj,B.obj,C.obj,D.obj,beta,E.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trsm
# ----
lib.ElTrsm_s.argtypes = [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElTrsm_s.restype = c_uint
lib.ElTrsm_d.argtypes = [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]
lib.ElTrsm_d.restype = c_uint
lib.ElTrsm_c.argtypes = [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]
lib.ElTrsm_c.restype = c_uint
lib.ElTrsm_z.argtypes = [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]
lib.ElTrsm_z.restype = c_uint
lib.ElTrsmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElTrsmDist_s.restype = c_uint
lib.ElTrsmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]
lib.ElTrsmDist_d.restype = c_uint
lib.ElTrsmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]
lib.ElTrsmDist_c.restype = c_uint
lib.ElTrsmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]
lib.ElTrsmDist_z.restype = c_uint
def Trsm(side,uplo,orient,diag,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrsm_s(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == dTag: lib.ElTrsm_d(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == cTag: lib.ElTrsm_c(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == zTag: lib.ElTrsm_z(side,uplo,orient,diag,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElTrsmDist_s(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == dTag: 
      lib.ElTrsmDist_d(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == cTag: 
      lib.ElTrsmDist_c(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == zTag: 
      lib.ElTrsmDist_z(side,uplo,orient,diag,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trstrm
# ------
lib.ElTrstrm_s.argtypes = [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElTrstrm_s.restype = c_uint
lib.ElTrstrm_d.argtypes = [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]
lib.ElTrstrm_d.restype = c_uint
lib.ElTrstrm_c.argtypes = [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]
lib.ElTrstrm_c.restype = c_uint
lib.ElTrstrm_z.argtypes = [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]
lib.ElTrstrm_z.restype = c_uint
lib.ElTrstrmDist_s.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,sType,c_void_p,c_void_p]
lib.ElTrstrmDist_s.restype = c_uint
lib.ElTrstrmDist_d.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,dType,c_void_p,c_void_p]
lib.ElTrstrmDist_d.restype = c_uint
lib.ElTrstrmDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,cType,c_void_p,c_void_p]
lib.ElTrstrmDist_c.restype = c_uint
lib.ElTrstrmDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_uint,zType,c_void_p,c_void_p]
lib.ElTrstrmDist_z.restype = c_uint
def Trstrm(side,uplo,orient,diag,alphaPre,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  alpha = TagToType(A.tag)(alphaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrstrm_s(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == dTag: lib.ElTrstrm_d(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == cTag: lib.ElTrstrm_c(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == zTag: lib.ElTrstrm_z(side,uplo,orient,diag,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElTrstrmDist_s(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == dTag: 
      lib.ElTrstrmDist_d(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == cTag: 
      lib.ElTrstrmDist_c(side,uplo,orient,diag,alpha,A.obj,B.obj)
    elif A.tag == zTag: 
      lib.ElTrstrmDist_z(side,uplo,orient,diag,alpha,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Trtrmm
# ------
lib.ElTrtrmm_s.argtypes = [c_uint,c_void_p]
lib.ElTrtrmm_s.restype = c_uint
lib.ElTrtrmm_d.argtypes = [c_uint,c_void_p]
lib.ElTrtrmm_d.restype = c_uint
lib.ElTrtrmm_c.argtypes = [c_uint,c_void_p]
lib.ElTrtrmm_c.restype = c_uint
lib.ElTrtrmm_z.argtypes = [c_uint,c_void_p]
lib.ElTrtrmm_z.restype = c_uint
lib.ElTrtrmmDist_s.argtypes = [c_uint,c_void_p]
lib.ElTrtrmmDist_s.restype = c_uint
lib.ElTrtrmmDist_d.argtypes = [c_uint,c_void_p]
lib.ElTrtrmmDist_d.restype = c_uint
lib.ElTrtrmmDist_c.argtypes = [c_uint,c_void_p]
lib.ElTrtrmmDist_c.restype = c_uint
lib.ElTrtrmmDist_z.argtypes = [c_uint,c_void_p]
lib.ElTrtrmmDist_z.restype = c_uint
def Trtrmm(uplo,A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTrtrmm_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElTrtrmm_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElTrtrmm_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElTrtrmm_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTrtrmmDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElTrtrmmDist_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElTrtrmmDist_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElTrtrmmDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Two-sided Trmm
# --------------
lib.ElTwoSidedTrmm_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmm_s.restype = c_uint
lib.ElTwoSidedTrmm_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmm_d.restype = c_uint
lib.ElTwoSidedTrmm_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmm_c.restype = c_uint
lib.ElTwoSidedTrmm_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmm_z.restype = c_uint
lib.ElTwoSidedTrmmDist_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmmDist_s.restype = c_uint
lib.ElTwoSidedTrmmDist_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmmDist_d.restype = c_uint
lib.ElTwoSidedTrmmDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmmDist_c.restype = c_uint
lib.ElTwoSidedTrmmDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrmmDist_z.restype = c_uint
def TwoSidedTrmm(uplo,diag,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoSidedTrmm_s(uplo,diag,A.obj,B.obj)
    elif A.tag == dTag: lib.ElTwoSidedTrmm_d(uplo,diag,A.obj,B.obj)
    elif A.tag == cTag: lib.ElTwoSidedTrmm_c(uplo,diag,A.obj,B.obj)
    elif A.tag == zTag: lib.ElTwoSidedTrmm_z(uplo,diag,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoSidedTrmmDist_s(uplo,diag,A.obj,B.obj)
    elif A.tag == dTag: lib.ElTwoSidedTrmmDist_d(uplo,diag,A.obj,B.obj)
    elif A.tag == cTag: lib.ElTwoSidedTrmmDist_c(uplo,diag,A.obj,B.obj)
    elif A.tag == zTag: lib.ElTwoSidedTrmmDist_z(uplo,diag,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Two-sided Trsm
# --------------
lib.ElTwoSidedTrsm_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsm_s.restype = c_uint
lib.ElTwoSidedTrsm_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsm_d.restype = c_uint
lib.ElTwoSidedTrsm_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsm_c.restype = c_uint
lib.ElTwoSidedTrsm_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsm_z.restype = c_uint
lib.ElTwoSidedTrsmDist_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsmDist_s.restype = c_uint
lib.ElTwoSidedTrsmDist_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsmDist_d.restype = c_uint
lib.ElTwoSidedTrsmDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsmDist_c.restype = c_uint
lib.ElTwoSidedTrsmDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElTwoSidedTrsmDist_z.restype = c_uint
def TwoSidedTrsm(uplo,diag,A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTwoSidedTrsm_s(uplo,diag,A.obj,B.obj)
    elif A.tag == dTag: lib.ElTwoSidedTrsm_d(uplo,diag,A.obj,B.obj)
    elif A.tag == cTag: lib.ElTwoSidedTrsm_c(uplo,diag,A.obj,B.obj)
    elif A.tag == zTag: lib.ElTwoSidedTrsm_z(uplo,diag,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTwoSidedTrsmDist_s(uplo,diag,A.obj,B.obj)
    elif A.tag == dTag: lib.ElTwoSidedTrsmDist_d(uplo,diag,A.obj,B.obj)
    elif A.tag == cTag: lib.ElTwoSidedTrsmDist_c(uplo,diag,A.obj,B.obj)
    elif A.tag == zTag: lib.ElTwoSidedTrsmDist_z(uplo,diag,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
