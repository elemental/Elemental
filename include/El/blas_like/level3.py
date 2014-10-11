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

# Symm
# ----
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
def Symm(side,uplo,alphaPre,A,B,betaPre,C):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymm_s(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSymm_d(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElSymm_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElSymm_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmDist_s(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSymmDist_d(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElSymmDist_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElSymmDist_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Hemm
# ----
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
def Hemm(side,uplo,alphaPre,A,B,betaPre,C):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymm_s(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSymm_d(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElHemm_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElHemm_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmDist_s(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSymmDist_d(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElHemmDist_c(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElHemmDist_z(side,uplo,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Syrk
# ----
lib.ElSyrk_s.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElSyrk_s.restype = c_uint
lib.ElSyrk_d.argtypes = [c_uint,c_uint,sType,c_void_p,dType,c_void_p]
lib.ElSyrk_d.restype = c_uint
lib.ElSyrk_c.argtypes = [c_uint,c_uint,sType,c_void_p,cType,c_void_p]
lib.ElSyrk_c.restype = c_uint
lib.ElSyrk_z.argtypes = [c_uint,c_uint,sType,c_void_p,zType,c_void_p]
lib.ElSyrk_z.restype = c_uint
lib.ElSyrkDist_s.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElSyrkDist_s.restype = c_uint
lib.ElSyrkDist_d.argtypes = [c_uint,c_uint,sType,c_void_p,dType,c_void_p]
lib.ElSyrkDist_d.restype = c_uint
lib.ElSyrkDist_c.argtypes = [c_uint,c_uint,sType,c_void_p,cType,c_void_p]
lib.ElSyrkDist_c.restype = c_uint
lib.ElSyrkDist_z.argtypes = [c_uint,c_uint,sType,c_void_p,zType,c_void_p]
lib.ElSyrkDist_z.restype = c_uint
def Syrk(uplo,orient,alphaPre,A,betaPre,C):
  if A.tag != C.tag: raise Exception('Datatypes of A and C must match')
  if type(A) is not type(C): raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyrk_s(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyrk_d(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElSyrk_c(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElSyrk_z(uplo,orient,alpha,A.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSyrkDist_s(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyrkDist_d(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElSyrkDist_c(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElSyrkDist_z(uplo,orient,alpha,A.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Herk
# ----
lib.ElHerk_c.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElHerk_c.restype = c_uint
lib.ElHerk_z.argtypes = [c_uint,c_uint,sType,c_void_p,dType,c_void_p]
lib.ElHerk_z.restype = c_uint
lib.ElHerkDist_c.argtypes = [c_uint,c_uint,sType,c_void_p,sType,c_void_p]
lib.ElHerkDist_c.restype = c_uint
lib.ElHerkDist_z.argtypes = [c_uint,c_uint,sType,c_void_p,dType,c_void_p]
lib.ElHerkDist_z.restype = c_uint
def Herk(uplo,orient,alpha,A,beta,C):
  if A.tag != C.tag: raise Exception('Datatypes of A and C must match')
  if type(A) is not type(C): raise Exception('Matrix types must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyrk_s(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyrk_d(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElHerk_c(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElHerk_z(uplo,orient,alpha,A.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSyrkDist_s(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyrkDist_d(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElHerkDist_c(uplo,orient,alpha,A.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElHerkDist_z(uplo,orient,alpha,A.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Syr2k
# -----
lib.ElSyr2k_s.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSyr2k_s.restype = c_uint
lib.ElSyr2k_d.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSyr2k_d.restype = c_uint
lib.ElSyr2k_c.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSyr2k_c.restype = c_uint
lib.ElSyr2k_z.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSyr2k_z.restype = c_uint
lib.ElSyr2kDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElSyr2kDist_s.restype = c_uint
lib.ElSyr2kDist_d.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,dType,c_void_p]
lib.ElSyr2kDist_d.restype = c_uint
lib.ElSyr2kDist_c.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,cType,c_void_p]
lib.ElSyr2kDist_c.restype = c_uint
lib.ElSyr2kDist_z.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,zType,c_void_p]
lib.ElSyr2kDist_z.restype = c_uint
def Syr2k(uplo,orient,alphaPre,A,B,betaPre,C):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  alpha = TagToType(A.tag)(alphaPre)
  beta = TagToType(A.tag)(betaPre)
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyr2k_s(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyr2k_d(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElSyr2k_c(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElSyr2k_z(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElSyr2kDist_s(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: 
      lib.ElSyr2kDist_d(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: 
      lib.ElSyr2kDist_c(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: 
      lib.ElSyr2kDist_z(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Her2k
# ----
lib.ElHer2k_c.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHer2k_c.restype = c_uint
lib.ElHer2k_z.argtypes = [c_uint,c_uint,sType,c_void_p,c_void_p,dType,c_void_p]
lib.ElHer2k_z.restype = c_uint
lib.ElHer2kDist_c.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,sType,c_void_p]
lib.ElHer2kDist_c.restype = c_uint
lib.ElHer2kDist_z.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,dType,c_void_p]
lib.ElHer2kDist_z.restype = c_uint
def Her2k(uplo,orient,alpha,A,B,beta,C):
  if A.tag != B.tag or B.tag != C.tag: 
    raise Exception('Datatypes of {A,B,C} must match')
  if type(A) is not type(B) or type(B) is not type(C): 
    raise Exception('Matrix types must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSyr2k_s(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: lib.ElSyr2k_d(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: lib.ElHer2k_c(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: lib.ElHer2k_z(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: 
      lib.ElSyr2kDist_s(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == dTag: 
      lib.ElSyr2kDist_d(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == cTag: 
      lib.ElHer2kDist_c(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    elif A.tag == zTag: 
      lib.ElHer2kDist_z(uplo,orient,alpha,A.obj,B.obj,beta,C.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
