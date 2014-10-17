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
from ctypes import CFUNCTYPE

# Hermitian function
# ==================
lib.ElRealHermitianFunction_s.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(sType,sType)]
lib.ElRealHermitianFunction_s.restype = c_uint
lib.ElRealHermitianFunction_d.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]
lib.ElRealHermitianFunction_d.restype = c_uint
lib.ElRealHermitianFunction_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(sType,sType)]
lib.ElRealHermitianFunction_c.restype = c_uint
lib.ElRealHermitianFunction_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]
lib.ElRealHermitianFunction_z.restype = c_uint
lib.ElRealHermitianFunctionDist_s.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(sType,sType)]
lib.ElRealHermitianFunctionDist_s.restype = c_uint
lib.ElRealHermitianFunctionDist_d.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]
lib.ElRealHermitianFunctionDist_d.restype = c_uint
lib.ElRealHermitianFunctionDist_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(sType,sType)]
lib.ElRealHermitianFunctionDist_c.restype = c_uint
lib.ElRealHermitianFunctionDist_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]
lib.ElRealHermitianFunctionDist_z.restype = c_uint

lib.ElComplexHermitianFunction_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(cType,sType)]
lib.ElComplexHermitianFunction_c.restype = c_uint
lib.ElComplexHermitianFunction_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(zType,dType)]
lib.ElComplexHermitianFunction_z.restype = c_uint
lib.ElComplexHermitianFunctionDist_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(cType,sType)]
lib.ElComplexHermitianFunctionDist_c.restype = c_uint
lib.ElComplexHermitianFunctionDist_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(zType,dType)]
lib.ElComplexHermitianFunctionDist_z.restype = c_uint
def HermitianFunction(uplo,A,func,isReal=True):
  if type(A) is Matrix:
    if isReal:
      cFunc = CFUNCTYPE(TagToType(Base(A.tag)),TagToType(Base(A.tag)))(func)
      if   A.tag == sTag: lib.ElRealHermitianFunction_s(uplo,A.obj,cFunc)
      elif A.tag == dTag: lib.ElRealHermitianFunction_d(uplo,A.obj,cFunc)
      elif A.tag == cTag: lib.ElRealHermitianFunction_c(uplo,A.obj,cFunc)
      elif A.tag == zTag: lib.ElRealHermitianFunction_z(uplo,A.obj,cFunc)
      else: raise Exception('Unsupported datatype')
    else:
      cFunc = CFUNCTYPE(TagToType(A.tag),TagToType(Base(A.tag)))(func)
      if   A.tag == cTag: lib.ElComplexHermitianFunction_c(uplo,A.obj,cFunc)
      elif A.tag == zTag: lib.ElComplexHermitianFunction_z(uplo,A.obj,cFunc)
      else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if isReal:
      cFunc = CFUNCTYPE(TagToType(Base(A.tag)),TagToType(Base(A.tag)))(func)
      if   A.tag == sTag: lib.ElRealHermitianFunctionDist_s(uplo,A.obj,cFunc)
      elif A.tag == dTag: lib.ElRealHermitianFunctionDist_d(uplo,A.obj,cFunc)
      elif A.tag == cTag: lib.ElRealHermitianFunctionDist_c(uplo,A.obj,cFunc)
      elif A.tag == zTag: lib.ElRealHermitianFunctionDist_z(uplo,A.obj,cFunc)
      else: raise Exception('Unsupported datatype')
    else:
      cFunc = CFUNCTYPE(TagToType(A.tag),TagToType(Base(A.tag)))(func)
      if   A.tag == cTag: lib.ElComplexHermitianFunctionDist_c(uplo,A.obj,cFunc)
      elif A.tag == zTag: lib.ElComplexHermitianFunctionDist_z(uplo,A.obj,cFunc)
      else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Inverse
# =======

# General
# -------
lib.ElInverse_s.argtypes = [c_void_p]
lib.ElInverse_s.restype = c_uint
lib.ElInverse_d.argtypes = [c_void_p]
lib.ElInverse_d.restype = c_uint
lib.ElInverse_c.argtypes = [c_void_p]
lib.ElInverse_c.restype = c_uint
lib.ElInverse_z.argtypes = [c_void_p]
lib.ElInverse_z.restype = c_uint
lib.ElInverseDist_s.argtypes = [c_void_p]
lib.ElInverseDist_s.restype = c_uint
lib.ElInverseDist_d.argtypes = [c_void_p]
lib.ElInverseDist_d.restype = c_uint
lib.ElInverseDist_c.argtypes = [c_void_p]
lib.ElInverseDist_c.restype = c_uint
lib.ElInverseDist_z.argtypes = [c_void_p]
lib.ElInverseDist_z.restype = c_uint
def Inverse(A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElInverse_s(A.obj)
    elif A.tag == dTag: lib.ElInverse_d(A.obj)
    elif A.tag == cTag: lib.ElInverse_c(A.obj)
    elif A.tag == zTag: lib.ElInverse_z(A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElInverseDist_s(A.obj)
    elif A.tag == dTag: lib.ElInverseDist_d(A.obj)
    elif A.tag == cTag: lib.ElInverseDist_c(A.obj)
    elif A.tag == zTag: lib.ElInverseDist_z(A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# After LU factorization with partial pivoting
# --------------------------------------------
lib.ElInverseAfterLUPartialPiv_s.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPiv_s.restype = c_uint
lib.ElInverseAfterLUPartialPiv_d.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPiv_d.restype = c_uint
lib.ElInverseAfterLUPartialPiv_c.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPiv_c.restype = c_uint
lib.ElInverseAfterLUPartialPiv_z.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPiv_z.restype = c_uint
lib.ElInverseAfterLUPartialPivDist_s.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPivDist_s.restype = c_uint
lib.ElInverseAfterLUPartialPivDist_d.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPivDist_d.restype = c_uint
lib.ElInverseAfterLUPartialPivDist_c.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPivDist_c.restype = c_uint
lib.ElInverseAfterLUPartialPivDist_z.argtypes = [c_void_p,c_void_p]
lib.ElInverseAfterLUPartialPivDist_z.restype = c_uint
def InverseAfterLUPartialPiv(A,p):
  if type(A) is not type(p):
    raise Exception('Types of A and p must match')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElInverseAfterLUPartialPiv_s(A.obj,p.obj)
    elif A.tag == dTag: lib.ElInverseAfterLUPartialPiv_d(A.obj,p.obj)
    elif A.tag == cTag: lib.ElInverseAfterLUPartialPiv_c(A.obj,p.obj)
    elif A.tag == zTag: lib.ElInverseAfterLUPartialPiv_z(A.obj,p.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElInverseAfterLUPartialPivDist_s(A.obj,p.obj)
    elif A.tag == dTag: lib.ElInverseAfterLUPartialPivDist_d(A.obj,p.obj)
    elif A.tag == cTag: lib.ElInverseAfterLUPartialPivDist_c(A.obj,p.obj)
    elif A.tag == zTag: lib.ElInverseAfterLUPartialPivDist_z(A.obj,p.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Invert a Hermitian Positive-Definite matrix
# -------------------------------------------
lib.ElHPDInverse_s.argtypes = [c_uint,c_void_p]
lib.ElHPDInverse_s.restype = c_uint
lib.ElHPDInverse_d.argtypes = [c_uint,c_void_p]
lib.ElHPDInverse_d.restype = c_uint
lib.ElHPDInverse_c.argtypes = [c_uint,c_void_p]
lib.ElHPDInverse_c.restype = c_uint
lib.ElHPDInverse_z.argtypes = [c_uint,c_void_p]
lib.ElHPDInverse_z.restype = c_uint
lib.ElHPDInverseDist_s.argtypes = [c_uint,c_void_p]
lib.ElHPDInverseDist_s.restype = c_uint
lib.ElHPDInverseDist_d.argtypes = [c_uint,c_void_p]
lib.ElHPDInverseDist_d.restype = c_uint
lib.ElHPDInverseDist_c.argtypes = [c_uint,c_void_p]
lib.ElHPDInverseDist_c.restype = c_uint
lib.ElHPDInverseDist_z.argtypes = [c_uint,c_void_p]
lib.ElHPDInverseDist_z.restype = c_uint
def HPDInverse(uplo,A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPDInverse_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHPDInverse_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHPDInverse_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHPDInverse_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPDInverseDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHPDInverseDist_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHPDInverseDist_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHPDInverseDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Invert a symmetric/Hermitian matrix
# -----------------------------------
lib.ElSymmetricInverse_s.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverse_s.restype = c_uint
lib.ElSymmetricInverse_d.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverse_d.restype = c_uint
lib.ElSymmetricInverse_c.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverse_c.restype = c_uint
lib.ElSymmetricInverse_z.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverse_z.restype = c_uint
lib.ElSymmetricInverseDist_s.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverseDist_s.restype = c_uint
lib.ElSymmetricInverseDist_d.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverseDist_d.restype = c_uint
lib.ElSymmetricInverseDist_c.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverseDist_c.restype = c_uint
lib.ElSymmetricInverseDist_z.argtypes = [c_uint,c_void_p]
lib.ElSymmetricInverseDist_z.restype = c_uint
lib.ElHermitianInverse_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianInverse_c.restype = c_uint
lib.ElHermitianInverse_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianInverse_z.restype = c_uint
lib.ElHermitianInverseDist_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianInverseDist_c.restype = c_uint
lib.ElHermitianInverseDist_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianInverseDist_z.restype = c_uint
def SymmetricInverse(uplo,A,conjugate=False):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricInverse_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElSymmetricInverse_d(uplo,A.obj)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianInverse_c(uplo,A.obj)
      else:         lib.ElSymmetricInverse_c(uplo,A.obj)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianInverse_z(uplo,A.obj)
      else:         lib.ElSymmetricInverse_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricInverseDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElSymmetricInverseDist_d(uplo,A.obj)
    elif A.tag == cTag: 
      if conjugate: lib.ElHermitianInverseDist_c(uplo,A.obj)
      else:         lib.ElSymmetricInverseDist_c(uplo,A.obj)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianInverseDist_z(uplo,A.obj)
      else:         lib.ElSymmetricInverseDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
def HermitianInverse(uplo,A):
  SymmetricInverse(uplo,A,True)

# Triangular
# ----------
lib.ElTriangularInverse_s.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverse_s.restype = c_uint
lib.ElTriangularInverse_d.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverse_d.restype = c_uint
lib.ElTriangularInverse_c.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverse_c.restype = c_uint
lib.ElTriangularInverse_z.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverse_z.restype = c_uint
lib.ElTriangularInverseDist_s.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverseDist_s.restype = c_uint
lib.ElTriangularInverseDist_d.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverseDist_d.restype = c_uint
lib.ElTriangularInverseDist_c.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverseDist_c.restype = c_uint
lib.ElTriangularInverseDist_z.argtypes = [c_uint,c_uint,c_void_p]
lib.ElTriangularInverseDist_z.restype = c_uint
def TriangularInverse(uplo,diag,A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTriangularInverse_s(uplo,diag,A.obj)
    elif A.tag == dTag: lib.ElTriangularInverse_d(uplo,diag,A.obj)
    elif A.tag == cTag: lib.ElTriangularInverse_c(uplo,diag,A.obj)
    elif A.tag == zTag: lib.ElTriangularInverse_z(uplo,diag,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTriangularInverseDist_s(uplo,diag,A.obj)
    elif A.tag == dTag: lib.ElTriangularInverseDist_d(uplo,diag,A.obj)
    elif A.tag == cTag: lib.ElTriangularInverseDist_c(uplo,diag,A.obj)
    elif A.tag == zTag: lib.ElTriangularInverseDist_z(uplo,diag,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Pseudoinverse
# =============
lib.ElPseudoinverse_s.argtypes = [c_void_p]
lib.ElPseudoinverse_s.restype = c_uint
lib.ElPseudoinverse_d.argtypes = [c_void_p]
lib.ElPseudoinverse_d.restype = c_uint
lib.ElPseudoinverse_c.argtypes = [c_void_p]
lib.ElPseudoinverse_c.restype = c_uint
lib.ElPseudoinverse_z.argtypes = [c_void_p]
lib.ElPseudoinverse_z.restype = c_uint
lib.ElPseudoinverseDist_s.argtypes = [c_void_p]
lib.ElPseudoinverseDist_s.restype = c_uint
lib.ElPseudoinverseDist_d.argtypes = [c_void_p]
lib.ElPseudoinverseDist_d.restype = c_uint
lib.ElPseudoinverseDist_c.argtypes = [c_void_p]
lib.ElPseudoinverseDist_c.restype = c_uint
lib.ElPseudoinverseDist_z.argtypes = [c_void_p]
lib.ElPseudoinverseDist_z.restype = c_uint
def Pseudoinverse(A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElPseudoinverse_s(A.obj)
    elif A.tag == dTag: lib.ElPseudoinverse_d(A.obj)
    elif A.tag == cTag: lib.ElPseudoinverse_c(A.obj)
    elif A.tag == zTag: lib.ElPseudoinverse_z(A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElPseudoinverseDist_s(A.obj)
    elif A.tag == dTag: lib.ElPseudoinverseDist_d(A.obj)
    elif A.tag == cTag: lib.ElPseudoinverseDist_c(A.obj)
    elif A.tag == zTag: lib.ElPseudoinverseDist_z(A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHermitianPseudoinverse_s.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverse_s.restype = c_uint
lib.ElHermitianPseudoinverse_d.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverse_d.restype = c_uint
lib.ElHermitianPseudoinverse_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverse_c.restype = c_uint
lib.ElHermitianPseudoinverse_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverse_z.restype = c_uint
lib.ElHermitianPseudoinverseDist_s.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverseDist_s.restype = c_uint
lib.ElHermitianPseudoinverseDist_d.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverseDist_d.restype = c_uint
lib.ElHermitianPseudoinverseDist_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverseDist_c.restype = c_uint
lib.ElHermitianPseudoinverseDist_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianPseudoinverseDist_z.restype = c_uint
def HermitianPseudoinverse(uplo,A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHermitianPseudoinverse_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHermitianPseudoinverse_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHermitianPseudoinverse_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHermitianPseudoinverse_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHermitianPseudoinverseDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHermitianPseudoinverseDist_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHermitianPseudoinverseDist_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHermitianPseudoinverseDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Sign
# ====
lib.ElSign_s.argtypes = [c_void_p]
lib.ElSign_s.restype = c_uint
lib.ElSign_d.argtypes = [c_void_p]
lib.ElSign_d.restype = c_uint
lib.ElSign_c.argtypes = [c_void_p]
lib.ElSign_c.restype = c_uint
lib.ElSign_z.argtypes = [c_void_p]
lib.ElSign_z.restype = c_uint
lib.ElSignDecomp_s.argtypes = [c_void_p,c_void_p]
lib.ElSignDecomp_s.restype = c_uint
lib.ElSignDecomp_d.argtypes = [c_void_p,c_void_p]
lib.ElSignDecomp_d.restype = c_uint
lib.ElSignDecomp_c.argtypes = [c_void_p,c_void_p]
lib.ElSignDecomp_c.restype = c_uint
lib.ElSignDecomp_z.argtypes = [c_void_p,c_void_p]
lib.ElSignDecomp_z.restype = c_uint
lib.ElSignDist_s.argtypes = [c_void_p]
lib.ElSignDist_s.restype = c_uint
lib.ElSignDist_d.argtypes = [c_void_p]
lib.ElSignDist_d.restype = c_uint
lib.ElSignDist_c.argtypes = [c_void_p]
lib.ElSignDist_c.restype = c_uint
lib.ElSignDist_z.argtypes = [c_void_p]
lib.ElSignDist_z.restype = c_uint
lib.ElSignDecompDist_s.argtypes = [c_void_p,c_void_p]
lib.ElSignDecompDist_s.restype = c_uint
lib.ElSignDecompDist_d.argtypes = [c_void_p,c_void_p]
lib.ElSignDecompDist_d.restype = c_uint
lib.ElSignDecompDist_c.argtypes = [c_void_p,c_void_p]
lib.ElSignDecompDist_c.restype = c_uint
lib.ElSignDecompDist_z.argtypes = [c_void_p,c_void_p]
lib.ElSignDecompDist_z.restype = c_uint
def Sign(A,fullDecomp=False):
  if type(A) is Matrix:
    if fullDecomp:
      N = Matrix(A.tag)
      if   A.tag == sTag: lib.ElSignDecomp_s(A.obj,N.obj)
      elif A.tag == dTag: lib.ElSignDecomp_d(A.obj,N.obj)
      elif A.tag == cTag: lib.ElSignDecomp_c(A.obj,N.obj)
      elif A.tag == zTag: lib.ElSignDecomp_z(A.obj,N.obj)
      else: raise Exception('Unsupported datatype')
      return N
    else:
      if   A.tag == sTag: lib.ElSign_s(A.obj)
      elif A.tag == dTag: lib.ElSign_d(A.obj)
      elif A.tag == cTag: lib.ElSign_c(A.obj)
      elif A.tag == zTag: lib.ElSign_z(A.obj)
      else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if fullDecomp:
      N = Matrix(A.tag)
      if   A.tag == sTag: lib.ElSignDecompDist_s(A.obj,N.obj)
      elif A.tag == dTag: lib.ElSignDecompDist_d(A.obj,N.obj)
      elif A.tag == cTag: lib.ElSignDecompDist_c(A.obj,N.obj)
      elif A.tag == zTag: lib.ElSignDecompDist_z(A.obj,N.obj)
      else: raise Exception('Unsupported datatype')
      return N
    else:
      if   A.tag == sTag: lib.ElSignDist_s(A.obj)
      elif A.tag == dTag: lib.ElSignDist_d(A.obj)
      elif A.tag == cTag: lib.ElSignDist_c(A.obj)
      elif A.tag == zTag: lib.ElSignDist_z(A.obj)
      else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHermitianSign_s.argtypes = [c_uint,c_void_p]
lib.ElHermitianSign_s.restype = c_uint
lib.ElHermitianSign_d.argtypes = [c_uint,c_void_p]
lib.ElHermitianSign_d.restype = c_uint
lib.ElHermitianSign_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianSign_c.restype = c_uint
lib.ElHermitianSign_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianSign_z.restype = c_uint
lib.ElHermitianSignDecomp_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecomp_s.restype = c_uint
lib.ElHermitianSignDecomp_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecomp_d.restype = c_uint
lib.ElHermitianSignDecomp_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecomp_c.restype = c_uint
lib.ElHermitianSignDecomp_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecomp_z.restype = c_uint
lib.ElHermitianSignDist_s.argtypes = [c_uint,c_void_p]
lib.ElHermitianSignDist_s.restype = c_uint
lib.ElHermitianSignDist_d.argtypes = [c_uint,c_void_p]
lib.ElHermitianSignDist_d.restype = c_uint
lib.ElHermitianSignDist_c.argtypes = [c_uint,c_void_p]
lib.ElHermitianSignDist_c.restype = c_uint
lib.ElHermitianSignDist_z.argtypes = [c_uint,c_void_p]
lib.ElHermitianSignDist_z.restype = c_uint
lib.ElHermitianSignDecompDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecompDist_s.restype = c_uint
lib.ElHermitianSignDecompDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecompDist_d.restype = c_uint
lib.ElHermitianSignDecompDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecompDist_c.restype = c_uint
lib.ElHermitianSignDecompDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElHermitianSignDecompDist_z.restype = c_uint
def HermitianSign(uplo,A,fullDecomp=False):
  if type(A) is Matrix:
    if fullDecomp:
      N = Matrix(A.tag)
      if   A.tag == sTag: lib.ElHermitianSignDecomp_s(uplo,A.obj,N.obj)
      elif A.tag == dTag: lib.ElHermitianSignDecomp_d(uplo,A.obj,N.obj)
      elif A.tag == cTag: lib.ElHermitianSignDecomp_c(uplo,A.obj,N.obj)
      elif A.tag == zTag: lib.ElHermitianSignDecomp_z(uplo,A.obj,N.obj)
      else: raise Exception('Unsupported datatype')
      return N
    else:
      if   A.tag == sTag: lib.ElHermitianSign_s(uplo,A.obj)
      elif A.tag == dTag: lib.ElHermitianSign_d(uplo,A.obj)
      elif A.tag == cTag: lib.ElHermitianSign_c(uplo,A.obj)
      elif A.tag == zTag: lib.ElHermitianSign_z(uplo,A.obj)
      else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if fullDecomp:
      N = Matrix(A.tag)
      if   A.tag == sTag: lib.ElHermitianSignDecompDist_s(uplo,A.obj,N.obj)
      elif A.tag == dTag: lib.ElHermitianSignDecompDist_d(uplo,A.obj,N.obj)
      elif A.tag == cTag: lib.ElHermitianSignDecompDist_c(uplo,A.obj,N.obj)
      elif A.tag == zTag: lib.ElHermitianSignDecompDist_z(uplo,A.obj,N.obj)
      else: raise Exception('Unsupported datatype')
      return N
    else:
      if   A.tag == sTag: lib.ElHermitianSignDist_s(uplo,A.obj)
      elif A.tag == dTag: lib.ElHermitianSignDist_d(uplo,A.obj)
      elif A.tag == cTag: lib.ElHermitianSignDist_c(uplo,A.obj)
      elif A.tag == zTag: lib.ElHermitianSignDist_z(uplo,A.obj)
      else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Square-root
# ===========
lib.ElSquareRoot_s.argtypes = [c_void_p]
lib.ElSquareRoot_s.restype = c_uint
lib.ElSquareRoot_d.argtypes = [c_void_p]
lib.ElSquareRoot_d.restype = c_uint
lib.ElSquareRoot_c.argtypes = [c_void_p]
lib.ElSquareRoot_c.restype = c_uint
lib.ElSquareRoot_z.argtypes = [c_void_p]
lib.ElSquareRoot_z.restype = c_uint
lib.ElSquareRootDist_s.argtypes = [c_void_p]
lib.ElSquareRootDist_s.restype = c_uint
lib.ElSquareRootDist_d.argtypes = [c_void_p]
lib.ElSquareRootDist_d.restype = c_uint
lib.ElSquareRootDist_c.argtypes = [c_void_p]
lib.ElSquareRootDist_c.restype = c_uint
lib.ElSquareRootDist_z.argtypes = [c_void_p]
lib.ElSquareRootDist_z.restype = c_uint
def SquareRoot(A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSquareRoot_s(A.obj)
    elif A.tag == dTag: lib.ElSquareRoot_d(A.obj)
    elif A.tag == cTag: lib.ElSquareRoot_c(A.obj)
    elif A.tag == zTag: lib.ElSquareRoot_z(A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSquareRootDist_s(A.obj)
    elif A.tag == dTag: lib.ElSquareRootDist_d(A.obj)
    elif A.tag == cTag: lib.ElSquareRootDist_c(A.obj)
    elif A.tag == zTag: lib.ElSquareRootDist_z(A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElHPSDSquareRoot_s.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRoot_s.restype = c_uint
lib.ElHPSDSquareRoot_d.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRoot_d.restype = c_uint
lib.ElHPSDSquareRoot_c.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRoot_c.restype = c_uint
lib.ElHPSDSquareRoot_z.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRoot_z.restype = c_uint
lib.ElHPSDSquareRootDist_s.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRootDist_s.restype = c_uint
lib.ElHPSDSquareRootDist_d.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRootDist_d.restype = c_uint
lib.ElHPSDSquareRootDist_c.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRootDist_c.restype = c_uint
lib.ElHPSDSquareRootDist_z.argtypes = [c_uint,c_void_p]
lib.ElHPSDSquareRootDist_z.restype = c_uint
def HPSDSquareRoot(uplo,A):
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPSDSquareRoot_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHPSDSquareRoot_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHPSDSquareRoot_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHPSDSquareRoot_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPSDSquareRootDist_s(uplo,A.obj)
    elif A.tag == dTag: lib.ElHPSDSquareRootDist_d(uplo,A.obj)
    elif A.tag == cTag: lib.ElHPSDSquareRootDist_c(uplo,A.obj)
    elif A.tag == zTag: lib.ElHPSDSquareRootDist_z(uplo,A.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')
