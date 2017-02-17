#
#  Copyright (c) 2009-2016, Jack Poulson
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
lib.ElRealHermitianFunction_d.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]
lib.ElRealHermitianFunction_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(sType,sType)]
lib.ElRealHermitianFunction_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]
lib.ElRealHermitianFunctionDist_s.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(sType,sType)]
lib.ElRealHermitianFunctionDist_d.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]
lib.ElRealHermitianFunctionDist_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(sType,sType)]
lib.ElRealHermitianFunctionDist_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(dType,dType)]

lib.ElComplexHermitianFunction_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(cType,sType)]
lib.ElComplexHermitianFunction_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(zType,dType)]
lib.ElComplexHermitianFunctionDist_c.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(cType,sType)]
lib.ElComplexHermitianFunctionDist_z.argtypes = \
  [c_uint,c_void_p,CFUNCTYPE(zType,dType)]

def HermitianFunction(uplo,A,func,isReal=True):
  if isReal:
    cFunc = CFUNCTYPE(TagToType(Base(A.tag)),TagToType(Base(A.tag)))(func)
  else:
    cFunc = CFUNCTYPE(TagToType(A.tag),TagToType(Base(A.tag)))(func)
  args = [uplo,A.obj,cFunc]
  if type(A) is Matrix:
    if isReal:
      if   A.tag == sTag: lib.ElRealHermitianFunction_s(*args)
      elif A.tag == dTag: lib.ElRealHermitianFunction_d(*args)
      elif A.tag == cTag: lib.ElRealHermitianFunction_c(*args)
      elif A.tag == zTag: lib.ElRealHermitianFunction_z(*args)
      else: DataExcept()
    else:
      if   A.tag == cTag: lib.ElComplexHermitianFunction_c(*args)
      elif A.tag == zTag: lib.ElComplexHermitianFunction_z(*args)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if isReal:
      if   A.tag == sTag: lib.ElRealHermitianFunctionDist_s(*args)
      elif A.tag == dTag: lib.ElRealHermitianFunctionDist_d(*args)
      elif A.tag == cTag: lib.ElRealHermitianFunctionDist_c(*args)
      elif A.tag == zTag: lib.ElRealHermitianFunctionDist_z(*args)
      else: DataExcept()
    else:
      if   A.tag == cTag: lib.ElComplexHermitianFunctionDist_c(*args)
      elif A.tag == zTag: lib.ElComplexHermitianFunctionDist_z(*args)
      else: DataExcept()
  else: TypeExcept()

# Inverse
# =======

# General
# -------
lib.ElInverse_s.argtypes = \
lib.ElInverse_d.argtypes = \
lib.ElInverse_c.argtypes = \
lib.ElInverse_z.argtypes = \
lib.ElInverseDist_s.argtypes = \
lib.ElInverseDist_d.argtypes = \
lib.ElInverseDist_c.argtypes = \
lib.ElInverseDist_z.argtypes = \
  [c_void_p]
def Inverse(A):
  args = [A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElInverse_s(*args)
    elif A.tag == dTag: lib.ElInverse_d(*args)
    elif A.tag == cTag: lib.ElInverse_c(*args)
    elif A.tag == zTag: lib.ElInverse_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElInverseDist_s(*args)
    elif A.tag == dTag: lib.ElInverseDist_d(*args)
    elif A.tag == cTag: lib.ElInverseDist_c(*args)
    elif A.tag == zTag: lib.ElInverseDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# After LU factorization with partial pivoting
# --------------------------------------------
lib.ElInverseAfterLUPartialPiv_s.argtypes = \
lib.ElInverseAfterLUPartialPiv_d.argtypes = \
lib.ElInverseAfterLUPartialPiv_c.argtypes = \
lib.ElInverseAfterLUPartialPiv_z.argtypes = \
lib.ElInverseAfterLUPartialPivDist_s.argtypes = \
lib.ElInverseAfterLUPartialPivDist_d.argtypes = \
lib.ElInverseAfterLUPartialPivDist_c.argtypes = \
lib.ElInverseAfterLUPartialPivDist_z.argtypes = \
  [c_void_p,c_void_p]
def InverseAfterLUPartialPiv(A,P):
  args = [A.obj,P.obj]
  if type(A) is Matrix:
    if type(P) is not Permutation:
      raise Exception('Expected P to be a Permutation')
    if   A.tag == sTag: lib.ElInverseAfterLUPartialPiv_s(*args)
    elif A.tag == dTag: lib.ElInverseAfterLUPartialPiv_d(*args)
    elif A.tag == cTag: lib.ElInverseAfterLUPartialPiv_c(*args)
    elif A.tag == zTag: lib.ElInverseAfterLUPartialPiv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(P) is not DistPermutation:
      raise Exception('Expected P to be a DistPermutation')
    if   A.tag == sTag: lib.ElInverseAfterLUPartialPivDist_s(*args)
    elif A.tag == dTag: lib.ElInverseAfterLUPartialPivDist_d(*args)
    elif A.tag == cTag: lib.ElInverseAfterLUPartialPivDist_c(*args)
    elif A.tag == zTag: lib.ElInverseAfterLUPartialPivDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Invert a Hermitian Positive-Definite matrix
# -------------------------------------------
lib.ElHPDInverse_s.argtypes = \
lib.ElHPDInverse_d.argtypes = \
lib.ElHPDInverse_c.argtypes = \
lib.ElHPDInverse_z.argtypes = \
lib.ElHPDInverseDist_s.argtypes = \
lib.ElHPDInverseDist_d.argtypes = \
lib.ElHPDInverseDist_c.argtypes = \
lib.ElHPDInverseDist_z.argtypes = \
  [c_uint,c_void_p]
def HPDInverse(uplo,A):
  args = [uplo,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPDInverse_s(*args)
    elif A.tag == dTag: lib.ElHPDInverse_d(*args)
    elif A.tag == cTag: lib.ElHPDInverse_c(*args)
    elif A.tag == zTag: lib.ElHPDInverse_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPDInverseDist_s(*args)
    elif A.tag == dTag: lib.ElHPDInverseDist_d(*args)
    elif A.tag == cTag: lib.ElHPDInverseDist_c(*args)
    elif A.tag == zTag: lib.ElHPDInverseDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Invert a symmetric/Hermitian matrix
# -----------------------------------
lib.ElSymmetricInverse_s.argtypes = \
lib.ElSymmetricInverse_d.argtypes = \
lib.ElSymmetricInverse_c.argtypes = \
lib.ElSymmetricInverse_z.argtypes = \
lib.ElSymmetricInverseDist_s.argtypes = \
lib.ElSymmetricInverseDist_d.argtypes = \
lib.ElSymmetricInverseDist_c.argtypes = \
lib.ElSymmetricInverseDist_z.argtypes = \
lib.ElHermitianInverse_c.argtypes = \
lib.ElHermitianInverse_z.argtypes = \
lib.ElHermitianInverseDist_c.argtypes = \
lib.ElHermitianInverseDist_z.argtypes = \
  [c_uint,c_void_p]
def SymmetricInverse(uplo,A,conjugate=False):
  args = [uplo,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricInverse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricInverse_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianInverse_c(*args)
      else:         lib.ElSymmetricInverse_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianInverse_z(*args)
      else:         lib.ElSymmetricInverse_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricInverseDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricInverseDist_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianInverseDist_c(*args)
      else:         lib.ElSymmetricInverseDist_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianInverseDist_z(*args)
      else:         lib.ElSymmetricInverseDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
def HermitianInverse(uplo,A):
  SymmetricInverse(uplo,A,True)

# Triangular
# ----------
lib.ElTriangularInverse_s.argtypes = \
lib.ElTriangularInverse_d.argtypes = \
lib.ElTriangularInverse_c.argtypes = \
lib.ElTriangularInverse_z.argtypes = \
lib.ElTriangularInverseDist_s.argtypes = \
lib.ElTriangularInverseDist_d.argtypes = \
lib.ElTriangularInverseDist_c.argtypes = \
lib.ElTriangularInverseDist_z.argtypes = \
  [c_uint,c_uint,c_void_p]
def TriangularInverse(uplo,diag,A):
  args = [uplo,diag,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElTriangularInverse_s(*args)
    elif A.tag == dTag: lib.ElTriangularInverse_d(*args)
    elif A.tag == cTag: lib.ElTriangularInverse_c(*args)
    elif A.tag == zTag: lib.ElTriangularInverse_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElTriangularInverseDist_s(*args)
    elif A.tag == dTag: lib.ElTriangularInverseDist_d(*args)
    elif A.tag == cTag: lib.ElTriangularInverseDist_c(*args)
    elif A.tag == zTag: lib.ElTriangularInverseDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Pseudoinverse
# =============
lib.ElPseudoinverse_s.argtypes = \
lib.ElPseudoinverse_d.argtypes = \
lib.ElPseudoinverse_c.argtypes = \
lib.ElPseudoinverse_z.argtypes = \
lib.ElPseudoinverseDist_s.argtypes = \
lib.ElPseudoinverseDist_d.argtypes = \
lib.ElPseudoinverseDist_c.argtypes = \
lib.ElPseudoinverseDist_z.argtypes = \
  [c_void_p]
def Pseudoinverse(A):
  args = [A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElPseudoinverse_s(*args)
    elif A.tag == dTag: lib.ElPseudoinverse_d(*args)
    elif A.tag == cTag: lib.ElPseudoinverse_c(*args)
    elif A.tag == zTag: lib.ElPseudoinverse_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElPseudoinverseDist_s(*args)
    elif A.tag == dTag: lib.ElPseudoinverseDist_d(*args)
    elif A.tag == cTag: lib.ElPseudoinverseDist_c(*args)
    elif A.tag == zTag: lib.ElPseudoinverseDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElHermitianPseudoinverse_s.argtypes = \
lib.ElHermitianPseudoinverse_d.argtypes = \
lib.ElHermitianPseudoinverse_c.argtypes = \
lib.ElHermitianPseudoinverse_z.argtypes = \
lib.ElHermitianPseudoinverseDist_s.argtypes = \
lib.ElHermitianPseudoinverseDist_d.argtypes = \
lib.ElHermitianPseudoinverseDist_c.argtypes = \
lib.ElHermitianPseudoinverseDist_z.argtypes = \
  [c_uint,c_void_p]
def HermitianPseudoinverse(uplo,A):
  args = [uplo,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHermitianPseudoinverse_s(*args)
    elif A.tag == dTag: lib.ElHermitianPseudoinverse_d(*args)
    elif A.tag == cTag: lib.ElHermitianPseudoinverse_c(*args)
    elif A.tag == zTag: lib.ElHermitianPseudoinverse_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHermitianPseudoinverseDist_s(*args)
    elif A.tag == dTag: lib.ElHermitianPseudoinverseDist_d(*args)
    elif A.tag == cTag: lib.ElHermitianPseudoinverseDist_c(*args)
    elif A.tag == zTag: lib.ElHermitianPseudoinverseDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Sign
# ====
lib.ElSign_s.argtypes = \
lib.ElSign_d.argtypes = \
lib.ElSign_c.argtypes = \
lib.ElSign_z.argtypes = \
lib.ElSignDist_s.argtypes = \
lib.ElSignDist_d.argtypes = \
lib.ElSignDist_c.argtypes = \
lib.ElSignDist_z.argtypes = \
  [c_void_p]

lib.ElSignDecomp_s.argtypes = \
lib.ElSignDecomp_d.argtypes = \
lib.ElSignDecomp_c.argtypes = \
lib.ElSignDecomp_z.argtypes = \
lib.ElSignDecompDist_s.argtypes = \
lib.ElSignDecompDist_d.argtypes = \
lib.ElSignDecompDist_c.argtypes = \
lib.ElSignDecompDist_z.argtypes = \
  [c_void_p,c_void_p]

def Sign(A,fullDecomp=False):
  if type(A) is Matrix:
    if fullDecomp:
      N = Matrix(A.tag)
      args = [A.obj,N.obj]
      if   A.tag == sTag: lib.ElSignDecomp_s(*args)
      elif A.tag == dTag: lib.ElSignDecomp_d(*args)
      elif A.tag == cTag: lib.ElSignDecomp_c(*args)
      elif A.tag == zTag: lib.ElSignDecomp_z(*args)
      else: DataExcept()
      return N
    else:
      args = [A.obj]
      if   A.tag == sTag: lib.ElSign_s(*args)
      elif A.tag == dTag: lib.ElSign_d(*args)
      elif A.tag == cTag: lib.ElSign_c(*args)
      elif A.tag == zTag: lib.ElSign_z(*args)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if fullDecomp:
      N = Matrix(A.tag)
      args = [A.obj,N.obj]
      if   A.tag == sTag: lib.ElSignDecompDist_s(*args)
      elif A.tag == dTag: lib.ElSignDecompDist_d(*args)
      elif A.tag == cTag: lib.ElSignDecompDist_c(*args)
      elif A.tag == zTag: lib.ElSignDecompDist_z(*args)
      else: DataExcept()
      return N
    else:
      args = [A.obj]
      if   A.tag == sTag: lib.ElSignDist_s(*args)
      elif A.tag == dTag: lib.ElSignDist_d(*args)
      elif A.tag == cTag: lib.ElSignDist_c(*args)
      elif A.tag == zTag: lib.ElSignDist_z(*args)
      else: DataExcept()
  else: TypeExcept()

lib.ElHermitianSign_s.argtypes = \
lib.ElHermitianSign_d.argtypes = \
lib.ElHermitianSign_c.argtypes = \
lib.ElHermitianSign_z.argtypes = \
lib.ElHermitianSignDist_s.argtypes = \
lib.ElHermitianSignDist_d.argtypes = \
lib.ElHermitianSignDist_c.argtypes = \
lib.ElHermitianSignDist_z.argtypes = \
  [c_uint,c_void_p]

lib.ElHermitianSignDecomp_s.argtypes = \
lib.ElHermitianSignDecomp_d.argtypes = \
lib.ElHermitianSignDecomp_c.argtypes = \
lib.ElHermitianSignDecomp_z.argtypes = \
lib.ElHermitianSignDecompDist_s.argtypes = \
lib.ElHermitianSignDecompDist_d.argtypes = \
lib.ElHermitianSignDecompDist_c.argtypes = \
lib.ElHermitianSignDecompDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p]

def HermitianSign(uplo,A,fullDecomp=False):
  if type(A) is Matrix:
    if fullDecomp:
      N = Matrix(A.tag)
      args = [uplo,A.obj,N.obj]
      if   A.tag == sTag: lib.ElHermitianSignDecomp_s(*args)
      elif A.tag == dTag: lib.ElHermitianSignDecomp_d(*args)
      elif A.tag == cTag: lib.ElHermitianSignDecomp_c(*args)
      elif A.tag == zTag: lib.ElHermitianSignDecomp_z(*args)
      else: DataExcept()
      return N
    else:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElHermitianSign_s(*args)
      elif A.tag == dTag: lib.ElHermitianSign_d(*args)
      elif A.tag == cTag: lib.ElHermitianSign_c(*args)
      elif A.tag == zTag: lib.ElHermitianSign_z(*args)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if fullDecomp:
      N = Matrix(A.tag)
      args = [uplo,A.obj,N.obj]
      if   A.tag == sTag: lib.ElHermitianSignDecompDist_s(*args)
      elif A.tag == dTag: lib.ElHermitianSignDecompDist_d(*args)
      elif A.tag == cTag: lib.ElHermitianSignDecompDist_c(*args)
      elif A.tag == zTag: lib.ElHermitianSignDecompDist_z(*args)
      else: DataExcept()
      return N
    else:
      args = [uplo,A.obj]
      if   A.tag == sTag: lib.ElHermitianSignDist_s(*args)
      elif A.tag == dTag: lib.ElHermitianSignDist_d(*args)
      elif A.tag == cTag: lib.ElHermitianSignDist_c(*args)
      elif A.tag == zTag: lib.ElHermitianSignDist_z(*args)
      else: DataExcept()
  else: TypeExcept()

# Square-root
# ===========
lib.ElSquareRoot_s.argtypes = \
lib.ElSquareRoot_d.argtypes = \
lib.ElSquareRoot_c.argtypes = \
lib.ElSquareRoot_z.argtypes = \
lib.ElSquareRootDist_s.argtypes = \
lib.ElSquareRootDist_d.argtypes = \
lib.ElSquareRootDist_c.argtypes = \
lib.ElSquareRootDist_z.argtypes = \
  [c_void_p]

def SquareRoot(A):
  args = [A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSquareRoot_s(*args)
    elif A.tag == dTag: lib.ElSquareRoot_d(*args)
    elif A.tag == cTag: lib.ElSquareRoot_c(*args)
    elif A.tag == zTag: lib.ElSquareRoot_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSquareRootDist_s(*args)
    elif A.tag == dTag: lib.ElSquareRootDist_d(*args)
    elif A.tag == cTag: lib.ElSquareRootDist_c(*args)
    elif A.tag == zTag: lib.ElSquareRootDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElHPSDSquareRoot_s.argtypes = \
lib.ElHPSDSquareRoot_d.argtypes = \
lib.ElHPSDSquareRoot_c.argtypes = \
lib.ElHPSDSquareRoot_z.argtypes = \
lib.ElHPSDSquareRootDist_s.argtypes = \
lib.ElHPSDSquareRootDist_d.argtypes = \
lib.ElHPSDSquareRootDist_c.argtypes = \
lib.ElHPSDSquareRootDist_z.argtypes = \
  [c_uint,c_void_p]

def HPSDSquareRoot(uplo,A):
  args = [uplo,A.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPSDSquareRoot_s(*args)
    elif A.tag == dTag: lib.ElHPSDSquareRoot_d(*args)
    elif A.tag == cTag: lib.ElHPSDSquareRoot_c(*args)
    elif A.tag == zTag: lib.ElHPSDSquareRoot_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPSDSquareRootDist_s(*args)
    elif A.tag == dTag: lib.ElHPSDSquareRootDist_d(*args)
    elif A.tag == cTag: lib.ElHPSDSquareRootDist_c(*args)
    elif A.tag == zTag: lib.ElHPSDSquareRootDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
