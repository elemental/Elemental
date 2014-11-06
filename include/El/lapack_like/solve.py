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

# Gaussian elimination
# ====================
lib.ElGaussianElimination_s.argtypes = [c_void_p,c_void_p]
lib.ElGaussianElimination_s.restype = c_uint
lib.ElGaussianElimination_d.argtypes = [c_void_p,c_void_p]
lib.ElGaussianElimination_d.restype = c_uint
lib.ElGaussianElimination_c.argtypes = [c_void_p,c_void_p]
lib.ElGaussianElimination_c.restype = c_uint
lib.ElGaussianElimination_z.argtypes = [c_void_p,c_void_p]
lib.ElGaussianElimination_z.restype = c_uint
lib.ElGaussianEliminationDist_s.argtypes = [c_void_p,c_void_p]
lib.ElGaussianEliminationDist_s.restype = c_uint
lib.ElGaussianEliminationDist_d.argtypes = [c_void_p,c_void_p]
lib.ElGaussianEliminationDist_d.restype = c_uint
lib.ElGaussianEliminationDist_c.argtypes = [c_void_p,c_void_p]
lib.ElGaussianEliminationDist_c.restype = c_uint
lib.ElGaussianEliminationDist_z.argtypes = [c_void_p,c_void_p]
lib.ElGaussianEliminationDist_z.restype = c_uint
def GaussianElimination(A,B):
  args = [A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGaussianElimination_s(*args)
    elif A.tag == dTag: lib.ElGaussianElimination_d(*args)
    elif A.tag == cTag: lib.ElGaussianElimination_c(*args)
    elif A.tag == zTag: lib.ElGaussianElimination_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGaussianEliminationDist_s(*args)
    elif A.tag == dTag: lib.ElGaussianEliminationDist_d(*args)
    elif A.tag == cTag: lib.ElGaussianEliminationDist_c(*args)
    elif A.tag == zTag: lib.ElGaussianEliminationDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# General Linear Model
# ====================
lib.ElGLM_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLM_s.restype = c_uint
lib.ElGLM_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLM_d.restype = c_uint
lib.ElGLM_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLM_c.restype = c_uint
lib.ElGLM_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLM_z.restype = c_uint
lib.ElGLMDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLMDist_s.restype = c_uint
lib.ElGLMDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLMDist_d.restype = c_uint
lib.ElGLMDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLMDist_c.restype = c_uint
lib.ElGLMDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLMDist_z.restype = c_uint
def GLM(A,B,D,Y):
  if type(A) is not type(B) or type(B) is not type(D) or type(D) is not type(Y):
    raise Exception('Matrix types of {A,B,D,Y} must match')
  if A.tag != B.tag or B.tag != D.tag or D.tag != Y.tag:
    raise Exception('Datatypes of {A,B,D,Y} must match')
  args = [A.obj,B.obj,D.obj,Y.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGLM_s(*args)
    elif A.tag == dTag: lib.ElGLM_d(*args)
    elif A.tag == cTag: lib.ElGLM_c(*args)
    elif A.tag == zTag: lib.ElGLM_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGLMDist_s(*args)
    elif A.tag == dTag: lib.ElGLMDist_d(*args)
    elif A.tag == cTag: lib.ElGLMDist_c(*args)
    elif A.tag == zTag: lib.ElGLMDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Symmetric/Hermitian solves
# ==========================
lib.ElSymmetricSolve_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolve_s.restype = c_uint
lib.ElSymmetricSolve_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolve_d.restype = c_uint
lib.ElSymmetricSolve_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolve_c.restype = c_uint
lib.ElSymmetricSolve_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolve_z.restype = c_uint
lib.ElSymmetricSolveDist_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolveDist_s.restype = c_uint
lib.ElSymmetricSolveDist_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolveDist_d.restype = c_uint
lib.ElSymmetricSolveDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolveDist_c.restype = c_uint
lib.ElSymmetricSolveDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolveDist_z.restype = c_uint
lib.ElSymmetricSolveDistSparse_s.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveDistSparse_s.restype = c_uint
lib.ElSymmetricSolveDistSparse_d.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveDistSparse_d.restype = c_uint
lib.ElSymmetricSolveDistSparse_c.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveDistSparse_c.restype = c_uint
lib.ElSymmetricSolveDistSparse_z.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveDistSparse_z.restype = c_uint
lib.ElHermitianSolveDistSparse_c.argtypes = [c_void_p,c_void_p]
lib.ElHermitianSolveDistSparse_c.restype = c_uint
lib.ElHermitianSolveDistSparse_z.argtypes = [c_void_p,c_void_p]
lib.ElHermitianSolveDistSparse_z.restype = c_uint
def SymmetricSolve(A,B,conjugate=False,uplo=LOWER,orient=NORMAL):
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('Expected RHS to be a Matrix')
    args = [uplo,orient,A.obj,B.obj]
    if   A.tag == sTag: lib.ElSymmetricSolve_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSolve_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianSolve_c(*args)
      else:         lib.ElSymmetricSolve_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianSolve_z(*args)
      else:         lib.ElSymmetricSolve_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(B) is not DistMatrix:
      raise Exception('Expected RHS to be a DistMatrix')
    args = [uplo,orient,A.obj,B.obj]
    if   A.tag == sTag: lib.ElSymmetricSolveDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSolveDist_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianSolveDist_c(*args)
      else:         lib.ElSymmetricSolveDist_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianSolveDist_z(*args)
      else:         lib.ElSymmetricSolveDist_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if orient != NORMAL:
      raise Exception('Non-normal sparse orientations not yet supported')
    if type(B) is not DistMultiVec:
      raise Exception('Expected RHS to be a DistMultiVec')
    args = [A.obj,B.obj]
    if   A.tag == sTag: lib.ElSymmetricSolveDistSparse_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSolveDistSparse_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianSolveDistSparse_c(*args)
      else:         lib.ElSymmetricSolveDistSparse_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianSolveDistSparse_z(*args)
      else:         lib.ElSymmetricSolveDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

def HermitianSolve(uplo,orient,A,B):
  SymmetricSolve(uplo,orient,A,B,True)

# Hermitian positive-definite solve
# =================================
lib.ElHPDSolve_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolve_s.restype = c_uint
lib.ElHPDSolve_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolve_d.restype = c_uint
lib.ElHPDSolve_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolve_c.restype = c_uint
lib.ElHPDSolve_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolve_z.restype = c_uint
lib.ElHPDSolveDist_s.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolveDist_s.restype = c_uint
lib.ElHPDSolveDist_d.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolveDist_d.restype = c_uint
lib.ElHPDSolveDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolveDist_c.restype = c_uint
lib.ElHPDSolveDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElHPDSolveDist_z.restype = c_uint
def HPDSolve(uplo,orient,A,B):
  if type(A) is not type(B):
    raise Exception('Matrix types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  args = [uplo,orient,A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPDSolve_s(*args)
    elif A.tag == dTag: lib.ElHPDSolve_d(*args)
    elif A.tag == cTag: lib.ElHPDSolve_c(*args)
    elif A.tag == zTag: lib.ElHPDSolve_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPDSolveDist_s(*args)
    elif A.tag == dTag: lib.ElHPDSolveDist_d(*args)
    elif A.tag == cTag: lib.ElHPDSolveDist_c(*args)
    elif A.tag == zTag: lib.ElHPDSolveDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Least squares
# =============
lib.ElLeastSquares_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquares_s.restype = c_uint
lib.ElLeastSquares_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquares_d.restype = c_uint
lib.ElLeastSquares_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquares_c.restype = c_uint
lib.ElLeastSquares_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquares_z.restype = c_uint
lib.ElLeastSquaresDist_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDist_s.restype = c_uint
lib.ElLeastSquaresDist_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDist_d.restype = c_uint
lib.ElLeastSquaresDist_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDist_c.restype = c_uint
lib.ElLeastSquaresDist_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDist_z.restype = c_uint
lib.ElLeastSquaresDistSparse_s.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_s.restype = c_uint
lib.ElLeastSquaresDistSparse_d.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_d.restype = c_uint
lib.ElLeastSquaresDistSparse_c.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_c.restype = c_uint
lib.ElLeastSquaresDistSparse_z.argtypes = [c_uint,c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_z.restype = c_uint
def LeastSquares(A,B,orient=NORMAL):
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('RHS was expected to be a Matrix')
    X = Matrix(A.tag)
    args = [orient,A.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElLeastSquares_s(*args)
    elif A.tag == dTag: lib.ElLeastSquares_d(*args)
    elif A.tag == cTag: lib.ElLeastSquares_c(*args)
    elif A.tag == zTag: lib.ElLeastSquares_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    if type(B) is not DistMatrix:
      raise Exception('RHS was expected to be a DistMatrix')
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [orient,A.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElLeastSquaresDist_s(*args)
    elif A.tag == dTag: lib.ElLeastSquaresDist_d(*args)
    elif A.tag == cTag: lib.ElLeastSquaresDist_c(*args)
    elif A.tag == zTag: lib.ElLeastSquaresDist_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistSparseMatrix:
    if type(B) is not DistMultiVec:
      raise Exception('RHS was expected to be a DistMultiVec')
    X = DistMultiVec(A.tag,A.Comm())
    args = [orient,A.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElLeastSquaresDistSparse_s(*args)
    elif A.tag == dTag: lib.ElLeastSquaresDistSparse_d(*args)
    elif A.tag == cTag: lib.ElLeastSquaresDistSparse_c(*args)
    elif A.tag == zTag: lib.ElLeastSquaresDistSparse_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

# Equality-constrained least squares
# ==================================
lib.ElLSE_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSE_s.restype = c_uint
lib.ElLSE_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSE_d.restype = c_uint
lib.ElLSE_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSE_c.restype = c_uint
lib.ElLSE_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSE_z.restype = c_uint
lib.ElLSEDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSEDist_s.restype = c_uint
lib.ElLSEDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSEDist_d.restype = c_uint
lib.ElLSEDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSEDist_c.restype = c_uint
lib.ElLSEDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSEDist_z.restype = c_uint
def LSE(A,B,C,D):
  if type(A) is not type(B) or type(B) is not type(C) or type(C) is not type(D):
    raise Exception('Matrix types of {A,B,C,D} must match')
  if A.tag != B.tag or B.tag != C.tag or C.tag != D.tag:
    raise Exception('Datatypes of {A,B,C,D} must match')
  if type(A) is Matrix:
    X = Matrix(A.tag)
    args = [A.obj,B.obj,C.obj,D.obj,X.obj]
    if   A.tag == sTag: lib.ElLSE_s(*args)
    elif A.tag == dTag: lib.ElLSE_d(*args)
    elif A.tag == cTag: lib.ElLSE_c(*args)
    elif A.tag == zTag: lib.ElLSE_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,B.obj,C.obj,D.obj,X.obj]
    if   A.tag == sTag: lib.ElLSEDist_s(*args)
    elif A.tag == dTag: lib.ElLSEDist_d(*args)
    elif A.tag == cTag: lib.ElLSEDist_c(*args)
    elif A.tag == zTag: lib.ElLSEDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

# Multishift Hessenberg solve
# ===========================
lib.ElMultiShiftHessSolve_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolve_s.restype = c_uint
lib.ElMultiShiftHessSolve_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolve_d.restype = c_uint
lib.ElMultiShiftHessSolve_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolve_c.restype = c_uint
lib.ElMultiShiftHessSolve_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolve_z.restype = c_uint
lib.ElMultiShiftHessSolveDist_s.argtypes = \
  [c_uint,c_uint,sType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolveDist_s.restype = c_uint
lib.ElMultiShiftHessSolveDist_d.argtypes = \
  [c_uint,c_uint,dType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolveDist_d.restype = c_uint
lib.ElMultiShiftHessSolveDist_c.argtypes = \
  [c_uint,c_uint,cType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolveDist_c.restype = c_uint
lib.ElMultiShiftHessSolveDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,c_void_p]
lib.ElMultiShiftHessSolveDist_z.restype = c_uint
def MultiShiftHessSolve(uplo,orient,alphaPre,H,shifts,X):
  if type(H) is not type(shifts) or type(shifts) is not type(X):
    raise Exception('Matrix types of {H,shifts,X} must match')
  if H.tag != shifts.tag or shifts.tag != X.tag:
    raise Exception('Datatypes of {H,shifts,X} must match')
  alpha = TagToType(H.tag)(alphaPre)
  if type(H) is Matrix:
    X = Matrix(H.tag)
    args = [uplo,orient,alpha,H.obj,shifts.obj,X.obj]
    if   A.tag == sTag: lib.ElMultiShiftHessSolve_s(*args)
    elif A.tag == dTag: lib.ElMultiShiftHessSolve_d(*args)
    elif A.tag == cTag: lib.ElMultiShiftHessSolve_c(*args)
    elif A.tag == zTag: lib.ElMultiShiftHessSolve_z(*args)
    else: DataExcept()
    return X
  elif type(H) is DistMatrix:
    X = DistMatrix(H.tag,MC,MR,H.Grid())
    args = [uplo,orient,alpha,H.obj,shifts.obj,X.obj]
    if   A.tag == sTag: lib.ElMultiShiftHessSolveDist_s(*args)
    elif A.tag == dTag: lib.ElMultiShiftHessSolveDist_d(*args)
    elif A.tag == cTag: lib.ElMultiShiftHessSolveDist_c(*args)
    elif A.tag == zTag: lib.ElMultiShiftHessSolveDist_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

# Ridge regression
# ================
(RIDGE_CHOLESKY,RIDGE_QR,RIDGE_SVD)=(0,1,2)
lib.ElRidge_s.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidge_s.restype = c_uint
lib.ElRidge_d.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidge_d.restype = c_uint
lib.ElRidge_c.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidge_c.restype = c_uint
lib.ElRidge_z.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidge_z.restype = c_uint
lib.ElRidgeDist_s.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidgeDist_s.restype = c_uint
lib.ElRidgeDist_d.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidgeDist_d.restype = c_uint
lib.ElRidgeDist_c.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidgeDist_c.restype = c_uint
lib.ElRidgeDist_z.argtypes = [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidgeDist_z.restype = c_uint
lib.ElRidgeDistSparse_s.argtypes = [c_void_p,c_void_p,sType,c_void_p]
lib.ElRidgeDistSparse_s.restype = c_uint
lib.ElRidgeDistSparse_d.argtypes = [c_void_p,c_void_p,sType,c_void_p]
lib.ElRidgeDistSparse_d.restype = c_uint
lib.ElRidgeDistSparse_c.argtypes = [c_void_p,c_void_p,sType,c_void_p]
lib.ElRidgeDistSparse_c.restype = c_uint
lib.ElRidgeDistSparse_z.argtypes = [c_void_p,c_void_p,sType,c_void_p]
lib.ElRidgeDistSparse_z.restype = c_uint
def Ridge(A,B,alpha,alg=RIDGE_CHOLESKY):
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('Expected RHS to be a Matrix')
    X = Matrix(A.tag)
    args = [A.obj,B.obj,alpha,X.obj,alg]
    if   A.tag == sTag: lib.ElRidge_s(*args)
    elif A.tag == dTag: lib.ElRidge_d(*args)
    elif A.tag == cTag: lib.ElRidge_c(*args)
    elif A.tag == zTag: lib.ElRidge_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    if type(B) is not DistMatrix:
      raise Exception('Expected RHS to be a DistMatrix')
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,B.obj,alpha,X.obj,alg]
    if   A.tag == sTag: lib.ElRidgeDist_s(*args)
    elif A.tag == dTag: lib.ElRidgeDist_d(*args)
    elif A.tag == cTag: lib.ElRidgeDist_c(*args)
    elif A.tag == zTag: lib.ElRidgeDist_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistSparseMatrix:
    if type(B) is not DistMultiVec:
      raise Exception('Expected RHS to be a DistMultiVec')
    X = DistMultiVec(A.tag,A.Comm())
    args = [A.obj,B.obj,alpha,X.obj]
    if   A.tag == sTag: lib.ElRidgeDistSparse_s(*args)
    elif A.tag == dTag: lib.ElRidgeDistSparse_d(*args)
    elif A.tag == cTag: lib.ElRidgeDistSparse_c(*args)
    elif A.tag == zTag: lib.ElRidgeDistSparse_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

# Tikhonov regularization
# =======================
(TIKHONOV_CHOLESKY,TIKHONOV_QR)=(0,1)

lib.ElTikhonov_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonov_s.restype = c_uint
lib.ElTikhonov_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonov_d.restype = c_uint
lib.ElTikhonov_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonov_c.restype = c_uint
lib.ElTikhonov_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonov_z.restype = c_uint
lib.ElTikhonovDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonovDist_s.restype = c_uint
lib.ElTikhonovDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonovDist_d.restype = c_uint
lib.ElTikhonovDist_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonovDist_c.restype = c_uint
lib.ElTikhonovDist_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonovDist_z.restype = c_uint
lib.ElTikhonovDistSparse_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElTikhonovDistSparse_s.restype = c_uint
lib.ElTikhonovDistSparse_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElTikhonovDistSparse_d.restype = c_uint
lib.ElTikhonovDistSparse_c.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElTikhonovDistSparse_c.restype = c_uint
lib.ElTikhonovDistSparse_z.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElTikhonovDistSparse_z.restype = c_uint

def Tikhonov(A,B,Gamma,alg=TIKHONOV_CHOLESKY):
  if type(A) is not type(Gamma):
    raise Exception('Matrix types of A and Gamma must match')
  if A.tag != B.tag or B.tag != Gamma.tag:
    raise Exception('Datatypes of {A,B,Gamma} must match')
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('Expected RHS to be a Matrix')
    X = Matrix(A.tag)
    args = [A.obj,B.obj,Gamma.obj,X.obj,alg]
    if   A.tag == sTag: lib.ElTikhonov_s(*args)
    elif A.tag == dTag: lib.ElTikhonov_d(*args)
    elif A.tag == cTag: lib.ElTikhonov_c(*args)
    elif A.tag == zTag: lib.ElTikhonov_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    if type(B) is not DistMatrix:
      raise Exception('Expected RHS to be a DistMatrix')
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,B.obj,Gamma.obj,X.obj,alg]
    if   A.tag == sTag: lib.ElTikhonovDist_s(*args)
    elif A.tag == dTag: lib.ElTikhonovDist_d(*args)
    elif A.tag == cTag: lib.ElTikhonovDist_c(*args)
    elif A.tag == zTag: lib.ElTikhonovDist_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistSparseMatrix:
    if type(B) is not DistMultiVec:
      raise Exception('Expected RHS to be a DistMultiVec')
    X = DistMultiVec(A.tag,A.Comm())
    args = [A.obj,B.obj,Gamma.obj,X.obj]
    if   A.tag == sTag: lib.ElTikhonovDistSparse_s(*args)
    elif A.tag == dTag: lib.ElTikhonovDistSparse_d(*args)
    elif A.tag == cTag: lib.ElTikhonovDistSparse_c(*args)
    elif A.tag == zTag: lib.ElTikhonovDistSparse_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()
