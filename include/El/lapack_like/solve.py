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
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGaussianElimination_s(A.tag,B.tag)
    elif A.tag == dTag: lib.ElGaussianElimination_d(A.tag,B.tag)
    elif A.tag == cTag: lib.ElGaussianElimination_c(A.tag,B.tag)
    elif A.tag == zTag: lib.ElGaussianElimination_z(A.tag,B.tag)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGaussianEliminationDist_s(A.tag,B.tag)
    elif A.tag == dTag: lib.ElGaussianEliminationDist_d(A.tag,B.tag)
    elif A.tag == cTag: lib.ElGaussianEliminationDist_c(A.tag,B.tag)
    elif A.tag == zTag: lib.ElGaussianEliminationDist_z(A.tag,B.tag)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

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
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElGLM_s(A.obj,B.obj,D.obj,Y.obj)
    elif A.tag == dTag: lib.ElGLM_d(A.obj,B.obj,D.obj,Y.obj)
    elif A.tag == cTag: lib.ElGLM_c(A.obj,B.obj,D.obj,Y.obj)
    elif A.tag == zTag: lib.ElGLM_z(A.obj,B.obj,D.obj,Y.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElGLMDist_s(A.obj,B.obj,D.obj,Y.obj)
    elif A.tag == dTag: lib.ElGLMDist_d(A.obj,B.obj,D.obj,Y.obj)
    elif A.tag == cTag: lib.ElGLMDist_c(A.obj,B.obj,D.obj,Y.obj)
    elif A.tag == zTag: lib.ElGLMDist_z(A.obj,B.obj,D.obj,Y.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

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
def SymmetricSolve(uplo,orient,A,B,conjugate=False):
  if type(A) is not type(B):
    raise Exception('Types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSymmetricSolve_s(uplo,orient,A.tag,B.tag)
    elif A.tag == dTag: lib.ElSymmetricSolve_d(uplo,orient,A.tag,B.tag)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianSolve_c(uplo,orient,A.tag,B.tag)
      else:         lib.ElSymmetricSolve_c(uplo,orient,A.tag,B.tag)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianSolve_z(uplo,orient,A.tag,B.tag)
      else:         lib.ElSymmetricSolve_z(uplo,orient,A.tag,B.tag)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSymmetricSolveDist_s(uplo,orient,A.tag,B.tag)
    elif A.tag == dTag: lib.ElSymmetricSolveDist_d(uplo,orient,A.tag,B.tag)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianSolveDist_c(uplo,orient,A.tag,B.tag)
      else:         lib.ElSymmetricSolveDist_c(uplo,orient,A.tag,B.tag)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianSolveDist_z(uplo,orient,A.tag,B.tag)
      else:         lib.ElSymmetricSolveDist_z(uplo,orient,A.tag,B.tag)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

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
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHPDSolve_s(uplo,orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElHPDSolve_d(uplo,orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElHPDSolve_c(uplo,orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElHPDSolve_z(uplo,orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHPDSolveDist_s(uplo,orient,A.obj,B.obj)
    elif A.tag == dTag: lib.ElHPDSolveDist_d(uplo,orient,A.obj,B.obj)
    elif A.tag == cTag: lib.ElHPDSolveDist_c(uplo,orient,A.obj,B.obj)
    elif A.tag == zTag: lib.ElHPDSolveDist_z(uplo,orient,A.obj,B.obj)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

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
def LeastSquares(orient,A,B):
  if type(A) is not type(B):
    raise Exception('Matrix types of A and B must match')
  if A.tag != B.tag:
    raise ExceptioN('Datatypes of A and B must match')
  if type(A) is Matrix:
    X = Matrix(A.tag)
    if   A.tag == sTag: lib.ElLeastSquares_s(orient,A.obj,B.obj,X.obj)
    elif A.tag == dTag: lib.ElLeastSquares_d(orient,A.obj,B.obj,X.obj)
    elif A.tag == cTag: lib.ElLeastSquares_c(orient,A.obj,B.obj,X.obj)
    elif A.tag == zTag: lib.ElLeastSquares_z(orient,A.obj,B.obj,X.obj)
    else: raise Exception('Unsupported datatype') 
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    if   A.tag == sTag: lib.ElLeastSquaresDist_s(orient,A.obj,B.obj,X.obj)
    elif A.tag == dTag: lib.ElLeastSquaresDist_d(orient,A.obj,B.obj,X.obj)
    elif A.tag == cTag: lib.ElLeastSquaresDist_c(orient,A.obj,B.obj,X.obj)
    elif A.tag == zTag: lib.ElLeastSquaresDist_z(orient,A.obj,B.obj,X.obj)
    else: raise Exception('Unsupported datatype') 
    return X
  else: raise Exception('Unsupported matrix type')

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
    if   A.tag == sTag: lib.ElLSE_s(A.obj,B.obj,C.obj,D.obj,X.obj)
    elif A.tag == dTag: lib.ElLSE_d(A.obj,B.obj,C.obj,D.obj,X.obj)
    elif A.tag == cTag: lib.ElLSE_c(A.obj,B.obj,C.obj,D.obj,X.obj)
    elif A.tag == zTag: lib.ElLSE_z(A.obj,B.obj,C.obj,D.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    if   A.tag == sTag: lib.ElLSEDist_s(A.obj,B.obj,C.obj,D.obj,X.obj)
    elif A.tag == dTag: lib.ElLSEDist_d(A.obj,B.obj,C.obj,D.obj,X.obj)
    elif A.tag == cTag: lib.ElLSEDist_c(A.obj,B.obj,C.obj,D.obj,X.obj)
    elif A.tag == zTag: lib.ElLSEDist_z(A.obj,B.obj,C.obj,D.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  else: raise Exception('Unsupported matrix type')

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
    if   A.tag == sTag: 
      lib.ElMultiShiftHessSolve_s(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    elif A.tag == dTag:
      lib.ElMultiShiftHessSolve_d(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    elif A.tag == cTag:
      lib.ElMultiShiftHessSolve_c(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    elif A.tag == zTag:
      lib.ElMultiShiftHessSolve_z(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  elif type(H) is DistMatrix:
    X = DistMatrix(H.tag,MC,MR,H.Grid())
    if   A.tag == sTag: 
      lib.ElMultiShiftHessSolveDist_s(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    elif A.tag == dTag:
      lib.ElMultiShiftHessSolveDist_d(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    elif A.tag == cTag:
      lib.ElMultiShiftHessSolveDist_c(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    elif A.tag == zTag:
      lib.ElMultiShiftHessSolveDist_z(uplo,orient,alpha,H.obj,shifts.obj,X.obj)
    else: raise Exception('Unsupported datatype')
    return X
  else: raise Exception('Unsupported matrix type')

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
def Ridge(A,B,alpha,alg=RIDGE_CHOLESKY):
  if type(A) is not type(B):
    raise Exception('Matrix types of A and B must match')
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElRidge_s(A.obj,B.obj,alpha,X.obj,alg)
    elif A.tag == dTag: lib.ElRidge_d(A.obj,B.obj,alpha,X.obj,alg)
    elif A.tag == cTag: lib.ElRidge_c(A.obj,B.obj,alpha,X.obj,alg)
    elif A.tag == zTag: lib.ElRidge_z(A.obj,B.obj,alpha,X.obj,alg)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElRidgeDist_s(A.obj,B.obj,alpha,X.obj,alg)
    elif A.tag == dTag: lib.ElRidgeDist_d(A.obj,B.obj,alpha,X.obj,alg)
    elif A.tag == cTag: lib.ElRidgeDist_c(A.obj,B.obj,alpha,X.obj,alg)
    elif A.tag == zTag: lib.ElRidgeDist_z(A.obj,B.obj,alpha,X.obj,alg)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

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
def Tikhonov(A,B,Gamma,alg=TIKHONOV_CHOLESKY):
  if type(A) is not type(B) or type(B) is not type(Gamma):
    raise Exception('Matrix types of {A,B,Gamma} must match')
  if A.tag != B.tag or B.tag != Gamma.tag:
    raise Exception('Datatypes of {A,B,Gamma} must match')
  if type(A) is Matrix:
    X = Matrix(A.tag)
    if   A.tag == sTag: lib.ElTikhonov_s(A.obj,B.obj,Gamma.obj,X.obj,alg)
    elif A.tag == dTag: lib.ElTikhonov_d(A.obj,B.obj,Gamma.obj,X.obj,alg)
    elif A.tag == cTag: lib.ElTikhonov_c(A.obj,B.obj,Gamma.obj,X.obj,alg)
    elif A.tag == zTag: lib.ElTikhonov_z(A.obj,B.obj,Gamma.obj,X.obj,alg)
    else: raise Exception('Unsupported datatype')
    return X
  elif type(A) is DistMatrix:
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    if   A.tag == sTag: lib.ElTikhonovDist_s(A.obj,B.obj,Gamma.obj,X.obj,alg)
    elif A.tag == dTag: lib.ElTikhonovDist_d(A.obj,B.obj,Gamma.obj,X.obj,alg)
    elif A.tag == cTag: lib.ElTikhonovDist_c(A.obj,B.obj,Gamma.obj,X.obj,alg)
    elif A.tag == zTag: lib.ElTikhonovDist_z(A.obj,B.obj,Gamma.obj,X.obj,alg)
    else: raise Exception('Unsupported datatype')
    return X
  else: raise Exception('Unsupported matrix type')
