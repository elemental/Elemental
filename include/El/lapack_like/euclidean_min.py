#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
import ctypes

# Euclidean minimization
# ======================

# General Linear Model
# --------------------
lib.ElGLM_s.argtypes = \
lib.ElGLM_d.argtypes = \
lib.ElGLM_c.argtypes = \
lib.ElGLM_z.argtypes = \
lib.ElGLMDist_s.argtypes = \
lib.ElGLMDist_d.argtypes = \
lib.ElGLMDist_c.argtypes = \
lib.ElGLMDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElGLM_s.restype = \
lib.ElGLM_d.restype = \
lib.ElGLM_c.restype = \
lib.ElGLM_z.restype = \
lib.ElGLMDist_s.restype = \
lib.ElGLMDist_d.restype = \
lib.ElGLMDist_c.restype = \
lib.ElGLMDist_z.restype = \
  c_uint

def GLM(A,B,D):
  if type(A) is not type(B) or type(B) is not type(D):
    raise Exception('Matrix types of {A,B,D} must match')
  if A.tag != B.tag or B.tag != D.tag:
    raise Exception('Datatypes of {A,B,D} must match')

  if type(A) is Matrix:
    Y = Matrix(A.tag)
    args = [A.obj,B.obj,D.obj,Y.obj]
    if   A.tag == sTag: lib.ElGLM_s(*args)
    elif A.tag == dTag: lib.ElGLM_d(*args)
    elif A.tag == cTag: lib.ElGLM_c(*args)
    elif A.tag == zTag: lib.ElGLM_z(*args)
    else: DataExcept()
    return Y
  elif type(A) is DistMatrix:
    Y = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,B.obj,D.obj,Y.obj]
    if   A.tag == sTag: lib.ElGLMDist_s(*args)
    elif A.tag == dTag: lib.ElGLMDist_d(*args)
    elif A.tag == cTag: lib.ElGLMDist_c(*args)
    elif A.tag == zTag: lib.ElGLMDist_z(*args)
    else: DataExcept()
    return Y
  else: TypeExcept()

# Least squares
# -------------
lib.ElLeastSquares_s.argtypes = \
lib.ElLeastSquares_d.argtypes = \
lib.ElLeastSquares_c.argtypes = \
lib.ElLeastSquares_z.argtypes = \
lib.ElLeastSquaresDist_s.argtypes = \
lib.ElLeastSquaresDist_d.argtypes = \
lib.ElLeastSquaresDist_c.argtypes = \
lib.ElLeastSquaresDist_z.argtypes = \
lib.ElLeastSquaresSparse_s.argtypes = \
lib.ElLeastSquaresSparse_d.argtypes = \
lib.ElLeastSquaresSparse_c.argtypes = \
lib.ElLeastSquaresSparse_z.argtypes = \
lib.ElLeastSquaresDistSparse_s.argtypes = \
lib.ElLeastSquaresDistSparse_d.argtypes = \
lib.ElLeastSquaresDistSparse_c.argtypes = \
lib.ElLeastSquaresDistSparse_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p]

lib.ElLeastSquares_s.restype = \
lib.ElLeastSquares_d.restype = \
lib.ElLeastSquares_c.restype = \
lib.ElLeastSquares_z.restype = \
lib.ElLeastSquaresDist_s.restype = \
lib.ElLeastSquaresDist_d.restype = \
lib.ElLeastSquaresDist_c.restype = \
lib.ElLeastSquaresDist_z.restype = \
lib.ElLeastSquaresSparse_s.restype = \
lib.ElLeastSquaresSparse_d.restype = \
lib.ElLeastSquaresSparse_c.restype = \
lib.ElLeastSquaresSparse_z.restype = \
lib.ElLeastSquaresDistSparse_s.restype = \
lib.ElLeastSquaresDistSparse_d.restype = \
lib.ElLeastSquaresDistSparse_c.restype = \
lib.ElLeastSquaresDistSparse_z.restype = \
  c_uint

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
  elif type(A) is SparseMatrix:
    if type(B) is not Matrix:
      raise Exception('RHS was expected to be a Matrix')
    X = Matrix(A.tag)
    args = [orient,A.obj,B.obj,X.obj]
    if   A.tag == sTag: lib.ElLeastSquaresSparse_s(*args)
    elif A.tag == dTag: lib.ElLeastSquaresSparse_d(*args)
    elif A.tag == cTag: lib.ElLeastSquaresSparse_c(*args)
    elif A.tag == zTag: lib.ElLeastSquaresSparse_z(*args)
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
# ----------------------------------
lib.ElLSE_s.argtypes = \
lib.ElLSE_d.argtypes = \
lib.ElLSE_c.argtypes = \
lib.ElLSE_z.argtypes = \
lib.ElLSEDist_s.argtypes = \
lib.ElLSEDist_d.argtypes = \
lib.ElLSEDist_c.argtypes = \
lib.ElLSEDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElLSE_s.restype = \
lib.ElLSE_d.restype = \
lib.ElLSE_c.restype = \
lib.ElLSE_z.restype = \
lib.ElLSEDist_s.restype = \
lib.ElLSEDist_d.restype = \
lib.ElLSEDist_c.restype = \
lib.ElLSEDist_z.restype = \
  c_uint

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

# Ridge regression
# ----------------
(RIDGE_CHOLESKY,RIDGE_QR,RIDGE_SVD)=(0,1,2)

lib.ElRidge_s.argtypes = \
lib.ElRidge_c.argtypes = \
lib.ElRidgeDist_s.argtypes = \
lib.ElRidgeDist_c.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidgeDistSparse_s.argtypes = \
lib.ElRidgeDistSparse_c.argtypes = \
  [c_void_p,c_void_p,sType,c_void_p]

lib.ElRidge_d.argtypes = \
lib.ElRidge_z.argtypes = \
lib.ElRidgeDist_d.argtypes = \
lib.ElRidgeDist_z.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p,c_uint]
lib.ElRidgeDistSparse_d.argtypes = \
lib.ElRidgeDistSparse_z.argtypes = \
  [c_void_p,c_void_p,dType,c_void_p]

lib.ElRidge_s.restype = \
lib.ElRidge_d.restype = \
lib.ElRidge_c.restype = \
lib.ElRidge_z.restype = \
lib.ElRidgeDist_s.restype = \
lib.ElRidgeDist_d.restype = \
lib.ElRidgeDist_c.restype = \
lib.ElRidgeDist_z.restype = \
lib.ElRidgeDistSparse_s.restype = \
lib.ElRidgeDistSparse_d.restype = \
lib.ElRidgeDistSparse_c.restype = \
lib.ElRidgeDistSparse_z.restype = \
  c_uint

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
# -----------------------
(TIKHONOV_CHOLESKY,TIKHONOV_QR)=(0,1)

lib.ElTikhonov_s.argtypes = \
lib.ElTikhonov_d.argtypes = \
lib.ElTikhonov_c.argtypes = \
lib.ElTikhonov_z.argtypes = \
lib.ElTikhonovDist_s.argtypes = \
lib.ElTikhonovDist_d.argtypes = \
lib.ElTikhonovDist_c.argtypes = \
lib.ElTikhonovDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonovDistSparse_s.argtypes = \
lib.ElTikhonovDistSparse_d.argtypes = \
lib.ElTikhonovDistSparse_c.argtypes = \
lib.ElTikhonovDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p]

lib.ElTikhonov_s.restype = c_uint
lib.ElTikhonov_d.restype = c_uint
lib.ElTikhonov_c.restype = c_uint
lib.ElTikhonov_z.restype = c_uint
lib.ElTikhonovDist_s.restype = c_uint
lib.ElTikhonovDist_d.restype = c_uint
lib.ElTikhonovDist_c.restype = c_uint
lib.ElTikhonovDist_z.restype = c_uint
lib.ElTikhonovDistSparse_s.restype = c_uint
lib.ElTikhonovDistSparse_d.restype = c_uint
lib.ElTikhonovDistSparse_c.restype = c_uint
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
