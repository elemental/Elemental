#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *
from factor import *
import ctypes

# Least squares
# =============
class LeastSquaresCtrl_s(ctypes.Structure):
  _fields_ = [("scaleTwoNorm",bType),("basisSize",iType),("alpha",sType),
              ("solveCtrl",RegSolveCtrl_s),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElLeastSquaresCtrlDefault_s(pointer(self))
class LeastSquaresCtrl_d(ctypes.Structure):
  _fields_ = [("scaleTwoNorm",bType),("basisSize",iType),("alpha",dType),
              ("solveCtrl",RegSolveCtrl_d),
              ("equilibrate",bType),("progress",bType),("time",bType)]
  def __init__(self):
    lib.ElLeastSquaresCtrlDefault_d(pointer(self))

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

lib.ElLeastSquaresXSparse_s.argtypes = \
lib.ElLeastSquaresXSparse_c.argtypes = \
lib.ElLeastSquaresXDistSparse_s.argtypes = \
lib.ElLeastSquaresXDistSparse_c.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,LeastSquaresCtrl_s]

lib.ElLeastSquaresXSparse_d.argtypes = \
lib.ElLeastSquaresXSparse_z.argtypes = \
lib.ElLeastSquaresXDistSparse_d.argtypes = \
lib.ElLeastSquaresXDistSparse_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,LeastSquaresCtrl_d]

def LeastSquares(A,B,ctrl=None,orient=NORMAL):
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
    argsCtrl = [orient,A.obj,B.obj,X.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElLeastSquaresSparse_s(*args)
      else:          lib.ElLeastSquaresXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElLeastSquaresSparse_d(*args)
      else:          lib.ElLeastSquaresXSparse_d(*argsCtrl)
    elif A.tag == cTag: 
      if ctrl==None: lib.ElLeastSquaresSparse_c(*args)
      else:          lib.ElLeastSquaresXSparse_c(*argsCtrl)
    elif A.tag == zTag: 
      if ctrl==None: lib.ElLeastSquaresSparse_z(*args)
      else:          lib.ElLeastSquaresXSparse_z(*argsCtrl)
    else: DataExcept()
    return X
  elif type(A) is DistSparseMatrix:
    if type(B) is not DistMultiVec:
      raise Exception('RHS was expected to be a DistMultiVec')
    X = DistMultiVec(A.tag,A.Comm())
    args = [orient,A.obj,B.obj,X.obj]
    argsCtrl = [orient,A.obj,B.obj,X.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElLeastSquaresDistSparse_s(*args)
      else:          lib.ElLeastSquaresXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElLeastSquaresDistSparse_d(*args)
      else:          lib.ElLeastSquaresXDistSparse_d(*argsCtrl)
    elif A.tag == cTag: 
      if ctrl==None: lib.ElLeastSquaresDistSparse_c(*args)
      else:          lib.ElLeastSquaresXDistSparse_c(*argsCtrl)
    elif A.tag == zTag: 
      if ctrl==None: lib.ElLeastSquaresDistSparse_z(*args)
      else:          lib.ElLeastSquaresXDistSparse_z(*argsCtrl)
    else: DataExcept()
    return X
  else: TypeExcept()

# Ridge regression
# ================
(RIDGE_CHOLESKY,RIDGE_QR,RIDGE_SVD)=(0,1,2)

lib.ElRidge_s.argtypes = \
lib.ElRidge_c.argtypes = \
lib.ElRidgeDist_s.argtypes = \
lib.ElRidgeDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,sType,c_void_p,c_uint]
lib.ElRidgeSparse_s.argtypes = \
lib.ElRidgeSparse_c.argtypes = \
lib.ElRidgeDistSparse_s.argtypes = \
lib.ElRidgeDistSparse_c.argtypes = \
  [c_uint,c_void_p,c_void_p,sType,c_void_p]

lib.ElRidge_d.argtypes = \
lib.ElRidge_z.argtypes = \
lib.ElRidgeDist_d.argtypes = \
lib.ElRidgeDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,dType,c_void_p,c_uint]
lib.ElRidgeSparse_d.argtypes = \
lib.ElRidgeSparse_z.argtypes = \
lib.ElRidgeDistSparse_d.argtypes = \
lib.ElRidgeDistSparse_z.argtypes = \
  [c_uint,c_void_p,c_void_p,dType,c_void_p]

def Ridge(A,B,gamma,orient=NORMAL,alg=RIDGE_CHOLESKY):
  if A.tag != B.tag:
    raise Exception('Datatypes of A and B must match')
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('Expected RHS to be a Matrix')
    X = Matrix(A.tag)
    args = [orient,A.obj,B.obj,gamma,X.obj,alg]
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
    args = [orient,A.obj,B.obj,gamma,X.obj,alg]
    if   A.tag == sTag: lib.ElRidgeDist_s(*args)
    elif A.tag == dTag: lib.ElRidgeDist_d(*args)
    elif A.tag == cTag: lib.ElRidgeDist_c(*args)
    elif A.tag == zTag: lib.ElRidgeDist_z(*args)
    else: DataExcept()
    return X
  elif type(A) is SparseMatrix:
    if type(B) is not Matrix:
      raise Exception('Expected RHS to be a Matrix')
    X = Matrix(A.tag)
    args = [orient,A.obj,B.obj,gamma,X.obj]
    if   A.tag == sTag: lib.ElRidgeSparse_s(*args)
    elif A.tag == dTag: lib.ElRidgeSparse_d(*args)
    elif A.tag == cTag: lib.ElRidgeSparse_c(*args)
    elif A.tag == zTag: lib.ElRidgeSparse_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistSparseMatrix:
    if type(B) is not DistMultiVec:
      raise Exception('Expected RHS to be a DistMultiVec')
    X = DistMultiVec(A.tag,A.Comm())
    args = [orient,A.obj,B.obj,gamma,X.obj]
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

lib.ElTikhonov_s.argtypes = \
lib.ElTikhonov_d.argtypes = \
lib.ElTikhonov_c.argtypes = \
lib.ElTikhonov_z.argtypes = \
lib.ElTikhonovDist_s.argtypes = \
lib.ElTikhonovDist_d.argtypes = \
lib.ElTikhonovDist_c.argtypes = \
lib.ElTikhonovDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p,c_uint]
lib.ElTikhonovSparse_s.argtypes = \
lib.ElTikhonovSparse_d.argtypes = \
lib.ElTikhonovSparse_c.argtypes = \
lib.ElTikhonovSparse_z.argtypes = \
lib.ElTikhonovDistSparse_s.argtypes = \
lib.ElTikhonovDistSparse_d.argtypes = \
lib.ElTikhonovDistSparse_c.argtypes = \
lib.ElTikhonovDistSparse_z.argtypes = \
  [c_uint,c_void_p,c_void_p,c_void_p,c_void_p]

def Tikhonov(A,B,G,orient=NORMAL,alg=TIKHONOV_CHOLESKY):
  if type(A) is not type(G):
    raise Exception('Matrix types of A and G must match')
  if A.tag != B.tag or B.tag != G.tag:
    raise Exception('Datatypes of {A,B,G} must match')
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('Expected RHS to be a Matrix')
    X = Matrix(A.tag)
    args = [orient,A.obj,B.obj,G.obj,X.obj,alg]
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
    args = [orient,A.obj,B.obj,G.obj,X.obj,alg]
    if   A.tag == sTag: lib.ElTikhonovDist_s(*args)
    elif A.tag == dTag: lib.ElTikhonovDist_d(*args)
    elif A.tag == cTag: lib.ElTikhonovDist_c(*args)
    elif A.tag == zTag: lib.ElTikhonovDist_z(*args)
    else: DataExcept()
    return X
  elif type(A) is SparseMatrix:
    if type(B) is not Matrix:
      raise Exception('Expected RHS to be a Matrix')
    X = Matrix(A.tag)
    args = [orient,A.obj,B.obj,G.obj,X.obj]
    if   A.tag == sTag: lib.ElTikhonovSparse_s(*args)
    elif A.tag == dTag: lib.ElTikhonovSparse_d(*args)
    elif A.tag == cTag: lib.ElTikhonovSparse_c(*args)
    elif A.tag == zTag: lib.ElTikhonovSparse_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistSparseMatrix:
    if type(B) is not DistMultiVec:
      raise Exception('Expected RHS to be a DistMultiVec')
    X = DistMultiVec(A.tag,A.Comm())
    args = [orient,A.obj,B.obj,G.obj,X.obj]
    if   A.tag == sTag: lib.ElTikhonovDistSparse_s(*args)
    elif A.tag == dTag: lib.ElTikhonovDistSparse_d(*args)
    elif A.tag == cTag: lib.ElTikhonovDistSparse_c(*args)
    elif A.tag == zTag: lib.ElTikhonovDistSparse_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()

# Equality-constrained least squares
# ==================================
lib.ElLSE_s.argtypes = \
lib.ElLSE_d.argtypes = \
lib.ElLSE_c.argtypes = \
lib.ElLSE_z.argtypes = \
lib.ElLSEDist_s.argtypes = \
lib.ElLSEDist_d.argtypes = \
lib.ElLSEDist_c.argtypes = \
lib.ElLSEDist_z.argtypes = \
lib.ElLSESparse_s.argtypes = \
lib.ElLSESparse_d.argtypes = \
lib.ElLSESparse_c.argtypes = \
lib.ElLSESparse_z.argtypes = \
lib.ElLSEDistSparse_s.argtypes = \
lib.ElLSEDistSparse_d.argtypes = \
lib.ElLSEDistSparse_c.argtypes = \
lib.ElLSEDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLSEXSparse_s.argtypes = \
lib.ElLSEXSparse_c.argtypes = \
lib.ElLSEXDistSparse_s.argtypes = \
lib.ElLSEXDistSparse_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,LeastSquaresCtrl_s]
lib.ElLSEXSparse_d.argtypes = \
lib.ElLSEXSparse_z.argtypes = \
lib.ElLSEXDistSparse_d.argtypes = \
lib.ElLSEXDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,LeastSquaresCtrl_d]

def LSE(A,B,C,D,ctrl=None):
  if type(A) is not type(B):
    raise Exception('Matrix types of A and B must match')
  if type(C) is not type(D):
    raise Exception('Matrix types of C and D must match')
  if A.tag != B.tag or B.tag != C.tag or C.tag != D.tag:
    raise Exception('Datatypes of {A,B,C,D} must match')
  if type(A) is Matrix:
    if type(C) is not Matrix:
      raise Exception('Expected C to be a Matrix')
    X = Matrix(A.tag)
    args = [A.obj,B.obj,C.obj,D.obj,X.obj]
    if   A.tag == sTag: lib.ElLSE_s(*args)
    elif A.tag == dTag: lib.ElLSE_d(*args)
    elif A.tag == cTag: lib.ElLSE_c(*args)
    elif A.tag == zTag: lib.ElLSE_z(*args)
    else: DataExcept()
    return X
  elif type(A) is DistMatrix:
    if type(C) is not DistMatrix:
      raise Exception('Expected C to be a DistMatrix')
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,B.obj,C.obj,D.obj,X.obj]
    if   A.tag == sTag: lib.ElLSEDist_s(*args)
    elif A.tag == dTag: lib.ElLSEDist_d(*args)
    elif A.tag == cTag: lib.ElLSEDist_c(*args)
    elif A.tag == zTag: lib.ElLSEDist_z(*args)
    else: DataExcept()
    return X
  elif type(A) is SparseMatrix:
    if type(C) is not Matrix:
      raise Exception('Expected C to be a Matrix')
    X = Matrix(A.tag)
    args = [A.obj,B.obj,C.obj,D.obj,X.obj]
    argsCtrl = [A.obj,B.obj,C.obj,D.obj,X.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElLSESparse_s(*args)
      else:          lib.ElLSEXSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElLSESparse_d(*args)
      else:          lib.ElLSEXSparse_d(*argsCtrl)
    elif A.tag == cTag: 
      if ctrl==None: lib.ElLSESparse_c(*args)
      else:          lib.ElLSEXSparse_c(*argsCtrl)
    elif A.tag == zTag: 
      if ctrl==None: lib.ElLSESparse_z(*args)
      else:          lib.ElLSEXSparse_z(*argsCtrl)
    else: DataExcept()
    return X
  elif type(A) is DistSparseMatrix:
    if type(C) is not DistMultiVec:
      raise Exception('Expected C to be a DistMultiVec')
    X = DistMultiVec(A.tag,A.Comm())
    args = [A.obj,B.obj,C.obj,D.obj,X.obj]
    argsCtrl = [A.obj,B.obj,C.obj,D.obj,X.obj,ctrl]
    if   A.tag == sTag: 
      if ctrl==None: lib.ElLSEDistSparse_s(*args)
      else:          lib.ElLSEXDistSparse_s(*argsCtrl)
    elif A.tag == dTag: 
      if ctrl==None: lib.ElLSEDistSparse_d(*args)
      else:          lib.ElLSEXDistSparse_d(*argsCtrl)
    elif A.tag == cTag: 
      if ctrl==None: lib.ElLSEDistSparse_c(*args)
      else:          lib.ElLSEXDistSparse_c(*argsCtrl)
    elif A.tag == zTag: 
      if ctrl==None: lib.ElLSEDistSparse_z(*args)
      else:          lib.ElLSEXDistSparse_z(*argsCtrl)
    else: DataExcept()
    return X
  else: TypeExcept()

# General Linear Model
# ====================
lib.ElGLM_s.argtypes = \
lib.ElGLM_d.argtypes = \
lib.ElGLM_c.argtypes = \
lib.ElGLM_z.argtypes = \
lib.ElGLMDist_s.argtypes = \
lib.ElGLMDist_d.argtypes = \
lib.ElGLMDist_c.argtypes = \
lib.ElGLMDist_z.argtypes = \
lib.ElGLMSparse_s.argtypes = \
lib.ElGLMSparse_d.argtypes = \
lib.ElGLMSparse_c.argtypes = \
lib.ElGLMSparse_z.argtypes = \
lib.ElGLMDistSparse_s.argtypes = \
lib.ElGLMDistSparse_d.argtypes = \
lib.ElGLMDistSparse_c.argtypes = \
lib.ElGLMDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElGLMSparse_s.argtypes = \
lib.ElGLMSparse_c.argtypes = \
lib.ElGLMDistSparse_s.argtypes = \
lib.ElGLMDistSparse_c.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,LeastSquaresCtrl_s]
lib.ElGLMSparse_d.argtypes = \
lib.ElGLMSparse_z.argtypes = \
lib.ElGLMDistSparse_d.argtypes = \
lib.ElGLMDistSparse_z.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,LeastSquaresCtrl_d]

def GLM(A,B,D,ctrl=None):
  if type(A) is not type(B):
    raise Exception('Expected types of A and B to match')
  if A.tag != B.tag or B.tag != D.tag:
    raise Exception('Datatypes of {A,B,D} must match')

  if type(A) is Matrix:
    if type(D) is not Matrix:
      raise Exception('Expected D to be a Matrix')
    X = Matrix(A.tag)
    Y = Matrix(A.tag)
    args = [A.obj,B.obj,D.obj,X.obj,Y.obj]
    if   A.tag == sTag: lib.ElGLM_s(*args)
    elif A.tag == dTag: lib.ElGLM_d(*args)
    elif A.tag == cTag: lib.ElGLM_c(*args)
    elif A.tag == zTag: lib.ElGLM_z(*args)
    else: DataExcept()
    return X, Y
  elif type(A) is DistMatrix:
    if type(D) is not DistMatrix:
      raise Exception('Expected D to be a DistMatrix')
    X = DistMatrix(A.tag,MC,MR,A.Grid())
    Y = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,B.obj,D.obj,X.obj,Y.obj]
    if   A.tag == sTag: lib.ElGLMDist_s(*args)
    elif A.tag == dTag: lib.ElGLMDist_d(*args)
    elif A.tag == cTag: lib.ElGLMDist_c(*args)
    elif A.tag == zTag: lib.ElGLMDist_z(*args)
    else: DataExcept()
    return X, Y
  elif type(A) is SparseMatrix:
    if type(D) is not Matrix:
      raise Exception('Expected D to be a Matrix')
    X = Matrix(A.tag)
    Y = Matrix(A.tag)
    args = [A.obj,B.obj,D.obj,X.obj,Y.obj]
    argsCtrl = [A.obj,B.obj,D.obj,X.obj,Y.obj,ctrl]
    if   A.tag == sTag:
      if ctrl==None: lib.ElGLMSparse_s(*args)
      else:          lib.ElGLMXSparse_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl==None: lib.ElGLMSparse_d(*args)
      else:          lib.ElGLMXSparse_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl==None: lib.ElGLMSparse_c(*args)
      else:          lib.ElGLMXSparse_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl==None: lib.ElGLMSparse_z(*args)
      else:          lib.ElGLMXSparse_z(*argsCtrl)
    else: DataExcept()
    return X, Y
  elif type(A) is DistSparseMatrix:
    if type(D) is not DistMultiVec:
      raise Exception('Expected D to be a DistMultiVec')
    X = DistMultiVec(A.tag)
    Y = DistMultiVec(A.tag)
    args = [A.obj,B.obj,D.obj,X.obj,Y.obj]
    argsCtrl = [A.obj,B.obj,D.obj,X.obj,Y.obj,ctrl]
    if   A.tag == sTag:
      if ctrl==None: lib.ElGLMDistSparse_s(*args)
      else:          lib.ElGLMXDistSparse_s(*argsCtrl)
    elif A.tag == dTag:
      if ctrl==None: lib.ElGLMDistSparse_d(*args)
      else:          lib.ElGLMXDistSparse_d(*argsCtrl)
    elif A.tag == cTag:
      if ctrl==None: lib.ElGLMDistSparse_c(*args)
      else:          lib.ElGLMXDistSparse_c(*argsCtrl)
    elif A.tag == zTag:
      if ctrl==None: lib.ElGLMDistSparse_z(*args)
      else:          lib.ElGLMXDistSparse_z(*argsCtrl)
    else: DataExcept()
    return X, Y
  else: TypeExcept()

