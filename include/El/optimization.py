#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

from ctypes import CFUNCTYPE

# Basis pursuit
# =============
lib.ElBasisPursuit_s.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuit_s.restype = c_uint
lib.ElBasisPursuit_d.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuit_d.restype = c_uint
lib.ElBasisPursuit_c.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuit_c.restype = c_uint
lib.ElBasisPursuit_z.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuit_z.restype = c_uint
lib.ElBasisPursuitDist_s.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuitDist_s.restype = c_uint
lib.ElBasisPursuitDist_d.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuitDist_d.restype = c_uint
lib.ElBasisPursuitDist_c.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuitDist_c.restype = c_uint
lib.ElBasisPursuitDist_z.argtypes = [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElBasisPursuitDist_z.restype = c_uint
def BasisPursuit(A,b):
  if type(A) is not type(b):
    raise Exception('Types of A and b must match')
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  numIts = iType()
  if type(A) is Matrix:
    z = Matrix(A.tag)
    args = [A.obj,b.obj,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElBasisPursuit_s(*args)
    elif A.tag == dTag: lib.ElBasisPursuit_d(*args)
    elif A.tag == cTag: lib.ElBasisPursuit_c(*args)
    elif A.tag == zTag: lib.ElBasisPursuit_z(*args)
    else: DataExcept()
    return z, numIts
  elif type(A) is DistMatrix:
    z = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElBasisPursuitDist_s(*args)
    elif A.tag == dTag: lib.ElBasisPursuitDist_d(*args)
    elif A.tag == cTag: lib.ElBasisPursuitDist_c(*args)
    elif A.tag == zTag: lib.ElBasisPursuitDist_z(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Least Absolute Shrinkage and Selection Operator
# ===============================================
lib.ElLasso_s.argtypes = [c_void_p,c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElLasso_s.restype = c_uint
lib.ElLasso_d.argtypes = [c_void_p,c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElLasso_d.restype = c_uint
lib.ElLasso_c.argtypes = [c_void_p,c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElLasso_c.restype = c_uint
lib.ElLasso_z.argtypes = [c_void_p,c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElLasso_z.restype = c_uint
lib.ElLassoDist_s.argtypes = [c_void_p,c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElLassoDist_s.restype = c_uint
lib.ElLassoDist_d.argtypes = [c_void_p,c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElLassoDist_d.restype = c_uint
lib.ElLassoDist_c.argtypes = [c_void_p,c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElLassoDist_c.restype = c_uint
lib.ElLassoDist_z.argtypes = [c_void_p,c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElLassoDist_z.restype = c_uint
def Lasso(A,b,lamb):
  if type(A) is not type(b):
    raise Exception('Types of A and b must match')
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  numIts = iType()
  if type(A) is Matrix:
    z = Matrix(A.tag)
    args = [A.obj,b.obj,lamb,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElLasso_s(*args)
    elif A.tag == dTag: lib.ElLasso_d(*args)
    elif A.tag == cTag: lib.ElLasso_c(*args)
    elif A.tag == zTag: lib.ElLasso_z(*args)
    else: DataExcept()
    return z, numIts
  elif type(A) is DistMatrix:
    z = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,lamb,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElLassoDist_s(*args)
    elif A.tag == dTag: lib.ElLassoDist_d(*args)
    elif A.tag == cTag: lib.ElLassoDist_c(*args)
    elif A.tag == zTag: lib.ElLassoDist_z(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Linear program
# ==============
lib.ElLinearProgram_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgram_s.restype = c_uint
lib.ElLinearProgram_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgram_d.restype = c_uint
lib.ElLinearProgramDist_s.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramDist_s.restype = c_uint
lib.ElLinearProgramDist_d.argtypes = [c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramDist_d.restype = c_uint
def LinearProgram(A,b,c):
  if type(A) is not type(b) or type(b) is not type(c):
    raise Exception('Types of {A,b,c} must match')
  if A.tag != b.tag or b.tag != c.tag:
    raise Exception('Datatypes of {A,b,c} must match')
  if type(A) is Matrix:
    x = Matrix(A.tag)
    args = [A.obj,b.obj,c.obj,x.obj]
    if   A.tag == sTag: lib.ElLinearProgram_s(*args)
    elif A.tag == dTag: lib.ElLinearProgram_d(*args)
    else: DataExcept()
    return x
  elif type(A) is DistMatrix:
    x = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,c.obj,x.obj]
    if   A.tag == sTag: lib.ElLinearProgramDist_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramDist_d(*args)
    else: DataExcept()
    return x
  else: TypeExcept()

lib.ElLinearProgramIPF_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPF_s.restype = c_uint
lib.ElLinearProgramIPF_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPF_d.restype = c_uint
lib.ElLinearProgramIPFDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPFDist_s.restype = c_uint
lib.ElLinearProgramIPFDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPFDist_d.restype = c_uint
lib.ElLinearProgramIPFSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPFSparse_s.restype = c_uint
lib.ElLinearProgramIPFSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPFSparse_d.restype = c_uint
lib.ElLinearProgramIPFDistSparse_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPFDistSparse_s.restype = c_uint
lib.ElLinearProgramIPFDistSparse_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramIPFDistSparse_d.restype = c_uint
def LinearProgramIPF(A,b,c,s,x,l):
  args = [A.obj,b.obj,c.obj,s.obj,x.obj,l.obj]
  if type(A) is Matrix:
    if type(b) is not Matrix or type(c) is not Matrix or \
       type(x) is not Matrix or type(l) is not Matrix or \
       type(s) is not Matrix:
      raise Exception('Expected {b,c,x,l,s} to be of type Matrix')
    if   A.tag == sTag: lib.ElLinearProgramIPF_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramIPF_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix or type(c) is not DistMatrix or \
       type(x) is not DistMatrix or type(l) is not DistMatrix or \
       type(s) is not DistMatrix:
      raise Exception('Expected {b,c,x,l,s} to be of type Matrix')
    if   A.tag == sTag: lib.ElLinearProgramIPFDist_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramIPFDist_d(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(b) is not Matrix or type(c) is not Matrix or \
       type(x) is not Matrix or type(l) is not Matrix or \
       type(s) is not Matrix:
      raise Exception('Expected {b,c,x,l,s} to be of type Matrix')
    if   A.tag == sTag: lib.ElLinearProgramIPFSparse_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramIPFSparse_d(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(b) is not DistMultiVec or type(c) is not DistMultiVec or \
       type(x) is not DistMultiVec or type(l) is not DistMultiVec or \
       type(s) is not DistMultiVec:
      raise Exception('Expected {b,c,x,l,s} to be of type DistMultiVec')
    if   A.tag == sTag: lib.ElLinearProgramIPFDistSparse_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramIPFDistSparse_d(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElLinearProgramMPC_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramMPC_s.restype = c_uint
lib.ElLinearProgramMPC_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramMPC_d.restype = c_uint
lib.ElLinearProgramMPCDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramMPCDist_s.restype = c_uint
lib.ElLinearProgramMPCDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,c_void_p,c_void_p]
lib.ElLinearProgramMPCDist_d.restype = c_uint
def LinearProgramMPC(A,b,c,s,x,l):
  args = [A.obj,b.obj,c.obj,s.obj,x.obj,l.obj]
  if type(A) is Matrix:
    if type(b) is not Matrix or type(c) is not Matrix or \
       type(x) is not Matrix or type(l) is not Matrix or \
       type(s) is not Matrix:
      raise Exception('Expected {b,c,x,l,s} to be of type Matrix')
    if   A.tag == sTag: lib.ElLinearProgramMPC_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramMPC_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(b) is not DistMatrix or type(c) is not DistMatrix or \
       type(x) is not DistMatrix or type(l) is not DistMatrix or \
       type(s) is not DistMatrix:
      raise Exception('Expected {b,c,x,l,s} to be of type Matrix')
    if   A.tag == sTag: lib.ElLinearProgramMPCDist_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramMPCDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElLinearProgramADMM_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElLinearProgramADMM_s.restype = c_uint
lib.ElLinearProgramADMM_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElLinearProgramADMM_d.restype = c_uint
lib.ElLinearProgramADMMDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElLinearProgramADMMDist_s.restype = c_uint
lib.ElLinearProgramADMMDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElLinearProgramADMMDist_d.restype = c_uint
def LinearProgramADMM(A,b,c):
  if type(A) is not type(b) or type(b) is not type(c):
    raise Exception('Types of {A,b,c} must match')
  if A.tag != b.tag or b.tag != c.tag:
    raise Exception('Datatypes of {A,b,c} must match')
  numIts = iType()
  if type(A) is Matrix:
    z = Matrix(A.tag)
    args = [A.obj,b.obj,c.obj,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElLinearProgram_s(*args)
    elif A.tag == dTag: lib.ElLinearProgram_d(*args)
    else: DataExcept()
    return z, numIts
  elif type(A) is DistMatrix:
    z = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,b.obj,c.obj,z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElLinearProgramDist_s(*args)
    elif A.tag == dTag: lib.ElLinearProgramDist_d(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Logistic regression
# ===================
(NO_PENALTY,L1_PENALTY,L2_PENALTY)=(0,1,2)

lib.ElLogisticRegression_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,sType,c_uint,POINTER(iType)]
lib.ElLogisticRegression_s.restype = c_uint
lib.ElLogisticRegression_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,dType,c_uint,POINTER(iType)]
lib.ElLogisticRegression_d.restype = c_uint
lib.ElLogisticRegressionDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,sType,c_uint,POINTER(iType)]
lib.ElLogisticRegressionDist_s.restype = c_uint
lib.ElLogisticRegressionDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,dType,c_uint,POINTER(iType)]
lib.ElLogisticRegressionDist_d.restype = c_uint
def LogisticRegression(G,q,gamma,penalty=L1_PENALTY):
  if type(G) is not type(q):
    raise Exception('Types of G and q must match')
  if G.tag != q.tag:
    raise Exception('Datatypes of G and q must match')
  numIts = iType()
  if type(G) is Matrix:
    z = Matrix(G.tag)
    args = [G.obj,q.obj,z.obj,gamma,penalty,pointer(numIts)]
    if   G.tag == sTag: lib.ElLogisticRegression_s(*args)
    elif G.tag == dTag: lib.ElLogisticRegression_d(*args)
    else: DataExcept()
    return z, numIts
  elif type(G) is DistMatrix:
    z = DistMatrix(G.tag,MC,MR,G.Grid())
    args = [G.obj,q.obj,z.obj,gamma,penalty,pointer(numIts)]
    if   G.tag == sTag: lib.ElLogisticRegressionDist_s(*args)
    elif G.tag == dTag: lib.ElLogisticRegressionDist_d(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Model fit
# =========
lib.ElModelFit_s.argtypes = \
  [CFUNCTYPE(None,c_void_p,sType),CFUNCTYPE(None,c_void_p,sType),
   c_void_p,c_void_p,c_void_p,sType,POINTER(iType)]
lib.ElModelFit_s.restype = c_uint
lib.ElModelFit_d.argtypes = \
  [CFUNCTYPE(None,c_void_p,dType),CFUNCTYPE(None,c_void_p,dType),
   c_void_p,c_void_p,c_void_p,dType,POINTER(iType)]
lib.ElModelFit_d.restype = c_uint
lib.ElModelFitDist_s.argtypes = \
  [CFUNCTYPE(None,c_void_p,sType),CFUNCTYPE(None,c_void_p,sType),
   c_void_p,c_void_p,c_void_p,sType,POINTER(iType)]
lib.ElModelFitDist_s.restype = c_uint
lib.ElModelFitDist_d.argtypes = \
  [CFUNCTYPE(None,c_void_p,dType),CFUNCTYPE(None,c_void_p,dType),
   c_void_p,c_void_p,c_void_p,dType,POINTER(iType)]
lib.ElModelFitDist_d.restype = c_uint
def ModelFit(lossProx,regProx,A,b,rho=1.2):
  if type(A) is not type(b):
    raise Exception('Types of A and b must match')
  if A.tag != b.tag:
    raise Exception('Datatypes of A and b must match')
  numIts = iType()
  cLoss = CFUNCTYPE(None,c_void_p,TagToType(A.tag))(lossProx)
  cReg = CFUNCTYPE(None,c_void_p,TagToType(A.tag))(regProx)
  if type(A) is Matrix:
    w = Matrix(A.tag)
    args = [cLoss,cReg,A.obj,b.obj,w.obj,rho,pointer(numIts)]
    if   A.tag == sTag: lib.ElModelFit_s(*args)
    elif A.tag == dTag: lib.ElModelFit_d(*args)
    else: DataExcept()
    return w, numIts
  elif type(A) is DistMatrix:
    w = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [cLoss,cReg,A.obj,b.obj,w.obj,rho,pointer(numIts)]
    if   A.tag == sTag: lib.ElModelFitDist_s(*args)
    elif A.tag == dTag: lib.ElModelFitDist_d(*args)
    else: DataExcept()
    return w, numIts
  else: TypeExcept()

# Non-negative matrix factorization
# =================================
lib.ElNMF_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNMF_s.restype = c_uint
lib.ElNMF_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNMF_d.restype = c_uint
lib.ElNMFDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNMFDist_s.restype = c_uint
lib.ElNMFDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElNMFDist_d.restype = c_uint
def NMF(A):
  args = [A.obj,X.obj,Y.obj]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElNMF_s(*args)
    elif A.tag == dTag: lib.ElNMF_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElNMFDist_s(*args)
    elif A.tag == dTag: lib.ElNMFDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Non-negative least squares
# ==========================
lib.ElNonNegativeLeastSquares_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElNonNegativeLeastSquares_s.restype = c_uint
lib.ElNonNegativeLeastSquares_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElNonNegativeLeastSquares_d.restype = c_uint
lib.ElNonNegativeLeastSquaresDist_s.argtypes = \
  [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElNonNegativeLeastSquaresDist_s.restype = c_uint
lib.ElNonNegativeLeastSquaresDist_d.argtypes = \
  [c_void_p,c_void_p,c_void_p,POINTER(iType)]
lib.ElNonNegativeLeastSquaresDist_d.restype = c_uint
def NonNegativeLeastSquares(A,Y):
  if type(A) is not type(Y):
    raise Exception('Types of A and Y must match')
  if A.tag != Y.tag:
    raise Exception('Datatypes of A and Y must match')
  numIts = iType()
  if type(A) is Matrix:
    Z = Matrix(A.tag)
    args = [A.obj,Y.obj,Z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElNonNegativeLeastSquares_s(*args)
    elif A.tag == dTag: lib.ElNonNegativeLeastSquares_d(*args)
    else: DataExcept()
    return Z, numIts
  elif type(A) is DistMatrix:
    Z = DistMatrix(A.tag,MC,MR,A.Grid())
    args = [A.obj,Y.obj,Z.obj,pointer(numIts)]
    if   A.tag == sTag: lib.ElNonNegativeLeastSquaresDist_s(*args)
    elif A.tag == dTag: lib.ElNonNegativeLeastSquaresDist_d(*args)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Quadratic program
# =================
lib.ElQuadraticProgram_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,POINTER(iType)]
lib.ElQuadraticProgram_s.restype = c_uint
lib.ElQuadraticProgram_d.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,POINTER(iType)]
lib.ElQuadraticProgram_d.restype = c_uint
lib.ElQuadraticProgramDist_s.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,POINTER(iType)]
lib.ElQuadraticProgramDist_s.restype = c_uint
lib.ElQuadraticProgramDist_d.argtypes = \
  [c_void_p,c_void_p,sType,sType,c_void_p,POINTER(iType)]
lib.ElQuadraticProgramDist_d.restype = c_uint
def QuadraticProgram(P,S,lb,ub):
  if type(P) is not type(S):
    raise Exception('Types of P and S must match')
  if P.tag != S.tag:
    raise Exception('Datatypes of P and S must match')
  numIts = iType()
  if type(P) is Matrix:
    Z = Matrix(P.tag)
    args = [P.obj,S.obj,lb,ub,Z.obj,pointer(numIts)]
    if   P.tag == sTag: lib.ElQuadraticProgram_s(*args)
    elif P.tag == dTag: lib.ElQuadraticProgram_d(*args)
    else: DataExcept()
    return Z, numIts
  elif type(P) is DistMatrix:
    Z = DistMatrix(P.tag,MC,MR,P.Grid())
    args = [P.obj,S.obj,lb,ub,Z.obj,pointer(numIts)]
    if   P.tag == sTag: lib.ElQuadraticProgramDist_s(*args)
    elif P.tag == dTag: lib.ElQuadraticProgramDist_d(*args)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Robust Principal Component Analysis
# ===================================
lib.ElRPCA_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCA_s.restype = c_uint
lib.ElRPCA_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCA_d.restype = c_uint
lib.ElRPCA_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCA_c.restype = c_uint
lib.ElRPCA_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCA_z.restype = c_uint
lib.ElRPCADist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCADist_s.restype = c_uint
lib.ElRPCADist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCADist_d.restype = c_uint
lib.ElRPCADist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCADist_c.restype = c_uint
lib.ElRPCADist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElRPCADist_z.restype = c_uint
def RPCA(M):
  if type(M) is Matrix:
    L = Matrix(M.tag)
    S = Matrix(M.tag)
    args = [M.obj,L.obj,S.obj]
    if   M.tag == sTag: lib.ElRPCA_s(*args)
    elif M.tag == dTag: lib.ElRPCA_d(*args)
    elif M.tag == cTag: lib.ElRPCA_c(*args)
    elif M.tag == zTag: lib.ElRPCA_z(*args)
    return L, S
  elif type(M) is DistMatrix:
    L = DistMatrix(M.tag,MC,MR,M.Grid())
    S = DistMatrix(M.tag,MC,MR,M.Grid())
    args = [M.obj,L.obj,S.obj]
    if   M.tag == sTag: lib.ElRPCADist_s(*args)
    elif M.tag == dTag: lib.ElRPCADist_d(*args)
    elif M.tag == cTag: lib.ElRPCADist_c(*args)
    elif M.tag == zTag: lib.ElRPCADist_z(*args)
    return L, S
  else: TypeExcept()

# Sparse inverse covariance selection
# ===================================
lib.ElSparseInvCov_s.argtypes = [c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElSparseInvCov_s.restype = c_uint
lib.ElSparseInvCov_d.argtypes = [c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElSparseInvCov_d.restype = c_uint
lib.ElSparseInvCov_c.argtypes = [c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElSparseInvCov_c.restype = c_uint
lib.ElSparseInvCov_z.argtypes = [c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElSparseInvCov_z.restype = c_uint
lib.ElSparseInvCovDist_s.argtypes = [c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElSparseInvCovDist_s.restype = c_uint
lib.ElSparseInvCovDist_d.argtypes = [c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElSparseInvCovDist_d.restype = c_uint
lib.ElSparseInvCovDist_c.argtypes = [c_void_p,sType,c_void_p,POINTER(iType)]
lib.ElSparseInvCovDist_c.restype = c_uint
lib.ElSparseInvCovDist_z.argtypes = [c_void_p,dType,c_void_p,POINTER(iType)]
lib.ElSparseInvCovDist_z.restype = c_uint
def SparseInvCov(D,lamb):
  numIts = iType()
  if type(D) is Matrix:
    Z = Matrix(D.tag)
    args = [D.obj,lamb,Z.obj,pointer(numIts)]
    if   D.tag == sTag: lib.ElSparseInvCov_s(*args)
    elif D.tag == dTag: lib.ElSparseInvCov_d(*args)
    elif D.tag == cTag: lib.ElSparseInvCov_c(*args)
    elif D.tag == zTag: lib.ElSparseInvCov_z(*args)
    else: DataExcept()
    return Z, numIts
  elif type(D) is DistMatrix:
    Z = DistMatrix(D.tag,MC,MR,D.Grid())
    args = [D.obj,lamb,Z.obj,pointer(numIts)]
    if   D.tag == sTag: lib.ElSparseInvCovDist_s(*args)
    elif D.tag == dTag: lib.ElSparseInvCovDist_d(*args)
    elif D.tag == cTag: lib.ElSparseInvCovDist_c(*args)
    elif D.tag == zTag: lib.ElSparseInvCovDist_z(*args)
    else: DataExcept()
    return Z, numIts
  else: TypeExcept()

# Support Vector Machine
# ======================
lib.ElSVM_s.argtypes = [c_void_p,c_void_p,c_void_p,sType,POINTER(iType)]
lib.ElSVM_s.restype = c_uint
lib.ElSVM_d.argtypes = [c_void_p,c_void_p,c_void_p,dType,POINTER(iType)]
lib.ElSVM_d.restype = c_uint
lib.ElSVMDist_s.argtypes = [c_void_p,c_void_p,c_void_p,sType,POINTER(iType)]
lib.ElSVMDist_s.restype = c_uint
lib.ElSVMDist_d.argtypes = [c_void_p,c_void_p,c_void_p,dType,POINTER(iType)]
lib.ElSVMDist_d.restype = c_uint
def SVM(G,q,gamma):
  if type(G) is not type(q):
    raise Exception('Types of G and q must match')
  if G.tag != q.tag:
    raise Exception('Datatypes of G and q must match')
  numIts = iType()
  if type(G) is Matrix:
    z = Matrix(G.tag)
    args = [G.obj,q.obj,z.obj,gamma,pointer(numIts)]
    if   G.tag == sTag: lib.ElSVM_s(*args)
    elif G.tag == dTag: lib.ElSVM_d(*args)
    else: DataExcept()
    return z, numIts
  elif type(G) is DistMatrix:
    z = DistMatrix(G.tag,MC,MR,G.Grid())
    args = [G.obj,q.obj,z.obj,gamma,pointer(numIts)]
    if   G.tag == sTag: lib.ElSVMDist_s(*args)
    elif G.tag == dTag: lib.ElSVMDist_d(*args)
    else: DataExcept()
    return z, numIts
  else: TypeExcept()

# Utilities
# =========

# Clipping
# --------
lib.ElLowerClip_s.argtypes = [c_void_p,sType]
lib.ElLowerClip_s.restype = c_uint
lib.ElLowerClip_d.argtypes = [c_void_p,dType]
lib.ElLowerClip_d.restype = c_uint
lib.ElLowerClipDist_s.argtypes = [c_void_p,sType]
lib.ElLowerClipDist_s.restype = c_uint
lib.ElLowerClipDist_d.argtypes = [c_void_p,dType]
lib.ElLowerClipDist_d.restype = c_uint
def LowerClip(X,lowerBound=0):
  args = [X.obj,lowerBound]
  if type(X) is Matrix:
    if   X.tag == sTag: lib.ElLowerClip_s(*args)
    elif X.tag == dTag: lib.ElLowerClip_d(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == sTag: lib.ElLowerClipDist_s(*args)
    elif X.tag == dTag: lib.ElLowerClipDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElUpperClip_s.argtypes = [c_void_p,sType]
lib.ElUpperClip_s.restype = c_uint
lib.ElUpperClip_d.argtypes = [c_void_p,dType]
lib.ElUpperClip_d.restype = c_uint
lib.ElUpperClipDist_s.argtypes = [c_void_p,sType]
lib.ElUpperClipDist_s.restype = c_uint
lib.ElUpperClipDist_d.argtypes = [c_void_p,dType]
lib.ElUpperClipDist_d.restype = c_uint
def UpperClip(X,upperBound=1):
  args = [X.obj,upperBound]
  if type(X) is Matrix:
    if   X.tag == sTag: lib.ElUpperClip_s(*args)
    elif X.tag == dTag: lib.ElUpperClip_d(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == sTag: lib.ElUpperClipDist_s(*args)
    elif X.tag == dTag: lib.ElUpperClipDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElClip_s.argtypes = [c_void_p,sType,sType]
lib.ElClip_s.restype = c_uint
lib.ElClip_d.argtypes = [c_void_p,dType,dType]
lib.ElClip_d.restype = c_uint
lib.ElClipDist_s.argtypes = [c_void_p,sType,sType]
lib.ElClipDist_s.restype = c_uint
lib.ElClipDist_d.argtypes = [c_void_p,dType,dType]
lib.ElClipDist_d.restype = c_uint
def Clip(X,lowerBound=0,upperBound=1):
  args = [X.obj,lowerBound,upperBound]
  if type(X) is Matrix:
    if   X.tag == sTag: lib.ElClip_s(*args)
    elif X.tag == dTag: lib.ElClip_d(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == sTag: lib.ElClipDist_s(*args)
    elif X.tag == dTag: lib.ElClipDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Coherence
# ---------
lib.ElCoherence_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElCoherence_s.restype = c_uint
lib.ElCoherence_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElCoherence_d.restype = c_uint
lib.ElCoherence_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElCoherence_c.restype = c_uint
lib.ElCoherence_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElCoherence_z.restype = c_uint
lib.ElCoherenceDist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElCoherenceDist_s.restype = c_uint
lib.ElCoherenceDist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElCoherenceDist_d.restype = c_uint
lib.ElCoherenceDist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElCoherenceDist_c.restype = c_uint
lib.ElCoherenceDist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElCoherenceDist_z.restype = c_uint
def Coherence(A):
  value = TagToType(Base(A.tag))()
  args = [A.obj,pointer(value)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElCoherence_s(*args)
    elif A.tag == dTag: lib.ElCoherence_d(*args)
    elif A.tag == cTag: lib.ElCoherence_c(*args)
    elif A.tag == zTag: lib.ElCoherence_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElCoherenceDist_s(*args)
    elif A.tag == dTag: lib.ElCoherenceDist_d(*args)
    elif A.tag == cTag: lib.ElCoherenceDist_c(*args)
    elif A.tag == zTag: lib.ElCoherenceDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return value

# Covariance
# ----------
lib.ElCovariance_s.argtypes = [c_void_p,c_void_p]
lib.ElCovariance_s.restype = c_uint
lib.ElCovariance_d.argtypes = [c_void_p,c_void_p]
lib.ElCovariance_d.restype = c_uint
lib.ElCovariance_c.argtypes = [c_void_p,c_void_p]
lib.ElCovariance_c.restype = c_uint
lib.ElCovariance_z.argtypes = [c_void_p,c_void_p]
lib.ElCovariance_z.restype = c_uint
lib.ElCovarianceDist_s.argtypes = [c_void_p,c_void_p]
lib.ElCovarianceDist_s.restype = c_uint
lib.ElCovarianceDist_d.argtypes = [c_void_p,c_void_p]
lib.ElCovarianceDist_d.restype = c_uint
lib.ElCovarianceDist_c.argtypes = [c_void_p,c_void_p]
lib.ElCovarianceDist_c.restype = c_uint
lib.ElCovarianceDist_z.argtypes = [c_void_p,c_void_p]
lib.ElCovarianceDist_z.restype = c_uint
def Covariance(D):
  if type(D) is Matrix:
    S = Matrix(D.tag)
    args = [D.obj,S.obj]
    if   D.tag == sTag: lib.ElCovariance_s(*args)
    elif D.tag == dTag: lib.ElCovariance_d(*args)
    elif D.tag == cTag: lib.ElCovariance_c(*args)
    elif D.tag == zTag: lib.ElCovariance_z(*args)
    else: DataExcept()
    return S
  elif type(D) is DistMatrix:
    S = DistMatrix(D.tag,MC,MR,D.Grid())
    args = [D.obj,S.obj]
    if   D.tag == sTag: lib.ElCovarianceDist_s(*args)
    elif D.tag == dTag: lib.ElCovarianceDist_d(*args)
    elif D.tag == cTag: lib.ElCovarianceDist_c(*args)
    elif D.tag == zTag: lib.ElCovarianceDist_z(*args)
    else: DataExcept()
    return S
  else: TypeExcept()

# Frobenius-norm proximal map
# ---------------------------
lib.ElFrobeniusProx_s.argtypes = [c_void_p,sType]
lib.ElFrobeniusProx_s.restype = c_uint
lib.ElFrobeniusProx_d.argtypes = [c_void_p,dType]
lib.ElFrobeniusProx_d.restype = c_uint
lib.ElFrobeniusProx_c.argtypes = [c_void_p,sType]
lib.ElFrobeniusProx_c.restype = c_uint
lib.ElFrobeniusProx_z.argtypes = [c_void_p,dType]
lib.ElFrobeniusProx_z.restype = c_uint
lib.ElFrobeniusProxDist_s.argtypes = [c_void_p,sType]
lib.ElFrobeniusProxDist_s.restype = c_uint
lib.ElFrobeniusProxDist_d.argtypes = [c_void_p,dType]
lib.ElFrobeniusProxDist_d.restype = c_uint
lib.ElFrobeniusProxDist_c.argtypes = [c_void_p,sType]
lib.ElFrobeniusProxDist_c.restype = c_uint
lib.ElFrobeniusProxDist_z.argtypes = [c_void_p,dType]
lib.ElFrobeniusProxDist_z.restype = c_uint
def FrobeniusProx(A,rho):
  args = [A.obj,rho]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElFrobeniusProx_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusProx_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusProx_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusProx_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElFrobeniusProxDist_s(*args)
    elif A.tag == dTag: lib.ElFrobeniusProxDist_d(*args)
    elif A.tag == cTag: lib.ElFrobeniusProxDist_c(*args)
    elif A.tag == zTag: lib.ElFrobeniusProxDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hinge-loss proximal map
# -----------------------
lib.ElHingeLossProx_s.argtypes = [c_void_p,sType]
lib.ElHingeLossProx_s.restype = c_uint
lib.ElHingeLossProx_d.argtypes = [c_void_p,dType]
lib.ElHingeLossProx_d.restype = c_uint
lib.ElHingeLossProxDist_s.argtypes = [c_void_p,sType]
lib.ElHingeLossProxDist_s.restype = c_uint
lib.ElHingeLossProxDist_d.argtypes = [c_void_p,dType]
lib.ElHingeLossProxDist_d.restype = c_uint
def HingeLossProx(A,rho):
  args = [A.obj,rho]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElHingeLossProx_s(*args)
    elif A.tag == dTag: lib.ElHingeLossProx_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElHingeLossProxDist_s(*args)
    elif A.tag == dTag: lib.ElHingeLossProxDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Log barrier
# -----------
lib.ElLogBarrier_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElLogBarrier_s.restype = c_uint
lib.ElLogBarrier_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElLogBarrier_d.restype = c_uint
lib.ElLogBarrier_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElLogBarrier_c.restype = c_uint
lib.ElLogBarrier_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElLogBarrier_z.restype = c_uint
lib.ElLogBarrierDist_s.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElLogBarrierDist_s.restype = c_uint
lib.ElLogBarrierDist_d.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElLogBarrierDist_d.restype = c_uint
lib.ElLogBarrierDist_c.argtypes = [c_uint,c_void_p,POINTER(sType)]
lib.ElLogBarrierDist_c.restype = c_uint
lib.ElLogBarrierDist_z.argtypes = [c_uint,c_void_p,POINTER(dType)]
lib.ElLogBarrierDist_z.restype = c_uint
def LogBarrier(uplo,A):
  barrier = TagToType(Base(A.tag))()
  args = [uplo,A.obj,pointer(barrier)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLogBarrier_s(*args)
    elif A.tag == dTag: lib.ElLogBarrier_d(*args)
    elif A.tag == cTag: lib.ElLogBarrier_c(*args)
    elif A.tag == zTag: lib.ElLogBarrier_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLogBarrierDist_s(*args)
    elif A.tag == dTag: lib.ElLogBarrierDist_d(*args)
    elif A.tag == cTag: lib.ElLogBarrierDist_c(*args)
    elif A.tag == zTag: lib.ElLogBarrierDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return barrier

# Log-det divergence
# ------------------
lib.ElLogDetDiv_s.argtypes = [c_uint,c_void_p,c_void_p,POINTER(sType)]
lib.ElLogDetDiv_s.restype = c_uint
lib.ElLogDetDiv_d.argtypes = [c_uint,c_void_p,c_void_p,POINTER(dType)]
lib.ElLogDetDiv_d.restype = c_uint
lib.ElLogDetDiv_c.argtypes = [c_uint,c_void_p,c_void_p,POINTER(sType)]
lib.ElLogDetDiv_c.restype = c_uint
lib.ElLogDetDiv_z.argtypes = [c_uint,c_void_p,c_void_p,POINTER(dType)]
lib.ElLogDetDiv_z.restype = c_uint
lib.ElLogDetDivDist_s.argtypes = [c_uint,c_void_p,c_void_p,POINTER(sType)]
lib.ElLogDetDivDist_s.restype = c_uint
lib.ElLogDetDivDist_d.argtypes = [c_uint,c_void_p,c_void_p,POINTER(dType)]
lib.ElLogDetDivDist_d.restype = c_uint
lib.ElLogDetDivDist_c.argtypes = [c_uint,c_void_p,c_void_p,POINTER(sType)]
lib.ElLogDetDivDist_c.restype = c_uint
lib.ElLogDetDivDist_z.argtypes = [c_uint,c_void_p,c_void_p,POINTER(dType)]
lib.ElLogDetDivDist_z.restype = c_uint
def LogDetDiv(uplo,A,B):
  div = TagToType(Base(A.tag))()
  args = [uplo,A.obj,B.obj,pointer(div)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLogDetDiv_s(*args)
    elif A.tag == dTag: lib.ElLogDetDiv_d(*args)
    elif A.tag == cTag: lib.ElLogDetDiv_c(*args)
    elif A.tag == zTag: lib.ElLogDetDiv_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLogDetDivDist_s(*args)
    elif A.tag == dTag: lib.ElLogDetDivDist_d(*args)
    elif A.tag == cTag: lib.ElLogDetDivDist_c(*args)
    elif A.tag == zTag: lib.ElLogDetDivDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return div

# Logistic proximal map
# ---------------------
lib.ElLogisticProx_s.argtypes = [c_void_p,sType]
lib.ElLogisticProx_s.restype = c_uint
lib.ElLogisticProx_d.argtypes = [c_void_p,dType]
lib.ElLogisticProx_d.restype = c_uint
lib.ElLogisticProxDist_s.argtypes = [c_void_p,sType]
lib.ElLogisticProxDist_s.restype = c_uint
lib.ElLogisticProxDist_d.argtypes = [c_void_p,dType]
lib.ElLogisticProxDist_d.restype = c_uint
def LogisticProx(A,rho):
  args = [A.obj,rho]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElLogisticProx_s(*args)
    elif A.tag == dTag: lib.ElLogisticProx_d(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElLogisticProxDist_s(*args)
    elif A.tag == dTag: lib.ElLogisticProxDist_d(*args)
    else: DataExcept()
  else: TypeExcept()

# Singular-value soft-thresholding
# --------------------------------
lib.ElSVT_s.argtypes = [c_void_p,sType,bType]
lib.ElSVT_s.restype = c_uint
lib.ElSVT_d.argtypes = [c_void_p,dType,bType]
lib.ElSVT_d.restype = c_uint
lib.ElSVT_c.argtypes = [c_void_p,sType,bType]
lib.ElSVT_c.restype = c_uint
lib.ElSVT_z.argtypes = [c_void_p,dType,bType]
lib.ElSVT_z.restype = c_uint
lib.ElSVTDist_s.argtypes = [c_void_p,sType,bType]
lib.ElSVTDist_s.restype = c_uint
lib.ElSVTDist_d.argtypes = [c_void_p,dType,bType]
lib.ElSVTDist_d.restype = c_uint
lib.ElSVTDist_c.argtypes = [c_void_p,sType,bType]
lib.ElSVTDist_c.restype = c_uint
lib.ElSVTDist_z.argtypes = [c_void_p,dType,bType]
lib.ElSVTDist_z.restype = c_uint
def SVT(A,rho,relative=False):
  args = [A.obj,rho,relative]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSVT_s(*args)
    elif A.tag == dTag: lib.ElSVT_d(*args)
    elif A.tag == cTag: lib.ElSVT_c(*args)
    elif A.tag == zTag: lib.ElSVT_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSVTDist_s(*args)
    elif A.tag == dTag: lib.ElSVTDist_d(*args)
    elif A.tag == cTag: lib.ElSVTDist_c(*args)
    elif A.tag == zTag: lib.ElSVTDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Soft-thresholding
# -----------------
lib.ElSoftThreshold_s.argtypes = [c_void_p,sType,bType]
lib.ElSoftThreshold_s.restype = c_uint
lib.ElSoftThreshold_d.argtypes = [c_void_p,dType,bType]
lib.ElSoftThreshold_d.restype = c_uint
lib.ElSoftThreshold_c.argtypes = [c_void_p,sType,bType]
lib.ElSoftThreshold_c.restype = c_uint
lib.ElSoftThreshold_z.argtypes = [c_void_p,dType,bType]
lib.ElSoftThreshold_z.restype = c_uint
lib.ElSoftThresholdDist_s.argtypes = [c_void_p,sType,bType]
lib.ElSoftThresholdDist_s.restype = c_uint
lib.ElSoftThresholdDist_d.argtypes = [c_void_p,dType,bType]
lib.ElSoftThresholdDist_d.restype = c_uint
lib.ElSoftThresholdDist_c.argtypes = [c_void_p,sType,bType]
lib.ElSoftThresholdDist_c.restype = c_uint
lib.ElSoftThresholdDist_z.argtypes = [c_void_p,dType,bType]
lib.ElSoftThresholdDist_z.restype = c_uint
def SoftThreshold(A,rho,relative=False):
  args = [A.obj,rho,relative]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElSoftThreshold_s(*args)
    elif A.tag == dTag: lib.ElSoftThreshold_d(*args)
    elif A.tag == cTag: lib.ElSoftThreshold_c(*args)
    elif A.tag == zTag: lib.ElSoftThreshold_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElSoftThresholdDist_s(*args)
    elif A.tag == dTag: lib.ElSoftThresholdDist_d(*args)
    elif A.tag == cTag: lib.ElSoftThresholdDist_c(*args)
    elif A.tag == zTag: lib.ElSoftThresholdDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
