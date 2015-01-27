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

# Linear solvers
# ==============

# Linear solve
# ------------
lib.ElLinearSolve_s.argtypes = \
lib.ElLinearSolve_d.argtypes = \
lib.ElLinearSolve_c.argtypes = \
lib.ElLinearSolve_z.argtypes = \
lib.ElLinearSolveDist_s.argtypes = \
lib.ElLinearSolveDist_d.argtypes = \
lib.ElLinearSolveDist_c.argtypes = \
lib.ElLinearSolveDist_z.argtypes = \
lib.ElLinearSolveSparse_s.argtypes = \
lib.ElLinearSolveSparse_d.argtypes = \
lib.ElLinearSolveSparse_c.argtypes = \
lib.ElLinearSolveSparse_z.argtypes = \
lib.ElLinearSolveDistSparse_s.argtypes = \
lib.ElLinearSolveDistSparse_d.argtypes = \
lib.ElLinearSolveDistSparse_c.argtypes = \
lib.ElLinearSolveDistSparse_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElLinearSolve_s.restype = \
lib.ElLinearSolve_d.restype = \
lib.ElLinearSolve_c.restype = \
lib.ElLinearSolve_z.restype = \
lib.ElLinearSolveDist_s.restype = \
lib.ElLinearSolveDist_d.restype = \
lib.ElLinearSolveDist_c.restype = \
lib.ElLinearSolveDist_z.restype = \
lib.ElLinearSolveSparse_s.restype = \
lib.ElLinearSolveSparse_d.restype = \
lib.ElLinearSolveSparse_c.restype = \
lib.ElLinearSolveSparse_z.restype = \
lib.ElLinearSolveDistSparse_s.restype = \
lib.ElLinearSolveDistSparse_d.restype = \
lib.ElLinearSolveDistSparse_c.restype = \
lib.ElLinearSolveDistSparse_z.restype = \
  c_uint

def LinearSolve(A,B):
  args = [A.obj,B.obj]
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('Expected B to be a Matrix')
    if   A.tag == sTag: lib.ElLinearSolve_s(*args)
    elif A.tag == dTag: lib.ElLinearSolve_d(*args)
    elif A.tag == cTag: lib.ElLinearSolve_c(*args)
    elif A.tag == zTag: lib.ElLinearSolve_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(B) is not DistMatrix:
      raise Exception('Expected B to be a DistMatrix')
    if   A.tag == sTag: lib.ElLinearSolveDist_s(*args)
    elif A.tag == dTag: lib.ElLinearSolveDist_d(*args)
    elif A.tag == cTag: lib.ElLinearSolveDist_c(*args)
    elif A.tag == zTag: lib.ElLinearSolveDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if type(B) is not Matrix:
      raise Exception('Expected B to be a Matrix')
    if   A.tag == sTag: lib.ElLinearSolveSparse_s(*args)
    elif A.tag == dTag: lib.ElLinearSolveSparse_d(*args)
    elif A.tag == cTag: lib.ElLinearSolveSparse_c(*args)
    elif A.tag == zTag: lib.ElLinearSolveSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if type(B) is not DistMultiVec:
      raise Exception('Expected B to be a DistMultiVec')
    if   A.tag == sTag: lib.ElLinearSolveDistSparse_s(*args)
    elif A.tag == dTag: lib.ElLinearSolveDistSparse_d(*args)
    elif A.tag == cTag: lib.ElLinearSolveDistSparse_c(*args)
    elif A.tag == zTag: lib.ElLinearSolveDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Symmetric/Hermitian solves
# --------------------------
lib.ElSymmetricSolve_s.argtypes = \
lib.ElSymmetricSolve_d.argtypes = \
lib.ElSymmetricSolve_c.argtypes = \
lib.ElSymmetricSolve_z.argtypes = \
lib.ElSymmetricSolveDist_s.argtypes = \
lib.ElSymmetricSolveDist_d.argtypes = \
lib.ElSymmetricSolveDist_c.argtypes = \
lib.ElSymmetricSolveDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p]
lib.ElSymmetricSolveDistSparse_s.argtypes = \
lib.ElSymmetricSolveDistSparse_d.argtypes = \
lib.ElSymmetricSolveDistSparse_c.argtypes = \
lib.ElSymmetricSolveDistSparse_z.argtypes = \
lib.ElHermitianSolveDistSparse_c.argtypes = \
lib.ElHermitianSolveDistSparse_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElSymmetricSolve_s.restype = \
lib.ElSymmetricSolve_d.restype = \
lib.ElSymmetricSolve_c.restype = \
lib.ElSymmetricSolve_z.restype = \
lib.ElSymmetricSolveDist_s.restype = \
lib.ElSymmetricSolveDist_d.restype = \
lib.ElSymmetricSolveDist_c.restype = \
lib.ElSymmetricSolveDist_z.restype = \
lib.ElSymmetricSolveDistSparse_s.restype = \
lib.ElSymmetricSolveDistSparse_d.restype = \
lib.ElSymmetricSolveDistSparse_c.restype = \
lib.ElSymmetricSolveDistSparse_z.restype = \
lib.ElHermitianSolveDistSparse_c.restype = \
lib.ElHermitianSolveDistSparse_z.restype = \
  c_uint

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

def HermitianSolve(A,B,uplo,orient):
  SymmetricSolve(A,B,True,uplo,orient)

# Hermitian positive-definite solve
# ---------------------------------
lib.ElHPDSolve_s.argtypes = \
lib.ElHPDSolve_d.argtypes = \
lib.ElHPDSolve_c.argtypes = \
lib.ElHPDSolve_z.argtypes = \
lib.ElHPDSolveDist_s.argtypes = \
lib.ElHPDSolveDist_d.argtypes = \
lib.ElHPDSolveDist_c.argtypes = \
lib.ElHPDSolveDist_z.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p]

lib.ElHPDSolve_s.restype = \
lib.ElHPDSolve_d.restype = \
lib.ElHPDSolve_c.restype = \
lib.ElHPDSolve_z.restype = \
lib.ElHPDSolveDist_s.restype = \
lib.ElHPDSolveDist_d.restype = \
lib.ElHPDSolveDist_c.restype = \
lib.ElHPDSolveDist_z.restype = \
  c_uint

def HPDSolve(A,B,uplo=LOWER,orient=NORMAL):
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

# Multishift Hessenberg solve
# ---------------------------
lib.ElMultiShiftHessSolve_s.argtypes = \
lib.ElMultiShiftHessSolve_d.argtypes = \
lib.ElMultiShiftHessSolve_c.argtypes = \
lib.ElMultiShiftHessSolve_z.argtypes = \
lib.ElMultiShiftHessSolveDist_s.argtypes = \
lib.ElMultiShiftHessSolveDist_d.argtypes = \
lib.ElMultiShiftHessSolveDist_c.argtypes = \
lib.ElMultiShiftHessSolveDist_z.argtypes = \
  [c_uint,c_uint,zType,c_void_p,c_void_p,c_void_p]

lib.ElMultiShiftHessSolve_s.restype = \
lib.ElMultiShiftHessSolve_d.restype = \
lib.ElMultiShiftHessSolve_c.restype = \
lib.ElMultiShiftHessSolve_z.restype = \
lib.ElMultiShiftHessSolveDist_s.restype = \
lib.ElMultiShiftHessSolveDist_d.restype = \
lib.ElMultiShiftHessSolveDist_c.restype = \
lib.ElMultiShiftHessSolveDist_z.restype = \
  c_uint

def MultiShiftHessSolve(H,shifts,X,alphaPre=1,uplo=LOWER,orient=NORMAL):
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
