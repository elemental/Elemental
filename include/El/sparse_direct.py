#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

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
# TODO: Combine with the preexisting SparseSolve routine
def SymmetricSolveSparse(A,X,conjugate=False):
  if type(A) is DistSparseMatrix:
    if type(X) is not DistMultiVec:
      raise Exception('X expected to be a DistMultiVec')
    args = [A.obj,X.obj]
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
def HermitianSolveSparse(A,X):
  SymmetricSolveSparse(A,X,True)

lib.ElLeastSquaresDistSparse_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_s.restype = c_uint
lib.ElLeastSquaresDistSparse_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_d.restype = c_uint
lib.ElLeastSquaresDistSparse_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_c.restype = c_uint
lib.ElLeastSquaresDistSparse_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElLeastSquaresDistSparse_z.restype = c_uint
def LeastSquares(A,Y):
  if type(A) is DistSparseMatrix:
    if type(Y) is not DistMultiVec:
      raise Exception("Expected Y to be a DistMultiVec")
    if A.tag != Y.tag:
      raise Exception("Expected datatypes of A and Y to match")
    X = DistMultiVec(A.tag)
    args = [A.obj,Y.obj,X.obj]
    if   A.tag == sTag: lib.ElLeastSquaresDistSparse_s(*args)
    elif A.tag == dTag: lib.ElLeastSquaresDistSparse_d(*args)
    elif A.tag == cTag: lib.ElLeastSquaresDistSparse_c(*args)
    elif A.tag == zTag: lib.ElLeastSquaresDistSparse_z(*args)
    else: DataExcept()
    return X
  else: TypeExcept()
