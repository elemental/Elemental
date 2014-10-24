#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from El.core import *

lib.ElSymmetricSolveSparseDist_s.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveSparseDist_s.restype = c_uint
lib.ElSymmetricSolveSparseDist_d.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveSparseDist_d.restype = c_uint
lib.ElSymmetricSolveSparseDist_c.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveSparseDist_c.restype = c_uint
lib.ElSymmetricSolveSparseDist_z.argtypes = [c_void_p,c_void_p]
lib.ElSymmetricSolveSparseDist_z.restype = c_uint
lib.ElHermitianSolveSparseDist_c.argtypes = [c_void_p,c_void_p]
lib.ElHermitianSolveSparseDist_c.restype = c_uint
lib.ElHermitianSolveSparseDist_z.argtypes = [c_void_p,c_void_p]
lib.ElHermitianSolveSparseDist_z.restype = c_uint
# TODO: Combine with the preexisting SparseSolve routine
def SymmetricSolveSparse(A,X,conjugate=False):
  if type(A) is DistSparseMatrix:
    if type(X) is not DistMultiVec:
      TypeExcept()
    args = [A.obj,X.obj]
    if   A.tag == sTag: lib.ElSymmetricSolveSparseDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSolveSparseDist_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElHermitianSolveSparseDist_c(*args)
      else:         lib.ElSymmetricSolveSparseDist_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElHermitianSolveSparseDist_z(*args)
      else:         lib.ElSymmetricSolveSparseDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
def HermitianSolveSparse(A,X):
  SymmetricSolveSparse(A,X,True)
