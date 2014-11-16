#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from ..core import *

from ctypes import CFUNCTYPE

# BLAS 1
# ======

# Axpy
# ----
lib.ElAxpy_i.argtypes = [iType,c_void_p,c_void_p]
lib.ElAxpy_i.restype = c_uint
lib.ElAxpy_s.argtypes = [sType,c_void_p,c_void_p]
lib.ElAxpy_s.restype = c_uint
lib.ElAxpy_d.argtypes = [dType,c_void_p,c_void_p]
lib.ElAxpy_d.restype = c_uint
lib.ElAxpy_c.argtypes = [cType,c_void_p,c_void_p]
lib.ElAxpy_c.restype = c_uint
lib.ElAxpy_z.argtypes = [zType,c_void_p,c_void_p]
lib.ElAxpy_z.restype = c_uint
lib.ElAxpyDist_i.argtypes = [iType,c_void_p,c_void_p]
lib.ElAxpyDist_i.restype = c_uint
lib.ElAxpyDist_s.argtypes = [sType,c_void_p,c_void_p]
lib.ElAxpyDist_s.restype = c_uint
lib.ElAxpyDist_d.argtypes = [dType,c_void_p,c_void_p]
lib.ElAxpyDist_d.restype = c_uint
lib.ElAxpyDist_c.argtypes = [cType,c_void_p,c_void_p]
lib.ElAxpyDist_c.restype = c_uint
lib.ElAxpyDist_z.argtypes = [zType,c_void_p,c_void_p]
lib.ElAxpyDist_z.restype = c_uint
lib.ElAxpySparse_i.argtypes = [iType,c_void_p,c_void_p]
lib.ElAxpySparse_i.restype = c_uint
lib.ElAxpySparse_s.argtypes = [sType,c_void_p,c_void_p]
lib.ElAxpySparse_s.restype = c_uint
lib.ElAxpySparse_d.argtypes = [dType,c_void_p,c_void_p]
lib.ElAxpySparse_d.restype = c_uint
lib.ElAxpySparse_c.argtypes = [cType,c_void_p,c_void_p]
lib.ElAxpySparse_c.restype = c_uint
lib.ElAxpySparse_z.argtypes = [zType,c_void_p,c_void_p]
lib.ElAxpySparse_z.restype = c_uint
lib.ElAxpyDistSparse_i.argtypes = [iType,c_void_p,c_void_p]
lib.ElAxpyDistSparse_i.restype = c_uint
lib.ElAxpyDistSparse_s.argtypes = [sType,c_void_p,c_void_p]
lib.ElAxpyDistSparse_s.restype = c_uint
lib.ElAxpyDistSparse_d.argtypes = [dType,c_void_p,c_void_p]
lib.ElAxpyDistSparse_d.restype = c_uint
lib.ElAxpyDistSparse_c.argtypes = [cType,c_void_p,c_void_p]
lib.ElAxpyDistSparse_c.restype = c_uint
lib.ElAxpyDistSparse_z.argtypes = [zType,c_void_p,c_void_p]
lib.ElAxpyDistSparse_z.restype = c_uint
def Axpy(alphaPre,X,Y):
  if type(X) is not type(Y): raise Exception('Types of X and Y must match')
  if X.tag != Y.tag: raise Exception('Datatypes of X and Y must match')
  alpha = TagToType(X.tag)(alphaPre)
  args = [alpha,X.obj,Y.obj]
  if type(X) is Matrix:
    if   X.tag == iTag: lib.ElAxpy_i(*args)
    elif X.tag == sTag: lib.ElAxpy_s(*args)
    elif X.tag == dTag: lib.ElAxpy_d(*args)
    elif X.tag == cTag: lib.ElAxpy_c(*args)
    elif X.tag == zTag: lib.ElAxpy_z(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == iTag: lib.ElAxpyDist_i(*args)
    elif X.tag == sTag: lib.ElAxpyDist_s(*args)
    elif X.tag == dTag: lib.ElAxpyDist_d(*args)
    elif X.tag == cTag: lib.ElAxpyDist_c(*args)
    elif X.tag == zTag: lib.ElAxpyDist_z(*args)
    else: DataExcept()
  elif type(X) is SparseMatrix:
    if   X.tag == iTag: lib.ElAxpySparse_i(*args)
    elif X.tag == sTag: lib.ElAxpySparse_s(*args)
    elif X.tag == dTag: lib.ElAxpySparse_d(*args)
    elif X.tag == cTag: lib.ElAxpySparse_c(*args)
    elif X.tag == zTag: lib.ElAxpySparse_z(*args)
    else: DataExcept()
  elif type(X) is DistSparseMatrix:
    if   X.tag == iTag: lib.ElAxpyDistSparse_i(*args)
    elif X.tag == sTag: lib.ElAxpyDistSparse_s(*args)
    elif X.tag == dTag: lib.ElAxpyDistSparse_d(*args)
    elif X.tag == cTag: lib.ElAxpyDistSparse_c(*args)
    elif X.tag == zTag: lib.ElAxpyDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# AxpyTrapezoid
# -------------
lib.ElAxpyTrapezoid_i.argtypes = [c_uint,iType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoid_i.restype = c_uint
lib.ElAxpyTrapezoid_s.argtypes = [c_uint,sType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoid_s.restype = c_uint
lib.ElAxpyTrapezoid_d.argtypes = [c_uint,dType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoid_d.restype = c_uint
lib.ElAxpyTrapezoid_c.argtypes = [c_uint,cType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoid_c.restype = c_uint
lib.ElAxpyTrapezoid_z.argtypes = [c_uint,zType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoid_z.restype = c_uint
lib.ElAxpyTrapezoidDist_i.argtypes = [c_uint,iType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDist_i.restype = c_uint
lib.ElAxpyTrapezoidDist_s.argtypes = [c_uint,sType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDist_s.restype = c_uint
lib.ElAxpyTrapezoidDist_d.argtypes = [c_uint,dType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDist_d.restype = c_uint
lib.ElAxpyTrapezoidDist_c.argtypes = [c_uint,cType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDist_c.restype = c_uint
lib.ElAxpyTrapezoidDist_z.argtypes = [c_uint,zType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDist_z.restype = c_uint
lib.ElAxpyTrapezoidSparse_i.argtypes = [c_uint,iType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidSparse_i.restype = c_uint
lib.ElAxpyTrapezoidSparse_s.argtypes = [c_uint,sType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidSparse_s.restype = c_uint
lib.ElAxpyTrapezoidSparse_d.argtypes = [c_uint,dType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidSparse_d.restype = c_uint
lib.ElAxpyTrapezoidSparse_c.argtypes = [c_uint,cType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidSparse_c.restype = c_uint
lib.ElAxpyTrapezoidSparse_z.argtypes = [c_uint,zType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidSparse_z.restype = c_uint
lib.ElAxpyTrapezoidDistSparse_i.argtypes = [c_uint,iType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDistSparse_i.restype = c_uint
lib.ElAxpyTrapezoidDistSparse_s.argtypes = [c_uint,sType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDistSparse_s.restype = c_uint
lib.ElAxpyTrapezoidDistSparse_d.argtypes = [c_uint,dType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDistSparse_d.restype = c_uint
lib.ElAxpyTrapezoidDistSparse_c.argtypes = [c_uint,cType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDistSparse_c.restype = c_uint
lib.ElAxpyTrapezoidDistSparse_z.argtypes = [c_uint,zType,c_void_p,c_void_p,iType]
lib.ElAxpyTrapezoidDistSparse_z.restype = c_uint
def AxpyTriangle(uplo,alphaPre,X,Y,offset=0):
  if type(X) is not type(Y): raise Exception('Types of X and Y must match')
  if X.tag != Y.tag: raise Exception('Datatypes of X and Y must match')
  alpha = TagToType(X.tag)(alphaPre)
  args = [uplo,alpha,X.obj,Y.obj,offset]
  if type(X) is Matrix:
    if   X.tag == iTag: lib.ElAxpyTriangle_i(*args)
    elif X.tag == sTag: lib.ElAxpyTriangle_s(*args)
    elif X.tag == dTag: lib.ElAxpyTriangle_d(*args)
    elif X.tag == cTag: lib.ElAxpyTriangle_c(*args)
    elif X.tag == zTag: lib.ElAxpyTriangle_z(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == iTag: lib.ElAxpyTriangleDist_i(*args)
    elif X.tag == sTag: lib.ElAxpyTriangleDist_s(*args)
    elif X.tag == dTag: lib.ElAxpyTriangleDist_d(*args)
    elif X.tag == cTag: lib.ElAxpyTriangleDist_c(*args)
    elif X.tag == zTag: lib.ElAxpyTriangleDist_z(*args)
    else: DataExcept()
  elif type(X) is SparseMatrix:
    if   X.tag == iTag: lib.ElAxpyTriangleSparse_i(*args)
    elif X.tag == sTag: lib.ElAxpyTriangleSparse_s(*args)
    elif X.tag == dTag: lib.ElAxpyTriangleSparse_d(*args)
    elif X.tag == cTag: lib.ElAxpyTriangleSparse_c(*args)
    elif X.tag == zTag: lib.ElAxpyTriangleSparse_z(*args)
    else: DataExcept()
  elif type(X) is DistSparseMatrix:
    if   X.tag == iTag: lib.ElAxpyTriangleDistSparse_i(*args)
    elif X.tag == sTag: lib.ElAxpyTriangleDistSparse_s(*args)
    elif X.tag == dTag: lib.ElAxpyTriangleDistSparse_d(*args)
    elif X.tag == cTag: lib.ElAxpyTriangleDistSparse_c(*args)
    elif X.tag == zTag: lib.ElAxpyTriangleDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Column norms
# ------------
lib.ElColumnNormsDistMultiVec_s.argtypes = [c_void_p,c_void_p]
lib.ElColumnNormsDistMultiVec_s.restype = c_uint
lib.ElColumnNormsDistMultiVec_d.argtypes = [c_void_p,c_void_p]
lib.ElColumnNormsDistMultiVec_d.restype = c_uint
lib.ElColumnNormsDistMultiVec_c.argtypes = [c_void_p,c_void_p]
lib.ElColumnNormsDistMultiVec_c.restype = c_uint
lib.ElColumnNormsDistMultiVec_z.argtypes = [c_void_p,c_void_p]
lib.ElColumnNormsDistMultiVec_z.restype = c_uint
def ColumnNorms(A):
  if type(A) is DistMultiVec:
    norms = Matrix(TagToType(Base(A.tag)))
    args = [A.obj,norms.obj]
    if   A.tag == sTag: lib.ElColumnNormsDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElColumnNormsDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElColumnNormsDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElColumnNormsDistMultiVec_z(*args)
    else: DataExcept()
    return norms
  else: TypeExcept()

# Conjugate
# ---------
lib.ElConjugate_c.argtypes = [c_void_p]
lib.ElConjugate_c.restype = c_uint
lib.ElConjugate_z.argtypes = [c_void_p]
lib.ElConjugate_z.restype = c_uint
lib.ElConjugateDist_c.argtypes = [c_void_p]
lib.ElConjugateDist_c.restype = c_uint
lib.ElConjugateDist_z.argtypes = [c_void_p]
lib.ElConjugateDist_z.restype = c_uint
def Conjugate(A):
  args = [A.obj]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElConjugate_c(*args)
    elif A.tag == zTag: lib.ElConjugate_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElConjugateDist_c(*args)
    elif A.tag == zTag: lib.ElConjugateDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Copy
# ----
lib.ElCopy_s.argtypes = [c_void_p,c_void_p]
lib.ElCopy_s.restype = c_uint
lib.ElCopy_d.argtypes = [c_void_p,c_void_p]
lib.ElCopy_d.restype = c_uint
lib.ElCopy_c.argtypes = [c_void_p,c_void_p]
lib.ElCopy_c.restype = c_uint
lib.ElCopy_z.argtypes = [c_void_p,c_void_p]
lib.ElCopy_z.restype = c_uint
lib.ElCopyDist_s.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_s.restype = c_uint
lib.ElCopyDist_d.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_d.restype = c_uint
lib.ElCopyDist_c.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_c.restype = c_uint
lib.ElCopyDist_z.argtypes = [c_void_p,c_void_p]
lib.ElCopyDist_z.restype = c_uint
lib.ElCopyGraph.argtypes = [c_void_p,c_void_p]
lib.ElCopyGraph.restype = c_uint
lib.ElCopyDistGraph.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistGraph.restype = c_uint
lib.ElCopySparse_i.argtypes = [c_void_p,c_void_p]
lib.ElCopySparse_i.restype = c_uint
lib.ElCopySparse_s.argtypes = [c_void_p,c_void_p]
lib.ElCopySparse_s.restype = c_uint
lib.ElCopySparse_d.argtypes = [c_void_p,c_void_p]
lib.ElCopySparse_d.restype = c_uint
lib.ElCopySparse_c.argtypes = [c_void_p,c_void_p]
lib.ElCopySparse_c.restype = c_uint
lib.ElCopySparse_z.argtypes = [c_void_p,c_void_p]
lib.ElCopySparse_z.restype = c_uint
lib.ElCopySparseToDense_i.argtypes = [c_void_p,c_void_p]
lib.ElCopySparseToDense_i.restype = c_uint
lib.ElCopySparseToDense_s.argtypes = [c_void_p,c_void_p]
lib.ElCopySparseToDense_s.restype = c_uint
lib.ElCopySparseToDense_d.argtypes = [c_void_p,c_void_p]
lib.ElCopySparseToDense_d.restype = c_uint
lib.ElCopySparseToDense_c.argtypes = [c_void_p,c_void_p]
lib.ElCopySparseToDense_c.restype = c_uint
lib.ElCopySparseToDense_z.argtypes = [c_void_p,c_void_p]
lib.ElCopySparseToDense_z.restype = c_uint
lib.ElCopyDistSparse_i.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparse_i.restype = c_uint
lib.ElCopyDistSparse_s.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparse_s.restype = c_uint
lib.ElCopyDistSparse_d.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparse_d.restype = c_uint
lib.ElCopyDistSparse_c.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparse_c.restype = c_uint
lib.ElCopyDistSparse_z.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparse_z.restype = c_uint
lib.ElCopyDistSparseToDense_i.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparseToDense_i.restype = c_uint
lib.ElCopyDistSparseToDense_s.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparseToDense_s.restype = c_uint
lib.ElCopyDistSparseToDense_d.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparseToDense_d.restype = c_uint
lib.ElCopyDistSparseToDense_c.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparseToDense_c.restype = c_uint
lib.ElCopyDistSparseToDense_z.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistSparseToDense_z.restype = c_uint
lib.ElCopyDistMultiVec_i.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistMultiVec_i.restype = c_uint
lib.ElCopyDistMultiVec_s.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistMultiVec_s.restype = c_uint
lib.ElCopyDistMultiVec_d.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistMultiVec_d.restype = c_uint
lib.ElCopyDistMultiVec_c.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistMultiVec_c.restype = c_uint
lib.ElCopyDistMultiVec_z.argtypes = [c_void_p,c_void_p]
lib.ElCopyDistMultiVec_z.restype = c_uint
def Copy(A,B):
  if A.tag != B.tag:
    raise Exception('Copying between datatypes is not yet supported in Python')
  args = [A.obj,B.obj]
  if type(A) is Matrix:
    if type(B) is not Matrix:
      raise Exception('Expected B to be a Matrix')
    if   B.tag == iTag: lib.ElCopy_i(*args)
    elif B.tag == sTag: lib.ElCopy_s(*args)
    elif B.tag == dTag: lib.ElCopy_d(*args)
    elif B.tag == cTag: lib.ElCopy_c(*args)
    elif B.tag == zTag: lib.ElCopy_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if type(B) is not DistMatrix:
      raise Exception('Expected B to be a DistMatrix')
    if   B.tag == iTag: lib.ElCopyDist_i(*args)
    elif B.tag == sTag: lib.ElCopyDist_s(*args)
    elif B.tag == dTag: lib.ElCopyDist_d(*args)
    elif B.tag == cTag: lib.ElCopyDist_c(*args)
    elif B.tag == zTag: lib.ElCopyDist_z(*args)
    else: DataExcept()
  elif type(A) is Graph:
    if type(B) is not Graph:
      raise Exception('Expected B to be a Graph')
    lib.ElCopyGraph(*args)
  elif type(A) is DistGraph:
    if type(B) is not DistGraph:
      raise Exception('Expected B to be a DistGraph')
    lib.ElCopyDistGraph(*args)
  elif type(A) is SparseMatrix:
    if type(B) is SparseMatrix:
      if   A.tag == iTag: lib.ElCopySparse_i(*args)
      elif A.tag == sTag: lib.ElCopySparse_s(*args)
      elif A.tag == dTag: lib.ElCopySparse_d(*args)
      elif A.tag == cTag: lib.ElCopySparse_c(*args)
      elif A.tag == zTag: lib.ElCopySparse_z(*args)
      else: DataExcept()
    elif type(B) is Matrix:
      if   A.tag == iTag: lib.ElCopySparseToDense_i(*args)
      elif A.tag == sTag: lib.ElCopySparseToDense_s(*args)
      elif A.tag == dTag: lib.ElCopySparseToDense_d(*args)
      elif A.tag == cTag: lib.ElCopySparseToDense_c(*args)
      elif A.tag == zTag: lib.ElCopySparseToDense_z(*args)
      else: DataExcept()
    else:
      raise Exception('Expected B to be a (Sparse)Matrix')
  elif type(A) is DistSparseMatrix:
    if type(B) is DistSparseMatrix:
      if   A.tag == iTag: lib.ElCopyDistSparse_i(*args)
      elif A.tag == sTag: lib.ElCopyDistSparse_s(*args)
      elif A.tag == dTag: lib.ElCopyDistSparse_d(*args)
      elif A.tag == cTag: lib.ElCopyDistSparse_c(*args)
      elif A.tag == zTag: lib.ElCopyDistSparse_z(*args)
      else: DataExcept()
    elif type(B) is DistMatrix:
      if   A.tag == iTag: lib.ElCopyDistSparseToDense_i(*args)
      elif A.tag == sTag: lib.ElCopyDistSparseToDense_s(*args)
      elif A.tag == dTag: lib.ElCopyDistSparseToDense_d(*args)
      elif A.tag == cTag: lib.ElCopyDistSparseToDense_c(*args)
      elif A.tag == zTag: lib.ElCopyDistSparseToDense_z(*args)
      else: DataExcept()
    else:
      raise Exception('Expected B to be a Dist(Sparse)Matrix')
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElCopyDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElCopyDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElCopyDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElCopyDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElCopyDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElCopyGraphFromRoot.argtypes = [c_void_p,c_void_p]
lib.ElCopyGraphFromRoot.restype = c_uint
lib.ElCopyMultiVecFromRoot_i.argtypes = [c_void_p,c_void_p]
lib.ElCopyMultiVecFromRoot_i.restype = c_uint
lib.ElCopyMultiVecFromRoot_s.argtypes = [c_void_p,c_void_p]
lib.ElCopyMultiVecFromRoot_s.restype = c_uint
lib.ElCopyMultiVecFromRoot_d.argtypes = [c_void_p,c_void_p]
lib.ElCopyMultiVecFromRoot_d.restype = c_uint
lib.ElCopyMultiVecFromRoot_c.argtypes = [c_void_p,c_void_p]
lib.ElCopyMultiVecFromRoot_c.restype = c_uint
lib.ElCopyMultiVecFromRoot_z.argtypes = [c_void_p,c_void_p]
lib.ElCopyMultiVecFromRoot_z.restype = c_uint
def CopyFromRoot(ADist,ASeq):
  args = [ADist.obj,ASeq.obj]
  if type(ADist) is DistGraph:
    if type(ASeq) is not Graph:
      raise Exception("Expected the result to be a Graph")
    lib.ElCopyGraphFromRoot(*args)
  elif type(ADist) is DistSparseMatrix:
    if type(ASeq) is not SparseMatrix:
      raise Exception("Expected the result to be a SparseMatrix")
    if ASeq.tag != ADist.tag:
      raise Exception("Expected the result to be of the same type")
    if   ADist.tag == iTag: lib.ElCopySparseMatrixFromRoot_i(*args)
    elif ADist.tag == sTag: lib.ElCopySparseMatrixFromRoot_s(*args)
    elif ADist.tag == dTag: lib.ElCopySparseMatrixFromRoot_d(*args)
    elif ADist.tag == cTag: lib.ElCopySparseMatrixFromRoot_c(*args)
    elif ADist.tag == zTag: lib.ElCopySparseMatrixFromRoot_z(*args)
    else: DataExcept()
  elif type(ADist) is DistMultiVec:
    if type(ASeq) is not Matrix:
      raise Exception("Expected the result to be a Matrix")
    if ASeq.tag != ADist.tag:
      raise Exception("Expected the result to be of the same type")
    if   ADist.tag == iTag: lib.ElCopyMultiVecFromRoot_i(*args)
    elif ADist.tag == sTag: lib.ElCopyMultiVecFromRoot_s(*args)
    elif ADist.tag == dTag: lib.ElCopyMultiVecFromRoot_d(*args)
    elif ADist.tag == cTag: lib.ElCopyMultiVecFromRoot_c(*args)
    elif ADist.tag == zTag: lib.ElCopyMultiVecFromRoot_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElCopyMultiVecFromNonRoot_i.argtypes = [c_void_p,iType]
lib.ElCopyMultiVecFromNonRoot_i.restype = c_uint
lib.ElCopyMultiVecFromNonRoot_s.argtypes = [c_void_p,iType]
lib.ElCopyMultiVecFromNonRoot_s.restype = c_uint
lib.ElCopyMultiVecFromNonRoot_d.argtypes = [c_void_p,iType]
lib.ElCopyMultiVecFromNonRoot_d.restype = c_uint
lib.ElCopyMultiVecFromNonRoot_c.argtypes = [c_void_p,iType]
lib.ElCopyMultiVecFromNonRoot_c.restype = c_uint
lib.ElCopyMultiVecFromNonRoot_z.argtypes = [c_void_p,iType]
lib.ElCopyMultiVecFromNonRoot_z.restype = c_uint
lib.ElCopyGraphFromNonRoot.argtypes = [c_void_p,iType]
lib.ElCopyGraphFromNonRoot.restype = c_uint
def CopyFromNonRoot(ADist,root=0):
  args = [ADist.obj,root]
  if type(ADist) is DistGraph:
    lib.ElCopyGraphFromNonRoot(*args)
  elif type(ADist) is DistSparseMatrix:
    if   ADist.tag == iTag: lib.ElCopySparseMatrixFromNonRoot_i(*args)
    elif ADist.tag == sTag: lib.ElCopySparseMatrixFromNonRoot_s(*args)
    elif ADist.tag == dTag: lib.ElCopySparseMatrixFromNonRoot_d(*args)
    elif ADist.tag == cTag: lib.ElCopySparseMatrixFromNonRoot_c(*args)
    elif ADist.tag == zTag: lib.ElCopySparseMatrixFromNonRoot_z(*args)
    else: DataExcept()
  elif type(ADist) is DistMultiVec:
    if   ADist.tag == iTag: lib.ElCopyMultiVecFromNonRoot_i(*args)
    elif ADist.tag == sTag: lib.ElCopyMultiVecFromNonRoot_s(*args)
    elif ADist.tag == dTag: lib.ElCopyMultiVecFromNonRoot_d(*args)
    elif ADist.tag == cTag: lib.ElCopyMultiVecFromNonRoot_c(*args)
    elif ADist.tag == zTag: lib.ElCopyMultiVecFromNonRoot_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Diagonal scale
# --------------
lib.ElDiagonalScale_i.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalScale_i.restype = c_uint
lib.ElDiagonalScale_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalScale_s.restype = c_uint
lib.ElDiagonalScale_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalScale_d.restype = c_uint
lib.ElDiagonalScale_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalScale_c.restype = c_uint
lib.ElDiagonalScale_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalScale_z.restype = c_uint
lib.ElDiagonalScaleDist_i.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalScaleDist_i.restype = c_uint
lib.ElDiagonalScaleDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalScaleDist_s.restype = c_uint
lib.ElDiagonalScaleDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalScaleDist_d.restype = c_uint
lib.ElDiagonalScaleDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalScaleDist_c.restype = c_uint
lib.ElDiagonalScaleDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalScaleDist_z.restype = c_uint
def DiagonalScale(side,orient,d,X):
  if type(d) is not type(X): raise Exception('Matrix types must match')
  if d.tag != X.tag: raise Exception('Matrix datatypes must match')
  args = [side,d.obj,X.obj]
  argsCpx = [side,orient,d.obj,X.obj]
  if type(X) is Matrix:
    if   X.tag == iTag: lib.ElDiagonalScale_i(*args)
    elif X.tag == sTag: lib.ElDiagonalScale_s(*args)
    elif X.tag == dTag: lib.ElDiagonalScale_d(*args)
    elif X.tag == cTag: lib.ElDiagonalScale_c(*argsCpx)
    elif X.tag == zTag: lib.ElDiagonalScale_z(*argsCpx)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == iTag: lib.ElDiagonalScaleDist_i(*args)
    elif X.tag == sTag: lib.ElDiagonalScaleDist_s(*args)
    elif X.tag == dTag: lib.ElDiagonalScaleDist_d(*args)
    elif X.tag == cTag: lib.ElDiagonalScaleDist_c(*argsCpx)
    elif X.tag == zTag: lib.ElDiagonalScaleDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# Diagonal scale trapezoid
# ------------------------
lib.ElDiagonalScaleTrapezoid_i.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoid_i.restype = c_uint
lib.ElDiagonalScaleTrapezoid_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoid_s.restype = c_uint
lib.ElDiagonalScaleTrapezoid_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoid_d.restype = c_uint
lib.ElDiagonalScaleTrapezoid_c.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoid_c.restype = c_uint
lib.ElDiagonalScaleTrapezoid_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoid_z.restype = c_uint
lib.ElDiagonalScaleTrapezoidDist_i.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoidDist_i.restype = c_uint
lib.ElDiagonalScaleTrapezoidDist_s.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoidDist_s.restype = c_uint
lib.ElDiagonalScaleTrapezoidDist_d.argtypes = \
  [c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoidDist_d.restype = c_uint
lib.ElDiagonalScaleTrapezoidDist_c.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoidDist_c.restype = c_uint
lib.ElDiagonalScaleTrapezoidDist_z.argtypes = \
  [c_uint,c_uint,c_uint,c_void_p,c_void_p,iType]
lib.ElDiagonalScaleTrapezoidDist_z.restype = c_uint
def DiagonalScaleTrapezoid(side,uplo,orient,d,X,offset=0):
  if type(d) is not type(X): raise Exception('Matrix types must match')
  if d.tag != X.tag: raise Exception('Matrix datatypes must match')
  args = [side,uplo,d.obj,X.obj,offset]
  argsCpx = [side,uplo,orient,d.obj,X.obj,offset]
  if type(X) is Matrix:
    if   X.tag == iTag: lib.ElDiagonalScaleTrapezoid_i(*args)
    elif X.tag == sTag: lib.ElDiagonalScaleTrapezoid_s(*args)
    elif X.tag == dTag: lib.ElDiagonalScaleTrapezoid_d(*args)
    elif X.tag == cTag: lib.ElDiagonalScaleTrapezoid_c(*argsCpx)
    elif X.tag == zTag: lib.ElDiagonalScaleTrapezoid_z(*argsCpx)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == iTag: lib.ElDiagonalScaleTrapezoidDist_i(*args)
    elif X.tag == sTag: lib.ElDiagonalScaleTrapezoidDist_s(*args)
    elif X.tag == dTag: lib.ElDiagonalScaleTrapezoidDist_d(*args)
    elif X.tag == cTag: lib.ElDiagonalScaleTrapezoidDist_c(*argsCpx)
    elif X.tag == zTag: lib.ElDiagonalScaleTrapezoidDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# Diagonal solve
# --------------
lib.ElDiagonalSolve_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolve_s.restype = c_uint
lib.ElDiagonalSolve_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolve_d.restype = c_uint
lib.ElDiagonalSolve_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolve_c.restype = c_uint
lib.ElDiagonalSolve_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolve_z.restype = c_uint
lib.ElDiagonalSolveDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolveDist_s.restype = c_uint
lib.ElDiagonalSolveDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolveDist_d.restype = c_uint
lib.ElDiagonalSolveDist_c.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolveDist_c.restype = c_uint
lib.ElDiagonalSolveDist_z.argtypes = [c_uint,c_uint,c_void_p,c_void_p]
lib.ElDiagonalSolveDist_z.restype = c_uint
def DiagonalSolve(side,orient,d,X):
  if type(d) is not type(X): raise Exception('Matrix types must match')
  if d.tag != X.tag: raise Exception('Matrix datatypes must match')
  args = [side,d.obj,X.obj]
  argsCpx = [side,orient,d.obj,X.obj]
  if type(X) is Matrix:
    if   X.tag == sTag: lib.ElDiagonalSolve_s(*args)
    elif X.tag == dTag: lib.ElDiagonalSolve_d(*args)
    elif X.tag == cTag: lib.ElDiagonalSolve_c(*argsCpx)
    elif X.tag == zTag: lib.ElDiagonalSolve_z(*argsCpx)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   X.tag == sTag: lib.ElDiagonalSolveDist_s(*args)
    elif X.tag == dTag: lib.ElDiagonalSolveDist_d(*args)
    elif X.tag == cTag: lib.ElDiagonalSolveDist_c(*argsCpx)
    elif X.tag == zTag: lib.ElDiagonalSolveDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# Dot
# ---
lib.ElDot_i.argtypes = [c_void_p,c_void_p,POINTER(iType)]
lib.ElDot_i.restype = c_uint
lib.ElDot_s.argtypes = [c_void_p,c_void_p,POINTER(sType)]
lib.ElDot_s.restype = c_uint
lib.ElDot_d.argtypes = [c_void_p,c_void_p,POINTER(dType)]
lib.ElDot_d.restype = c_uint
lib.ElDot_c.argtypes = [c_void_p,c_void_p,POINTER(cType)]
lib.ElDot_c.restype = c_uint
lib.ElDot_z.argtypes = [c_void_p,c_void_p,POINTER(zType)]
lib.ElDot_z.restype = c_uint
lib.ElDotDist_i.argtypes = [c_void_p,c_void_p,POINTER(iType)]
lib.ElDotDist_i.restype = c_uint
lib.ElDotDist_s.argtypes = [c_void_p,c_void_p,POINTER(sType)]
lib.ElDotDist_s.restype = c_uint
lib.ElDotDist_d.argtypes = [c_void_p,c_void_p,POINTER(dType)]
lib.ElDotDist_d.restype = c_uint
lib.ElDotDist_c.argtypes = [c_void_p,c_void_p,POINTER(cType)]
lib.ElDotDist_c.restype = c_uint
lib.ElDotDist_z.argtypes = [c_void_p,c_void_p,POINTER(zType)]
lib.ElDotDist_z.restype = c_uint
def Dot(A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  prod = TagToType(A.tag)()
  args = [A.obj,B.obj,pointer(prod)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElDot_i(*args) 
    elif A.tag == sTag: lib.ElDot_s(*args)
    elif A.tag == dTag: lib.ElDot_d(*args)
    elif A.tag == cTag: lib.ElDot_c(*args)
    elif A.tag == zTag: lib.ElDot_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElDotDist_i(*args) 
    elif A.tag == sTag: lib.ElDotDist_s(*args)
    elif A.tag == dTag: lib.ElDotDist_d(*args)
    elif A.tag == cTag: lib.ElDotDist_c(*args)
    elif A.tag == zTag: lib.ElDotDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return prod

# Dotu
# ----
lib.ElDot_c.argtypes = [c_void_p,c_void_p,POINTER(cType)]
lib.ElDot_c.restype = c_uint
lib.ElDot_z.argtypes = [c_void_p,c_void_p,POINTER(zType)]
lib.ElDot_z.restype = c_uint
lib.ElDotDist_c.argtypes = [c_void_p,c_void_p,POINTER(cType)]
lib.ElDotDist_c.restype = c_uint
lib.ElDotDist_z.argtypes = [c_void_p,c_void_p,POINTER(zType)]
lib.ElDotDist_z.restype = c_uint
def Dotu(A,B):
  if type(A) is not type(B): raise Exception('Types of A and B must match')
  if A.tag != B.tag: raise Exception('Datatypes of A and B must match')
  prod = TagToType(A.tag)()
  args = [A.obj,B.obj,pointer(prod)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElDot_i(*args) 
    elif A.tag == sTag: lib.ElDot_s(*args)
    elif A.tag == dTag: lib.ElDot_d(*args)
    elif A.tag == cTag: lib.ElDotu_c(*args)
    elif A.tag == zTag: lib.ElDotu_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElDotDist_i(*args) 
    elif A.tag == sTag: lib.ElDotDist_s(*args)
    elif A.tag == dTag: lib.ElDotDist_d(*args)
    elif A.tag == cTag: lib.ElDotuDist_c(*args)
    elif A.tag == zTag: lib.ElDotuDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return prod

# Entrywise fill
# --------------
lib.ElEntrywiseFill_i.argtypes = [c_void_p,CFUNCTYPE(iType)]
lib.ElEntrywiseFill_i.restype = c_uint
lib.ElEntrywiseFill_s.argtypes = [c_void_p,CFUNCTYPE(sType)]
lib.ElEntrywiseFill_s.restype = c_uint
lib.ElEntrywiseFill_d.argtypes = [c_void_p,CFUNCTYPE(dType)]
lib.ElEntrywiseFill_d.restype = c_uint
lib.ElEntrywiseFill_c.argtypes = [c_void_p,CFUNCTYPE(cType)]
lib.ElEntrywiseFill_c.restype = c_uint
lib.ElEntrywiseFill_z.argtypes = [c_void_p,CFUNCTYPE(zType)]
lib.ElEntrywiseFill_z.restype = c_uint
lib.ElEntrywiseFillDist_i.argtypes = [c_void_p,CFUNCTYPE(iType)]
lib.ElEntrywiseFillDist_i.restype = c_uint
lib.ElEntrywiseFillDist_s.argtypes = [c_void_p,CFUNCTYPE(sType)]
lib.ElEntrywiseFillDist_s.restype = c_uint
lib.ElEntrywiseFillDist_d.argtypes = [c_void_p,CFUNCTYPE(dType)]
lib.ElEntrywiseFillDist_d.restype = c_uint
lib.ElEntrywiseFillDist_c.argtypes = [c_void_p,CFUNCTYPE(cType)]
lib.ElEntrywiseFillDist_c.restype = c_uint
lib.ElEntrywiseFillDist_z.argtypes = [c_void_p,CFUNCTYPE(zType)]
lib.ElEntrywiseFillDist_z.restype = c_uint
def EntrywiseFill(A,fill):
  cFill = CFUNCTYPE(TagToType(A.tag))(fill)
  args = [A.obj,cFill]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElEntrywiseFill_i(*args)
    elif A.tag == sTag: lib.ElEntrywiseFill_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseFill_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseFill_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseFill_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElEntrywiseFillDist_i(*args)
    elif A.tag == sTag: lib.ElEntrywiseFillDist_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseFillDist_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseFillDist_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseFillDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Entrywise map
# -------------
lib.ElEntrywiseMap_i.argtypes = [c_void_p,CFUNCTYPE(iType,iType)]
lib.ElEntrywiseMap_i.restype = c_uint
lib.ElEntrywiseMap_s.argtypes = [c_void_p,CFUNCTYPE(sType,sType)]
lib.ElEntrywiseMap_s.restype = c_uint
lib.ElEntrywiseMap_d.argtypes = [c_void_p,CFUNCTYPE(dType,dType)]
lib.ElEntrywiseMap_d.restype = c_uint
lib.ElEntrywiseMap_c.argtypes = [c_void_p,CFUNCTYPE(cType,cType)]
lib.ElEntrywiseMap_c.restype = c_uint
lib.ElEntrywiseMap_z.argtypes = [c_void_p,CFUNCTYPE(zType,zType)]
lib.ElEntrywiseMap_z.restype = c_uint
lib.ElEntrywiseMapDist_i.argtypes = [c_void_p,CFUNCTYPE(iType,iType)]
lib.ElEntrywiseMapDist_i.restype = c_uint
lib.ElEntrywiseMapDist_s.argtypes = [c_void_p,CFUNCTYPE(sType,sType)]
lib.ElEntrywiseMapDist_s.restype = c_uint
lib.ElEntrywiseMapDist_d.argtypes = [c_void_p,CFUNCTYPE(dType,dType)]
lib.ElEntrywiseMapDist_d.restype = c_uint
lib.ElEntrywiseMapDist_c.argtypes = [c_void_p,CFUNCTYPE(cType,cType)]
lib.ElEntrywiseMapDist_c.restype = c_uint
lib.ElEntrywiseMapDist_z.argtypes = [c_void_p,CFUNCTYPE(zType,zType)]
lib.ElEntrywiseMapDist_z.restype = c_uint
lib.ElEntrywiseMapSparse_i.argtypes = [c_void_p,CFUNCTYPE(iType,iType)]
lib.ElEntrywiseMapSparse_i.restype = c_uint
lib.ElEntrywiseMapSparse_s.argtypes = [c_void_p,CFUNCTYPE(sType,sType)]
lib.ElEntrywiseMapSparse_s.restype = c_uint
lib.ElEntrywiseMapSparse_d.argtypes = [c_void_p,CFUNCTYPE(dType,dType)]
lib.ElEntrywiseMapSparse_d.restype = c_uint
lib.ElEntrywiseMapSparse_c.argtypes = [c_void_p,CFUNCTYPE(cType,cType)]
lib.ElEntrywiseMapSparse_c.restype = c_uint
lib.ElEntrywiseMapSparse_z.argtypes = [c_void_p,CFUNCTYPE(zType,zType)]
lib.ElEntrywiseMapSparse_z.restype = c_uint
lib.ElEntrywiseMapDistSparse_i.argtypes = [c_void_p,CFUNCTYPE(iType,iType)]
lib.ElEntrywiseMapDistSparse_i.restype = c_uint
lib.ElEntrywiseMapDistSparse_s.argtypes = [c_void_p,CFUNCTYPE(sType,sType)]
lib.ElEntrywiseMapDistSparse_s.restype = c_uint
lib.ElEntrywiseMapDistSparse_d.argtypes = [c_void_p,CFUNCTYPE(dType,dType)]
lib.ElEntrywiseMapDistSparse_d.restype = c_uint
lib.ElEntrywiseMapDistSparse_c.argtypes = [c_void_p,CFUNCTYPE(cType,cType)]
lib.ElEntrywiseMapDistSparse_c.restype = c_uint
lib.ElEntrywiseMapDistSparse_z.argtypes = [c_void_p,CFUNCTYPE(zType,zType)]
lib.ElEntrywiseMapDistSparse_z.restype = c_uint

def EntrywiseMap(A,mapFunc):
  cMap = CFUNCTYPE(TagToType(A.tag),TagToType(A.tag))(mapFunc)
  args = [A.obj,cMap]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElEntrywiseMap_i(*args)
    elif A.tag == sTag: lib.ElEntrywiseMap_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseMap_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseMap_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseMap_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElEntrywiseMapDist_i(*args)
    elif A.tag == sTag: lib.ElEntrywiseMapDist_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseMapDist_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseMapDist_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseMapDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElEntrywiseMapSparse_i(*args)
    elif A.tag == sTag: lib.ElEntrywiseMapSparse_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseMapSparse_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseMapSparse_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseMapSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElEntrywiseMapDistSparse_i(*args)
    elif A.tag == sTag: lib.ElEntrywiseMapDistSparse_s(*args)
    elif A.tag == dTag: lib.ElEntrywiseMapDistSparse_d(*args)
    elif A.tag == cTag: lib.ElEntrywiseMapDistSparse_c(*args)
    elif A.tag == zTag: lib.ElEntrywiseMapDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Fill
# ----
lib.ElFill_i.argtypes = [c_void_p,iType]
lib.ElFill_i.restype = c_uint
lib.ElFill_s.argtypes = [c_void_p,sType]
lib.ElFill_s.restype = c_uint
lib.ElFill_d.argtypes = [c_void_p,dType]
lib.ElFill_d.restype = c_uint
lib.ElFill_c.argtypes = [c_void_p,cType]
lib.ElFill_c.restype = c_uint
lib.ElFill_z.argtypes = [c_void_p,zType]
lib.ElFill_z.restype = c_uint
lib.ElFillDist_i.argtypes = [c_void_p,iType]
lib.ElFillDist_i.restype = c_uint
lib.ElFillDist_s.argtypes = [c_void_p,sType]
lib.ElFillDist_s.restype = c_uint
lib.ElFillDist_d.argtypes = [c_void_p,dType]
lib.ElFillDist_d.restype = c_uint
lib.ElFillDist_c.argtypes = [c_void_p,cType]
lib.ElFillDist_c.restype = c_uint
lib.ElFillDist_z.argtypes = [c_void_p,zType]
lib.ElFillDist_z.restype = c_uint
def Fill(A,alphaPre):
  alpha = TagToType(A.tag)(alphaPre)
  args = [A.obj,alpha]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElFill_i(*args)
    elif A.tag == sTag: lib.ElFill_s(*args)
    elif A.tag == dTag: lib.ElFill_d(*args)
    elif A.tag == cTag: lib.ElFill_c(*args)
    elif A.tag == zTag: lib.ElFill_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElFillDist_i(*args)
    elif A.tag == sTag: lib.ElFillDist_s(*args)
    elif A.tag == dTag: lib.ElFillDist_d(*args)
    elif A.tag == cTag: lib.ElFillDist_c(*args)
    elif A.tag == zTag: lib.ElFillDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hadamard
# --------
lib.ElHadamard_i.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamard_i.restype = c_uint
lib.ElHadamard_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamard_s.restype = c_uint
lib.ElHadamard_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamard_d.restype = c_uint
lib.ElHadamard_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamard_c.restype = c_uint
lib.ElHadamard_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamard_z.restype = c_uint
lib.ElHadamardDist_i.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamardDist_i.restype = c_uint
lib.ElHadamardDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamardDist_s.restype = c_uint
lib.ElHadamardDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamardDist_d.restype = c_uint
lib.ElHadamardDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamardDist_c.restype = c_uint
lib.ElHadamardDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElHadamardDist_z.restype = c_uint
def Hadamard(A,B,C):
  if type(A) is not type(B) or type(B) is not type(C):
    raise Exception('Matrix types must match')
  if A.tag != B.tag or B.tag != C.tag:
    raise Exception('Matrix datatypes must match')
  args = [A.obj,B.obj,C.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElHadamard_i(*args)
    elif A.tag == sTag: lib.ElHadamard_s(*args)
    elif A.tag == dTag: lib.ElHadamard_d(*args)
    elif A.tag == cTag: lib.ElHadamard_c(*args)
    elif A.tag == zTag: lib.ElHadamard_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElHadamardDist_i(*args)
    elif A.tag == sTag: lib.ElHadamardDist_s(*args)
    elif A.tag == dTag: lib.ElHadamardDist_d(*args)
    elif A.tag == cTag: lib.ElHadamardDist_c(*args)
    elif A.tag == zTag: lib.ElHadamardDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Hilbert-Schmidt
# ---------------
lib.ElHilbertSchmidt_i.argtypes = [c_void_p,c_void_p,POINTER(iType)]
lib.ElHilbertSchmidt_i.restype = c_uint
lib.ElHilbertSchmidt_s.argtypes = [c_void_p,c_void_p,POINTER(sType)]
lib.ElHilbertSchmidt_s.restype = c_uint
lib.ElHilbertSchmidt_d.argtypes = [c_void_p,c_void_p,POINTER(dType)]
lib.ElHilbertSchmidt_d.restype = c_uint
lib.ElHilbertSchmidt_c.argtypes = [c_void_p,c_void_p,POINTER(cType)]
lib.ElHilbertSchmidt_c.restype = c_uint
lib.ElHilbertSchmidt_z.argtypes = [c_void_p,c_void_p,POINTER(zType)]
lib.ElHilbertSchmidt_z.restype = c_uint
lib.ElHilbertSchmidtDist_i.argtypes = [c_void_p,c_void_p,POINTER(iType)]
lib.ElHilbertSchmidtDist_i.restype = c_uint
lib.ElHilbertSchmidtDist_s.argtypes = [c_void_p,c_void_p,POINTER(sType)]
lib.ElHilbertSchmidtDist_s.restype = c_uint
lib.ElHilbertSchmidtDist_d.argtypes = [c_void_p,c_void_p,POINTER(dType)]
lib.ElHilbertSchmidtDist_d.restype = c_uint
lib.ElHilbertSchmidtDist_c.argtypes = [c_void_p,c_void_p,POINTER(cType)]
lib.ElHilbertSchmidtDist_c.restype = c_uint
lib.ElHilbertSchmidtDist_z.argtypes = [c_void_p,c_void_p,POINTER(zType)]
lib.ElHilbertSchmidtDist_z.restype = c_uint
def HilbertSchmidt(A,B):
  if type(A) is type(B): raise Exception('Matrix types must match')
  if A.tag != B.tag: raise Exception('Datatypes must match')
  prod = TagToType(A.tag)()
  args = [A.obj,B.obj,pointer(prod)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElHilbertSchmidt_i(*args)
    elif A.tag == sTag: lib.ElHilbertSchmidt_s(*args)
    elif A.tag == dTag: lib.ElHilbertSchmidt_d(*args)
    elif A.tag == cTag: lib.ElHilbertSchmidt_c(*args)
    elif A.tag == zTag: lib.ElHilbertSchmidt_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElHilbertSchmidtDist_i(*args)
    elif A.tag == sTag: lib.ElHilbertSchmidtDist_s(*args)
    elif A.tag == dTag: lib.ElHilbertSchmidtDist_d(*args)
    elif A.tag == cTag: lib.ElHilbertSchmidtDist_c(*args)
    elif A.tag == zTag: lib.ElHilbertSchmidtDist_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return prod

# Index dependent fill
# --------------------
lib.ElIndexDependentFill_i.argtypes = [c_void_p,CFUNCTYPE(iType,iType,iType)]
lib.ElIndexDependentFill_i.restype = c_uint
lib.ElIndexDependentFill_s.argtypes = [c_void_p,CFUNCTYPE(sType,iType,iType)]
lib.ElIndexDependentFill_s.restype = c_uint
lib.ElIndexDependentFill_d.argtypes = [c_void_p,CFUNCTYPE(dType,iType,iType)]
lib.ElIndexDependentFill_d.restype = c_uint
lib.ElIndexDependentFill_c.argtypes = [c_void_p,CFUNCTYPE(cType,iType,iType)]
lib.ElIndexDependentFill_c.restype = c_uint
lib.ElIndexDependentFill_z.argtypes = [c_void_p,CFUNCTYPE(zType,iType,iType)]
lib.ElIndexDependentFill_z.restype = c_uint
lib.ElIndexDependentFillDist_i.argtypes = \
  [c_void_p,CFUNCTYPE(iType,iType,iType)]
lib.ElIndexDependentFillDist_i.restype = c_uint
lib.ElIndexDependentFillDist_s.argtypes = \
  [c_void_p,CFUNCTYPE(sType,iType,iType)]
lib.ElIndexDependentFillDist_s.restype = c_uint
lib.ElIndexDependentFillDist_d.argtypes = \
  [c_void_p,CFUNCTYPE(dType,iType,iType)]
lib.ElIndexDependentFillDist_d.restype = c_uint
lib.ElIndexDependentFillDist_c.argtypes = \
  [c_void_p,CFUNCTYPE(cType,iType,iType)]
lib.ElIndexDependentFillDist_c.restype = c_uint
lib.ElIndexDependentFillDist_z.argtypes = \
  [c_void_p,CFUNCTYPE(zType,iType,iType)]
lib.ElIndexDependentFillDist_z.restype = c_uint
def IndexDependentFill(A,fill):
  cFill = CFUNCTYPE(TagToType(A.tag),iType,iType)(fill)
  args = [A.obj,cFill]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElIndexDependentFill_i(*args)
    elif A.tag == sTag: lib.ElIndexDependentFill_s(*args)
    elif A.tag == dTag: lib.ElIndexDependentFill_d(*args)
    elif A.tag == cTag: lib.ElIndexDependentFill_c(*args)
    elif A.tag == zTag: lib.ElIndexDependentFill_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElIndexDependentFillDist_i(*args)
    elif A.tag == sTag: lib.ElIndexDependentFillDist_s(*args)
    elif A.tag == dTag: lib.ElIndexDependentFillDist_d(*args)
    elif A.tag == cTag: lib.ElIndexDependentFillDist_c(*args)
    elif A.tag == zTag: lib.ElIndexDependentFillDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Index dependent map
# -------------------
lib.ElIndexDependentMap_i.argtypes = \
  [c_void_p,CFUNCTYPE(iType,iType,iType,iType)]
lib.ElIndexDependentMap_i.restype = c_uint
lib.ElIndexDependentMap_s.argtypes = \
  [c_void_p,CFUNCTYPE(sType,iType,iType,sType)]
lib.ElIndexDependentMap_s.restype = c_uint
lib.ElIndexDependentMap_d.argtypes = \
  [c_void_p,CFUNCTYPE(dType,iType,iType,dType)]
lib.ElIndexDependentMap_d.restype = c_uint
lib.ElIndexDependentMap_c.argtypes = \
  [c_void_p,CFUNCTYPE(cType,iType,iType,cType)]
lib.ElIndexDependentMap_c.restype = c_uint
lib.ElIndexDependentMap_z.argtypes = \
  [c_void_p,CFUNCTYPE(zType,iType,iType,zType)]
lib.ElIndexDependentMap_z.restype = c_uint
lib.ElIndexDependentMapDist_i.argtypes = \
  [c_void_p,CFUNCTYPE(iType,iType,iType,iType)]
lib.ElIndexDependentMapDist_i.restype = c_uint
lib.ElIndexDependentMapDist_s.argtypes = \
  [c_void_p,CFUNCTYPE(sType,iType,iType,sType)]
lib.ElIndexDependentMapDist_s.restype = c_uint
lib.ElIndexDependentMapDist_d.argtypes = \
  [c_void_p,CFUNCTYPE(dType,iType,iType,dType)]
lib.ElIndexDependentMapDist_d.restype = c_uint
lib.ElIndexDependentMapDist_c.argtypes = \
  [c_void_p,CFUNCTYPE(cType,iType,iType,cType)]
lib.ElIndexDependentMapDist_c.restype = c_uint
lib.ElIndexDependentMapDist_z.argtypes = \
  [c_void_p,CFUNCTYPE(zType,iType,iType,zType)]
lib.ElIndexDependentMapDist_z.restype = c_uint
def IndexDependentMap(A,mapFunc):
  typeA = TagToType(A)
  cMap = CFUNCTYPE(typeA,iType,iType,typeA)(mapFunc)
  args = [A.obj,cMap]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElIndexDependentMap_i(*args)
    elif A.tag == sTag: lib.ElIndexDependentMap_s(*args)
    elif A.tag == dTag: lib.ElIndexDependentMap_d(*args)
    elif A.tag == cTag: lib.ElIndexDependentMap_c(*args)
    elif A.tag == zTag: lib.ElIndexDependentMap_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElIndexDependentMapDist_i(*args)
    elif A.tag == sTag: lib.ElIndexDependentMapDist_s(*args)
    elif A.tag == dTag: lib.ElIndexDependentMapDist_d(*args)
    elif A.tag == cTag: lib.ElIndexDependentMapDist_c(*args)
    elif A.tag == zTag: lib.ElIndexDependentMapDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Make symmetric/Hermitian
# ------------------------
lib.ElMakeSymmetric_i.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetric_i.restype = c_uint
lib.ElMakeSymmetric_s.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetric_s.restype = c_uint
lib.ElMakeSymmetric_d.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetric_d.restype = c_uint
lib.ElMakeSymmetric_c.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetric_c.restype = c_uint
lib.ElMakeSymmetric_z.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetric_z.restype = c_uint
lib.ElMakeSymmetricDist_i.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDist_i.restype = c_uint
lib.ElMakeSymmetricDist_s.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDist_s.restype = c_uint
lib.ElMakeSymmetricDist_d.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDist_d.restype = c_uint
lib.ElMakeSymmetricDist_c.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDist_c.restype = c_uint
lib.ElMakeSymmetricDist_z.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDist_z.restype = c_uint
lib.ElMakeSymmetricSparse_i.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricSparse_i.restype = c_uint
lib.ElMakeSymmetricSparse_s.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricSparse_s.restype = c_uint
lib.ElMakeSymmetricSparse_d.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricSparse_d.restype = c_uint
lib.ElMakeSymmetricSparse_c.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricSparse_c.restype = c_uint
lib.ElMakeSymmetricSparse_z.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricSparse_z.restype = c_uint
lib.ElMakeSymmetricDistSparse_i.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDistSparse_i.restype = c_uint
lib.ElMakeSymmetricDistSparse_s.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDistSparse_s.restype = c_uint
lib.ElMakeSymmetricDistSparse_d.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDistSparse_d.restype = c_uint
lib.ElMakeSymmetricDistSparse_c.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDistSparse_c.restype = c_uint
lib.ElMakeSymmetricDistSparse_z.argtypes = [c_uint,c_void_p]
lib.ElMakeSymmetricDistSparse_z.restype = c_uint

lib.ElMakeHermitian_c.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitian_c.restype = c_uint
lib.ElMakeHermitian_z.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitian_z.restype = c_uint
lib.ElMakeHermitianDist_c.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitianDist_c.restype = c_uint
lib.ElMakeHermitianDist_z.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitianDist_z.restype = c_uint
lib.ElMakeHermitianSparse_c.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitianSparse_c.restype = c_uint
lib.ElMakeHermitianSparse_z.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitianSparse_z.restype = c_uint
lib.ElMakeHermitianDistSparse_c.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitianDistSparse_c.restype = c_uint
lib.ElMakeHermitianDistSparse_z.argtypes = [c_uint,c_void_p]
lib.ElMakeHermitianDistSparse_z.restype = c_uint

def MakeSymmetric(uplo,A,conj=False):
  args = [uplo,A.obj]
  if type(A) is Matrix: 
    if   A.tag == iTag: lib.ElMakeSymmetric_i(*args)
    elif A.tag == sTag: lib.ElMakeSymmetric_s(*args)
    elif A.tag == dTag: lib.ElMakeSymmetric_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElMakeHermitian_c(*args)
      else:    lib.ElMakeSymmetric_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElMakeHermitian_z(*args)
      else:    lib.ElMakeSymmetric_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMakeSymmetricDist_i(*args)
    elif A.tag == sTag: lib.ElMakeSymmetricDist_s(*args)
    elif A.tag == dTag: lib.ElMakeSymmetricDist_d(*args)
    elif A.tag == cTag: 
      if conj: lib.ElMakeHermitianDist_c(*args)
      else:    lib.ElMakeSymmetricDist_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElMakeHermitianDist_z(*args)
      else:    lib.ElMakeSymmetricDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElMakeSymmetricSparse_i(*args)
    elif A.tag == sTag: lib.ElMakeSymmetricSparse_s(*args)
    elif A.tag == dTag: lib.ElMakeSymmetricSparse_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElMakeHermitianSparse_c(*args)
      else:         lib.ElMakeHermitianSparse_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElMakeHermitianSparse_z(*args)
      else:         lib.ElMakeHermitianSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElMakeSymmetricDistSparse_i(*args)
    elif A.tag == sTag: lib.ElMakeSymmetricDistSparse_s(*args)
    elif A.tag == dTag: lib.ElMakeSymmetricDistSparse_d(*args)
    elif A.tag == cTag:
      if conjugate: lib.ElMakeHermitianDistSparse_c(*args)
      else:         lib.ElMakeHermitianDistSparse_c(*args)
    elif A.tag == zTag:
      if conjugate: lib.ElMakeHermitianDistSparse_z(*args)
      else:         lib.ElMakeHermitianDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

def MakeHermitian(uplo,A):
  MakeSymmetric(uplo,A,True)

# Make real
# ---------
lib.ElMakeReal_c.argtypes = [c_uint,c_void_p]
lib.ElMakeReal_c.restype = c_uint
lib.ElMakeReal_z.argtypes = [c_uint,c_void_p]
lib.ElMakeReal_z.restype = c_uint
lib.ElMakeRealDist_c.argtypes = [c_uint,c_void_p]
lib.ElMakeRealDist_c.restype = c_uint
lib.ElMakeRealDist_z.argtypes = [c_uint,c_void_p]
lib.ElMakeRealDist_z.restype = c_uint
def MakeReal(uplo,A):
  args = [uplo,A.obj]
  if type(A) is Matrix:
    if   A.tag == cTag: lib.ElMakeReal_c(*args)
    elif A.tag == zTag: lib.ElMakeReal_z(*args)
  elif type(A) is DistMatrix:
    if   A.tag == cTag: lib.ElMakeRealDist_c(*args)
    elif A.tag == zTag: lib.ElMakeRealDist_z(*args)
  else: TypeExcept()

# Make trapezoidal
# ----------------
lib.ElMakeTrapezoidal_i.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidal_i.restype = c_uint
lib.ElMakeTrapezoidal_s.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidal_s.restype = c_uint
lib.ElMakeTrapezoidal_d.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidal_d.restype = c_uint
lib.ElMakeTrapezoidal_c.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidal_c.restype = c_uint
lib.ElMakeTrapezoidal_z.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidal_z.restype = c_uint
lib.ElMakeTrapezoidalDist_i.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalDist_i.restype = c_uint
lib.ElMakeTrapezoidalDist_s.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalDist_s.restype = c_uint
lib.ElMakeTrapezoidalDist_d.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalDist_d.restype = c_uint
lib.ElMakeTrapezoidalDist_c.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalDist_c.restype = c_uint
lib.ElMakeTrapezoidalDist_z.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalDist_z.restype = c_uint
lib.ElMakeTrapezoidalSparse_i.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalSparse_i.restype = c_uint
lib.ElMakeTrapezoidalSparse_s.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalSparse_s.restype = c_uint
lib.ElMakeTrapezoidalSparse_d.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalSparse_d.restype = c_uint
lib.ElMakeTrapezoidalSparse_c.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalSparse_c.restype = c_uint
lib.ElMakeTrapezoidalSparse_z.argtypes = [c_uint,c_void_p,iType]
lib.ElMakeTrapezoidalSparse_z.restype = c_uint
def MakeTrapezoidal(uplo,A,offset=0):
  args = [uplo,A.obj,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMakeTrapezoidal_i(*args)
    elif A.tag == sTag: lib.ElMakeTrapezoidal_s(*args)
    elif A.tag == dTag: lib.ElMakeTrapezoidal_d(*args)
    elif A.tag == cTag: lib.ElMakeTrapezoidal_c(*args)
    elif A.tag == zTag: lib.ElMakeTrapezoidal_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMakeTrapezoidalDist_i(*args)
    elif A.tag == sTag: lib.ElMakeTrapezoidalDist_s(*args)
    elif A.tag == dTag: lib.ElMakeTrapezoidalDist_d(*args)
    elif A.tag == cTag: lib.ElMakeTrapezoidalDist_c(*args)
    elif A.tag == zTag: lib.ElMakeTrapezoidalDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElMakeTrapezoidalSparse_i(*args)
    elif A.tag == sTag: lib.ElMakeTrapezoidalSparse_s(*args)
    elif A.tag == dTag: lib.ElMakeTrapezoidalSparse_d(*args)
    elif A.tag == cTag: lib.ElMakeTrapezoidalSparse_c(*args)
    elif A.tag == zTag: lib.ElMakeTrapezoidalSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElMakeTrapezoidalDistSparse_i(*args)
    elif A.tag == sTag: lib.ElMakeTrapezoidalDistSparse_s(*args)
    elif A.tag == dTag: lib.ElMakeTrapezoidalDistSparse_d(*args)
    elif A.tag == cTag: lib.ElMakeTrapezoidalDistSparse_c(*args)
    elif A.tag == zTag: lib.ElMakeTrapezoidalDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Max
# ---
class ValueInt_i(ctypes.Structure):
  _fields_ = [("value",iType),("index",iType)]
class ValueInt_s(ctypes.Structure):
  _fields_ = [("value",sType),("index",iType)]
class ValueInt_d(ctypes.Structure):
  _fields_ = [("value",dType),("index",iType)]
class ValueInt_c(ctypes.Structure):
  _fields_ = [("value",cType),("index",iType)]
class ValueInt_z(ctypes.Structure):
  _fields_ = [("value",zType),("index",iType)]
def TagToValueInt(tag):
  if   tag == iTag: return ValueInt_i()
  elif tag == sTag: return ValueInt_s()
  elif tag == dTag: return ValueInt_d()
  elif tag == cTag: return ValueInt_c()
  elif tag == zTag: return ValueInt_z()
  else: DataExcept()

class ValueIntPair_i(ctypes.Structure):
  _fields_ = [("value",iType),("indices",(iType*2))]
class ValueIntPair_s(ctypes.Structure):
  _fields_ = [("value",sType),("indices",(iType*2))]
class ValueIntPair_d(ctypes.Structure):
  _fields_ = [("value",dType),("indices",(iType*2))]
class ValueIntPair_c(ctypes.Structure):
  _fields_ = [("value",cType),("indices",(iType*2))]
class ValueIntPair_z(ctypes.Structure):
  _fields_ = [("value",zType),("indices",(iType*2))]
def TagToValueIntPair(tag):
  if   tag == iTag: return ValueIntPair_i()
  elif tag == sTag: return ValueIntPair_s()
  elif tag == dTag: return ValueIntPair_d()
  elif tag == cTag: return ValueIntPair_c()
  elif tag == zTag: return ValueIntPair_z()
  else: DataExcept()

lib.ElMax_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMax_i.restype = c_uint
lib.ElMax_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMax_s.restype = c_uint
lib.ElMax_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMax_d.restype = c_uint
lib.ElMaxDist_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMaxDist_i.restype = c_uint
lib.ElMaxDist_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMaxDist_s.restype = c_uint
lib.ElMaxDist_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMaxDist_d.restype = c_uint
def Max(A):
  pair = TagToValueIntPair(A.tag)
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMax_i(*args) 
    elif A.tag == sTag: lib.ElMax_s(*args) 
    elif A.tag == dTag: lib.ElMax_d(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMaxDist_i(*args) 
    elif A.tag == sTag: lib.ElMaxDist_s(*args) 
    elif A.tag == dTag: lib.ElMaxDist_d(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElSymmetricMax_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMax_i.restype = c_uint
lib.ElSymmetricMax_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMax_s.restype = c_uint
lib.ElSymmetricMax_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMax_d.restype = c_uint
lib.ElSymmetricMaxDist_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMaxDist_i.restype = c_uint
lib.ElSymmetricMaxDist_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMaxDist_s.restype = c_uint
lib.ElSymmetricMaxDist_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMaxDist_d.restype = c_uint
def SymmetricMax(uplo,A):
  pair = TagToValueIntPair(A.tag)
  args = [uplo,A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSymmetricMax_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMax_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMax_d(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSymmetricMaxDist_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMaxDist_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMaxDist_d(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElVectorMax_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMax_i.restype = c_uint
lib.ElVectorMax_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMax_s.restype = c_uint
lib.ElVectorMax_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMax_d.restype = c_uint
lib.ElVectorMaxDist_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMaxDist_i.restype = c_uint
lib.ElVectorMaxDist_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMaxDist_s.restype = c_uint
lib.ElVectorMaxDist_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMaxDist_d.restype = c_uint
def VectorMax(A):
  pair = TagToValueInt(A.tag)
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElVectorMax_i(*args) 
    elif A.tag == sTag: lib.ElVectorMax_s(*args) 
    elif A.tag == dTag: lib.ElVectorMax_d(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElVectorMaxDist_i(*args) 
    elif A.tag == sTag: lib.ElVectorMaxDist_s(*args) 
    elif A.tag == dTag: lib.ElVectorMaxDist_d(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.index

# MaxAbs
# ------
lib.ElMaxAbs_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMaxAbs_i.restype = c_uint
lib.ElMaxAbs_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMaxAbs_s.restype = c_uint
lib.ElMaxAbs_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMaxAbs_d.restype = c_uint
lib.ElMaxAbs_c.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMaxAbs_c.restype = c_uint
lib.ElMaxAbs_z.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMaxAbs_z.restype = c_uint
lib.ElMaxAbsDist_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMaxAbsDist_i.restype = c_uint
lib.ElMaxAbsDist_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMaxAbsDist_s.restype = c_uint
lib.ElMaxAbsDist_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMaxAbsDist_d.restype = c_uint
lib.ElMaxAbsDist_c.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMaxAbsDist_c.restype = c_uint
lib.ElMaxAbsDist_z.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMaxAbsDist_z.restype = c_uint
def MaxAbs(A):
  pair = TagToValueIntPair(Base(A.tag))
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMaxAbs_i(*args) 
    elif A.tag == sTag: lib.ElMaxAbs_s(*args) 
    elif A.tag == dTag: lib.ElMaxAbs_d(*args) 
    elif A.tag == cTag: lib.ElMaxAbs_c(*args) 
    elif A.tag == zTag: lib.ElMaxAbs_z(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMaxAbsDist_i(*args) 
    elif A.tag == sTag: lib.ElMaxAbsDist_s(*args) 
    elif A.tag == dTag: lib.ElMaxAbsDist_d(*args) 
    elif A.tag == cTag: lib.ElMaxAbsDist_c(*args) 
    elif A.tag == zTag: lib.ElMaxAbsDist_z(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElSymmetricMaxAbs_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMaxAbs_i.restype = c_uint
lib.ElSymmetricMaxAbs_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMaxAbs_s.restype = c_uint
lib.ElSymmetricMaxAbs_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMaxAbs_d.restype = c_uint
lib.ElSymmetricMaxAbs_c.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMaxAbs_c.restype = c_uint
lib.ElSymmetricMaxAbs_z.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMaxAbs_z.restype = c_uint
lib.ElSymmetricMaxAbsDist_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMaxAbsDist_i.restype = c_uint
lib.ElSymmetricMaxAbsDist_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMaxAbsDist_s.restype = c_uint
lib.ElSymmetricMaxAbsDist_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMaxAbsDist_d.restype = c_uint
lib.ElSymmetricMaxAbsDist_c.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMaxAbsDist_c.restype = c_uint
lib.ElSymmetricMaxAbsDist_z.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMaxAbsDist_z.restype = c_uint
def SymmetricMaxAbs(uplo,A):
  pair = TagToValueIntPair(Base(A.tag))
  args = [uplo,A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSymmetricMaxAbs_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMaxAbs_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMaxAbs_d(*args) 
    elif A.tag == cTag: lib.ElSymmetricMaxAbs_c(*args) 
    elif A.tag == zTag: lib.ElSymmetricMaxAbs_z(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSymmetricMaxAbsDist_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMaxAbsDist_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMaxAbsDist_d(*args) 
    elif A.tag == cTag: lib.ElSymmetricMaxAbsDist_c(*args) 
    elif A.tag == zTag: lib.ElSymmetricMaxAbsDist_z(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElVectorMaxAbs_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMaxAbs_i.restype = c_uint
lib.ElVectorMaxAbs_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMaxAbs_s.restype = c_uint
lib.ElVectorMaxAbs_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMaxAbs_d.restype = c_uint
lib.ElVectorMaxAbs_c.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMaxAbs_c.restype = c_uint
lib.ElVectorMaxAbs_z.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMaxAbs_z.restype = c_uint
lib.ElVectorMaxAbsDist_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMaxAbsDist_i.restype = c_uint
lib.ElVectorMaxAbsDist_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMaxAbsDist_s.restype = c_uint
lib.ElVectorMaxAbsDist_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMaxAbsDist_d.restype = c_uint
lib.ElVectorMaxAbsDist_c.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMaxAbsDist_c.restype = c_uint
lib.ElVectorMaxAbsDist_z.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMaxAbsDist_z.restype = c_uint
def VectorMaxAbs(A):
  pair = TagToValueInt(Base(A.tag))
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElVectorMaxAbs_i(*args) 
    elif A.tag == sTag: lib.ElVectorMaxAbs_s(*args) 
    elif A.tag == dTag: lib.ElVectorMaxAbs_d(*args) 
    elif A.tag == cTag: lib.ElVectorMaxAbs_c(*args) 
    elif A.tag == zTag: lib.ElVectorMaxAbs_z(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElVectorMaxAbsDist_i(*args) 
    elif A.tag == sTag: lib.ElVectorMaxAbsDist_s(*args) 
    elif A.tag == dTag: lib.ElVectorMaxAbsDist_d(*args) 
    elif A.tag == cTag: lib.ElVectorMaxAbsDist_c(*args) 
    elif A.tag == zTag: lib.ElVectorMaxAbsDist_z(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.index

# Min
# ---
lib.ElMin_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMin_i.restype = c_uint
lib.ElMin_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMin_s.restype = c_uint
lib.ElMin_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMin_d.restype = c_uint
lib.ElMinDist_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMinDist_i.restype = c_uint
lib.ElMinDist_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMinDist_s.restype = c_uint
lib.ElMinDist_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMinDist_d.restype = c_uint
def Min(A):
  pair = TagToValueIntPair(A.tag)
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMin_i(*args) 
    elif A.tag == sTag: lib.ElMin_s(*args) 
    elif A.tag == dTag: lib.ElMin_d(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMinDist_i(*args) 
    elif A.tag == sTag: lib.ElMinDist_s(*args) 
    elif A.tag == dTag: lib.ElMinDist_d(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElSymmetricMin_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMin_i.restype = c_uint
lib.ElSymmetricMin_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMin_s.restype = c_uint
lib.ElSymmetricMin_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMin_d.restype = c_uint
lib.ElSymmetricMinDist_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMinDist_i.restype = c_uint
lib.ElSymmetricMinDist_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMinDist_s.restype = c_uint
lib.ElSymmetricMinDist_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMinDist_d.restype = c_uint
def SymmetricMin(uplo,A):
  pair = TagToValueIntPair(A.tag)
  args = [uplo,A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSymmetricMin_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMin_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMin_d(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSymmetricMinDist_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMinDist_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMinDist_d(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElVectorMin_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMin_i.restype = c_uint
lib.ElVectorMin_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMin_s.restype = c_uint
lib.ElVectorMin_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMin_d.restype = c_uint
lib.ElVectorMinDist_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMinDist_i.restype = c_uint
lib.ElVectorMinDist_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMinDist_s.restype = c_uint
lib.ElVectorMinDist_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMinDist_d.restype = c_uint
def VectorMin(A):
  pair = TagToValueInt(A.tag)
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElVectorMin_i(*args) 
    elif A.tag == sTag: lib.ElVectorMin_s(*args) 
    elif A.tag == dTag: lib.ElVectorMin_d(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElVectorMinDist_i(*args) 
    elif A.tag == sTag: lib.ElVectorMinDist_s(*args) 
    elif A.tag == dTag: lib.ElVectorMinDist_d(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.index

# MinAbs
# ------
lib.ElMinAbs_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMinAbs_i.restype = c_uint
lib.ElMinAbs_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMinAbs_s.restype = c_uint
lib.ElMinAbs_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMinAbs_d.restype = c_uint
lib.ElMinAbs_c.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMinAbs_c.restype = c_uint
lib.ElMinAbs_z.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMinAbs_z.restype = c_uint
lib.ElMinAbsDist_i.argtypes = [c_void_p,POINTER(ValueIntPair_i)]
lib.ElMinAbsDist_i.restype = c_uint
lib.ElMinAbsDist_s.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMinAbsDist_s.restype = c_uint
lib.ElMinAbsDist_d.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMinAbsDist_d.restype = c_uint
lib.ElMinAbsDist_c.argtypes = [c_void_p,POINTER(ValueIntPair_s)]
lib.ElMinAbsDist_c.restype = c_uint
lib.ElMinAbsDist_z.argtypes = [c_void_p,POINTER(ValueIntPair_d)]
lib.ElMinAbsDist_z.restype = c_uint
def MinAbs(A):
  pair = TagToValueIntPair(Base(A.tag))
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElMinAbs_i(*args) 
    elif A.tag == sTag: lib.ElMinAbs_s(*args) 
    elif A.tag == dTag: lib.ElMinAbs_d(*args) 
    elif A.tag == cTag: lib.ElMinAbs_c(*args) 
    elif A.tag == zTag: lib.ElMinAbs_z(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElMinAbsDist_i(*args) 
    elif A.tag == sTag: lib.ElMinAbsDist_s(*args) 
    elif A.tag == dTag: lib.ElMinAbsDist_d(*args) 
    elif A.tag == cTag: lib.ElMinAbsDist_c(*args) 
    elif A.tag == zTag: lib.ElMinAbsDist_z(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElSymmetricMinAbs_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMinAbs_i.restype = c_uint
lib.ElSymmetricMinAbs_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMinAbs_s.restype = c_uint
lib.ElSymmetricMinAbs_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMinAbs_d.restype = c_uint
lib.ElSymmetricMinAbs_c.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMinAbs_c.restype = c_uint
lib.ElSymmetricMinAbs_z.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMinAbs_z.restype = c_uint
lib.ElSymmetricMinAbsDist_i.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_i)]
lib.ElSymmetricMinAbsDist_i.restype = c_uint
lib.ElSymmetricMinAbsDist_s.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMinAbsDist_s.restype = c_uint
lib.ElSymmetricMinAbsDist_d.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMinAbsDist_d.restype = c_uint
lib.ElSymmetricMinAbsDist_c.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_s)]
lib.ElSymmetricMinAbsDist_c.restype = c_uint
lib.ElSymmetricMinAbsDist_z.argtypes = [c_uint,c_void_p,POINTER(ValueIntPair_d)]
lib.ElSymmetricMinAbsDist_z.restype = c_uint
def SymmetricMinAbs(uplo,A):
  pair = TagToValueIntPair(Base(A.tag))
  args = [uplo,A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSymmetricMinAbs_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMinAbs_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMinAbs_d(*args) 
    elif A.tag == cTag: lib.ElSymmetricMinAbs_c(*args) 
    elif A.tag == zTag: lib.ElSymmetricMinAbs_z(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSymmetricMinAbsDist_i(*args) 
    elif A.tag == sTag: lib.ElSymmetricMinAbsDist_s(*args) 
    elif A.tag == dTag: lib.ElSymmetricMinAbsDist_d(*args) 
    elif A.tag == cTag: lib.ElSymmetricMinAbsDist_c(*args) 
    elif A.tag == zTag: lib.ElSymmetricMinAbsDist_z(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.indices[0], pair.indices[1]

lib.ElVectorMinAbs_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMinAbs_i.restype = c_uint
lib.ElVectorMinAbs_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMinAbs_s.restype = c_uint
lib.ElVectorMinAbs_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMinAbs_d.restype = c_uint
lib.ElVectorMinAbs_c.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMinAbs_c.restype = c_uint
lib.ElVectorMinAbs_z.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMinAbs_z.restype = c_uint
lib.ElVectorMinAbsDist_i.argtypes = [c_void_p,POINTER(ValueInt_i)]
lib.ElVectorMinAbsDist_i.restype = c_uint
lib.ElVectorMinAbsDist_s.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMinAbsDist_s.restype = c_uint
lib.ElVectorMinAbsDist_d.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMinAbsDist_d.restype = c_uint
lib.ElVectorMinAbsDist_c.argtypes = [c_void_p,POINTER(ValueInt_s)]
lib.ElVectorMinAbsDist_c.restype = c_uint
lib.ElVectorMinAbsDist_z.argtypes = [c_void_p,POINTER(ValueInt_d)]
lib.ElVectorMinAbsDist_z.restype = c_uint
def VectorMinAbs(A):
  pair = TagToValueInt(Base(A.tag))
  args = [A.obj,pointer(pair)]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElVectorMinAbs_i(*args) 
    elif A.tag == sTag: lib.ElVectorMinAbs_s(*args) 
    elif A.tag == dTag: lib.ElVectorMinAbs_d(*args) 
    elif A.tag == cTag: lib.ElVectorMinAbs_c(*args) 
    elif A.tag == zTag: lib.ElVectorMinAbs_z(*args) 
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElVectorMinAbsDist_i(*args) 
    elif A.tag == sTag: lib.ElVectorMinAbsDist_s(*args) 
    elif A.tag == dTag: lib.ElVectorMinAbsDist_d(*args) 
    elif A.tag == cTag: lib.ElVectorMinAbsDist_c(*args) 
    elif A.tag == zTag: lib.ElVectorMinAbsDist_z(*args) 
    else: DataExcept()
  else: TypeExcept()
  return pair.value, pair.index

# Nrm2
# ----
lib.ElNrm2_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElNrm2_s.restype = c_uint
lib.ElNrm2_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElNrm2_d.restype = c_uint
lib.ElNrm2_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElNrm2_c.restype = c_uint
lib.ElNrm2_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElNrm2_z.restype = c_uint
lib.ElNrm2Dist_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElNrm2Dist_s.restype = c_uint
lib.ElNrm2Dist_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElNrm2Dist_d.restype = c_uint
lib.ElNrm2Dist_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElNrm2Dist_c.restype = c_uint
lib.ElNrm2Dist_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElNrm2Dist_z.restype = c_uint
lib.ElNrm2DistMultiVec_s.argtypes = [c_void_p,POINTER(sType)]
lib.ElNrm2DistMultiVec_s.restype = c_uint
lib.ElNrm2DistMultiVec_d.argtypes = [c_void_p,POINTER(dType)]
lib.ElNrm2DistMultiVec_d.restype = c_uint
lib.ElNrm2DistMultiVec_c.argtypes = [c_void_p,POINTER(sType)]
lib.ElNrm2DistMultiVec_c.restype = c_uint
lib.ElNrm2DistMultiVec_z.argtypes = [c_void_p,POINTER(dType)]
lib.ElNrm2DistMultiVec_z.restype = c_uint
def Nrm2(A):
  gamma = TagToType(Base(A.tag))()
  args = [A.obj,pointer(gamma)]
  if type(A) is Matrix:
    if   A.tag == sTag: lib.ElNrm2_s(*args)
    elif A.tag == dTag: lib.ElNrm2_d(*args)
    elif A.tag == cTag: lib.ElNrm2_c(*args)
    elif A.tag == zTag: lib.ElNrm2_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == sTag: lib.ElNrm2Dist_s(*args)
    elif A.tag == dTag: lib.ElNrm2Dist_d(*args)
    elif A.tag == cTag: lib.ElNrm2Dist_c(*args)
    elif A.tag == zTag: lib.ElNrm2Dist_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == sTag: lib.ElNrm2DistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElNrm2DistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElNrm2DistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElNrm2DistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()
  return gamma.value

# Scale
# -----
lib.ElScale_i.argtypes = [iType,c_void_p]
lib.ElScale_i.restype = c_uint
lib.ElScale_s.argtypes = [sType,c_void_p]
lib.ElScale_s.restype = c_uint
lib.ElScale_d.argtypes = [dType,c_void_p]
lib.ElScale_d.restype = c_uint
lib.ElScale_c.argtypes = [cType,c_void_p]
lib.ElScale_c.restype = c_uint
lib.ElScale_z.argtypes = [zType,c_void_p]
lib.ElScale_z.restype = c_uint
lib.ElScaleDist_i.argtypes = [iType,c_void_p]
lib.ElScaleDist_i.restype = c_uint
lib.ElScaleDist_s.argtypes = [sType,c_void_p]
lib.ElScaleDist_s.restype = c_uint
lib.ElScaleDist_d.argtypes = [dType,c_void_p]
lib.ElScaleDist_d.restype = c_uint
lib.ElScaleDist_c.argtypes = [cType,c_void_p]
lib.ElScaleDist_c.restype = c_uint
lib.ElScaleDist_z.argtypes = [zType,c_void_p]
lib.ElScaleDist_z.restype = c_uint
lib.ElScaleSparse_i.argtypes = [iType,c_void_p]
lib.ElScaleSparse_i.restype = c_uint
lib.ElScaleSparse_s.argtypes = [sType,c_void_p]
lib.ElScaleSparse_s.restype = c_uint
lib.ElScaleSparse_d.argtypes = [dType,c_void_p]
lib.ElScaleSparse_d.restype = c_uint
lib.ElScaleSparse_c.argtypes = [cType,c_void_p]
lib.ElScaleSparse_c.restype = c_uint
lib.ElScaleSparse_z.argtypes = [zType,c_void_p]
lib.ElScaleSparse_z.restype = c_uint
lib.ElScaleDistSparse_i.argtypes = [iType,c_void_p]
lib.ElScaleDistSparse_i.restype = c_uint
lib.ElScaleDistSparse_s.argtypes = [sType,c_void_p]
lib.ElScaleDistSparse_s.restype = c_uint
lib.ElScaleDistSparse_d.argtypes = [dType,c_void_p]
lib.ElScaleDistSparse_d.restype = c_uint
lib.ElScaleDistSparse_c.argtypes = [cType,c_void_p]
lib.ElScaleDistSparse_c.restype = c_uint
lib.ElScaleDistSparse_z.argtypes = [zType,c_void_p]
lib.ElScaleDistSparse_z.restype = c_uint
lib.ElScaleDistMultiVec_i.argtypes = [iType,c_void_p]
lib.ElScaleDistMultiVec_i.restype = c_uint
lib.ElScaleDistMultiVec_s.argtypes = [sType,c_void_p]
lib.ElScaleDistMultiVec_s.restype = c_uint
lib.ElScaleDistMultiVec_d.argtypes = [dType,c_void_p]
lib.ElScaleDistMultiVec_d.restype = c_uint
lib.ElScaleDistMultiVec_c.argtypes = [cType,c_void_p]
lib.ElScaleDistMultiVec_c.restype = c_uint
lib.ElScaleDistMultiVec_z.argtypes = [zType,c_void_p]
lib.ElScaleDistMultiVec_z.restype = c_uint
def Scale(alphaPre,A):
  alpha = TagToType(A.tag)(alphaPre)
  args = [alpha,A.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElScale_i(*args)
    elif A.tag == sTag: lib.ElScale_s(*args)
    elif A.tag == dTag: lib.ElScale_d(*args)
    elif A.tag == cTag: lib.ElScale_c(*args)
    elif A.tag == zTag: lib.ElScale_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElScaleDist_i(*args)
    elif A.tag == sTag: lib.ElScaleDist_s(*args)
    elif A.tag == dTag: lib.ElScaleDist_d(*args)
    elif A.tag == cTag: lib.ElScaleDist_c(*args)
    elif A.tag == zTag: lib.ElScaleDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElScaleSparse_i(*args)
    elif A.tag == sTag: lib.ElScaleSparse_s(*args)
    elif A.tag == dTag: lib.ElScaleSparse_d(*args)
    elif A.tag == cTag: lib.ElScaleSparse_c(*args)
    elif A.tag == zTag: lib.ElScaleSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElScaleDistSparse_i(*args)
    elif A.tag == sTag: lib.ElScaleDistSparse_s(*args)
    elif A.tag == dTag: lib.ElScaleDistSparse_d(*args)
    elif A.tag == cTag: lib.ElScaleDistSparse_c(*args)
    elif A.tag == zTag: lib.ElScaleDistSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElScaleDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElScaleDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElScaleDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElScaleDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElScaleDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Scale trapezoid
# ---------------
lib.ElScaleTrapezoid_i.argtypes = [iType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoid_i.restype = c_uint
lib.ElScaleTrapezoid_s.argtypes = [sType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoid_s.restype = c_uint
lib.ElScaleTrapezoid_d.argtypes = [dType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoid_d.restype = c_uint
lib.ElScaleTrapezoid_c.argtypes = [cType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoid_c.restype = c_uint
lib.ElScaleTrapezoid_z.argtypes = [zType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoid_z.restype = c_uint
lib.ElScaleTrapezoidDist_i.argtypes = [iType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDist_i.restype = c_uint
lib.ElScaleTrapezoidDist_s.argtypes = [sType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDist_s.restype = c_uint
lib.ElScaleTrapezoidDist_d.argtypes = [dType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDist_d.restype = c_uint
lib.ElScaleTrapezoidDist_c.argtypes = [cType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDist_c.restype = c_uint
lib.ElScaleTrapezoidDist_z.argtypes = [zType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDist_z.restype = c_uint
lib.ElScaleTrapezoidSparse_i.argtypes = [iType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidSparse_i.restype = c_uint
lib.ElScaleTrapezoidSparse_s.argtypes = [sType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidSparse_s.restype = c_uint
lib.ElScaleTrapezoidSparse_d.argtypes = [dType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidSparse_d.restype = c_uint
lib.ElScaleTrapezoidSparse_c.argtypes = [cType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidSparse_c.restype = c_uint
lib.ElScaleTrapezoidSparse_z.argtypes = [zType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidSparse_z.restype = c_uint
lib.ElScaleTrapezoidDistSparse_i.argtypes = [iType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDistSparse_i.restype = c_uint
lib.ElScaleTrapezoidDistSparse_s.argtypes = [sType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDistSparse_s.restype = c_uint
lib.ElScaleTrapezoidDistSparse_d.argtypes = [dType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDistSparse_d.restype = c_uint
lib.ElScaleTrapezoidDistSparse_c.argtypes = [cType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDistSparse_c.restype = c_uint
lib.ElScaleTrapezoidDistSparse_z.argtypes = [zType,c_uint,c_void_p,iType]
lib.ElScaleTrapezoidDistSparse_z.restype = c_uint
def ScaleTrapezoid(alphaPre,uplo,A,offset=0):
  alpha = TagToType(A.tag)(alphaPre)
  args = [alpha,uplo,A.obj,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElScaleTrapezoid_i(*args)
    elif A.tag == sTag: lib.ElScaleTrapezoid_s(*args)
    elif A.tag == dTag: lib.ElScaleTrapezoid_d(*args)
    elif A.tag == cTag: lib.ElScaleTrapezoid_c(*args)
    elif A.tag == zTag: lib.ElScaleTrapezoid_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElScaleTrapezoidDist_i(*args)
    elif A.tag == sTag: lib.ElScaleTrapezoidDist_s(*args)
    elif A.tag == dTag: lib.ElScaleTrapezoidDist_d(*args)
    elif A.tag == cTag: lib.ElScaleTrapezoidDist_c(*args)
    elif A.tag == zTag: lib.ElScaleTrapezoidDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElScaleTrapezoidSparse_i(*args)
    elif A.tag == sTag: lib.ElScaleTrapezoidSparse_s(*args)
    elif A.tag == dTag: lib.ElScaleTrapezoidSparse_d(*args)
    elif A.tag == cTag: lib.ElScaleTrapezoidSparse_c(*args)
    elif A.tag == zTag: lib.ElScaleTrapezoidSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElScaleTrapezoidDistSparse_i(*args)
    elif A.tag == sTag: lib.ElScaleTrapezoidDistSparse_s(*args)
    elif A.tag == dTag: lib.ElScaleTrapezoidDistSparse_d(*args)
    elif A.tag == cTag: lib.ElScaleTrapezoidDistSparse_c(*args)
    elif A.tag == zTag: lib.ElScaleTrapezoidDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Set diagonal
# ------------
lib.ElSetDiagonal_i.argtypes = [c_void_p,iType,iType]
lib.ElSetDiagonal_i.restype = c_uint
lib.ElSetDiagonal_s.argtypes = [c_void_p,sType,iType]
lib.ElSetDiagonal_s.restype = c_uint
lib.ElSetDiagonal_d.argtypes = [c_void_p,dType,iType]
lib.ElSetDiagonal_d.restype = c_uint
lib.ElSetDiagonal_c.argtypes = [c_void_p,cType,iType]
lib.ElSetDiagonal_c.restype = c_uint
lib.ElSetDiagonal_z.argtypes = [c_void_p,zType,iType]
lib.ElSetDiagonal_z.restype = c_uint
lib.ElSetDiagonalDist_i.argtypes = [c_void_p,iType,iType]
lib.ElSetDiagonalDist_i.restype = c_uint
lib.ElSetDiagonalDist_s.argtypes = [c_void_p,sType,iType]
lib.ElSetDiagonalDist_s.restype = c_uint
lib.ElSetDiagonalDist_d.argtypes = [c_void_p,dType,iType]
lib.ElSetDiagonalDist_d.restype = c_uint
lib.ElSetDiagonalDist_c.argtypes = [c_void_p,cType,iType]
lib.ElSetDiagonalDist_c.restype = c_uint
lib.ElSetDiagonalDist_z.argtypes = [c_void_p,zType,iType]
lib.ElSetDiagonalDist_z.restype = c_uint
def SetDiagonal(A,alphaPre,offset=0):
  alpha = TagToType(A.tag)(alphaPre)
  args = [A.obj,alpha,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSetDiagonal_i(*args)
    elif A.tag == sTag: lib.ElSetDiagonal_s(*args)
    elif A.tag == dTag: lib.ElSetDiagonal_d(*args)
    elif A.tag == cTag: lib.ElSetDiagonal_c(*args)
    elif A.tag == zTag: lib.ElSetDiagonal_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSetDiagonalDist_i(*args)
    elif A.tag == sTag: lib.ElSetDiagonalDist_s(*args)
    elif A.tag == dTag: lib.ElSetDiagonalDist_d(*args)
    elif A.tag == cTag: lib.ElSetDiagonalDist_c(*args)
    elif A.tag == zTag: lib.ElSetDiagonalDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Swap
# ----
lib.ElSwap_i.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwap_i.restype = c_uint
lib.ElSwap_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwap_s.restype = c_uint
lib.ElSwap_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwap_d.restype = c_uint
lib.ElSwap_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwap_c.restype = c_uint
lib.ElSwap_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwap_z.restype = c_uint
lib.ElSwapDist_i.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwapDist_i.restype = c_uint
lib.ElSwapDist_s.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwapDist_s.restype = c_uint
lib.ElSwapDist_d.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwapDist_d.restype = c_uint
lib.ElSwapDist_c.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwapDist_c.restype = c_uint
lib.ElSwapDist_z.argtypes = [c_uint,c_void_p,c_void_p]
lib.ElSwapDist_z.restype = c_uint
def Swap(orient,X,Y):
  if type(X) is not type(Y): raise Exception('Matrix types must match')
  if X.tag != Y.tag: raise Exception('Matrix datatypes must match')
  args = [orient,X.obj,Y.obj]
  if type(X) is Matrix:
    if   A.tag == iTag: lib.ElSwap_i(*args)
    elif A.tag == sTag: lib.ElSwap_s(*args)
    elif A.tag == dTag: lib.ElSwap_d(*args)
    elif A.tag == cTag: lib.ElSwap_c(*args)
    elif A.tag == zTag: lib.ElSwap_z(*args)
    else: DataExcept()
  elif type(X) is DistMatrix:
    if   A.tag == iTag: lib.ElSwapDist_i(*args)
    elif A.tag == sTag: lib.ElSwapDist_s(*args)
    elif A.tag == dTag: lib.ElSwapDist_d(*args)
    elif A.tag == cTag: lib.ElSwapDist_c(*args)
    elif A.tag == zTag: lib.ElSwapDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElRowSwap_i.argtypes = [c_void_p,iType,iType]
lib.ElRowSwap_i.restype = c_uint
lib.ElRowSwap_s.argtypes = [c_void_p,iType,iType]
lib.ElRowSwap_s.restype = c_uint
lib.ElRowSwap_d.argtypes = [c_void_p,iType,iType]
lib.ElRowSwap_d.restype = c_uint
lib.ElRowSwap_c.argtypes = [c_void_p,iType,iType]
lib.ElRowSwap_c.restype = c_uint
lib.ElRowSwap_z.argtypes = [c_void_p,iType,iType]
lib.ElRowSwap_z.restype = c_uint
lib.ElRowSwapDist_i.argtypes = [c_void_p,iType,iType]
lib.ElRowSwapDist_i.restype = c_uint
lib.ElRowSwapDist_s.argtypes = [c_void_p,iType,iType]
lib.ElRowSwapDist_s.restype = c_uint
lib.ElRowSwapDist_d.argtypes = [c_void_p,iType,iType]
lib.ElRowSwapDist_d.restype = c_uint
lib.ElRowSwapDist_c.argtypes = [c_void_p,iType,iType]
lib.ElRowSwapDist_c.restype = c_uint
lib.ElRowSwapDist_z.argtypes = [c_void_p,iType,iType]
lib.ElRowSwapDist_z.restype = c_uint
def RowSwap(A,iTo,iFrom):
  args = [A.obj,iTo,iFrom]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRowSwap_i(*args)
    elif A.tag == sTag: lib.ElRowSwap_s(*args)
    elif A.tag == dTag: lib.ElRowSwap_d(*args)
    elif A.tag == cTag: lib.ElRowSwap_c(*args)
    elif A.tag == zTag: lib.ElRowSwap_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElRowSwapDist_i(*args)
    elif A.tag == sTag: lib.ElRowSwapDist_s(*args)
    elif A.tag == dTag: lib.ElRowSwapDist_d(*args)
    elif A.tag == cTag: lib.ElRowSwapDist_c(*args)
    elif A.tag == zTag: lib.ElRowSwapDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElColSwap_i.argtypes = [c_void_p,iType,iType]
lib.ElColSwap_i.restype = c_uint
lib.ElColSwap_s.argtypes = [c_void_p,iType,iType]
lib.ElColSwap_s.restype = c_uint
lib.ElColSwap_d.argtypes = [c_void_p,iType,iType]
lib.ElColSwap_d.restype = c_uint
lib.ElColSwap_c.argtypes = [c_void_p,iType,iType]
lib.ElColSwap_c.restype = c_uint
lib.ElColSwap_z.argtypes = [c_void_p,iType,iType]
lib.ElColSwap_z.restype = c_uint
lib.ElColSwapDist_i.argtypes = [c_void_p,iType,iType]
lib.ElColSwapDist_i.restype = c_uint
lib.ElColSwapDist_s.argtypes = [c_void_p,iType,iType]
lib.ElColSwapDist_s.restype = c_uint
lib.ElColSwapDist_d.argtypes = [c_void_p,iType,iType]
lib.ElColSwapDist_d.restype = c_uint
lib.ElColSwapDist_c.argtypes = [c_void_p,iType,iType]
lib.ElColSwapDist_c.restype = c_uint
lib.ElColSwapDist_z.argtypes = [c_void_p,iType,iType]
lib.ElColSwapDist_z.restype = c_uint
def ColSwap(A,jTo,jFrom):
  args = [A.obj,jTo,jFrom]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElColSwap_i(*args)
    elif A.tag == sTag: lib.ElColSwap_s(*args)
    elif A.tag == dTag: lib.ElColSwap_d(*args)
    elif A.tag == cTag: lib.ElColSwap_c(*args)
    elif A.tag == zTag: lib.ElColSwap_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElColSwapDist_i(*args)
    elif A.tag == sTag: lib.ElColSwapDist_s(*args)
    elif A.tag == dTag: lib.ElColSwapDist_d(*args)
    elif A.tag == cTag: lib.ElColSwapDist_c(*args)
    elif A.tag == zTag: lib.ElColSwapDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

lib.ElSymmetricSwap_i.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwap_i.restype = c_uint
lib.ElSymmetricSwap_s.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwap_s.restype = c_uint
lib.ElSymmetricSwap_d.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwap_d.restype = c_uint
lib.ElSymmetricSwap_c.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwap_c.restype = c_uint
lib.ElSymmetricSwap_z.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwap_z.restype = c_uint
lib.ElSymmetricSwapDist_i.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwapDist_i.restype = c_uint
lib.ElSymmetricSwapDist_s.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwapDist_s.restype = c_uint
lib.ElSymmetricSwapDist_d.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwapDist_d.restype = c_uint
lib.ElSymmetricSwapDist_c.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwapDist_c.restype = c_uint
lib.ElSymmetricSwapDist_z.argtypes = [c_uint,c_void_p,iType,iType]
lib.ElSymmetricSwapDist_z.restype = c_uint
def SymmetricSwap(uplo,A,jTo,jFrom):
  args = [uplo,A.obj,jTo,jFrom]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElSymmetricSwap_i(*args)
    elif A.tag == sTag: lib.ElSymmetricSwap_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSwap_d(*args)
    elif A.tag == cTag: lib.ElSymmetricSwap_c(*args)
    elif A.tag == zTag: lib.ElSymmetricSwap_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElSymmetricSwapDist_i(*args)
    elif A.tag == sTag: lib.ElSymmetricSwapDist_s(*args)
    elif A.tag == dTag: lib.ElSymmetricSwapDist_d(*args)
    elif A.tag == cTag: lib.ElSymmetricSwapDist_c(*args)
    elif A.tag == zTag: lib.ElSymmetricSwapDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Transpose/Adjoint
# -----------------
lib.ElTranspose_i.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_i.restype = c_uint
lib.ElTranspose_s.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_s.restype = c_uint
lib.ElTranspose_d.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_d.restype = c_uint
lib.ElTranspose_c.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_c.restype = c_uint
lib.ElTranspose_z.argtypes = [c_void_p,c_void_p]
lib.ElTranspose_z.restype = c_uint
lib.ElTransposeDist_i.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_i.restype = c_uint
lib.ElTransposeDist_s.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_s.restype = c_uint
lib.ElTransposeDist_d.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_d.restype = c_uint
lib.ElTransposeDist_c.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_c.restype = c_uint
lib.ElTransposeDist_z.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDist_z.restype = c_uint
lib.ElTransposeSparse_i.argtypes = [c_void_p,c_void_p]
lib.ElTransposeSparse_i.restype = c_uint
lib.ElTransposeSparse_s.argtypes = [c_void_p,c_void_p]
lib.ElTransposeSparse_s.restype = c_uint
lib.ElTransposeSparse_d.argtypes = [c_void_p,c_void_p]
lib.ElTransposeSparse_d.restype = c_uint
lib.ElTransposeSparse_c.argtypes = [c_void_p,c_void_p]
lib.ElTransposeSparse_c.restype = c_uint
lib.ElTransposeSparse_z.argtypes = [c_void_p,c_void_p]
lib.ElTransposeSparse_z.restype = c_uint
lib.ElTransposeDistSparse_i.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDistSparse_i.restype = c_uint
lib.ElTransposeDistSparse_s.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDistSparse_s.restype = c_uint
lib.ElTransposeDistSparse_d.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDistSparse_d.restype = c_uint
lib.ElTransposeDistSparse_c.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDistSparse_c.restype = c_uint
lib.ElTransposeDistSparse_z.argtypes = [c_void_p,c_void_p]
lib.ElTransposeDistSparse_z.restype = c_uint

lib.ElAdjoint_c.argtypes = [c_void_p,c_void_p]
lib.ElAdjoint_c.restype = c_uint
lib.ElAdjoint_z.argtypes = [c_void_p,c_void_p]
lib.ElAdjoint_z.restype = c_uint
lib.ElAdjointDist_c.argtypes = [c_void_p,c_void_p]
lib.ElAdjointDist_c.restype = c_uint
lib.ElAdjointDist_z.argtypes = [c_void_p,c_void_p]
lib.ElAdjointDist_z.restype = c_uint
lib.ElAdjointSparse_c.argtypes = [c_void_p,c_void_p]
lib.ElAdjointSparse_c.restype = c_uint
lib.ElAdjointSparse_z.argtypes = [c_void_p,c_void_p]
lib.ElAdjointSparse_z.restype = c_uint
lib.ElAdjointDistSparse_c.argtypes = [c_void_p,c_void_p]
lib.ElAdjointDistSparse_c.restype = c_uint
lib.ElAdjointDistSparse_z.argtypes = [c_void_p,c_void_p]
lib.ElAdjointDistSparse_z.restype = c_uint

def Transpose(A,B,conj=False):
  if A.tag != B.tag:
    raise Exception('Transposing between datatypes not yet supported in Python')
  if type(A) is not type(B): raise Exception('Matrix types must match')
  args = [A.obj,B.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElTranspose_i(*args)
    elif A.tag == sTag: lib.ElTranspose_s(*args)
    elif A.tag == dTag: lib.ElTranspose_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElAdjoint_c(*args)
      else:    lib.ElTranspose_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElAdjoint_z(*args)
      else:    lib.ElTranspose_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElTransposeDist_i(*args)
    elif A.tag == sTag: lib.ElTransposeDist_s(*args)
    elif A.tag == dTag: lib.ElTransposeDist_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElAdjointDist_c(*args)
      else:    lib.ElTransposeDist_c(*args)
    elif A.tag == zTag: 
      if conj: lib.ElAdjointDist_z(*args)
      else:    lib.ElTransposeDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElTransposeSparse_i(*args)
    elif A.tag == sTag: lib.ElTransposeSparse_s(*args)
    elif A.tag == dTag: lib.ElTransposeSparse_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElAdjointSparse_c(*args)
      else:    lib.ElTransposeSparse_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElAdjointSparse_z(*args)
      else:    lib.ElTransposeSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElTransposeDistSparse_i(*args)
    elif A.tag == sTag: lib.ElTransposeDistSparse_s(*args)
    elif A.tag == dTag: lib.ElTransposeDistSparse_d(*args)
    elif A.tag == cTag:
      if conj: lib.ElAdjointDistSparse_c(*args)
      else:    lib.ElTransposeDistSparse_c(*args)
    elif A.tag == zTag:
      if conj: lib.ElAdjointDistSparse_z(*args)
      else:    lib.ElTransposeDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

def Adjoint(A,B):
  Transpose(A,B,True)

# Real part
# ---------
# TODO: Version which returns the result instead?
lib.ElRealPart_i.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_i.restype = c_uint
lib.ElRealPart_s.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_s.restype = c_uint
lib.ElRealPart_d.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_d.restype = c_uint
lib.ElRealPart_c.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_c.restype = c_uint
lib.ElRealPart_z.argtypes = [c_void_p,c_void_p]
lib.ElRealPart_z.restype = c_uint
lib.ElRealPartDist_i.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_i.restype = c_uint
lib.ElRealPartDist_s.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_s.restype = c_uint
lib.ElRealPartDist_d.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_d.restype = c_uint
lib.ElRealPartDist_c.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_c.restype = c_uint
lib.ElRealPartDist_z.argtypes = [c_void_p,c_void_p]
lib.ElRealPartDist_z.restype = c_uint
def RealPart(A,AReal):
  if AReal.tag != Base(A.tag):
    raise Exception('AReal must have the base datatype of A')
  if type(A) is not type(AReal): raise Exception('Matrix types must match')
  args = [A.obj,AReal.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElRealPart_i(*args)
    elif A.tag == sTag: lib.ElRealPart_s(*args)
    elif A.tag == dTag: lib.ElRealPart_d(*args)
    elif A.tag == cTag: lib.ElRealPart_c(*args)
    elif A.tag == zTag: lib.ElRealPart_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElRealPartDist_i(*args)
    elif A.tag == sTag: lib.ElRealPartDist_s(*args)
    elif A.tag == dTag: lib.ElRealPartDist_d(*args)
    elif A.tag == cTag: lib.ElRealPartDist_c(*args)
    elif A.tag == zTag: lib.ElRealPartDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Imaginary part
# --------------
# TODO: Version which returns the result instead?
lib.ElImagPart_i.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_i.restype = c_uint
lib.ElImagPart_s.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_s.restype = c_uint
lib.ElImagPart_d.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_d.restype = c_uint
lib.ElImagPart_c.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_c.restype = c_uint
lib.ElImagPart_z.argtypes = [c_void_p,c_void_p]
lib.ElImagPart_z.restype = c_uint
lib.ElImagPartDist_i.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_i.restype = c_uint
lib.ElImagPartDist_s.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_s.restype = c_uint
lib.ElImagPartDist_d.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_d.restype = c_uint
lib.ElImagPartDist_c.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_c.restype = c_uint
lib.ElImagPartDist_z.argtypes = [c_void_p,c_void_p]
lib.ElImagPartDist_z.restype = c_uint
def ImagPart(A,AImag):
  if AImag.tag != Base(A.tag):
    raise Exception('AImag must have the base datatype of A')
  if type(A) is not type(AImag): raise Exception('Matrix types must match')
  args = [A.obj,AImag.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElImagPart_i(*args)
    elif A.tag == sTag: lib.ElImagPart_s(*args)
    elif A.tag == dTag: lib.ElImagPart_d(*args)
    elif A.tag == cTag: lib.ElImagPart_c(*args)
    elif A.tag == zTag: lib.ElImagPart_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElImagPartDist_i(*args)
    elif A.tag == sTag: lib.ElImagPartDist_s(*args)
    elif A.tag == dTag: lib.ElImagPartDist_d(*args)
    elif A.tag == cTag: lib.ElImagPartDist_c(*args)
    elif A.tag == zTag: lib.ElImagPartDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Update diagonal
# ---------------
lib.ElUpdateDiagonal_i.argtypes = [c_void_p,iType,iType]
lib.ElUpdateDiagonal_i.restype = c_uint
lib.ElUpdateDiagonal_s.argtypes = [c_void_p,sType,iType]
lib.ElUpdateDiagonal_s.restype = c_uint
lib.ElUpdateDiagonal_d.argtypes = [c_void_p,dType,iType]
lib.ElUpdateDiagonal_d.restype = c_uint
lib.ElUpdateDiagonal_c.argtypes = [c_void_p,cType,iType]
lib.ElUpdateDiagonal_c.restype = c_uint
lib.ElUpdateDiagonal_z.argtypes = [c_void_p,zType,iType]
lib.ElUpdateDiagonal_z.restype = c_uint
lib.ElUpdateDiagonalDist_i.argtypes = [c_void_p,iType,iType]
lib.ElUpdateDiagonalDist_i.restype = c_uint
lib.ElUpdateDiagonalDist_s.argtypes = [c_void_p,sType,iType]
lib.ElUpdateDiagonalDist_s.restype = c_uint
lib.ElUpdateDiagonalDist_d.argtypes = [c_void_p,dType,iType]
lib.ElUpdateDiagonalDist_d.restype = c_uint
lib.ElUpdateDiagonalDist_c.argtypes = [c_void_p,cType,iType]
lib.ElUpdateDiagonalDist_c.restype = c_uint
lib.ElUpdateDiagonalDist_z.argtypes = [c_void_p,zType,iType]
lib.ElUpdateDiagonalDist_z.restype = c_uint
lib.ElUpdateDiagonalSparse_i.argtypes = [c_void_p,iType,iType]
lib.ElUpdateDiagonalSparse_i.restype = c_uint
lib.ElUpdateDiagonalSparse_s.argtypes = [c_void_p,sType,iType]
lib.ElUpdateDiagonalSparse_s.restype = c_uint
lib.ElUpdateDiagonalSparse_d.argtypes = [c_void_p,dType,iType]
lib.ElUpdateDiagonalSparse_d.restype = c_uint
lib.ElUpdateDiagonalSparse_c.argtypes = [c_void_p,cType,iType]
lib.ElUpdateDiagonalSparse_c.restype = c_uint
lib.ElUpdateDiagonalSparse_z.argtypes = [c_void_p,zType,iType]
lib.ElUpdateDiagonalSparse_z.restype = c_uint
lib.ElUpdateDiagonalDistSparse_i.argtypes = [c_void_p,iType,iType]
lib.ElUpdateDiagonalDistSparse_i.restype = c_uint
lib.ElUpdateDiagonalDistSparse_s.argtypes = [c_void_p,sType,iType]
lib.ElUpdateDiagonalDistSparse_s.restype = c_uint
lib.ElUpdateDiagonalDistSparse_d.argtypes = [c_void_p,dType,iType]
lib.ElUpdateDiagonalDistSparse_d.restype = c_uint
lib.ElUpdateDiagonalDistSparse_c.argtypes = [c_void_p,cType,iType]
lib.ElUpdateDiagonalDistSparse_c.restype = c_uint
lib.ElUpdateDiagonalDistSparse_z.argtypes = [c_void_p,zType,iType]
lib.ElUpdateDiagonalDistSparse_z.restype = c_uint
def UpdateDiagonal(A,alphaPre,offset=0):
  alpha = TagToType(A.tag)(alphaPre)
  args = [A.obj,alpha,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElUpdateDiagonal_i(*args)
    elif A.tag == sTag: lib.ElUpdateDiagonal_s(*args)
    elif A.tag == dTag: lib.ElUpdateDiagonal_d(*args)
    elif A.tag == cTag: lib.ElUpdateDiagonal_c(*args)
    elif A.tag == zTag: lib.ElUpdateDiagonal_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElUpdateDiagonalDist_i(*args)
    elif A.tag == sTag: lib.ElUpdateDiagonalDist_s(*args)
    elif A.tag == dTag: lib.ElUpdateDiagonalDist_d(*args)
    elif A.tag == cTag: lib.ElUpdateDiagonalDist_c(*args)
    elif A.tag == zTag: lib.ElUpdateDiagonalDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElUpdateDiagonalSparse_i(*args)
    elif A.tag == sTag: lib.ElUpdateDiagonalSparse_s(*args)
    elif A.tag == dTag: lib.ElUpdateDiagonalSparse_d(*args)
    elif A.tag == cTag: lib.ElUpdateDiagonalSparse_c(*args)
    elif A.tag == zTag: lib.ElUpdateDiagonalSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElUpdateDiagonalDistSparse_i(*args)
    elif A.tag == sTag: lib.ElUpdateDiagonalDistSparse_s(*args)
    elif A.tag == dTag: lib.ElUpdateDiagonalDistSparse_d(*args)
    elif A.tag == cTag: lib.ElUpdateDiagonalDistSparse_c(*args)
    elif A.tag == zTag: lib.ElUpdateDiagonalDistSparse_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Zero
# ----
lib.ElZero_i.argtypes = [c_void_p]
lib.ElZero_i.restype = c_uint
lib.ElZero_s.argtypes = [c_void_p]
lib.ElZero_s.restype = c_uint
lib.ElZero_d.argtypes = [c_void_p]
lib.ElZero_d.restype = c_uint
lib.ElZero_c.argtypes = [c_void_p]
lib.ElZero_c.restype = c_uint
lib.ElZero_z.argtypes = [c_void_p]
lib.ElZero_z.restype = c_uint
lib.ElZeroDist_i.argtypes = [c_void_p]
lib.ElZeroDist_i.restype = c_uint
lib.ElZeroDist_s.argtypes = [c_void_p]
lib.ElZeroDist_s.restype = c_uint
lib.ElZeroDist_d.argtypes = [c_void_p]
lib.ElZeroDist_d.restype = c_uint
lib.ElZeroDist_c.argtypes = [c_void_p]
lib.ElZeroDist_c.restype = c_uint
lib.ElZeroDist_z.argtypes = [c_void_p]
lib.ElZeroDist_z.restype = c_uint
lib.ElZeroSparse_i.argtypes = [c_void_p]
lib.ElZeroSparse_i.restype = c_uint
lib.ElZeroSparse_s.argtypes = [c_void_p]
lib.ElZeroSparse_s.restype = c_uint
lib.ElZeroSparse_d.argtypes = [c_void_p]
lib.ElZeroSparse_d.restype = c_uint
lib.ElZeroSparse_c.argtypes = [c_void_p]
lib.ElZeroSparse_c.restype = c_uint
lib.ElZeroSparse_z.argtypes = [c_void_p]
lib.ElZeroSparse_z.restype = c_uint
lib.ElZeroDistSparse_i.argtypes = [c_void_p]
lib.ElZeroDistSparse_i.restype = c_uint
lib.ElZeroDistSparse_s.argtypes = [c_void_p]
lib.ElZeroDistSparse_s.restype = c_uint
lib.ElZeroDistSparse_d.argtypes = [c_void_p]
lib.ElZeroDistSparse_d.restype = c_uint
lib.ElZeroDistSparse_c.argtypes = [c_void_p]
lib.ElZeroDistSparse_c.restype = c_uint
lib.ElZeroDistSparse_z.argtypes = [c_void_p]
lib.ElZeroDistSparse_z.restype = c_uint
lib.ElZeroDistMultiVec_i.argtypes = [c_void_p]
lib.ElZeroDistMultiVec_i.restype = c_uint
lib.ElZeroDistMultiVec_s.argtypes = [c_void_p]
lib.ElZeroDistMultiVec_s.restype = c_uint
lib.ElZeroDistMultiVec_d.argtypes = [c_void_p]
lib.ElZeroDistMultiVec_d.restype = c_uint
lib.ElZeroDistMultiVec_c.argtypes = [c_void_p]
lib.ElZeroDistMultiVec_c.restype = c_uint
lib.ElZeroDistMultiVec_z.argtypes = [c_void_p]
lib.ElZeroDistMultiVec_z.restype = c_uint

def Zero(A):
  args = [A.obj]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElZero_i(*args)
    elif A.tag == sTag: lib.ElZero_s(*args)
    elif A.tag == dTag: lib.ElZero_d(*args)
    elif A.tag == cTag: lib.ElZero_c(*args)
    elif A.tag == zTag: lib.ElZero_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElZeroDist_i(*args)
    elif A.tag == sTag: lib.ElZeroDist_s(*args)
    elif A.tag == dTag: lib.ElZeroDist_d(*args)
    elif A.tag == cTag: lib.ElZeroDist_c(*args)
    elif A.tag == zTag: lib.ElZeroDist_z(*args)
    else: DataExcept()
  elif type(A) is SparseMatrix:
    if   A.tag == iTag: lib.ElZeroSparse_i(*args)
    elif A.tag == sTag: lib.ElZeroSparse_s(*args)
    elif A.tag == dTag: lib.ElZeroSparse_d(*args)
    elif A.tag == cTag: lib.ElZeroSparse_c(*args)
    elif A.tag == zTag: lib.ElZeroSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistSparseMatrix:
    if   A.tag == iTag: lib.ElZeroDistSparse_i(*args)
    elif A.tag == sTag: lib.ElZeroDistSparse_s(*args)
    elif A.tag == dTag: lib.ElZeroDistSparse_d(*args)
    elif A.tag == cTag: lib.ElZeroDistSparse_c(*args)
    elif A.tag == zTag: lib.ElZeroDistSparse_z(*args)
    else: DataExcept()
  elif type(A) is DistMultiVec:
    if   A.tag == iTag: lib.ElZeroDistMultiVec_i(*args)
    elif A.tag == sTag: lib.ElZeroDistMultiVec_s(*args)
    elif A.tag == dTag: lib.ElZeroDistMultiVec_d(*args)
    elif A.tag == cTag: lib.ElZeroDistMultiVec_c(*args)
    elif A.tag == zTag: lib.ElZeroDistMultiVec_z(*args)
    else: DataExcept()
  else: TypeExcept()
