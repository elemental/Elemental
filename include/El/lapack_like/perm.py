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

lib.ElPermutationMetaSet.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermutationMetaSet.restype = c_uint
lib.ElPermutationMetaClear.argtypes = [c_void_p]
lib.ElPermutationMetaClear.restype = c_uint
lib.ElPermutationMetaTotalSend.argtypes = [c_void_p,POINTER(c_int)]
lib.ElPermutationMetaTotalSend.restype = c_uint
lib.ElPermutationMetaTotalRecv.argtypes = [c_void_p,POINTER(c_int)]
lib.ElPermutationMetaTotalRecv.restype = c_uint
lib.ElPermutationMetaScaleUp.argtypes = [c_void_p,iType]
lib.ElPermutationMetaScaleUp.restype = c_uint
lib.ElPermutationMetaScaleDown.argtypes = [c_void_p,iType]
lib.ElPermutationMetaScaleDown.restype = c_uint
class PermutationMeta(ctypes.Structure):
  _fields_ = [("align",iType),("comm",MPI_Comm),
              ("sendCounts",POINTER(iType)),("sendDispls",POINTER(iType)),
              ("recvCounts",POINTER(iType)),("recvDispls",POINTER(iType)),
              ("recvIdx",POINTER(iType)),   ("recvRanks",POINTER(iType))]
  def __init__(self,p,pInv):
    if type(p) is not DistMatrix or type(pInv) is not DistMatrix:
      raise Exception('Types of p and pInv must be DistMatrix')
    if p.tag != iType or pInv.tag != iType:
      raise Exception('p and pInv must be integral')
    lib.ElPermutationMetaSet(p.obj,pInv,obj,pointer(self))
  def Set(self,p,pInv):
    if type(p) is not DistMatrix or type(pInv) is not DistMatrix:
      raise Exception('Types of p and pInv must be DistMatrix')
    if p.tag != iType or pInv.tag != iType:
      raise Exception('p and pInv must be integral')
    lib.ElPermutationMetaSet(p.obj,pInv,obj,pointer(self))
  def Clear(self,p,pInv):
    lib.ElPermutationMetaClear(pointer(self))
  def TotalSend(self):
    total = c_int()
    lib.ElPermutationMetaTotalSend(pointer(self),pointer(total))
    return total
  def TotalRecv(self):
    total = c_int()
    lib.ElPermutationMetaTotalRecv(pointer(self),pointer(total))
    return total
  def ScaleUp(self,length):
    lib.ElPermutationMetaScaleUp(pointer(self),length)
  def ScaleDown(self,length):
    lib.ElPermutationMetaScaleDown(pointer(self),length)

# Apply column pivots
# ===================
lib.ElApplyColPivots_i.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivots_i.restype = c_uint
lib.ElApplyColPivots_s.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivots_s.restype = c_uint
lib.ElApplyColPivots_d.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivots_d.restype = c_uint
lib.ElApplyColPivots_c.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivots_c.restype = c_uint
lib.ElApplyColPivots_z.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivots_z.restype = c_uint
lib.ElApplyColPivotsDist_i.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivotsDist_i.restype = c_uint
lib.ElApplyColPivotsDist_s.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivotsDist_s.restype = c_uint
lib.ElApplyColPivotsDist_d.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivotsDist_d.restype = c_uint
lib.ElApplyColPivotsDist_c.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivotsDist_c.restype = c_uint
lib.ElApplyColPivotsDist_z.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyColPivotsDist_z.restype = c_uint
def ApplyColPivots(A,pivots,offset=0):
  if type(A) is not type(pivots):
    raise Exception('Types of A and pivots must match')
  if pivots.tag != iTag:
    raise Exception('pivots must be integral')
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElApplyColPivots_i(A.obj,pivots.obj,offset)
    elif A.tag == sTag: lib.ElApplyColPivots_s(A.obj,pivots.obj,offset)
    elif A.tag == dTag: lib.ElApplyColPivots_d(A.obj,pivots.obj,offset)
    elif A.tag == cTag: lib.ElApplyColPivots_c(A.obj,pivots.obj,offset)
    elif A.tag == zTag: lib.ElApplyColPivots_z(A.obj,pivots.obj,offset)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElApplyColPivotsDist_i(A.obj,pivots.obj,offset)
    elif A.tag == sTag: lib.ElApplyColPivotsDist_s(A.obj,pivots.obj,offset)
    elif A.tag == dTag: lib.ElApplyColPivotsDist_d(A.obj,pivots.obj,offset)
    elif A.tag == cTag: lib.ElApplyColPivotsDist_c(A.obj,pivots.obj,offset)
    elif A.tag == zTag: lib.ElApplyColPivotsDist_z(A.obj,pivots.obj,offset)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Apply row pivots
# ================
lib.ElApplyRowPivots_i.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivots_i.restype = c_uint
lib.ElApplyRowPivots_s.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivots_s.restype = c_uint
lib.ElApplyRowPivots_d.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivots_d.restype = c_uint
lib.ElApplyRowPivots_c.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivots_c.restype = c_uint
lib.ElApplyRowPivots_z.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivots_z.restype = c_uint
lib.ElApplyRowPivotsDist_i.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivotsDist_i.restype = c_uint
lib.ElApplyRowPivotsDist_s.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivotsDist_s.restype = c_uint
lib.ElApplyRowPivotsDist_d.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivotsDist_d.restype = c_uint
lib.ElApplyRowPivotsDist_c.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivotsDist_c.restype = c_uint
lib.ElApplyRowPivotsDist_z.argtypes = [c_void_p,c_void_p,iType]
lib.ElApplyRowPivotsDist_z.restype = c_uint
def ApplyRowPivots(A,pivots,offset=0):
  if type(A) is not type(pivots):
    raise Exception('Types of A and pivots must match')
  if pivots.tag != iTag:
    raise Exception('pivots must be integral')
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElApplyRowPivots_i(A.obj,pivots.obj,offset)
    elif A.tag == sTag: lib.ElApplyRowPivots_s(A.obj,pivots.obj,offset)
    elif A.tag == dTag: lib.ElApplyRowPivots_d(A.obj,pivots.obj,offset)
    elif A.tag == cTag: lib.ElApplyRowPivots_c(A.obj,pivots.obj,offset)
    elif A.tag == zTag: lib.ElApplyRowPivots_z(A.obj,pivots.obj,offset)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElApplyRowPivotsDist_i(A.obj,pivots.obj,offset)
    elif A.tag == sTag: lib.ElApplyRowPivotsDist_s(A.obj,pivots.obj,offset)
    elif A.tag == dTag: lib.ElApplyRowPivotsDist_d(A.obj,pivots.obj,offset)
    elif A.tag == cTag: lib.ElApplyRowPivotsDist_c(A.obj,pivots.obj,offset)
    elif A.tag == zTag: lib.ElApplyRowPivotsDist_z(A.obj,pivots.obj,offset)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Apply symmetric pivots
# ======================
lib.ElApplySymmetricPivots_i.argtypes = [c_uint,c_void_p,c_void_p,iType]
lib.ElApplySymmetricPivots_i.restype = c_uint
lib.ElApplySymmetricPivots_s.argtypes = [c_uint,c_void_p,c_void_p,iType]
lib.ElApplySymmetricPivots_s.restype = c_uint
lib.ElApplySymmetricPivots_d.argtypes = [c_uint,c_void_p,c_void_p,iType]
lib.ElApplySymmetricPivots_d.restype = c_uint
lib.ElApplySymmetricPivots_c.argtypes = [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplySymmetricPivots_c.restype = c_uint
lib.ElApplySymmetricPivots_z.argtypes = [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplySymmetricPivots_z.restype = c_uint
lib.ElApplySymmetricPivotsDist_i.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplySymmetricPivotsDist_i.restype = c_uint
lib.ElApplySymmetricPivotsDist_s.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplySymmetricPivotsDist_s.restype = c_uint
lib.ElApplySymmetricPivotsDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplySymmetricPivotsDist_d.restype = c_uint
lib.ElApplySymmetricPivotsDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplySymmetricPivotsDist_c.restype = c_uint
lib.ElApplySymmetricPivotsDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplySymmetricPivotsDist_z.restype = c_uint
def ApplySymmetricPivots(uplo,A,p,conjugate=False,offset=0):
  if type(A) is not type(p):
    raise Exception('Types of A and p must match')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == iTag: 
      lib.ElApplySymmetricPivots_i(uplo,A.obj,pivots.obj,offset)
    elif A.tag == sTag: 
      lib.ElApplySymmetricPivots_s(uplo,A.obj,pivots.obj,offset)
    elif A.tag == dTag: 
      lib.ElApplySymmetricPivots_d(uplo,A.obj,pivots.obj,offset)
    elif A.tag == cTag: 
      lib.ElApplySymmetricPivots_c(uplo,A.obj,pivots.obj,conjugate,offset)
    elif A.tag == zTag: 
      lib.ElApplySymmetricPivots_z(uplo,A.obj,pivots.obj,conjugate,offset)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: 
      lib.ElApplySymmetricPivotsDist_i(uplo,A.obj,pivots.obj,offset)
    elif A.tag == sTag: 
      lib.ElApplySymmetricPivotsDist_s(uplo,A.obj,pivots.obj,offset)
    elif A.tag == dTag: 
      lib.ElApplySymmetricPivotsDist_d(uplo,A.obj,pivots.obj,offset)
    elif A.tag == cTag: 
      lib.ElApplySymmetricPivotsDist_c(uplo,A.obj,pivots.obj,conjugate,offset)
    elif A.tag == zTag: 
      lib.ElApplySymmetricPivotsDist_z(uplo,A.obj,pivots.obj,conjugate,offset)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

lib.ElApplyInverseSymmetricPivots_i.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplyInverseSymmetricPivots_i.restype = c_uint
lib.ElApplyInverseSymmetricPivots_s.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplyInverseSymmetricPivots_s.restype = c_uint
lib.ElApplyInverseSymmetricPivots_d.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplyInverseSymmetricPivots_d.restype = c_uint
lib.ElApplyInverseSymmetricPivots_c.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplyInverseSymmetricPivots_c.restype = c_uint
lib.ElApplyInverseSymmetricPivots_z.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplyInverseSymmetricPivots_z.restype = c_uint
lib.ElApplyInverseSymmetricPivotsDist_i.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplyInverseSymmetricPivotsDist_i.restype = c_uint
lib.ElApplyInverseSymmetricPivotsDist_s.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplyInverseSymmetricPivotsDist_s.restype = c_uint
lib.ElApplyInverseSymmetricPivotsDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]
lib.ElApplyInverseSymmetricPivotsDist_d.restype = c_uint
lib.ElApplyInverseSymmetricPivotsDist_c.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplyInverseSymmetricPivotsDist_c.restype = c_uint
lib.ElApplyInverseSymmetricPivotsDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]
lib.ElApplyInverseSymmetricPivotsDist_z.restype = c_uint
def ApplyInverseSymmetricPivots(uplo,A,p,conjugate=False,offset=0):
  if type(A) is not type(p):
    raise Exception('Types of A and p must match')
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(A) is Matrix:
    if   A.tag == iTag: 
      lib.ElApplyInverseSymmetricPivots_i(uplo,A.obj,pivots.obj,offset)
    elif A.tag == sTag: 
      lib.ElApplyInverseSymmetricPivots_s(uplo,A.obj,pivots.obj,offset)
    elif A.tag == dTag: 
      lib.ElApplyInverseSymmetricPivots_d(uplo,A.obj,pivots.obj,offset)
    elif A.tag == cTag: 
      lib.ElApplyInverseSymmetricPivots_c \
      (uplo,A.obj,pivots.obj,conjugate,offset)
    elif A.tag == zTag: 
      lib.ElApplyInverseSymmetricPivots_z \
      (uplo,A.obj,pivots.obj,conjugate,offset)
    else: raise Exception('Unsupported datatype')
  elif type(A) is DistMatrix:
    if   A.tag == iTag: 
      lib.ElApplyInverseSymmetricPivotsDist_i(uplo,A.obj,pivots.obj,offset)
    elif A.tag == sTag: 
      lib.ElApplyInverseSymmetricPivotsDist_s(uplo,A.obj,pivots.obj,offset)
    elif A.tag == dTag: 
      lib.ElApplyInverseSymmetricPivotsDist_d(uplo,A.obj,pivots.obj,offset)
    elif A.tag == cTag: 
      lib.ElApplyInverseSymmetricPivotsDist_c \
      (uplo,A.obj,pivots.obj,conjugate,offset)
    elif A.tag == zTag: 
      lib.ElApplyInverseSymmetricPivotsDist_z \
      (uplo,A.obj,pivots.obj,conjugate,offset)
    else: raise Exception('Unsupported datatype')
  else: raise Exception('Unsupported matrix type')

# Explicit permutation
# ====================
# TODO

# Invert permutation
# ==================
# TODO

# Permute columns
# ===============
# TODO

# Permute rows
# ============
# TODO

# Pivot parity
# ============
# TODO

# Convert a pivot sequence to a partial permutation vector
# ========================================================
# TODO

# Convert a pivot sequence to a permutation vector
# ================================================
# TODO
