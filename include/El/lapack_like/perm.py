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
  _fields_ = [("align",iType),("comm",mpi.Comm),
              ("sendCounts",POINTER(iType)),("sendDispls",POINTER(iType)),
              ("recvCounts",POINTER(iType)),("recvDispls",POINTER(iType)),
              ("numSendIdx",iType),         ("numRecvIdx",iType),
              ("sendIdx",POINTER(iType)),   ("sendRanks",POINTER(iType)),
              ("recvIdx",POINTER(iType)),   ("recvRanks",POINTER(iType))]
  def __init__(self,p,pInv):
    if type(p) is not DistMatrix or type(pInv) is not DistMatrix:
      raise Exception('Types of p and pInv must be DistMatrix')
    if p.tag != iTag or pInv.tag != iTag:
      raise Exception('p and pInv must be integral')
    lib.ElPermutationMetaSet(p.obj,pInv.obj,pointer(self))
  def Set(self,p,pInv):
    if type(p) is not DistMatrix or type(pInv) is not DistMatrix:
      raise Exception('Types of p and pInv must be DistMatrix')
    if p.tag != iTag or pInv.tag != iTag:
      raise Exception('p and pInv must be integral')
    lib.ElPermutationMetaSet(p.obj,pInv.obj,pointer(self))
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
  args = [A.obj,pivots.obj,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElApplyColPivots_i(*args)
    elif A.tag == sTag: lib.ElApplyColPivots_s(*args)
    elif A.tag == dTag: lib.ElApplyColPivots_d(*args)
    elif A.tag == cTag: lib.ElApplyColPivots_c(*args)
    elif A.tag == zTag: lib.ElApplyColPivots_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElApplyColPivotsDist_i(*args)
    elif A.tag == sTag: lib.ElApplyColPivotsDist_s(*args)
    elif A.tag == dTag: lib.ElApplyColPivotsDist_d(*args)
    elif A.tag == cTag: lib.ElApplyColPivotsDist_c(*args)
    elif A.tag == zTag: lib.ElApplyColPivotsDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
  args = [A.obj,pivots.obj,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElApplyRowPivots_i(*args)
    elif A.tag == sTag: lib.ElApplyRowPivots_s(*args)
    elif A.tag == dTag: lib.ElApplyRowPivots_d(*args)
    elif A.tag == cTag: lib.ElApplyRowPivots_c(*args)
    elif A.tag == zTag: lib.ElApplyRowPivots_z(*args)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElApplyRowPivotsDist_i(*args)
    elif A.tag == sTag: lib.ElApplyRowPivotsDist_s(*args)
    elif A.tag == dTag: lib.ElApplyRowPivotsDist_d(*args)
    elif A.tag == cTag: lib.ElApplyRowPivotsDist_c(*args)
    elif A.tag == zTag: lib.ElApplyRowPivotsDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

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
  args = [uplo,A.obj,pivots.obj,offset]
  argsCpx = [uplo,A.obj,pivots.obj,conjugate,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElApplySymmetricPivots_i(*args)
    elif A.tag == sTag: lib.ElApplySymmetricPivots_s(*args)
    elif A.tag == dTag: lib.ElApplySymmetricPivots_d(*args)
    elif A.tag == cTag: lib.ElApplySymmetricPivots_c(*argsCpx)
    elif A.tag == zTag: lib.ElApplySymmetricPivots_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElApplySymmetricPivotsDist_i(*args)
    elif A.tag == sTag: lib.ElApplySymmetricPivotsDist_s(*args)
    elif A.tag == dTag: lib.ElApplySymmetricPivotsDist_d(*args)
    elif A.tag == cTag: lib.ElApplySymmetricPivotsDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElApplySymmetricPivotsDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

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
  args = [uplo,A.obj,pivots.obj,offset]
  argsCpx = [uplo,A.obj,pivots.obj,conjugate,offset]
  if type(A) is Matrix:
    if   A.tag == iTag: lib.ElApplyInverseSymmetricPivots_i(*args)
    elif A.tag == sTag: lib.ElApplyInverseSymmetricPivots_s(*args)
    elif A.tag == dTag: lib.ElApplyInverseSymmetricPivots_d(*args)
    elif A.tag == cTag: lib.ElApplyInverseSymmetricPivots_c(*argsCpx)
    elif A.tag == zTag: lib.ElApplyInverseSymmetricPivots_z(*argsCpx)
    else: DataExcept()
  elif type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElApplyInverseSymmetricPivotsDist_i(*args)
    elif A.tag == sTag: lib.ElApplyInverseSymmetricPivotsDist_s(*args)
    elif A.tag == dTag: lib.ElApplyInverseSymmetricPivotsDist_d(*args)
    elif A.tag == cTag: lib.ElApplyInverseSymmetricPivotsDist_c(*argsCpx)
    elif A.tag == zTag: lib.ElApplyInverseSymmetricPivotsDist_z(*argsCpx)
    else: DataExcept()
  else: TypeExcept()

# Explicit permutation
# ====================
lib.ElExplicitPermutation.argtypes = [c_void_p,c_void_p]
lib.ElExplicitPermutation.restype = c_uint
lib.ElExplicitPermutationDist.argtypes = [c_void_p,c_void_p]
lib.ElExplicitPermutationDist.restype = c_uint
def ExplicitPermutation(p):
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(p) is Matrix:
    P = Matrix(iTag)
    lib.ElExplicitPermutation(p.obj,P.obj)
    return P
  elif type(p) is DistMatrix:
    P = DistMatrix(iTag,MC,MR,p.Grid())
    lib.ElExplicitPermutationDist(p.obj,P.obj)
    return P
  else: TypeExcept()

# Invert permutation
# ==================
lib.ElInvertPermutation.argtypes = [c_void_p,c_void_p]
lib.ElInvertPermutation.restype = c_uint
lib.ElInvertPermutationDist.argtypes = [c_void_p,c_void_p]
lib.ElInvertPermutationDist.restype = c_uint
def InvertPermutation(p):
  if p.tag != iTag:
    raise Exception('p must be integral')
  if type(p) is Matrix:
    pInv = Matrix(iTag)
    lib.ElInvertPermutation(p.obj,pInv.obj)
    return pInv
  elif type(p) is DistMatrix:
    pInv = DistMatrix(iTag,VC,STAR,p.Grid())
    lib.ElInvertPermutationDist(p.obj,pInv.obj)
    return pInv
  else: TypeExcept()

# Parity of a permutation
# =======================
lib.ElPermutationParity.argtypes = [c_void_p,POINTER(bType)]
lib.ElPermutationParity.restype = c_uint
lib.ElPermutationParityDist.argtypes = [c_void_p,POINTER(bType)]
lib.ElPermutationParityDist.restype = c_uint
def PermutationParity(p):
  if p.tag != iTag:
    raise Exception('p must be integral')
  parity = bType()
  if type(p) is Matrix:
    lib.ElPermutationParity(p.obj,pointer(parity))
  elif type(p) is DistMatrix:
    lib.ElPermutationParityDist(p.obj,pointer(parity))
  else: TypeExcept()
  return parity

# Permute columns
# ===============
lib.ElPermuteCols_i.argtypes = [c_void_p,c_void_p]
lib.ElPermuteCols_i.restype = c_uint
lib.ElPermuteCols_s.argtypes = [c_void_p,c_void_p]
lib.ElPermuteCols_s.restype = c_uint
lib.ElPermuteCols_d.argtypes = [c_void_p,c_void_p]
lib.ElPermuteCols_d.restype = c_uint
lib.ElPermuteCols_c.argtypes = [c_void_p,c_void_p]
lib.ElPermuteCols_c.restype = c_uint
lib.ElPermuteCols_z.argtypes = [c_void_p,c_void_p]
lib.ElPermuteCols_z.restype = c_uint
lib.ElPermuteColsDist_i.argtypes = [c_void_p,c_void_p]
lib.ElPermuteColsDist_i.restype = c_uint
lib.ElPermuteColsDist_s.argtypes = [c_void_p,c_void_p]
lib.ElPermuteColsDist_s.restype = c_uint
lib.ElPermuteColsDist_d.argtypes = [c_void_p,c_void_p]
lib.ElPermuteColsDist_d.restype = c_uint
lib.ElPermuteColsDist_c.argtypes = [c_void_p,c_void_p]
lib.ElPermuteColsDist_c.restype = c_uint
lib.ElPermuteColsDist_z.argtypes = [c_void_p,c_void_p]
lib.ElPermuteColsDist_z.restype = c_uint
lib.ElPermuteColsBoth_i.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBoth_i.restype = c_uint
lib.ElPermuteColsBoth_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBoth_s.restype = c_uint
lib.ElPermuteColsBoth_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBoth_d.restype = c_uint
lib.ElPermuteColsBoth_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBoth_c.restype = c_uint
lib.ElPermuteColsBoth_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBoth_z.restype = c_uint
lib.ElPermuteColsBothDist_i.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBothDist_i.restype = c_uint
lib.ElPermuteColsBothDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBothDist_s.restype = c_uint
lib.ElPermuteColsBothDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBothDist_d.restype = c_uint
lib.ElPermuteColsBothDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBothDist_c.restype = c_uint
lib.ElPermuteColsBothDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteColsBothDist_z.restype = c_uint
def PermuteCols(A,p,pInv=None):
  if type(A) is Matrix:
    if pInv == None:
      if   A.tag == iTag: lib.ElPermuteCols_i(A.obj,p.obj)
      elif A.tag == sTag: lib.ElPermuteCols_s(A.obj,p.obj)
      elif A.tag == dTag: lib.ElPermuteCols_d(A.obj,p.obj)
      elif A.tag == cTag: lib.ElPermuteCols_c(A.obj,p.obj)
      elif A.tag == zTag: lib.ElPermuteCols_z(A.obj,p.obj)
      else: DataExcept()
    else:
      if type(pInv) != Matrix or pInv.tag != iTag:
        raise Exception('pInv must be an integer Matrix')
      if   A.tag == iTag: lib.ElPermuteColsBoth_i(A.obj,p.obj,pInv.obj)
      elif A.tag == sTag: lib.ElPermuteColsBoth_s(A.obj,p.obj,pInv.obj)
      elif A.tag == dTag: lib.ElPermuteColsBoth_d(A.obj,p.obj,pInv.obj)
      elif A.tag == cTag: lib.ElPermuteColsBoth_c(A.obj,p.obj,pInv.obj)
      elif A.tag == zTag: lib.ElPermuteColsBoth_z(A.obj,p.obj,pInv.obj)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if pInv == None:
      if   A.tag == iTag: lib.ElPermuteColsDist_i(A.obj,p.obj)
      elif A.tag == sTag: lib.ElPermuteColsDist_s(A.obj,p.obj)
      elif A.tag == dTag: lib.ElPermuteColsDist_d(A.obj,p.obj)
      elif A.tag == cTag: lib.ElPermuteColsDist_c(A.obj,p.obj)
      elif A.tag == zTag: lib.ElPermuteColsDist_z(A.obj,p.obj)
      else: DataExcept()
    else:
      if type(pInv) != Matrix or pInv.tag != iTag:
        raise Exception('pInv must be an integer Matrix')
      if   A.tag == iTag: lib.ElPermuteColsBothDist_i(A.obj,p.obj,pInv.obj)
      elif A.tag == sTag: lib.ElPermuteColsBothDist_s(A.obj,p.obj,pInv.obj)
      elif A.tag == dTag: lib.ElPermuteColsBothDist_d(A.obj,p.obj,pInv.obj)
      elif A.tag == cTag: lib.ElPermuteColsBothDist_c(A.obj,p.obj,pInv.obj)
      elif A.tag == zTag: lib.ElPermuteColsBothDist_z(A.obj,p.obj,pInv.obj)
      else: DataExcept()
  else: TypeExcept()

lib.ElPermuteColsMetaDist_i.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteColsMetaDist_i.restype = c_uint
lib.ElPermuteColsMetaDist_s.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteColsMetaDist_s.restype = c_uint
lib.ElPermuteColsMetaDist_d.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteColsMetaDist_d.restype = c_uint
lib.ElPermuteColsMetaDist_c.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteColsMetaDist_c.restype = c_uint
lib.ElPermuteColsMetaDist_z.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteColsMetaDist_z.restype = c_uint
def PermuteColsMeta(A,meta):
  if type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElPermuteColsMetaDist_i(A.obj,pointer(meta))
    elif A.tag == sTag: lib.ElPermuteColsMetaDist_s(A.obj,pointer(meta))
    elif A.tag == dTag: lib.ElPermuteColsMetaDist_d(A.obj,pointer(meta))
    elif A.tag == cTag: lib.ElPermuteColsMetaDist_c(A.obj,pointer(meta))
    elif A.tag == zTag: lib.ElPermuteColsMetaDist_z(A.obj,pointer(meta))
    else: DataExcept()
  else: TypeExcept()

# Permute rows
# ============
lib.ElPermuteRows_i.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRows_i.restype = c_uint
lib.ElPermuteRows_s.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRows_s.restype = c_uint
lib.ElPermuteRows_d.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRows_d.restype = c_uint
lib.ElPermuteRows_c.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRows_c.restype = c_uint
lib.ElPermuteRows_z.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRows_z.restype = c_uint
lib.ElPermuteRowsDist_i.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRowsDist_i.restype = c_uint
lib.ElPermuteRowsDist_s.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRowsDist_s.restype = c_uint
lib.ElPermuteRowsDist_d.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRowsDist_d.restype = c_uint
lib.ElPermuteRowsDist_c.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRowsDist_c.restype = c_uint
lib.ElPermuteRowsDist_z.argtypes = [c_void_p,c_void_p]
lib.ElPermuteRowsDist_z.restype = c_uint
lib.ElPermuteRowsBoth_i.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBoth_i.restype = c_uint
lib.ElPermuteRowsBoth_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBoth_s.restype = c_uint
lib.ElPermuteRowsBoth_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBoth_d.restype = c_uint
lib.ElPermuteRowsBoth_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBoth_c.restype = c_uint
lib.ElPermuteRowsBoth_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBoth_z.restype = c_uint
lib.ElPermuteRowsBothDist_i.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBothDist_i.restype = c_uint
lib.ElPermuteRowsBothDist_s.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBothDist_s.restype = c_uint
lib.ElPermuteRowsBothDist_d.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBothDist_d.restype = c_uint
lib.ElPermuteRowsBothDist_c.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBothDist_c.restype = c_uint
lib.ElPermuteRowsBothDist_z.argtypes = [c_void_p,c_void_p,c_void_p]
lib.ElPermuteRowsBothDist_z.restype = c_uint
def PermuteRows(A,p,pInv=None):
  if type(A) is Matrix:
    if pInv == None:
      args = [A.obj,p.obj]
      if   A.tag == iTag: lib.ElPermuteRows_i(*args)
      elif A.tag == sTag: lib.ElPermuteRows_s(*args)
      elif A.tag == dTag: lib.ElPermuteRows_d(*args)
      elif A.tag == cTag: lib.ElPermuteRows_c(*args)
      elif A.tag == zTag: lib.ElPermuteRows_z(*args)
      else: DataExcept()
    else:
      if type(pInv) != Matrix or pInv.tag != iTag:
        raise Exception('pInv must be integral')
      args = [A.obj,p.obj,pInv.obj]
      if   A.tag == iTag: lib.ElPermuteRowsBoth_i(*args)
      elif A.tag == sTag: lib.ElPermuteRowsBoth_s(*args)
      elif A.tag == dTag: lib.ElPermuteRowsBoth_d(*args)
      elif A.tag == cTag: lib.ElPermuteRowsBoth_c(*args)
      elif A.tag == zTag: lib.ElPermuteRowsBoth_z(*args)
      else: DataExcept()
  elif type(A) is DistMatrix:
    if pInv == None:
      args = [A.obj,p.obj]
      if   A.tag == iTag: lib.ElPermuteRowsDist_i(*args)
      elif A.tag == sTag: lib.ElPermuteRowsDist_s(*args)
      elif A.tag == dTag: lib.ElPermuteRowsDist_d(*args)
      elif A.tag == cTag: lib.ElPermuteRowsDist_c(*args)
      elif A.tag == zTag: lib.ElPermuteRowsDist_z(*args)
      else: DataExcept()
    else:
      if type(pInv) != Matrix or pInv.tag != iTag:
        raise Exception('pInv must be integral')
      args = [A.obj,p.obj,pInv.obj]
      if   A.tag == iTag: lib.ElPermuteRowsBothDist_i(*args)
      elif A.tag == sTag: lib.ElPermuteRowsBothDist_s(*args)
      elif A.tag == dTag: lib.ElPermuteRowsBothDist_d(*args)
      elif A.tag == cTag: lib.ElPermuteRowsBothDist_c(*args)
      elif A.tag == zTag: lib.ElPermuteRowsBothDist_z(*args)
      else: DataExcept()
  else: TypeExcept()

lib.ElPermuteRowsMetaDist_i.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteRowsMetaDist_i.restype = c_uint
lib.ElPermuteRowsMetaDist_s.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteRowsMetaDist_s.restype = c_uint
lib.ElPermuteRowsMetaDist_d.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteRowsMetaDist_d.restype = c_uint
lib.ElPermuteRowsMetaDist_c.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteRowsMetaDist_c.restype = c_uint
lib.ElPermuteRowsMetaDist_z.argtypes = [c_void_p,POINTER(PermutationMeta)]
lib.ElPermuteRowsMetaDist_z.restype = c_uint
def PermuteRowsMeta(A,meta):
  args = [A.obj,pointer(meta)]
  if type(A) is DistMatrix:
    if   A.tag == iTag: lib.ElPermuteRowsMetaDist_i(*args)
    elif A.tag == sTag: lib.ElPermuteRowsMetaDist_s(*args)
    elif A.tag == dTag: lib.ElPermuteRowsMetaDist_d(*args)
    elif A.tag == cTag: lib.ElPermuteRowsMetaDist_c(*args)
    elif A.tag == zTag: lib.ElPermuteRowsMetaDist_z(*args)
    else: DataExcept()
  else: TypeExcept()

# Pivot parity
# ============
lib.ElPivotParity.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElPivotParity.restype = c_uint
lib.ElPivotParityDist.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElPivotParityDist.restype = c_uint
def PivotParity(p,offset=0):
  if p.tag != iTag:
    raise Exception('p must be integral')
  parity = bType()
  args = [p.obj,offset,pointer(parity)]
  if type(p) is Matrix:       lib.ElPivotParity(*args)
  elif type(p) is DistMatrix: lib.ElPivotParityDist(*args)
  else: TypeExcept()
  return parity

# Convert a pivot sequence to a partial permutation vector
# ========================================================
lib.ElPivotsToPartialPermutation.argtypes = [c_void_p,c_void_p,c_void_p,iType]
lib.ElPivotsToPartialPermutation.restype = c_uint
lib.ElPivotsToPartialPermutationDist.argtypes = \
  [c_void_p,c_void_p,c_void_p,iType]
lib.ElPivotsToPartialPermutationDist.restype = c_uint
def PivotsToPartialPermutation(pivots,offset=0):
  if pivots.tag != iTag:
    raise Exception('pivots must be integral')
  if type(pivots) is Matrix:
    p = Matrix(iTag)
    pInv = Matrix(iTag)
    lib.ElPivotsToPartialPermutation(pivots.obj,p.obj,pInv.obj,offset)
    return p, pInv
  elif type(pivots) is DistMatrix:
    p = DistMatrix(iTag,VC,STAR,pivots.Grid())
    pInv = DistMatrix(iTag,VC,STAR,pivots.Grid())
    lib.ElPivotsToPartialPermutationDist(pivots.obj,p.obj,pInv.obj,offset)
    return p, pInv
  else: TypeExcept()

# Convert a pivot sequence to a permutation vector
# ================================================
lib.ElPivotsToPermutation.argtypes = [c_void_p,c_void_p,iType]
lib.ElPivotsToPermutation.restype = c_uint
lib.ElPivotsToPermutationDist.argtypes = [c_void_p,c_void_p,iType]
lib.ElPivotsToPermutationDist.restype = c_uint
def PivotsToPermutation(pivots,offset=0):
  if pivots.tag != iTag:
    raise Exception('pivots must be integral')
  if type(pivots) is Matrix:
    p = Matrix(iTag)
    lib.ElPivotsToPermutation(pivots.obj,p.obj,offset)
    return p
  elif type(pivots) is DistMatrix:
    p = DistMatrix(iTag,VC,STAR,pivots.Grid())
    lib.ElPivotsToPermutationDist(pivots.obj,p.obj,offset)
    return p
  else: TypeExcept()

lib.ElPivotsToInversePermutation.argtypes = [c_void_p,c_void_p,iType]
lib.ElPivotsToInversePermutation.restype = c_uint
lib.ElPivotsToInversePermutationDist.argtypes = [c_void_p,c_void_p,iType]
lib.ElPivotsToInversePermutationDist.restype = c_uint
def PivotsToInversePermutation(pivots,offset=0):
  if pivots.tag != iTag:
    raise Exception('pivots must be integral')
  if type(pivots) is Matrix:
    pInv = Matrix(iTag)
    lib.ElPivotsToInversePermutation(pivots.obj,pInv.obj,offset)
    return pInv
  elif type(pivots) is DistMatrix:
    pInv = DistMatrix(iTag,VC,STAR,pivots.Grid())
    lib.ElPivotsToInversePermutationDist(pivots.obj,pInv.obj,offset)
    return pInv
  else: TypeExcept()

