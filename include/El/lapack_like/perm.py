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
lib.ElPermutationMetaClear.argtypes = [c_void_p]
lib.ElPermutationMetaTotalSend.argtypes = [c_void_p,POINTER(c_int)]
lib.ElPermutationMetaTotalRecv.argtypes = [c_void_p,POINTER(c_int)]
lib.ElPermutationMetaScaleUp.argtypes = [c_void_p,iType]
lib.ElPermutationMetaScaleDown.argtypes = [c_void_p,iType]
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
lib.ElApplyColPivots_i.argtypes = \
lib.ElApplyColPivots_s.argtypes = \
lib.ElApplyColPivots_d.argtypes = \
lib.ElApplyColPivots_c.argtypes = \
lib.ElApplyColPivots_z.argtypes = \
lib.ElApplyColPivotsDist_i.argtypes = \
lib.ElApplyColPivotsDist_s.argtypes = \
lib.ElApplyColPivotsDist_d.argtypes = \
lib.ElApplyColPivotsDist_c.argtypes = \
lib.ElApplyColPivotsDist_z.argtypes = \
  [c_void_p,c_void_p,iType]
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
lib.ElApplyRowPivots_i.argtypes = \
lib.ElApplyRowPivots_s.argtypes = \
lib.ElApplyRowPivots_d.argtypes = \
lib.ElApplyRowPivots_c.argtypes = \
lib.ElApplyRowPivots_z.argtypes = \
lib.ElApplyRowPivotsDist_i.argtypes = \
lib.ElApplyRowPivotsDist_s.argtypes = \
lib.ElApplyRowPivotsDist_d.argtypes = \
lib.ElApplyRowPivotsDist_c.argtypes = \
lib.ElApplyRowPivotsDist_z.argtypes = \
  [c_void_p,c_void_p,iType]
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
lib.ElApplySymmetricPivots_i.argtypes = \
lib.ElApplySymmetricPivots_s.argtypes = \
lib.ElApplySymmetricPivots_d.argtypes = \
lib.ElApplySymmetricPivotsDist_i.argtypes = \
lib.ElApplySymmetricPivotsDist_s.argtypes = \
lib.ElApplySymmetricPivotsDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]

lib.ElApplySymmetricPivots_c.argtypes = \
lib.ElApplySymmetricPivots_z.argtypes = \
lib.ElApplySymmetricPivotsDist_c.argtypes = \
lib.ElApplySymmetricPivotsDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]

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
lib.ElApplyInverseSymmetricPivots_s.argtypes = \
lib.ElApplyInverseSymmetricPivots_d.argtypes = \
lib.ElApplyInverseSymmetricPivotsDist_i.argtypes = \
lib.ElApplyInverseSymmetricPivotsDist_s.argtypes = \
lib.ElApplyInverseSymmetricPivotsDist_d.argtypes = \
  [c_uint,c_void_p,c_void_p,iType]

lib.ElApplyInverseSymmetricPivots_c.argtypes = \
lib.ElApplyInverseSymmetricPivots_z.argtypes = \
lib.ElApplyInverseSymmetricPivotsDist_c.argtypes = \
lib.ElApplyInverseSymmetricPivotsDist_z.argtypes = \
  [c_uint,c_void_p,c_void_p,bType,iType]

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
lib.ElExplicitPermutationDist.argtypes = [c_void_p,c_void_p]
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
lib.ElInvertPermutationDist.argtypes = [c_void_p,c_void_p]
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
lib.ElPermutationParityDist.argtypes = [c_void_p,POINTER(bType)]
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
lib.ElPermuteCols_i.argtypes = \
lib.ElPermuteCols_s.argtypes = \
lib.ElPermuteCols_d.argtypes = \
lib.ElPermuteCols_c.argtypes = \
lib.ElPermuteCols_z.argtypes = \
lib.ElPermuteColsDist_i.argtypes = \
lib.ElPermuteColsDist_s.argtypes = \
lib.ElPermuteColsDist_d.argtypes = \
lib.ElPermuteColsDist_c.argtypes = \
lib.ElPermuteColsDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElPermuteColsBoth_i.argtypes = \
lib.ElPermuteColsBoth_s.argtypes = \
lib.ElPermuteColsBoth_d.argtypes = \
lib.ElPermuteColsBoth_c.argtypes = \
lib.ElPermuteColsBoth_z.argtypes = \
lib.ElPermuteColsBothDist_i.argtypes = \
lib.ElPermuteColsBothDist_s.argtypes = \
lib.ElPermuteColsBothDist_d.argtypes = \
lib.ElPermuteColsBothDist_c.argtypes = \
lib.ElPermuteColsBothDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

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

lib.ElPermuteColsMetaDist_i.argtypes = \
lib.ElPermuteColsMetaDist_s.argtypes = \
lib.ElPermuteColsMetaDist_d.argtypes = \
lib.ElPermuteColsMetaDist_c.argtypes = \
lib.ElPermuteColsMetaDist_z.argtypes = \
  [c_void_p,POINTER(PermutationMeta)]

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
lib.ElPermuteRows_i.argtypes = \
lib.ElPermuteRows_s.argtypes = \
lib.ElPermuteRows_d.argtypes = \
lib.ElPermuteRows_c.argtypes = \
lib.ElPermuteRows_z.argtypes = \
lib.ElPermuteRowsDist_i.argtypes = \
lib.ElPermuteRowsDist_s.argtypes = \
lib.ElPermuteRowsDist_d.argtypes = \
lib.ElPermuteRowsDist_c.argtypes = \
lib.ElPermuteRowsDist_z.argtypes = \
  [c_void_p,c_void_p]

lib.ElPermuteRowsBoth_i.argtypes = \
lib.ElPermuteRowsBoth_s.argtypes = \
lib.ElPermuteRowsBoth_d.argtypes = \
lib.ElPermuteRowsBoth_c.argtypes = \
lib.ElPermuteRowsBoth_z.argtypes = \
lib.ElPermuteRowsBothDist_i.argtypes = \
lib.ElPermuteRowsBothDist_s.argtypes = \
lib.ElPermuteRowsBothDist_d.argtypes = \
lib.ElPermuteRowsBothDist_c.argtypes = \
lib.ElPermuteRowsBothDist_z.argtypes = \
  [c_void_p,c_void_p,c_void_p]

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

lib.ElPermuteRowsMetaDist_i.argtypes = \
lib.ElPermuteRowsMetaDist_s.argtypes = \
lib.ElPermuteRowsMetaDist_d.argtypes = \
lib.ElPermuteRowsMetaDist_c.argtypes = \
lib.ElPermuteRowsMetaDist_z.argtypes = \
  [c_void_p,POINTER(PermutationMeta)]

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
lib.ElPivotParityDist.argtypes = [c_void_p,iType,POINTER(bType)]
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
lib.ElPivotsToPartialPermutation.argtypes = \
lib.ElPivotsToPartialPermutationDist.argtypes = \
  [c_void_p,c_void_p,c_void_p,iType]
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
lib.ElPivotsToPermutation.argtypes = \
lib.ElPivotsToPermutationDist.argtypes = \
  [c_void_p,c_void_p,iType]
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

lib.ElPivotsToInversePermutation.argtypes = \
lib.ElPivotsToInversePermutationDist.argtypes = \
  [c_void_p,c_void_p,iType]
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
