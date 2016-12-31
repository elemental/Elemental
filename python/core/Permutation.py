#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License,
#  which can be found in the LICENSE file in the root directory, or at
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
from imports     import mpi
from Matrix import Matrix
from DistMatrix import DistMatrix
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

class Permutation(object):
  lib.ElPermutationCreate.argtypes = [POINTER(c_void_p)]
  def __init__(self):
    self.obj = c_void_p()
    lib.ElPermutationCreate(pointer(self.obj))

  lib.ElPermutationDestroy.argtypes = [c_void_p]
  def Destroy(self):
    lib.ElPermutationDestroy(self.obj)

  lib.ElPermutationEmpty.argtypes = [c_void_p]
  def Empty(self):
    lib.ElPermutationEmpty(self.obj)

  lib.ElPermutationMakeIdentity.argtypes = [c_void_p,iType]
  def MakeIdentity(self,size):
    lib.ElPermutationMakeIdentity(self.obj,size)

  lib.ElPermutationReserveSwaps.argtypes = [c_void_p,iType]
  def ReserveSwaps(self,maxSwaps):
    lib.ElPermutationReserveSwaps(self.obj,maxSwaps)

  lib.ElPermutationMakeArbitrary.argtypes = [c_void_p]
  def MakeArbitrary(self):
    lib.ElPermutationMakeArbitrary(self.obj)

  lib.ElPermutationSwap.argtypes = [c_void_p,iType,iType]
  def Swap(self,origin,dest):
    lib.ElPermutationSwap(self.obj,origin,dest)

  lib.ElPermutationSwapSequence.argtypes = [c_void_p,c_void_p,iType]
  def SwapSequence(self,PAppend,offset=0):
    lib.ElPermutationSwapSequence(self.obj,PAppend.obj,offset)

  lib.ElPermutationSetImage.argtypes = [c_void_p,iType,iType]
  def SetImage(self,origin,dest):
    lib.ElPermutationSetImage(self.obj,origin,dest)

  lib.ElPermutationHeight.argtypes = [c_void_p,POINTER(iType)]
  def Height(self):
    height = iType()
    lib.ElPermutationHeight(self.obj,pointer(height))
    return height.value

  lib.ElPermutationWidth.argtypes = [c_void_p,POINTER(iType)]
  def Width(self):
    width = iType()
    lib.ElPermutationWidth(self.obj,pointer(width))
    return width.value

  lib.ElPermutationParity.argtypes = [c_void_p,POINTER(bType)]
  def Parity(self):
    parity = bType()
    lib.ElPermutationParity(self.obj,pointer(parity))
    return parity.value

  lib.ElPermutationIsSwapSequence.argtypes = [c_void_p,POINTER(bType)]
  def IsSwapSequence(self):
    isSwap = bType()
    lib.ElPermutationIsSwapSequence(self.obj,pointer(isSwap))
    return isSwap.value

  lib.ElPermutationIsImplicitSwapSequence.argtypes = [c_void_p,POINTER(bType)]
  def IsImplicitSwapSequence(self):
    isImplicit = bType()
    lib.ElPermutationIsImplicitSwapSequence(self.obj,pointer(isImplicit))
    return isImplicit.value

  lib.ElPermutationPermuteCols_i.argtypes = \
  lib.ElPermutationPermuteCols_s.argtypes = \
  lib.ElPermutationPermuteCols_d.argtypes = \
  lib.ElPermutationPermuteCols_c.argtypes = \
  lib.ElPermutationPermuteCols_z.argtypes = \
    [c_void_p,c_void_p,iType]
  def PermuteCols(A,offset=0):
    if type(A) is not Matrix:
      raise Exception('Expected A to be a Matrix')
    args = [self.obj,A.obj,offset]
    if   A.tag == iTag: lib.ElPermutationPermuteCols_i(*args)
    elif A.tag == sTag: lib.ElPermutationPermuteCols_s(*args)
    elif A.tag == dTag: lib.ElPermutationPermuteCols_d(*args)
    elif A.tag == cTag: lib.ElPermutationPermuteCols_c(*args)
    elif A.tag == zTag: lib.ElPermutationPermuteCols_z(*args)
    else: DataExcept()

  lib.ElPermutationPermuteRows_i.argtypes = \
  lib.ElPermutationPermuteRows_s.argtypes = \
  lib.ElPermutationPermuteRows_d.argtypes = \
  lib.ElPermutationPermuteRows_c.argtypes = \
  lib.ElPermutationPermuteRows_z.argtypes = \
    [c_void_p,c_void_p,iType]
  def PermuteRows(A,offset=0):
    if type(A) is not Matrix:
      raise Exception('Expected A to be a Matrix')
    args = [self.obj,A.obj,offset]
    if   A.tag == iTag: lib.ElPermutationPermuteRows_i(*args)
    elif A.tag == sTag: lib.ElPermutationPermuteRows_s(*args)
    elif A.tag == dTag: lib.ElPermutationPermuteRows_d(*args)
    elif A.tag == cTag: lib.ElPermutationPermuteRows_c(*args)
    elif A.tag == zTag: lib.ElPermutationPermuteRows_z(*args)
    else: DataExcept()

  lib.ElPermutationPermuteSymmetrically_i.argtypes = \
  lib.ElPermutationPermuteSymmetrically_s.argtypes = \
  lib.ElPermutationPermuteSymmetrically_d.argtypes = \
    [c_void_p,c_uint,c_void_p,iType]
  lib.ElPermutationPermuteSymmetrically_c.argtypes = \
  lib.ElPermutationPermuteSymmetrically_z.argtypes = \
    [c_void_p,c_uint,c_void_p,bType,iType]
  def PermuteSymmetrically(A,uplo=LOWER,conjugate=True,offset=0):
    if type(A) is not Matrix:
      raise Exception('Expected A to be a Matrix')
    args = [self.obj,uplo,A.obj,offset]
    argsCpx = [self.obj,uplo,A.obj,conjugate,offset]
    if   A.tag == iTag: lib.ElPermutationPermuteSymmetrically_i(*args)
    elif A.tag == sTag: lib.ElPermutationPermuteSymmetrically_s(*args)
    elif A.tag == dTag: lib.ElPermutationPermuteSymmetrically_d(*args)
    elif A.tag == cTag: lib.ElPermutationPermuteSymmetrically_c(*argsCpx)
    elif A.tag == zTag: lib.ElPermutationPermuteSymmetrically_z(*argsCpx)
    else: DataExcept()

  lib.ElPermutationExplicitVector.argtypes = [c_void_p,c_void_p]
  def ExplicitVector(self):
    p = Matrix(iTag)
    lib.ElPermutationExplicitVector(self.obj,p.obj)
    return p

  lib.ElPermutationExplicitMatrix.argtypes = [c_void_p,c_void_p]
  def ExplicitMatrix(self):
    P = Matrix(iTag)
    lib.ElPermutationExplicitMatrix(self.obj,P.obj)
    return P

class DistPermutation(object):
  lib.ElDistPermutationCreate.argtypes = [POINTER(c_void_p),c_void_p]
  def __init__(self,grid):
    self.obj = c_void_p()
    lib.ElDistPermutationCreate(pointer(self.obj),grid.obj)

  lib.ElDistPermutationDestroy.argtypes = [c_void_p]
  def Destroy(self):
    lib.ElDistPermutationDestroy(self.obj)

  lib.ElDistPermutationEmpty.argtypes = [c_void_p]
  def Empty(self):
    lib.ElDistPermutationEmpty(self.obj)

  lib.ElDistPermutationMakeIdentity.argtypes = [c_void_p,iType]
  def MakeIdentity(self,size):
    lib.ElDistPermutationMakeIdentity(self.obj,size)

  lib.ElDistPermutationReserveSwaps.argtypes = [c_void_p,iType]
  def ReserveSwaps(self,maxSwaps):
    lib.ElDistPermutationReserveSwaps(self.obj,maxSwaps)

  lib.ElDistPermutationMakeArbitrary.argtypes = [c_void_p]
  def MakeArbitrary(self):
    lib.ElDistPermutationMakeArbitrary(self.obj)

  lib.ElDistPermutationSwap.argtypes = [c_void_p,iType,iType]
  def Swap(self,origin,dest):
    lib.ElDistPermutationSwap(self.obj,origin,dest)

  lib.ElDistPermutationSwapSequence.argtypes = [c_void_p,c_void_p,iType]
  def SwapSequence(self,PAppend,offset=0):
    lib.ElDistPermutationSwapSequence(self.obj,PAppend.obj,offset)

  lib.ElDistPermutationSetImage.argtypes = [c_void_p,iType,iType]
  def SetImage(self,origin,dest):
    lib.ElDistPermutationSetImage(self.obj,origin,dest)

  lib.ElDistPermutationHeight.argtypes = [c_void_p,POINTER(iType)]
  def Height(self):
    height = iType()
    lib.ElDistPermutationHeight(self.obj,pointer(height))
    return height.value

  lib.ElDistPermutationWidth.argtypes = [c_void_p,POINTER(iType)]
  def Width(self):
    width = iType()
    lib.ElDistPermutationWidth(self.obj,pointer(width))
    return width.value

  lib.ElDistPermutationParity.argtypes = [c_void_p,POINTER(bType)]
  def Parity(self):
    parity = bType()
    lib.ElDistPermutationParity(self.obj,pointer(parity))
    return parity.value

  lib.ElDistPermutationIsSwapSequence.argtypes = [c_void_p,POINTER(bType)]
  def IsSwapSequence(self):
    isSwap = bType()
    lib.ElDistPermutationIsSwapSequence(self.obj,pointer(isSwap))
    return isSwap.value

  lib.ElDistPermutationIsImplicitSwapSequence.argtypes = \
    [c_void_p,POINTER(bType)]
  def IsImplicitSwapSequence(self):
    isImplicit = bType()
    lib.ElDistPermutationIsImplicitSwapSequence(self.obj,pointer(isImplicit))
    return isImplicit.value

  lib.ElDistPermutationPermuteCols_i.argtypes = \
  lib.ElDistPermutationPermuteCols_s.argtypes = \
  lib.ElDistPermutationPermuteCols_d.argtypes = \
  lib.ElDistPermutationPermuteCols_c.argtypes = \
  lib.ElDistPermutationPermuteCols_z.argtypes = \
    [c_void_p,c_void_p,iType]
  def PermuteCols(A,offset=0):
    if type(A) is not DistMatrix:
      raise Exception('Expected A to be a DistMatrix')
    args = [self.obj,A.obj,offset]
    if   A.tag == iTag: lib.ElDistPermutationPermuteCols_i(*args)
    elif A.tag == sTag: lib.ElDistPermutationPermuteCols_s(*args)
    elif A.tag == dTag: lib.ElDistPermutationPermuteCols_d(*args)
    elif A.tag == cTag: lib.ElDistPermutationPermuteCols_c(*args)
    elif A.tag == zTag: lib.ElDistPermutationPermuteCols_z(*args)
    else: DataExcept()

  lib.ElDistPermutationPermuteRows_i.argtypes = \
  lib.ElDistPermutationPermuteRows_s.argtypes = \
  lib.ElDistPermutationPermuteRows_d.argtypes = \
  lib.ElDistPermutationPermuteRows_c.argtypes = \
  lib.ElDistPermutationPermuteRows_z.argtypes = \
    [c_void_p,c_void_p,iType]
  def PermuteRows(A,offset=0):
    if type(A) is not DistMatrix:
      raise Exception('Expected A to be a DistMatrix')
    args = [self.obj,A.obj,offset]
    if   A.tag == iTag: lib.ElDistPermutationPermuteRows_i(*args)
    elif A.tag == sTag: lib.ElDistPermutationPermuteRows_s(*args)
    elif A.tag == dTag: lib.ElDistPermutationPermuteRows_d(*args)
    elif A.tag == cTag: lib.ElDistPermutationPermuteRows_c(*args)
    elif A.tag == zTag: lib.ElDistPermutationPermuteRows_z(*args)
    else: DataExcept()

  lib.ElDistPermutationPermuteSymmetrically_i.argtypes = \
  lib.ElDistPermutationPermuteSymmetrically_s.argtypes = \
  lib.ElDistPermutationPermuteSymmetrically_d.argtypes = \
    [c_void_p,c_uint,c_void_p,iType]
  lib.ElDistPermutationPermuteSymmetrically_c.argtypes = \
  lib.ElDistPermutationPermuteSymmetrically_z.argtypes = \
    [c_void_p,c_uint,c_void_p,bType,iType]
  def PermuteSymmetrically(A,uplo=LOWER,conjugate=True,offset=0):
    if type(A) is not DistMatrix:
      raise Exception('Expected A to be a DistMatrix')
    args = [self.obj,uplo,A.obj,offset]
    argsCpx = [self.obj,uplo,A.obj,conjugate,offset]
    if   A.tag == iTag: lib.ElDistPermutationPermuteSymmetrically_i(*args)
    elif A.tag == sTag: lib.ElDistPermutationPermuteSymmetrically_s(*args)
    elif A.tag == dTag: lib.ElDistPermutationPermuteSymmetrically_d(*args)
    elif A.tag == cTag: lib.ElDistPermutationPermuteSymmetrically_c(*argsCpx)
    elif A.tag == zTag: lib.ElDistPermutationPermuteSymmetrically_z(*argsCpx)
    else: DataExcept()

  lib.ElDistPermutationExplicitVector.argtypes = [c_void_p,c_void_p]
  def ExplicitVector(self):
    p = DistMatrix(iTag,VC,STAR)
    lib.ElDistPermutationExplicitVector(self.obj,p.obj)
    return p

  lib.ElDistPermutationExplicitMatrix.argtypes = [c_void_p,c_void_p]
  def ExplicitMatrix(self):
    P = DistMatrix(iTag,VC,STAR)
    lib.ElDistPermutationExplicitMatrix(self.obj,P.obj)
    return P
