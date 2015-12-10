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

class Permutation(object):
  lib.ElPermutationCreate.argtypes = [POINTER(c_void_p)]
  def __init__(self):
    self.obj = c_void_p()
    lib.ElPermutationCreate(pointer(self.obj))

  lib.ElPermutationDestroy.argtypes = [c_void_p]
  def Destroy(self):
    lib.ElPermutationDestroy(self.obj)  

  # TODO: More member functions, e.g., PermuteRows and PermuteCols

class DistPermutation(object):
  lib.ElDistPermutationCreate.argtypes = [POINTER(c_void_p),c_void_p]
  def __init__(self,grid):
    self.obj = c_void_p()
    lib.ElDistPermutationCreate(pointer(self.obj),grid.obj)

  lib.ElDistPermutationDestroy.argtypes = [c_void_p]
  def Destroy(self):
    lib.ElDistPermutationDestroy(self.obj)

  # TODO: More member functions, e.g., PermuteRows and PermuteCols

