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

import Grid

import Matrix as M

class DistMultiVec(object):
  # Constructors and destructors
  # ============================
  lib.ElDistMultiVecCreate_i.argtypes = \
  lib.ElDistMultiVecCreate_s.argtypes = \
  lib.ElDistMultiVecCreate_d.argtypes = \
  lib.ElDistMultiVecCreate_c.argtypes = \
  lib.ElDistMultiVecCreate_z.argtypes = \
    [POINTER(c_void_p),c_void_p]
  def __init__(self,tag=dTag,grid=Grid.Grid.Default(),create=True):
    self.obj = c_void_p()
    self.tag = tag
    CheckTag(tag)
    if create:
      args = [pointer(self.obj),grid.obj]
      if   tag == iTag: lib.ElDistMultiVecCreate_i(*args)
      elif tag == sTag: lib.ElDistMultiVecCreate_s(*args)
      elif tag == dTag: lib.ElDistMultiVecCreate_d(*args)
      elif tag == cTag: lib.ElDistMultiVecCreate_c(*args)
      elif tag == zTag: lib.ElDistMultiVecCreate_z(*args)
      else: DataExcept()

  lib.ElDistMultiVecDestroy_i.argtypes = \
  lib.ElDistMultiVecDestroy_s.argtypes = \
  lib.ElDistMultiVecDestroy_d.argtypes = \
  lib.ElDistMultiVecDestroy_c.argtypes = \
  lib.ElDistMultiVecDestroy_z.argtypes = \
    [c_void_p]
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMultiVecDestroy_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecDestroy_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecDestroy_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecDestroy_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecDestroy_z(*args)
    else: DataExcept()

  # Assignment and reconfiguration
  # ==============================
  lib.ElDistMultiVecEmpty_i.argtypes = \
  lib.ElDistMultiVecEmpty_s.argtypes = \
  lib.ElDistMultiVecEmpty_d.argtypes = \
  lib.ElDistMultiVecEmpty_c.argtypes = \
  lib.ElDistMultiVecEmpty_z.argtypes = \
    [c_void_p]
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMultiVecEmpty_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecEmpty_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecEmpty_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecEmpty_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecEmpty_z(*args)
    else: DataExcept()

  lib.ElDistMultiVecResize_i.argtypes = \
  lib.ElDistMultiVecResize_s.argtypes = \
  lib.ElDistMultiVecResize_d.argtypes = \
  lib.ElDistMultiVecResize_c.argtypes = \
  lib.ElDistMultiVecResize_z.argtypes = \
    [c_void_p,iType,iType]
  def Resize(self,height,width):
    args = [self.obj,height,width]
    if   self.tag == iTag: lib.ElDistMultiVecResize_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecResize_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecResize_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecResize_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecResize_z(*args)
    else: DataExcept()

  lib.ElDistMultiVecSetGrid_i.argtypes = \
  lib.ElDistMultiVecSetGrid_s.argtypes = \
  lib.ElDistMultiVecSetGrid_d.argtypes = \
  lib.ElDistMultiVecSetGrid_c.argtypes = \
  lib.ElDistMultiVecSetGrid_z.argtypes = \
    [c_void_p,c_void_p]
  def SetGrid(self,grid):
    args = [self.obj,grid.obj]
    if   self.tag == iTag: lib.ElDistMultiVecSetGrid_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecSetGrid_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecSetGrid_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecSetGrid_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecSetGrid_z(*args)
    else: DataExcept()

  # Queries
  # =======
  lib.ElDistMultiVecHeight_i.argtypes = \
  lib.ElDistMultiVecHeight_s.argtypes = \
  lib.ElDistMultiVecHeight_d.argtypes = \
  lib.ElDistMultiVecHeight_c.argtypes = \
  lib.ElDistMultiVecHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Height(self):
    height = iType()
    args = [self.obj,pointer(height)]
    if   self.tag == iTag: lib.ElDistMultiVecHeight_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecHeight_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecHeight_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecHeight_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecHeight_z(*args)
    else: DataExcept()
    return height.value

  lib.ElDistMultiVecWidth_i.argtypes = \
  lib.ElDistMultiVecWidth_s.argtypes = \
  lib.ElDistMultiVecWidth_d.argtypes = \
  lib.ElDistMultiVecWidth_c.argtypes = \
  lib.ElDistMultiVecWidth_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Width(self):
    width = iType()
    args = [self.obj,pointer(width)]
    if   self.tag == iTag: lib.ElDistMultiVecWidth_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecWidth_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecWidth_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecWidth_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecWidth_z(*args)
    else: DataExcept()
    return width.value

  lib.ElDistMultiVecFirstLocalRow_i.argtypes = \
  lib.ElDistMultiVecFirstLocalRow_s.argtypes = \
  lib.ElDistMultiVecFirstLocalRow_d.argtypes = \
  lib.ElDistMultiVecFirstLocalRow_c.argtypes = \
  lib.ElDistMultiVecFirstLocalRow_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def FirstLocalRow(self):
    firstLocalRow = iType()
    args = [self.obj,pointer(firstLocalRow)]
    if   self.tag == iTag: lib.ElDistMultiVecFirstLocalRow_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecFirstLocalRow_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecFirstLocalRow_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecFirstLocalRow_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecFirstLocalRow_z(*args)
    else: DataExcept()
    return firstLocalRow.value

  lib.ElDistMultiVecLocalHeight_i.argtypes = \
  lib.ElDistMultiVecLocalHeight_s.argtypes = \
  lib.ElDistMultiVecLocalHeight_d.argtypes = \
  lib.ElDistMultiVecLocalHeight_c.argtypes = \
  lib.ElDistMultiVecLocalHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def LocalHeight(self):
    localHeight = iType()
    args = [self.obj,pointer(localHeight)]
    if   self.tag == iTag: lib.ElDistMultiVecLocalHeight_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecLocalHeight_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecLocalHeight_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecLocalHeight_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecLocalHeight_z(*args)
    else: DataExcept()
    return localHeight.value

  lib.ElDistMultiVecMatrix_i.argtypes = \
  lib.ElDistMultiVecMatrix_s.argtypes = \
  lib.ElDistMultiVecMatrix_d.argtypes = \
  lib.ElDistMultiVecMatrix_c.argtypes = \
  lib.ElDistMultiVecMatrix_z.argtypes = \
  lib.ElDistMultiVecLockedMatrix_i.argtypes = \
  lib.ElDistMultiVecLockedMatrix_s.argtypes = \
  lib.ElDistMultiVecLockedMatrix_d.argtypes = \
  lib.ElDistMultiVecLockedMatrix_c.argtypes = \
  lib.ElDistMultiVecLockedMatrix_z.argtypes = \
    [c_void_p,POINTER(c_void_p)]
  def Matrix(self,locked=False):
    A = M.Matrix(dTag,False)
    args = [self.obj,pointer(A.obj)]
    if locked:
      if   self.tag == iTag: lib.ElDistMultiVecLockedMatrix_i(*args)
      elif self.tag == sTag: lib.ElDistMultiVecLockedMatrix_s(*args)
      elif self.tag == dTag: lib.ElDistMultiVecLockedMatrix_d(*args)
      elif self.tag == cTag: lib.ElDistMultiVecLockedMatrix_c(*args)
      elif self.tag == zTag: lib.ElDistMultiVecLockedMatrix_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistMultiVecMatrix_i(*args)
      elif self.tag == sTag: lib.ElDistMultiVecMatrix_s(*args)
      elif self.tag == dTag: lib.ElDistMultiVecMatrix_d(*args)
      elif self.tag == cTag: lib.ElDistMultiVecMatrix_c(*args)
      elif self.tag == zTag: lib.ElDistMultiVecMatrix_z(*args)
      else: DataExcept()
    return A

  lib.ElDistMultiVecGrid_i.argtypes = \
  lib.ElDistMultiVecGrid_s.argtypes = \
  lib.ElDistMultiVecGrid_d.argtypes = \
  lib.ElDistMultiVecGrid_c.argtypes = \
  lib.ElDistMultiVecGrid_z.argtypes = \
    [c_void_p,POINTER(c_void_p)]
  def Grid(self):
    grid = Grid.Grid(create=False)
    args = [self.obj,pointer(grid.obj)]
    if   self.tag == iTag: lib.ElDistMultiVecGrid_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecGrid_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecGrid_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecGrid_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecGrid_z(*args)
    else: DataExcept()
    return grid

  lib.ElDistMultiVecBlocksize_i.argtypes = \
  lib.ElDistMultiVecBlocksize_s.argtypes = \
  lib.ElDistMultiVecBlocksize_d.argtypes = \
  lib.ElDistMultiVecBlocksize_c.argtypes = \
  lib.ElDistMultiVecBlocksize_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Blocksize(self):
    blocksize = iType()
    args = [self.obj,pointer(blocksize)]
    if   self.tag == iTag: lib.ElDistMultiVecBlocksize_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecBlocksize_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecBlocksize_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecBlocksize_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecBlocksize_z(*args)
    else: DataExcept()
    return blocksize.value

  lib.ElDistMultiVecRowOwner_i.argtypes = \
  lib.ElDistMultiVecRowOwner_s.argtypes = \
  lib.ElDistMultiVecRowOwner_d.argtypes = \
  lib.ElDistMultiVecRowOwner_c.argtypes = \
  lib.ElDistMultiVecRowOwner_z.argtypes = \
    [c_void_p,iType,POINTER(c_int)]
  def RowOwner(self,i):
    owner = c_int()
    args = [self.obj,i,pointer(owner)]
    if   self.tag == iTag: lib.ElDistMultiVecRowOwner_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecRowOwner_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecRowOwner_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecRowOwner_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecRowOwner_z(*args)
    else: DataExcept()
    return owner.value

  lib.ElDistMultiVecGlobalRow_i.argtypes = \
  lib.ElDistMultiVecGlobalRow_s.argtypes = \
  lib.ElDistMultiVecGlobalRow_d.argtypes = \
  lib.ElDistMultiVecGlobalRow_c.argtypes = \
  lib.ElDistMultiVecGlobalRow_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def GlobalRow(self,iLoc):
    i = iType()
    args = [self.obj,iLoc,pointer(i)]
    if   self.tag == iTag: lib.ElDistMultiVecGlobalRow_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecGlobalRow_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecGlobalRow_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecGlobalRow_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecGlobalRow_z(*args)
    else: DataExcept()
    return i.value

  # Entrywise manipulation
  # ======================
  lib.ElDistMultiVecGet_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
  lib.ElDistMultiVecGet_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElDistMultiVecGet_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
  lib.ElDistMultiVecGet_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
  lib.ElDistMultiVecGet_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
  def Get(self,i,j):
    value = TagToType(self.tag)()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElDistMultiVecGet_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecGet_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecGet_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecGet_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecGet_z(*args)
    else: DataExcept()
    return ScalarData(value)

  lib.ElDistMultiVecSet_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistMultiVecSet_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMultiVecSet_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistMultiVecSet_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistMultiVecSet_z.argtypes = [c_void_p,iType,iType,zType]
  def Set(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElDistMultiVecSet_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecSet_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecSet_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecSet_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecSet_z(*args)
    else: DataExcept()

  lib.ElDistMultiVecReserve_i.argtypes = \
  lib.ElDistMultiVecReserve_s.argtypes = \
  lib.ElDistMultiVecReserve_d.argtypes = \
  lib.ElDistMultiVecReserve_c.argtypes = \
  lib.ElDistMultiVecReserve_z.argtypes = \
    [c_void_p,iType]
  def Reserve(self,numRemoteEntries):
    args = [self.obj,numRemoteEntries]
    if   self.tag == iTag: lib.ElDistMultiVecReserve_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecReserve_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecReserve_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecReserve_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecReserve_z(*args)
    else: DataExcept()

  lib.ElDistMultiVecQueueUpdate_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistMultiVecQueueUpdate_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMultiVecQueueUpdate_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistMultiVecQueueUpdate_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistMultiVecQueueUpdate_z.argtypes = [c_void_p,iType,iType,zType]
  def QueueUpdate(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElDistMultiVecQueueUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecQueueUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecQueueUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecQueueUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecQueueUpdate_z(*args)
    else: DataExcept()

  lib.ElDistMultiVecProcessQueues_i.argtypes = \
  lib.ElDistMultiVecProcessQueues_s.argtypes = \
  lib.ElDistMultiVecProcessQueues_d.argtypes = \
  lib.ElDistMultiVecProcessQueues_c.argtypes = \
  lib.ElDistMultiVecProcessQueues_z.argtypes = \
    [c_void_p]
  def ProcessQueues(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMultiVecProcessQueues_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecProcessQueues_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecProcessQueues_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecProcessQueues_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecProcessQueues_z(*args)
    else: DataExcept()

  lib.ElDistMultiVecGetLocal_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
  lib.ElDistMultiVecGetLocal_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElDistMultiVecGetLocal_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
  lib.ElDistMultiVecGetLocal_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
  lib.ElDistMultiVecGetLocal_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
  def GetLocal(self,iLocal,j):
    value = TagToType(self.tag)()
    args = [self.obj,iLocal,j,pointer(value)]
    if   self.tag == iTag: lib.ElDistMultiVecGetLocal_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecGetLocal_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecGetLocal_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecGetLocal_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecGetLocal_z(*args)
    else: DataExcept()
    return ScalarData(value)

  lib.ElDistMultiVecSetLocal_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistMultiVecSetLocal_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMultiVecSetLocal_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistMultiVecSetLocal_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistMultiVecSetLocal_z.argtypes = [c_void_p,iType,iType,zType]
  def SetLocal(self,iLocal,j,value):
    args = [self.obj,iLocal,j,value]
    if   self.tag == iTag: lib.ElDistMultiVecSetLocal_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecSetLocal_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecSetLocal_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecSetLocal_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecSetLocal_z(*args)
    else: DataExcept()

  lib.ElDistMultiVecUpdateLocal_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistMultiVecUpdateLocal_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMultiVecUpdateLocal_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistMultiVecUpdateLocal_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistMultiVecUpdateLocal_z.argtypes = [c_void_p,iType,iType,zType]
  def UpdateLocal(self,iLocal,j,value):
    args = [self.obj,iLocal,j,value]
    if   self.tag == iTag: lib.ElDistMultiVecUpdateLocal_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecUpdateLocal_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecUpdateLocal_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecUpdateLocal_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecUpdateLocal_z(*args)
    else: DataExcept()

  lib.ElGetContigSubmatrixDistMultiVec_i.argtypes = \
  lib.ElGetContigSubmatrixDistMultiVec_s.argtypes = \
  lib.ElGetContigSubmatrixDistMultiVec_d.argtypes = \
  lib.ElGetContigSubmatrixDistMultiVec_c.argtypes = \
  lib.ElGetContigSubmatrixDistMultiVec_z.argtypes = \
    [c_void_p,IndexRange,IndexRange,c_void_p]
  def __getitem__(self,indTup):
    iInd, jInd = indTup
    if isinstance(iInd,slice):
      if iInd.start == None:
        iInd = slice(0,iInd.stop,iInd.step)
      if iInd.stop == None:
        iInd = slice(iInd.start,self.Height(),iInd.step)
    if isinstance(jInd,slice):
      if jInd.start == None:
        jInd = slice(0,jInd.stop,jInd.step)
      if jInd.stop == None:
        jInd = slice(jInd.start,self.Width(),jInd.step)
    iRan = IndexRange(iInd)
    jRan = IndexRange(jInd)
    ASub = DistMultiVec(self.tag,self.Grid())
    args = [self.obj,iRan,jRan,ASub.obj]
    if   self.tag == iTag: lib.ElGetContigSubmatrixDistMultiVec_i(*args)
    elif self.tag == sTag: lib.ElGetContigSubmatrixDistMultiVec_s(*args)
    elif self.tag == dTag: lib.ElGetContigSubmatrixDistMultiVec_d(*args)
    elif self.tag == cTag: lib.ElGetContigSubmatrixDistMultiVec_c(*args)
    elif self.tag == zTag: lib.ElGetContigSubmatrixDistMultiVec_z(*args)
    else: DataExcept()
    return ASub
