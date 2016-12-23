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

import DistGraph as DG

class DistSparseMatrix(object):
  # Constructors and destructors
  # ============================
  lib.ElDistSparseMatrixCreate_i.argtypes = \
  lib.ElDistSparseMatrixCreate_s.argtypes = \
  lib.ElDistSparseMatrixCreate_d.argtypes = \
  lib.ElDistSparseMatrixCreate_c.argtypes = \
  lib.ElDistSparseMatrixCreate_z.argtypes = \
    [POINTER(c_void_p),c_void_p]
  def __init__(self,tag=dTag,grid=Grid.Grid.Default(),create=True):
    self.obj = c_void_p()
    self.tag = tag
    CheckTag(tag)
    if create:
      args = [pointer(self.obj),grid.obj]
      if   tag == iTag: lib.ElDistSparseMatrixCreate_i(*args)
      elif tag == sTag: lib.ElDistSparseMatrixCreate_s(*args)
      elif tag == dTag: lib.ElDistSparseMatrixCreate_d(*args)
      elif tag == cTag: lib.ElDistSparseMatrixCreate_c(*args)
      elif tag == zTag: lib.ElDistSparseMatrixCreate_z(*args)
      else: DataExcept()

  lib.ElDistSparseMatrixDestroy_i.argtypes = \
  lib.ElDistSparseMatrixDestroy_s.argtypes = \
  lib.ElDistSparseMatrixDestroy_d.argtypes = \
  lib.ElDistSparseMatrixDestroy_c.argtypes = \
  lib.ElDistSparseMatrixDestroy_z.argtypes = \
    [c_void_p]
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixDestroy_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixDestroy_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixDestroy_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixDestroy_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixDestroy_z(*args)
    else: DataExcept()

  # Assignment and reconfiguration
  # ==============================
  lib.ElDistSparseMatrixEmpty_i.argtypes = \
  lib.ElDistSparseMatrixEmpty_s.argtypes = \
  lib.ElDistSparseMatrixEmpty_d.argtypes = \
  lib.ElDistSparseMatrixEmpty_c.argtypes = \
  lib.ElDistSparseMatrixEmpty_z.argtypes = \
    [c_void_p]
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixEmpty_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixResize_i.argtypes = \
  lib.ElDistSparseMatrixResize_s.argtypes = \
  lib.ElDistSparseMatrixResize_d.argtypes = \
  lib.ElDistSparseMatrixResize_c.argtypes = \
  lib.ElDistSparseMatrixResize_z.argtypes = \
    [c_void_p,iType,iType]
  def Resize(self,height,width):
    args = [self.obj,height,width]
    if   self.tag == iTag: lib.ElDistSparseMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixResize_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixSetGrid_i.argtypes = \
  lib.ElDistSparseMatrixSetGrid_s.argtypes = \
  lib.ElDistSparseMatrixSetGrid_d.argtypes = \
  lib.ElDistSparseMatrixSetGrid_c.argtypes = \
  lib.ElDistSparseMatrixSetGrid_z.argtypes = \
    [c_void_p,c_void_p]
  def SetGrid(self,grid):
    args = [self.obj,grid.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixSetGrid_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixSetGrid_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixSetGrid_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixSetGrid_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixSetGrid_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixReserve_i.argtypes = \
  lib.ElDistSparseMatrixReserve_s.argtypes = \
  lib.ElDistSparseMatrixReserve_d.argtypes = \
  lib.ElDistSparseMatrixReserve_c.argtypes = \
  lib.ElDistSparseMatrixReserve_z.argtypes = \
    [c_void_p,iType,iType]
  def Reserve(self,numLocalEntries,numRemoteEntries=0):
    args = [self.obj,numLocalEntries,numRemoteEntries]
    if   self.tag == iTag: lib.ElDistSparseMatrixReserve_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixReserve_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixReserve_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixReserve_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixReserve_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistSparseMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistSparseMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistSparseMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistSparseMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
  def Update(self,row,col,value):
    args = [self.obj,row,col,value]
    if   self.tag == iTag: lib.ElDistSparseMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixUpdate_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixUpdateLocal_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistSparseMatrixUpdateLocal_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistSparseMatrixUpdateLocal_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistSparseMatrixUpdateLocal_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistSparseMatrixUpdateLocal_z.argtypes = [c_void_p,iType,iType,zType]
  def UpdateLocal(self,localRow,col,value):
    args = [self.obj,localRow,col,value]
    if   self.tag == iTag: lib.ElDistSparseMatrixUpdateLocal_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixUpdateLocal_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixUpdateLocal_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixUpdateLocal_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixUpdateLocal_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixZero_i.argtypes = \
  lib.ElDistSparseMatrixZero_s.argtypes = \
  lib.ElDistSparseMatrixZero_d.argtypes = \
  lib.ElDistSparseMatrixZero_c.argtypes = \
  lib.ElDistSparseMatrixZero_z.argtypes = \
    [c_void_p,iType,iType]
  def Zero(self,row,col):
    args = [self.obj,row,col]
    if   self.tag == iTag: lib.ElDistSparseMatrixZero_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixZero_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixZero_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixZero_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixZero_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixZeroLocal_i.argtypes = \
  lib.ElDistSparseMatrixZeroLocal_s.argtypes = \
  lib.ElDistSparseMatrixZeroLocal_d.argtypes = \
  lib.ElDistSparseMatrixZeroLocal_c.argtypes = \
  lib.ElDistSparseMatrixZeroLocal_z.argtypes = \
    [c_void_p,iType,iType]
  def ZeroLocal(self,localRow,col):
    args = [self.obj,localRow,col]
    if   self.tag == iTag: lib.ElDistSparseMatrixZeroLocal_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixZeroLocal_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixZeroLocal_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixZeroLocal_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixZeroLocal_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixQueueUpdate_i.argtypes = \
    [c_void_p,iType,iType,iType,bType]
  lib.ElDistSparseMatrixQueueUpdate_s.argtypes = \
    [c_void_p,iType,iType,sType,bType]
  lib.ElDistSparseMatrixQueueUpdate_d.argtypes = \
    [c_void_p,iType,iType,dType,bType]
  lib.ElDistSparseMatrixQueueUpdate_c.argtypes = \
    [c_void_p,iType,iType,cType,bType]
  lib.ElDistSparseMatrixQueueUpdate_z.argtypes = \
    [c_void_p,iType,iType,zType,bType]
  def QueueUpdate(self,row,col,value,passive=False):
    args = [self.obj,row,col,value,passive]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueUpdate_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixQueueLocalUpdate_i.argtypes = \
    [c_void_p,iType,iType,iType]
  lib.ElDistSparseMatrixQueueLocalUpdate_s.argtypes = \
    [c_void_p,iType,iType,sType]
  lib.ElDistSparseMatrixQueueLocalUpdate_d.argtypes = \
    [c_void_p,iType,iType,dType]
  lib.ElDistSparseMatrixQueueLocalUpdate_c.argtypes = \
    [c_void_p,iType,iType,cType]
  lib.ElDistSparseMatrixQueueLocalUpdate_z.argtypes = \
    [c_void_p,iType,iType,zType]
  def QueueLocalUpdate(self,localRow,col,value):
    args = [self.obj,localRow,col,value]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueLocalUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueLocalUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueLocalUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueLocalUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueLocalUpdate_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixQueueZero_i.argtypes = \
  lib.ElDistSparseMatrixQueueZero_s.argtypes = \
  lib.ElDistSparseMatrixQueueZero_d.argtypes = \
  lib.ElDistSparseMatrixQueueZero_c.argtypes = \
  lib.ElDistSparseMatrixQueueZero_z.argtypes = \
    [c_void_p,iType,iType,bType]
  def QueueZero(self,row,col,passive=False):
    args = [self.obj,row,col,passive]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueZero_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueZero_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueZero_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueZero_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueZero_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixQueueLocalZero_i.argtypes = \
  lib.ElDistSparseMatrixQueueLocalZero_s.argtypes = \
  lib.ElDistSparseMatrixQueueLocalZero_d.argtypes = \
  lib.ElDistSparseMatrixQueueLocalZero_c.argtypes = \
  lib.ElDistSparseMatrixQueueLocalZero_z.argtypes = \
    [c_void_p,iType,iType]
  def QueueLocalZero(self,localRow,col):
    args = [self.obj,localRow,col]
    if   self.tag == iTag: lib.ElDistSparseMatrixQueueLocalZero_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixQueueLocalZero_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixQueueLocalZero_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixQueueLocalZero_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixQueueLocalZero_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixProcessQueues_i.argtypes = \
  lib.ElDistSparseMatrixProcessQueues_s.argtypes = \
  lib.ElDistSparseMatrixProcessQueues_d.argtypes = \
  lib.ElDistSparseMatrixProcessQueues_c.argtypes = \
  lib.ElDistSparseMatrixProcessQueues_z.argtypes = \
    [c_void_p]
  def ProcessQueues(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixProcessQueues_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixProcessQueues_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixProcessQueues_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixProcessQueues_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixProcessQueues_z(*args)
    else: DataExcept()

  lib.ElDistSparseMatrixProcessLocalQueues_i.argtypes = \
  lib.ElDistSparseMatrixProcessLocalQueues_s.argtypes = \
  lib.ElDistSparseMatrixProcessLocalQueues_d.argtypes = \
  lib.ElDistSparseMatrixProcessLocalQueues_c.argtypes = \
  lib.ElDistSparseMatrixProcessLocalQueues_z.argtypes = \
    [c_void_p]
  def ProcessLocalQueues(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistSparseMatrixProcessLocalQueues_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixProcessLocalQueues_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixProcessLocalQueues_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixProcessLocalQueues_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixProcessLocalQueues_z(*args)
    else: DataExcept()

  # Queries
  # =======
  lib.ElDistSparseMatrixHeight_i.argtypes = \
  lib.ElDistSparseMatrixHeight_s.argtypes = \
  lib.ElDistSparseMatrixHeight_d.argtypes = \
  lib.ElDistSparseMatrixHeight_c.argtypes = \
  lib.ElDistSparseMatrixHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Height(self):
    height = iType()
    args = [self.obj,pointer(height)]
    if   self.tag == iTag: lib.ElDistSparseMatrixHeight_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixHeight_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixHeight_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixHeight_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixHeight_z(*args)
    else: DataExcept()
    return height.value

  lib.ElDistSparseMatrixWidth_i.argtypes = \
  lib.ElDistSparseMatrixWidth_s.argtypes = \
  lib.ElDistSparseMatrixWidth_d.argtypes = \
  lib.ElDistSparseMatrixWidth_c.argtypes = \
  lib.ElDistSparseMatrixWidth_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Width(self):
    width = iType()
    args = [self.obj,pointer(width)]
    if   self.tag == iTag: lib.ElDistSparseMatrixWidth_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixWidth_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixWidth_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixWidth_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixWidth_z(*args)
    else: DataExcept()
    return width.value

  lib.ElDistSparseMatrixDistGraph_i.argtypes = \
  lib.ElDistSparseMatrixDistGraph_s.argtypes = \
  lib.ElDistSparseMatrixDistGraph_d.argtypes = \
  lib.ElDistSparseMatrixDistGraph_c.argtypes = \
  lib.ElDistSparseMatrixDistGraph_z.argtypes = \
  lib.ElDistSparseMatrixLockedDistGraph_i.argtypes = \
  lib.ElDistSparseMatrixLockedDistGraph_s.argtypes = \
  lib.ElDistSparseMatrixLockedDistGraph_d.argtypes = \
  lib.ElDistSparseMatrixLockedDistGraph_c.argtypes = \
  lib.ElDistSparseMatrixLockedDistGraph_z.argtypes = \
    [c_void_p,POINTER(c_void_p)]
  def DistGraph(self,locked=False):
    graph = DG.DistGraph(mpi.COMM_WORLD(),False)
    args = [self.obj,pointer(graph.obj)]
    if locked:
      if   self.tag == iTag: lib.ElDistSparseMatrixLockedDistGraph_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixLockedDistGraph_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixLockedDistGraph_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixLockedDistGraph_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixLockedDistGraph_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistSparseMatrixDistGraph_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixDistGraph_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixDistGraph_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixDistGraph_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixDistGraph_z(*args)
      else: DataExcept()
    return graph

  lib.ElDistSparseMatrixFirstLocalRow_i.argtypes = \
  lib.ElDistSparseMatrixFirstLocalRow_s.argtypes = \
  lib.ElDistSparseMatrixFirstLocalRow_d.argtypes = \
  lib.ElDistSparseMatrixFirstLocalRow_c.argtypes = \
  lib.ElDistSparseMatrixFirstLocalRow_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def FirstLocalRow(self):
    firstLocalRow = iType()
    args = [self.obj,pointer(firstLocalRow)]
    if   self.tag == iTag: lib.ElDistSparseMatrixFirstLocalRow_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixFirstLocalRow_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixFirstLocalRow_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixFirstLocalRow_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixFirstLocalRow_z(*args)
    else: DataExcept()
    return firstLocalRow.value

  lib.ElDistSparseMatrixLocalHeight_i.argtypes = \
  lib.ElDistSparseMatrixLocalHeight_s.argtypes = \
  lib.ElDistSparseMatrixLocalHeight_d.argtypes = \
  lib.ElDistSparseMatrixLocalHeight_c.argtypes = \
  lib.ElDistSparseMatrixLocalHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def LocalHeight(self):
    localHeight = iType()
    args = [self.obj,pointer(localHeight)]
    if   self.tag == iTag: lib.ElDistSparseMatrixLocalHeight_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixLocalHeight_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixLocalHeight_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixLocalHeight_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixLocalHeight_z(*args)
    else: DataExcept()
    return localHeight.value

  lib.ElDistSparseMatrixNumLocalEntries_i.argtypes = \
  lib.ElDistSparseMatrixNumLocalEntries_s.argtypes = \
  lib.ElDistSparseMatrixNumLocalEntries_d.argtypes = \
  lib.ElDistSparseMatrixNumLocalEntries_c.argtypes = \
  lib.ElDistSparseMatrixNumLocalEntries_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def NumLocalEntries(self):
    numLocalEntries = iType()
    args = [self.obj,pointer(numLocalEntries)]
    if   self.tag == iTag: lib.ElDistSparseMatrixNumLocalEntries_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixNumLocalEntries_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixNumLocalEntries_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixNumLocalEntries_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixNumLocalEntries_z(*args)
    else: DataExcept()
    return numLocalEntries.value

  lib.ElDistSparseMatrixCapacity_i.argtypes = \
  lib.ElDistSparseMatrixCapacity_s.argtypes = \
  lib.ElDistSparseMatrixCapacity_d.argtypes = \
  lib.ElDistSparseMatrixCapacity_c.argtypes = \
  lib.ElDistSparseMatrixCapacity_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Capacity(self):
    capacity = iType()
    args = [self.obj,pointer(capacity)]
    if   self.tag == iTag: lib.ElDistSparseMatrixCapacity_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixCapacity_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixCapacity_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixCapacity_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixCapacity_z(*args)
    else: DataExcept()
    return capacity.value

  lib.ElDistSparseMatrixLocallyConsistent_i.argtypes = \
  lib.ElDistSparseMatrixLocallyConsistent_s.argtypes = \
  lib.ElDistSparseMatrixLocallyConsistent_d.argtypes = \
  lib.ElDistSparseMatrixLocallyConsistent_c.argtypes = \
  lib.ElDistSparseMatrixLocallyConsistent_z.argtypes = \
    [c_void_p,POINTER(bType)]
  def LocallyConsistent(self):
    consistent = bType()
    args = [self.obj,pointer(consistent)]
    if   self.tag == iTag: lib.ElDistSparseMatrixLocallyConsistent_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixLocallyConsistent_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixLocallyConsistent_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixLocallyConsistent_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixLocallyConsistent_z(*args)
    else: DataExcept()
    return consistent.value

  lib.ElDistSparseMatrixGrid_i.argtypes = \
  lib.ElDistSparseMatrixGrid_s.argtypes = \
  lib.ElDistSparseMatrixGrid_d.argtypes = \
  lib.ElDistSparseMatrixGrid_c.argtypes = \
  lib.ElDistSparseMatrixGrid_z.argtypes = \
    [c_void_p,POINTER(c_void_p)]
  def Grid(self):
    grid = Grid.Grid(create=False)
    args = [self.obj,pointer(grid.obj)]
    if   self.tag == iTag: lib.ElDistSparseMatrixGrid_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixGrid_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixGrid_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixGrid_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixGrid_z(*args)
    else: DataExcept()
    return grid

  lib.ElDistSparseMatrixBlocksize_i.argtypes = \
  lib.ElDistSparseMatrixBlocksize_s.argtypes = \
  lib.ElDistSparseMatrixBlocksize_d.argtypes = \
  lib.ElDistSparseMatrixBlocksize_c.argtypes = \
  lib.ElDistSparseMatrixBlocksize_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Blocksize(self):
    blocksize = iType()
    args = [self.obj,pointer(blocksize)]
    if   self.tag == iTag: lib.ElDistSparseMatrixBlocksize_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixBlocksize_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixBlocksize_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixBlocksize_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixBlocksize_z(*args)
    else: DataExcept()
    return blocksize.value

  lib.ElDistSparseMatrixRowOwner_i.argtypes = \
  lib.ElDistSparseMatrixRowOwner_s.argtypes = \
  lib.ElDistSparseMatrixRowOwner_d.argtypes = \
  lib.ElDistSparseMatrixRowOwner_c.argtypes = \
  lib.ElDistSparseMatrixRowOwner_z.argtypes = \
    [c_void_p,iType,POINTER(c_int)]
  def RowOwner(self,i):
    owner = c_int()
    args = [self.obj,i,pointer(owner)]
    if   self.tag == iTag: lib.ElDistSparseMatrixRowOwner_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixRowOwner_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixRowOwner_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixRowOwner_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixRowOwner_z(*args)
    else: DataExcept()
    return owner.value

  lib.ElDistSparseMatrixGlobalRow_i.argtypes = \
  lib.ElDistSparseMatrixGlobalRow_s.argtypes = \
  lib.ElDistSparseMatrixGlobalRow_d.argtypes = \
  lib.ElDistSparseMatrixGlobalRow_c.argtypes = \
  lib.ElDistSparseMatrixGlobalRow_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def GlobalRow(self,iLoc):
    i = iType()
    args = [self.obj,iLoc,pointer(i)]
    if   self.tag == iTag: lib.ElDistSparseMatrixGlobalRow_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixGlobalRow_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixGlobalRow_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixGlobalRow_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixGlobalRow_z(*args)
    else: DataExcept()
    return i.value

  lib.ElDistSparseMatrixRow_i.argtypes = \
  lib.ElDistSparseMatrixRow_s.argtypes = \
  lib.ElDistSparseMatrixRow_d.argtypes = \
  lib.ElDistSparseMatrixRow_c.argtypes = \
  lib.ElDistSparseMatrixRow_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def Row(self,localInd):
    row = iType()
    args = [self.obj,localInd,pointer(row)]
    if   self.tag == iTag: lib.ElDistSparseMatrixRow_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixRow_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixRow_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixRow_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixRow_z(*args)
    else: DataExcept()
    return row.value

  lib.ElDistSparseMatrixCol_i.argtypes = \
  lib.ElDistSparseMatrixCol_s.argtypes = \
  lib.ElDistSparseMatrixCol_d.argtypes = \
  lib.ElDistSparseMatrixCol_c.argtypes = \
  lib.ElDistSparseMatrixCol_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def Col(self,localInd):
    col = iType()
    args = [self.obj,localInd,pointer(col)]
    if   self.tag == iTag: lib.ElDistSparseMatrixCol_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixCol_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixCol_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixCol_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixCol_z(*args)
    else: DataExcept()
    return col.value

  lib.ElDistSparseMatrixValue_i.argtypes = [c_void_p,iType,POINTER(iType)]
  lib.ElDistSparseMatrixValue_s.argtypes = [c_void_p,iType,POINTER(sType)]
  lib.ElDistSparseMatrixValue_d.argtypes = [c_void_p,iType,POINTER(dType)]
  lib.ElDistSparseMatrixValue_c.argtypes = [c_void_p,iType,POINTER(cType)]
  lib.ElDistSparseMatrixValue_z.argtypes = [c_void_p,iType,POINTER(zType)]
  def Value(self,localInd):
    value =  TagToType(self.tag)()
    args = [self.obj,localInd,pointer(value)]
    if   self.tag == iTag: lib.ElDistSparseMatrixValue_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixValue_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixValue_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixValue_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixValue_z(*args)
    else: DataExcept()
    return ScalarData(value)

  lib.ElDistSparseMatrixRowOffset_i.argtypes = \
  lib.ElDistSparseMatrixRowOffset_s.argtypes = \
  lib.ElDistSparseMatrixRowOffset_d.argtypes = \
  lib.ElDistSparseMatrixRowOffset_c.argtypes = \
  lib.ElDistSparseMatrixRowOffset_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def RowOffset(self,localRow):
    offset = iType()
    args = [self.obj,localRow,pointer(offset)]
    if   self.tag == iTag: lib.ElDistSparseMatrixRowOffset_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixRowOffset_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixRowOffset_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixRowOffset_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixRowOffset_z(*args)
    else: DataExcept()
    return offset.value

  lib.ElDistSparseMatrixOffset_i.argtypes = \
  lib.ElDistSparseMatrixOffset_s.argtypes = \
  lib.ElDistSparseMatrixOffset_d.argtypes = \
  lib.ElDistSparseMatrixOffset_c.argtypes = \
  lib.ElDistSparseMatrixOffset_z.argtypes = \
    [c_void_p,iType,iType,POINTER(iType)]
  def Offset(self,localRow,col):
    offset = iType()
    args = [self.obj,localRow,col,pointer(offset)]
    if   self.tag == iTag: lib.ElDistSparseMatrixOffset_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixOffset_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixOffset_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixOffset_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixOffset_z(*args)
    else: DataExcept()
    return offset.value

  lib.ElDistSparseMatrixNumConnections_i.argtypes = \
  lib.ElDistSparseMatrixNumConnections_s.argtypes = \
  lib.ElDistSparseMatrixNumConnections_d.argtypes = \
  lib.ElDistSparseMatrixNumConnections_c.argtypes = \
  lib.ElDistSparseMatrixNumConnections_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def NumConnections(self,localRow):
    numConnections = iType()
    args = [self.obj,localRow,pointer(numConnections)]
    if   self.tag == iTag: lib.ElDistSparseMatrixNumConnections_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixNumConnections_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixNumConnections_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixNumConnections_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixNumConnections_z(*args)
    else: DataExcept()
    return numConnections.value

  lib.ElDistSparseMatrixImbalance_i.argtypes = \
  lib.ElDistSparseMatrixImbalance_s.argtypes = \
  lib.ElDistSparseMatrixImbalance_d.argtypes = \
  lib.ElDistSparseMatrixImbalance_c.argtypes = \
  lib.ElDistSparseMatrixImbalance_z.argtypes = \
    [c_void_p,POINTER(dType)]
  def Imbalance(self):
    imbalance = dType()
    args = [self.obj,pointer(imbalance)]
    if   self.tag == iTag: lib.ElDistSparseMatrixImbalance_i(*args)
    elif self.tag == sTag: lib.ElDistSparseMatrixImbalance_s(*args)
    elif self.tag == dTag: lib.ElDistSparseMatrixImbalance_d(*args)
    elif self.tag == cTag: lib.ElDistSparseMatrixImbalance_c(*args)
    elif self.tag == zTag: lib.ElDistSparseMatrixImbalance_z(*args)
    else: DataExcept()
    return imbalance.value

  lib.ElDistSparseMatrixSourceBuffer_i.argtypes = \
  lib.ElDistSparseMatrixSourceBuffer_s.argtypes = \
  lib.ElDistSparseMatrixSourceBuffer_d.argtypes = \
  lib.ElDistSparseMatrixSourceBuffer_c.argtypes = \
  lib.ElDistSparseMatrixSourceBuffer_z.argtypes = \
  lib.ElDistSparseMatrixLockedSourceBuffer_i.argtypes = \
  lib.ElDistSparseMatrixLockedSourceBuffer_s.argtypes = \
  lib.ElDistSparseMatrixLockedSourceBuffer_d.argtypes = \
  lib.ElDistSparseMatrixLockedSourceBuffer_c.argtypes = \
  lib.ElDistSparseMatrixLockedSourceBuffer_z.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def SourceBuffer(self,locked=False):
    sourceBuf = POINTER(iType)()
    args = [self.obj,pointer(sourceBuf)]
    if locked:
      if   self.tag == iTag: lib.ElDistSparseMatrixLockedSourceBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixLockedSourceBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixLockedSourceBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixLockedSourceBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixLockedSourceBuffer_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistSparseMatrixSourceBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixSourceBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixSourceBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixSourceBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixSourceBuffer_z(*args)
      else: DataExcept()
    return sourceBuf

  lib.ElDistSparseMatrixTargetBuffer_i.argtypes = \
  lib.ElDistSparseMatrixTargetBuffer_s.argtypes = \
  lib.ElDistSparseMatrixTargetBuffer_d.argtypes = \
  lib.ElDistSparseMatrixTargetBuffer_c.argtypes = \
  lib.ElDistSparseMatrixTargetBuffer_z.argtypes = \
  lib.ElDistSparseMatrixLockedTargetBuffer_i.argtypes = \
  lib.ElDistSparseMatrixLockedTargetBuffer_s.argtypes = \
  lib.ElDistSparseMatrixLockedTargetBuffer_d.argtypes = \
  lib.ElDistSparseMatrixLockedTargetBuffer_c.argtypes = \
  lib.ElDistSparseMatrixLockedTargetBuffer_z.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def TargetBuffer(self,locked=False):
    targetBuf = POINTER(iType)()
    args = [self.obj,pointer(targetBuf)]
    if locked:
      if   self.tag == iTag: lib.ElDistSparseMatrixLockedTargetBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixLockedTargetBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixLockedTargetBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixLockedTargetBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixLockedTargetBuffer_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistSparseMatrixTargetBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixTargetBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixTargetBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixTargetBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixTargetBuffer_z(*args)
      else: DataExcept()
    return targetBuf

  lib.ElDistSparseMatrixValueBuffer_i.argtypes = \
  lib.ElDistSparseMatrixLockedValueBuffer_i.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  lib.ElDistSparseMatrixValueBuffer_s.argtypes = \
  lib.ElDistSparseMatrixLockedValueBuffer_s.argtypes = \
    [c_void_p,POINTER(POINTER(sType))]
  lib.ElDistSparseMatrixValueBuffer_d.argtypes = \
  lib.ElDistSparseMatrixLockedValueBuffer_d.argtypes = \
    [c_void_p,POINTER(POINTER(dType))]
  lib.ElDistSparseMatrixValueBuffer_c.argtypes = \
  lib.ElDistSparseMatrixLockedValueBuffer_c.argtypes = \
    [c_void_p,POINTER(POINTER(cType))]
  lib.ElDistSparseMatrixValueBuffer_z.argtypes = \
  lib.ElDistSparseMatrixLockedValueBuffer_z.argtypes = \
    [c_void_p,POINTER(POINTER(zType))]
  def ValueBuffer(self,locked=False):
    valueBuf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(valueBuf)]
    if locked:
      if   self.tag == iTag: lib.ElDistSparseMatrixLockedValueBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixLockedValueBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixLockedValueBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixLockedValueBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixLockedValueBuffer_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistSparseMatrixValueBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistSparseMatrixValueBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistSparseMatrixValueBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistSparseMatrixValueBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistSparseMatrixValueBuffer_z(*args)
      else: DataExcept()
    return valueBuf

  lib.ElGetContigSubmatrixDistSparse_i.argtypes = \
  lib.ElGetContigSubmatrixDistSparse_s.argtypes = \
  lib.ElGetContigSubmatrixDistSparse_d.argtypes = \
  lib.ElGetContigSubmatrixDistSparse_c.argtypes = \
  lib.ElGetContigSubmatrixDistSparse_z.argtypes = \
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
    ASub = DistSparseMatrix(self.tag,self.Grid())
    args = [self.obj,iRan,jRan,ASub.obj]
    if   self.tag == iTag: lib.ElGetContigSubmatrixDistSparse_i(*args)
    elif self.tag == sTag: lib.ElGetContigSubmatrixDistSparse_s(*args)
    elif self.tag == dTag: lib.ElGetContigSubmatrixDistSparse_d(*args)
    elif self.tag == cTag: lib.ElGetContigSubmatrixDistSparse_c(*args)
    elif self.tag == zTag: lib.ElGetContigSubmatrixDistSparse_z(*args)
    else: DataExcept()
    return ASub
