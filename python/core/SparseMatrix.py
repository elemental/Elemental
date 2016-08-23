#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import Graph as G

class SparseMatrix(object):
  # Constructors and destructors
  # ============================
  lib.ElSparseMatrixCreate_i.argtypes = \
  lib.ElSparseMatrixCreate_s.argtypes = \
  lib.ElSparseMatrixCreate_d.argtypes = \
  lib.ElSparseMatrixCreate_c.argtypes = \
  lib.ElSparseMatrixCreate_z.argtypes = \
    [POINTER(c_void_p)]
  def __init__(self,tag=dTag,create=True):
    self.obj = c_void_p()
    self.tag = tag
    CheckTag(tag)
    if create:
      args = [pointer(self.obj)]
      if   tag == iTag: lib.ElSparseMatrixCreate_i(*args)
      elif tag == sTag: lib.ElSparseMatrixCreate_s(*args)
      elif tag == dTag: lib.ElSparseMatrixCreate_d(*args)
      elif tag == cTag: lib.ElSparseMatrixCreate_c(*args)
      elif tag == zTag: lib.ElSparseMatrixCreate_z(*args)
      else: DataExcept()

  lib.ElSparseMatrixDestroy_i.argtypes = \
  lib.ElSparseMatrixDestroy_s.argtypes = \
  lib.ElSparseMatrixDestroy_d.argtypes = \
  lib.ElSparseMatrixDestroy_c.argtypes = \
  lib.ElSparseMatrixDestroy_z.argtypes = \
    [c_void_p]
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElSparseMatrixDestroy_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixDestroy_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixDestroy_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixDestroy_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixDestroy_z(*args)
    else: DataExcept()

  # Assignment and reconfiguration
  # ==============================
  lib.ElSparseMatrixEmpty_i.argtypes = \
  lib.ElSparseMatrixEmpty_s.argtypes = \
  lib.ElSparseMatrixEmpty_d.argtypes = \
  lib.ElSparseMatrixEmpty_c.argtypes = \
  lib.ElSparseMatrixEmpty_z.argtypes = \
    [c_void_p]
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElSparseMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixEmpty_z(*args)
    else: DataExcept()

  lib.ElSparseMatrixResize_i.argtypes = \
  lib.ElSparseMatrixResize_s.argtypes = \
  lib.ElSparseMatrixResize_d.argtypes = \
  lib.ElSparseMatrixResize_c.argtypes = \
  lib.ElSparseMatrixResize_z.argtypes = \
    [c_void_p,iType,iType]
  def Resize(self,height,width):
    args = [self.obj,height,width]
    if   self.tag == iTag: lib.ElSparseMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixResize_z(*args)
    else: DataExcept()

  lib.ElSparseMatrixReserve_i.argtypes = \
  lib.ElSparseMatrixReserve_s.argtypes = \
  lib.ElSparseMatrixReserve_d.argtypes = \
  lib.ElSparseMatrixReserve_c.argtypes = \
  lib.ElSparseMatrixReserve_z.argtypes = \
    [c_void_p,iType]
  def Reserve(self,numEntries):
    args = [self.obj,numEntries]
    if   self.tag == iTag: lib.ElSparseMatrixReserve_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixReserve_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixReserve_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixReserve_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixReserve_z(*args)
    else: DataExcept()

  lib.ElSparseMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElSparseMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElSparseMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElSparseMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElSparseMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
  def Update(self,row,col,value):
    args = [self.obj,row,col,value]
    if   self.tag == iTag: lib.ElSparseMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixUpdate_z(*args)
    else: DataExcept()

  lib.ElSparseMatrixZero_i.argtypes = \
  lib.ElSparseMatrixZero_s.argtypes = \
  lib.ElSparseMatrixZero_d.argtypes = \
  lib.ElSparseMatrixZero_c.argtypes = \
  lib.ElSparseMatrixZero_z.argtypes = \
    [c_void_p,iType,iType]
  def Zero(self,row,col):
    args = [self.obj,row,col]
    if   self.tag == iTag: lib.ElSparseMatrixZero_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixZero_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixZero_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixZero_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixZero_z(*args)
    else: DataExcept()

  lib.ElSparseMatrixQueueUpdate_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElSparseMatrixQueueUpdate_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElSparseMatrixQueueUpdate_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElSparseMatrixQueueUpdate_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElSparseMatrixQueueUpdate_z.argtypes = [c_void_p,iType,iType,zType]
  def QueueUpdate(self,row,col,value):
    args = [self.obj,row,col,value]
    if   self.tag == iTag: lib.ElSparseMatrixQueueUpdate_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixQueueUpdate_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixQueueUpdate_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixQueueUpdate_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixQueueUpdate_z(*args)
    else: DataExcept()

  lib.ElSparseMatrixQueueZero_i.argtypes = \
  lib.ElSparseMatrixQueueZero_s.argtypes = \
  lib.ElSparseMatrixQueueZero_d.argtypes = \
  lib.ElSparseMatrixQueueZero_c.argtypes = \
  lib.ElSparseMatrixQueueZero_z.argtypes = \
    [c_void_p,iType,iType]
  def QueueZero(self,row,col):
    args = [self.obj,row,col]
    if   self.tag == iTag: lib.ElSparseMatrixQueueZero_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixQueueZero_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixQueueZero_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixQueueZero_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixQueueZero_z(*args)
    else: DataExcept()

  lib.ElSparseMatrixProcessQueues_i.argtypes = \
  lib.ElSparseMatrixProcessQueues_s.argtypes = \
  lib.ElSparseMatrixProcessQueues_d.argtypes = \
  lib.ElSparseMatrixProcessQueues_c.argtypes = \
  lib.ElSparseMatrixProcessQueues_z.argtypes = \
    [c_void_p]
  def ProcessQueues(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElSparseMatrixProcessQueues_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixProcessQueues_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixProcessQueues_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixProcessQueues_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixProcessQueues_z(*args)
    else: DataExcept()

  # Queries
  # =======
  lib.ElSparseMatrixHeight_i.argtypes = \
  lib.ElSparseMatrixHeight_s.argtypes = \
  lib.ElSparseMatrixHeight_d.argtypes = \
  lib.ElSparseMatrixHeight_c.argtypes = \
  lib.ElSparseMatrixHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Height(self):
    height = iType()
    args = [self.obj,pointer(height)]
    if   self.tag == iTag: lib.ElSparseMatrixHeight_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixHeight_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixHeight_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixHeight_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixHeight_z(*args)
    else: DataExcept()
    return height.value

  lib.ElSparseMatrixWidth_i.argtypes = \
  lib.ElSparseMatrixWidth_s.argtypes = \
  lib.ElSparseMatrixWidth_d.argtypes = \
  lib.ElSparseMatrixWidth_c.argtypes = \
  lib.ElSparseMatrixWidth_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Width(self):
    width = iType()
    args = [self.obj,pointer(width)]
    if   self.tag == iTag: lib.ElSparseMatrixWidth_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixWidth_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixWidth_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixWidth_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixWidth_z(*args)
    else: DataExcept()
    return width.value

  lib.ElSparseMatrixNumEntries_i.argtypes = \
  lib.ElSparseMatrixNumEntries_s.argtypes = \
  lib.ElSparseMatrixNumEntries_d.argtypes = \
  lib.ElSparseMatrixNumEntries_c.argtypes = \
  lib.ElSparseMatrixNumEntries_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def NumEntries(self):
    numEntries = iType()
    args = [self.obj,pointer(numEntries)]
    if   self.tag == iTag: lib.ElSparseMatrixNumEntries_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixNumEntries_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixNumEntries_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixNumEntries_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixNumEntries_z(*args)
    else: DataExcept()
    return numEntries.value

  lib.ElSparseMatrixCapacity_i.argtypes = \
  lib.ElSparseMatrixCapacity_s.argtypes = \
  lib.ElSparseMatrixCapacity_d.argtypes = \
  lib.ElSparseMatrixCapacity_c.argtypes = \
  lib.ElSparseMatrixCapacity_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Capacity(self):
    capacity = iType()
    args = [self.obj,pointer(capacity)]
    if   self.tag == iTag: lib.ElSparseMatrixCapacity_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixCapacity_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixCapacity_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixCapacity_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixCapacity_z(*args)
    else: DataExcept()
    return capacity.value

  lib.ElSparseMatrixConsistent_i.argtypes = \
  lib.ElSparseMatrixConsistent_s.argtypes = \
  lib.ElSparseMatrixConsistent_d.argtypes = \
  lib.ElSparseMatrixConsistent_c.argtypes = \
  lib.ElSparseMatrixConsistent_z.argtypes = \
    [c_void_p,POINTER(bType)]
  def Consistent(self):
    consistent = bType()
    args = [self.obj,pointer(consistent)]
    if   self.tag == iTag: lib.ElSparseMatrixConsistent_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixConsistent_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixConsistent_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixConsistent_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixConsistent_z(*args)
    else: DataExcept()
    return consistent.value

  lib.ElSparseMatrixGraph_i.argtypes = \
  lib.ElSparseMatrixGraph_s.argtypes = \
  lib.ElSparseMatrixGraph_d.argtypes = \
  lib.ElSparseMatrixGraph_c.argtypes = \
  lib.ElSparseMatrixGraph_z.argtypes = \
  lib.ElSparseMatrixLockedGraph_i.argtypes = \
  lib.ElSparseMatrixLockedGraph_s.argtypes = \
  lib.ElSparseMatrixLockedGraph_d.argtypes = \
  lib.ElSparseMatrixLockedGraph_c.argtypes = \
  lib.ElSparseMatrixLockedGraph_z.argtypes = \
    [c_void_p,POINTER(c_void_p)]
  def Graph(self,locked=False):
    graph = G.Graph(False)
    args = [self.obj,pointer(graph.obj)]
    if locked:
      if   self.tag == iTag: lib.ElSparseMatrixLockedGraph_i(*args)  
      elif self.tag == sTag: lib.ElSparseMatrixLockedGraph_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixLockedGraph_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixLockedGraph_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixLockedGraph_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElSparseMatrixGraph_i(*args)  
      elif self.tag == sTag: lib.ElSparseMatrixGraph_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixGraph_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixGraph_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixGraph_z(*args)
      else: DataExcept()
    return graph

  lib.ElSparseMatrixRow_i.argtypes = \
  lib.ElSparseMatrixRow_s.argtypes = \
  lib.ElSparseMatrixRow_d.argtypes = \
  lib.ElSparseMatrixRow_c.argtypes = \
  lib.ElSparseMatrixRow_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def Row(self,index):
    row = iType() 
    args = [self.obj,index,pointer(row)]
    if   self.tag == iTag: lib.ElSparseMatrixRow_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixRow_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixRow_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixRow_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixRow_z(*args)
    else: DataExcept()
    return row.value

  lib.ElSparseMatrixCol_i.argtypes = \
  lib.ElSparseMatrixCol_s.argtypes = \
  lib.ElSparseMatrixCol_d.argtypes = \
  lib.ElSparseMatrixCol_c.argtypes = \
  lib.ElSparseMatrixCol_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def Col(self,index):
    col = iType() 
    args = [self.obj,index,pointer(col)]
    if   self.tag == iTag: lib.ElSparseMatrixCol_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixCol_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixCol_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixCol_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixCol_z(*args)
    else: DataExcept()
    return col.value

  lib.ElSparseMatrixValue_i.argtypes = [c_void_p,iType,POINTER(iType)]
  lib.ElSparseMatrixValue_s.argtypes = [c_void_p,iType,POINTER(sType)]
  lib.ElSparseMatrixValue_d.argtypes = [c_void_p,iType,POINTER(dType)]
  lib.ElSparseMatrixValue_c.argtypes = [c_void_p,iType,POINTER(cType)]
  lib.ElSparseMatrixValue_z.argtypes = [c_void_p,iType,POINTER(zType)]
  def Value(self,index):
    value =  TagToType(self.tag)()
    args = [self.obj,index,pointer(value)]
    if   self.tag == iTag: lib.ElSparseMatrixValue_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixValue_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixValue_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixValue_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixValue_z(*args)
    else: DataExcept()
    return ScalarData(value)

  lib.ElSparseMatrixRowOffset_i.argtypes = \
  lib.ElSparseMatrixRowOffset_s.argtypes = \
  lib.ElSparseMatrixRowOffset_d.argtypes = \
  lib.ElSparseMatrixRowOffset_c.argtypes = \
  lib.ElSparseMatrixRowOffset_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def RowOffset(self,row):
    offset = iType()
    args = [self.obj,row,pointer(offset)]
    if   self.tag == iTag: lib.ElSparseMatrixRowOffset_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixRowOffset_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixRowOffset_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixRowOffset_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixRowOffset_z(*args)
    else: DataExcept()
    return offset.value

  lib.ElSparseMatrixOffset_i.argtypes = \
  lib.ElSparseMatrixOffset_s.argtypes = \
  lib.ElSparseMatrixOffset_d.argtypes = \
  lib.ElSparseMatrixOffset_c.argtypes = \
  lib.ElSparseMatrixOffset_z.argtypes = \
    [c_void_p,iType,iType,POINTER(iType)]
  def Offset(self,row,col):
    offset = iType()
    args = [self.obj,row,col,pointer(offset)]
    if   self.tag == iTag: lib.ElSparseMatrixOffset_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixOffset_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixOffset_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixOffset_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixOffset_z(*args)
    else: DataExcept()
    return offset.value

  lib.ElSparseMatrixNumConnections_i.argtypes = \
  lib.ElSparseMatrixNumConnections_s.argtypes = \
  lib.ElSparseMatrixNumConnections_d.argtypes = \
  lib.ElSparseMatrixNumConnections_c.argtypes = \
  lib.ElSparseMatrixNumConnections_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def NumConnections(self,row):
    numConnections = iType()
    args = [self.obj,row,pointer(numConnections)]
    if   self.tag == iTag: lib.ElSparseMatrixNumConnections_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixNumConnections_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixNumConnections_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixNumConnections_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixNumConnections_z(*args)
    else: DataExcept()
    return numConnections.value

  lib.ElSparseMatrixSourceBuffer_i.argtypes = \
  lib.ElSparseMatrixSourceBuffer_s.argtypes = \
  lib.ElSparseMatrixSourceBuffer_d.argtypes = \
  lib.ElSparseMatrixSourceBuffer_c.argtypes = \
  lib.ElSparseMatrixSourceBuffer_z.argtypes = \
  lib.ElSparseMatrixLockedSourceBuffer_i.argtypes = \
  lib.ElSparseMatrixLockedSourceBuffer_s.argtypes = \
  lib.ElSparseMatrixLockedSourceBuffer_d.argtypes = \
  lib.ElSparseMatrixLockedSourceBuffer_c.argtypes = \
  lib.ElSparseMatrixLockedSourceBuffer_z.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def SourceBuffer(self,locked=False):
    sourceBuf = POINTER(iType)()
    args = [self.obj,pointer(sourceBuf)]    
    if locked:
      if   self.tag == iTag: lib.ElSparseMatrixLockedSourceBuffer_i(*args)
      elif self.tag == sTag: lib.ElSparseMatrixLockedSourceBuffer_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixLockedSourceBuffer_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixLockedSourceBuffer_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixLockedSourceBuffer_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElSparseMatrixSourceBuffer_i(*args)
      elif self.tag == sTag: lib.ElSparseMatrixSourceBuffer_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixSourceBuffer_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixSourceBuffer_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixSourceBuffer_z(*args)
      else: DataExcept()
    return sourceBuf

  lib.ElSparseMatrixTargetBuffer_i.argtypes = \
  lib.ElSparseMatrixTargetBuffer_s.argtypes = \
  lib.ElSparseMatrixTargetBuffer_d.argtypes = \
  lib.ElSparseMatrixTargetBuffer_c.argtypes = \
  lib.ElSparseMatrixTargetBuffer_z.argtypes = \
  lib.ElSparseMatrixLockedTargetBuffer_i.argtypes = \
  lib.ElSparseMatrixLockedTargetBuffer_s.argtypes = \
  lib.ElSparseMatrixLockedTargetBuffer_d.argtypes = \
  lib.ElSparseMatrixLockedTargetBuffer_c.argtypes = \
  lib.ElSparseMatrixLockedTargetBuffer_z.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def TargetBuffer(self,locked=False):
    targetBuf = POINTER(iType)()
    args = [self.obj,pointer(targetBuf)]    
    if locked:
      if   self.tag == iTag: lib.ElSparseMatrixLockedTargetBuffer_i(*args)
      elif self.tag == sTag: lib.ElSparseMatrixLockedTargetBuffer_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixLockedTargetBuffer_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixLockedTargetBuffer_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixLockedTargetBuffer_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElSparseMatrixTargetBuffer_i(*args)
      elif self.tag == sTag: lib.ElSparseMatrixTargetBuffer_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixTargetBuffer_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixTargetBuffer_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixTargetBuffer_z(*args)
      else: DataExcept()
    return targetBuf

  lib.ElSparseMatrixValueBuffer_i.argtypes = \
  lib.ElSparseMatrixLockedValueBuffer_i.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  lib.ElSparseMatrixValueBuffer_s.argtypes = \
  lib.ElSparseMatrixLockedValueBuffer_s.argtypes = \
    [c_void_p,POINTER(POINTER(sType))]
  lib.ElSparseMatrixValueBuffer_d.argtypes = \
  lib.ElSparseMatrixLockedValueBuffer_d.argtypes = \
    [c_void_p,POINTER(POINTER(dType))]
  lib.ElSparseMatrixValueBuffer_c.argtypes = \
  lib.ElSparseMatrixLockedValueBuffer_c.argtypes = \
    [c_void_p,POINTER(POINTER(cType))]
  lib.ElSparseMatrixValueBuffer_z.argtypes = \
  lib.ElSparseMatrixLockedValueBuffer_z.argtypes = \
    [c_void_p,POINTER(POINTER(zType))]
  def ValueBuffer(self,locked=False):
    valueBuf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(valueBuf)]
    if locked:
      if   self.tag == iTag: lib.ElSparseMatrixLockedValueBuffer_i(*args)
      elif self.tag == sTag: lib.ElSparseMatrixLockedValueBuffer_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixLockedValueBuffer_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixLockedValueBuffer_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixLockedValueBuffer_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElSparseMatrixValueBuffer_i(*args)
      elif self.tag == sTag: lib.ElSparseMatrixValueBuffer_s(*args)
      elif self.tag == dTag: lib.ElSparseMatrixValueBuffer_d(*args)
      elif self.tag == cTag: lib.ElSparseMatrixValueBuffer_c(*args)
      elif self.tag == zTag: lib.ElSparseMatrixValueBuffer_z(*args)
      else: DataExcept()
    return valueBuf

  lib.ElGetContigSubmatrixSparse_i.argtypes = \
  lib.ElGetContigSubmatrixSparse_s.argtypes = \
  lib.ElGetContigSubmatrixSparse_d.argtypes = \
  lib.ElGetContigSubmatrixSparse_c.argtypes = \
  lib.ElGetContigSubmatrixSparse_z.argtypes = \
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
    ASub = SparseMatrix(self.tag)
    args = [self.obj,iRan,jRan,ASub.obj]
    if   self.tag == iTag: lib.ElGetContigSubmatrixSparse_i(*args)
    elif self.tag == sTag: lib.ElGetContigSubmatrixSparse_s(*args)
    elif self.tag == dTag: lib.ElGetContigSubmatrixSparse_d(*args)
    elif self.tag == cTag: lib.ElGetContigSubmatrixSparse_c(*args)
    elif self.tag == zTag: lib.ElGetContigSubmatrixSparse_z(*args)
    else: DataExcept()
    return ASub
