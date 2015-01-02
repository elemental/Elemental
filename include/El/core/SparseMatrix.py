#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import Graph as G

# SparseMatrix
# ============

lib.ElSparseMatrixCreate_i.argtypes = [POINTER(c_void_p)]
lib.ElSparseMatrixCreate_i.restype = c_uint
lib.ElSparseMatrixCreate_s.argtypes = [POINTER(c_void_p)]
lib.ElSparseMatrixCreate_s.restype = c_uint
lib.ElSparseMatrixCreate_d.argtypes = [POINTER(c_void_p)]
lib.ElSparseMatrixCreate_d.restype = c_uint
lib.ElSparseMatrixCreate_c.argtypes = [POINTER(c_void_p)]
lib.ElSparseMatrixCreate_c.restype = c_uint
lib.ElSparseMatrixCreate_z.argtypes = [POINTER(c_void_p)]
lib.ElSparseMatrixCreate_z.restype = c_uint

lib.ElSparseMatrixDestroy_i.argtypes = [c_void_p]
lib.ElSparseMatrixDestroy_i.restype = c_uint
lib.ElSparseMatrixDestroy_s.argtypes = [c_void_p]
lib.ElSparseMatrixDestroy_s.restype = c_uint
lib.ElSparseMatrixDestroy_d.argtypes = [c_void_p]
lib.ElSparseMatrixDestroy_d.restype = c_uint
lib.ElSparseMatrixDestroy_c.argtypes = [c_void_p]
lib.ElSparseMatrixDestroy_c.restype = c_uint
lib.ElSparseMatrixDestroy_z.argtypes = [c_void_p]
lib.ElSparseMatrixDestroy_z.restype = c_uint

lib.ElSparseMatrixEmpty_i.argtypes = [c_void_p]
lib.ElSparseMatrixEmpty_i.restype = c_uint
lib.ElSparseMatrixEmpty_s.argtypes = [c_void_p]
lib.ElSparseMatrixEmpty_s.restype = c_uint
lib.ElSparseMatrixEmpty_d.argtypes = [c_void_p]
lib.ElSparseMatrixEmpty_d.restype = c_uint
lib.ElSparseMatrixEmpty_c.argtypes = [c_void_p]
lib.ElSparseMatrixEmpty_c.restype = c_uint
lib.ElSparseMatrixEmpty_z.argtypes = [c_void_p]
lib.ElSparseMatrixEmpty_z.restype = c_uint

lib.ElSparseMatrixResize_i.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixResize_i.restype = c_uint
lib.ElSparseMatrixResize_s.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixResize_s.restype = c_uint
lib.ElSparseMatrixResize_d.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixResize_d.restype = c_uint
lib.ElSparseMatrixResize_c.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixResize_c.restype = c_uint
lib.ElSparseMatrixResize_z.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixResize_z.restype = c_uint

lib.ElSparseMatrixReserve_i.argtypes = [c_void_p,iType]
lib.ElSparseMatrixReserve_i.restype = c_uint
lib.ElSparseMatrixReserve_s.argtypes = [c_void_p,iType]
lib.ElSparseMatrixReserve_s.restype = c_uint
lib.ElSparseMatrixReserve_d.argtypes = [c_void_p,iType]
lib.ElSparseMatrixReserve_d.restype = c_uint
lib.ElSparseMatrixReserve_c.argtypes = [c_void_p,iType]
lib.ElSparseMatrixReserve_c.restype = c_uint
lib.ElSparseMatrixReserve_z.argtypes = [c_void_p,iType]
lib.ElSparseMatrixReserve_z.restype = c_uint

lib.ElSparseMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElSparseMatrixUpdate_i.restype = c_uint
lib.ElSparseMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElSparseMatrixUpdate_s.restype = c_uint
lib.ElSparseMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElSparseMatrixUpdate_d.restype = c_uint
lib.ElSparseMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElSparseMatrixUpdate_c.restype = c_uint
lib.ElSparseMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElSparseMatrixUpdate_z.restype = c_uint

lib.ElSparseMatrixZero_i.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixZero_i.restype = c_uint
lib.ElSparseMatrixZero_s.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixZero_s.restype = c_uint
lib.ElSparseMatrixZero_d.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixZero_d.restype = c_uint
lib.ElSparseMatrixZero_c.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixZero_c.restype = c_uint
lib.ElSparseMatrixZero_z.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixZero_z.restype = c_uint

lib.ElSparseMatrixQueueUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElSparseMatrixQueueUpdate_i.restype = c_uint
lib.ElSparseMatrixQueueUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElSparseMatrixQueueUpdate_s.restype = c_uint
lib.ElSparseMatrixQueueUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElSparseMatrixQueueUpdate_d.restype = c_uint
lib.ElSparseMatrixQueueUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElSparseMatrixQueueUpdate_c.restype = c_uint
lib.ElSparseMatrixQueueUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElSparseMatrixQueueUpdate_z.restype = c_uint

lib.ElSparseMatrixQueueZero_i.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixQueueZero_i.restype = c_uint
lib.ElSparseMatrixQueueZero_s.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixQueueZero_s.restype = c_uint
lib.ElSparseMatrixQueueZero_d.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixQueueZero_d.restype = c_uint
lib.ElSparseMatrixQueueZero_c.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixQueueZero_c.restype = c_uint
lib.ElSparseMatrixQueueZero_z.argtypes = [c_void_p,iType,iType]
lib.ElSparseMatrixQueueZero_z.restype = c_uint

lib.ElSparseMatrixMakeConsistent_i.argtypes = [c_void_p]
lib.ElSparseMatrixMakeConsistent_i.restype = c_uint
lib.ElSparseMatrixMakeConsistent_s.argtypes = [c_void_p]
lib.ElSparseMatrixMakeConsistent_s.restype = c_uint
lib.ElSparseMatrixMakeConsistent_d.argtypes = [c_void_p]
lib.ElSparseMatrixMakeConsistent_d.restype = c_uint
lib.ElSparseMatrixMakeConsistent_c.argtypes = [c_void_p]
lib.ElSparseMatrixMakeConsistent_c.restype = c_uint
lib.ElSparseMatrixMakeConsistent_z.argtypes = [c_void_p]
lib.ElSparseMatrixMakeConsistent_z.restype = c_uint

lib.ElSparseMatrixHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixHeight_i.restype = c_uint
lib.ElSparseMatrixHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixHeight_s.restype = c_uint
lib.ElSparseMatrixHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixHeight_d.restype = c_uint
lib.ElSparseMatrixHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixHeight_c.restype = c_uint
lib.ElSparseMatrixHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixHeight_z.restype = c_uint

lib.ElSparseMatrixWidth_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixWidth_i.restype = c_uint
lib.ElSparseMatrixWidth_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixWidth_s.restype = c_uint
lib.ElSparseMatrixWidth_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixWidth_d.restype = c_uint
lib.ElSparseMatrixWidth_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixWidth_c.restype = c_uint
lib.ElSparseMatrixWidth_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixWidth_z.restype = c_uint

lib.ElSparseMatrixNumEntries_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixNumEntries_i.restype = c_uint
lib.ElSparseMatrixNumEntries_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixNumEntries_s.restype = c_uint
lib.ElSparseMatrixNumEntries_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixNumEntries_d.restype = c_uint
lib.ElSparseMatrixNumEntries_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixNumEntries_c.restype = c_uint
lib.ElSparseMatrixNumEntries_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixNumEntries_z.restype = c_uint

lib.ElSparseMatrixCapacity_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixCapacity_i.restype = c_uint
lib.ElSparseMatrixCapacity_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixCapacity_s.restype = c_uint
lib.ElSparseMatrixCapacity_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixCapacity_d.restype = c_uint
lib.ElSparseMatrixCapacity_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixCapacity_c.restype = c_uint
lib.ElSparseMatrixCapacity_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElSparseMatrixCapacity_z.restype = c_uint

lib.ElSparseMatrixConsistent_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElSparseMatrixConsistent_i.restype = c_uint
lib.ElSparseMatrixConsistent_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElSparseMatrixConsistent_s.restype = c_uint
lib.ElSparseMatrixConsistent_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElSparseMatrixConsistent_d.restype = c_uint
lib.ElSparseMatrixConsistent_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElSparseMatrixConsistent_c.restype = c_uint
lib.ElSparseMatrixConsistent_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElSparseMatrixConsistent_z.restype = c_uint

lib.ElSparseMatrixGraph_i.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixGraph_i.restype = c_uint
lib.ElSparseMatrixGraph_s.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixGraph_s.restype = c_uint
lib.ElSparseMatrixGraph_d.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixGraph_d.restype = c_uint
lib.ElSparseMatrixGraph_c.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixGraph_c.restype = c_uint
lib.ElSparseMatrixGraph_z.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixGraph_z.restype = c_uint

lib.ElSparseMatrixLockedGraph_i.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixLockedGraph_i.restype = c_uint
lib.ElSparseMatrixLockedGraph_s.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixLockedGraph_s.restype = c_uint
lib.ElSparseMatrixLockedGraph_d.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixLockedGraph_d.restype = c_uint
lib.ElSparseMatrixLockedGraph_c.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixLockedGraph_c.restype = c_uint
lib.ElSparseMatrixLockedGraph_z.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElSparseMatrixLockedGraph_z.restype = c_uint

lib.ElSparseMatrixRow_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixRow_i.restype = c_uint
lib.ElSparseMatrixRow_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixRow_s.restype = c_uint
lib.ElSparseMatrixRow_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixRow_d.restype = c_uint
lib.ElSparseMatrixRow_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixRow_c.restype = c_uint
lib.ElSparseMatrixRow_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixRow_z.restype = c_uint

lib.ElSparseMatrixCol_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixCol_i.restype = c_uint
lib.ElSparseMatrixCol_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixCol_s.restype = c_uint
lib.ElSparseMatrixCol_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixCol_d.restype = c_uint
lib.ElSparseMatrixCol_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixCol_c.restype = c_uint
lib.ElSparseMatrixCol_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixCol_z.restype = c_uint

lib.ElSparseMatrixValue_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixValue_i.restype = c_uint
lib.ElSparseMatrixValue_s.argtypes = [c_void_p,iType,POINTER(sType)]
lib.ElSparseMatrixValue_s.restype = c_uint
lib.ElSparseMatrixValue_d.argtypes = [c_void_p,iType,POINTER(dType)]
lib.ElSparseMatrixValue_d.restype = c_uint
lib.ElSparseMatrixValue_c.argtypes = [c_void_p,iType,POINTER(cType)]
lib.ElSparseMatrixValue_c.restype = c_uint
lib.ElSparseMatrixValue_z.argtypes = [c_void_p,iType,POINTER(zType)]
lib.ElSparseMatrixValue_z.restype = c_uint

lib.ElSparseMatrixEntryOffset_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixEntryOffset_i.restype = c_uint
lib.ElSparseMatrixEntryOffset_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixEntryOffset_s.restype = c_uint
lib.ElSparseMatrixEntryOffset_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixEntryOffset_d.restype = c_uint
lib.ElSparseMatrixEntryOffset_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixEntryOffset_c.restype = c_uint
lib.ElSparseMatrixEntryOffset_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixEntryOffset_z.restype = c_uint

lib.ElSparseMatrixNumConnections_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixNumConnections_i.restype = c_uint
lib.ElSparseMatrixNumConnections_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixNumConnections_s.restype = c_uint
lib.ElSparseMatrixNumConnections_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixNumConnections_d.restype = c_uint
lib.ElSparseMatrixNumConnections_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixNumConnections_c.restype = c_uint
lib.ElSparseMatrixNumConnections_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElSparseMatrixNumConnections_z.restype = c_uint

lib.ElSparseMatrixSourceBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixSourceBuffer_i.restype = c_uint
lib.ElSparseMatrixSourceBuffer_s.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixSourceBuffer_s.restype = c_uint
lib.ElSparseMatrixSourceBuffer_d.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixSourceBuffer_d.restype = c_uint
lib.ElSparseMatrixSourceBuffer_c.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixSourceBuffer_c.restype = c_uint
lib.ElSparseMatrixSourceBuffer_z.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixSourceBuffer_z.restype = c_uint

lib.ElSparseMatrixLockedSourceBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedSourceBuffer_i.restype = c_uint
lib.ElSparseMatrixLockedSourceBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedSourceBuffer_s.restype = c_uint
lib.ElSparseMatrixLockedSourceBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedSourceBuffer_d.restype = c_uint
lib.ElSparseMatrixLockedSourceBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedSourceBuffer_c.restype = c_uint
lib.ElSparseMatrixLockedSourceBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedSourceBuffer_z.restype = c_uint

lib.ElSparseMatrixTargetBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixTargetBuffer_i.restype = c_uint
lib.ElSparseMatrixTargetBuffer_s.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixTargetBuffer_s.restype = c_uint
lib.ElSparseMatrixTargetBuffer_d.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixTargetBuffer_d.restype = c_uint
lib.ElSparseMatrixTargetBuffer_c.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixTargetBuffer_c.restype = c_uint
lib.ElSparseMatrixTargetBuffer_z.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixTargetBuffer_z.restype = c_uint

lib.ElSparseMatrixLockedTargetBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedTargetBuffer_i.restype = c_uint
lib.ElSparseMatrixLockedTargetBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedTargetBuffer_s.restype = c_uint
lib.ElSparseMatrixLockedTargetBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedTargetBuffer_d.restype = c_uint
lib.ElSparseMatrixLockedTargetBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedTargetBuffer_c.restype = c_uint
lib.ElSparseMatrixLockedTargetBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedTargetBuffer_z.restype = c_uint

lib.ElSparseMatrixValueBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixValueBuffer_i.restype = c_uint
lib.ElSparseMatrixValueBuffer_s.argtypes = [c_void_p,POINTER(POINTER(sType))]
lib.ElSparseMatrixValueBuffer_s.restype = c_uint
lib.ElSparseMatrixValueBuffer_d.argtypes = [c_void_p,POINTER(POINTER(dType))]
lib.ElSparseMatrixValueBuffer_d.restype = c_uint
lib.ElSparseMatrixValueBuffer_c.argtypes = [c_void_p,POINTER(POINTER(cType))]
lib.ElSparseMatrixValueBuffer_c.restype = c_uint
lib.ElSparseMatrixValueBuffer_z.argtypes = [c_void_p,POINTER(POINTER(zType))]
lib.ElSparseMatrixValueBuffer_z.restype = c_uint

lib.ElSparseMatrixLockedValueBuffer_i.argtypes = \
  [c_void_p,POINTER(POINTER(iType))]
lib.ElSparseMatrixLockedValueBuffer_i.restype = c_uint
lib.ElSparseMatrixLockedValueBuffer_s.argtypes = \
  [c_void_p,POINTER(POINTER(sType))]
lib.ElSparseMatrixLockedValueBuffer_s.restype = c_uint
lib.ElSparseMatrixLockedValueBuffer_d.argtypes = \
  [c_void_p,POINTER(POINTER(dType))]
lib.ElSparseMatrixLockedValueBuffer_d.restype = c_uint
lib.ElSparseMatrixLockedValueBuffer_c.argtypes = \
  [c_void_p,POINTER(POINTER(cType))]
lib.ElSparseMatrixLockedValueBuffer_c.restype = c_uint
lib.ElSparseMatrixLockedValueBuffer_z.argtypes = \
  [c_void_p,POINTER(POINTER(zType))]
lib.ElSparseMatrixLockedValueBuffer_z.restype = c_uint

class SparseMatrix(object):
  # Constructors and destructors
  # ============================
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
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElSparseMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixEmpty_z(*args)
    else: DataExcept()
  def Resize(self,height,width):
    args = [self.obj,height,width]
    if   self.tag == iTag: lib.ElSparseMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixResize_z(*args)
    else: DataExcept()
  def Reserve(self,numEntries):
    args = [self.obj,numEntries]
    if   self.tag == iTag: lib.ElSparseMatrixReserve_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixReserve_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixReserve_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixReserve_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixReserve_z(*args)
    else: DataExcept()
  def Update(self,row,col,value):
    args = [self.obj,row,col,value]
    if   self.tag == iTag: lib.ElSparseMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixUpdate_z(*args)
    else: DataExcept()
  def Zero(self,row,col):
    args = [self.obj,row,col]
    if   self.tag == iTag: lib.ElSparseMatrixZero_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixZero_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixZero_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixZero_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixZero_z(*args)
    else: DataExcept()
  def QueueUpdate(self,row,col,value):
    args = [self.obj,row,col,value]
    if   self.tag == iTag: lib.ElSparseMatrixQueueUpdate_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixQueueUpdate_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixQueueUpdate_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixQueueUpdate_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixQueueUpdate_z(*args)
    else: DataExcept()
  def QueueZero(self,row,col):
    args = [self.obj,row,col]
    if   self.tag == iTag: lib.ElSparseMatrixQueueZero_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixQueueZero_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixQueueZero_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixQueueZero_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixQueueZero_z(*args)
    else: DataExcept()
  def MakeConsistent(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElSparseMatrixMakeConsistent_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixMakeConsistent_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixMakeConsistent_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixMakeConsistent_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixMakeConsistent_z(*args)
    else: DataExcept()
  # Queries
  # =======
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
  def Graph(self):
    graph = G.Graph(False)
    args = [self.obj,pointer(graph.obj)]
    if   self.tag == iTag: lib.ElSparseMatrixGraph_i(*args)  
    elif self.tag == sTag: lib.ElSparseMatrixGraph_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixGraph_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixGraph_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixGraph_z(*args)
    else: DataExcept()
    return graph
  def LockedGraph(self):
    graph = G.Graph(False)
    args = [self.obj,pointer(graph.obj)]
    if   self.tag == iTag: lib.ElSparseMatrixLockedGraph_i(*args)  
    elif self.tag == sTag: lib.ElSparseMatrixLockedGraph_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixLockedGraph_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixLockedGraph_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixLockedGraph_z(*args)
    else: DataExcept()
    return graph
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
  def Value(self,index):
    value =  TagToType(self.tag)()
    args = [self.obj,index,pointer(value)]
    if   self.tag == iTag: lib.ElSparseMatrixValue_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixValue_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixValue_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixValue_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixValue_z(*args)
    else: DataExcept()
    return value.value
  def EntryOffset(self,row):
    offset = iType()
    args = [self.obj,row,pointer(offset)]
    if   self.tag == iTag: lib.ElSparseMatrixEntryOffset_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixEntryOffset_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixEntryOffset_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixEntryOffset_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixEntryOffset_z(*args)
    else: DataExcept()
    return offset.value
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
  def SourceBuffer(self):
    sourceBuf = POINTER(iType)()
    args = [self.obj,pointer(sourceBuf)]    
    if   self.tag == iTag: lib.ElSparseMatrixSourceBuffer_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixSourceBuffer_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixSourceBuffer_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixSourceBuffer_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixSourceBuffer_z(*args)
    else: DataExcept()
    return sourceBuf
  def LockedSourceBuffer(self):
    sourceBuf = POINTER(iType)()
    args = [self.obj,pointer(sourceBuf)]    
    if   self.tag == iTag: lib.ElSparseMatrixLockedSourceBuffer_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixLockedSourceBuffer_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixLockedSourceBuffer_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixLockedSourceBuffer_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixLockedSourceBuffer_z(*args)
    else: DataExcept()
    return sourceBuf
  def TargetBuffer(self):
    targetBuf = POINTER(iType)()
    args = [self.obj,pointer(targetBuf)]    
    if   self.tag == iTag: lib.ElSparseMatrixTargetBuffer_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixTargetBuffer_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixTargetBuffer_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixTargetBuffer_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixTargetBuffer_z(*args)
    else: DataExcept()
    return targetBuf
  def LockedTargetBuffer(self):
    targetBuf = POINTER(iType)()
    args = [self.obj,pointer(targetBuf)]    
    if   self.tag == iTag: lib.ElSparseMatrixLockedTargetBuffer_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixLockedTargetBuffer_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixLockedTargetBuffer_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixLockedTargetBuffer_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixLockedTargetBuffer_z(*args)
    else: DataExcept()
    return targetBuf
  def ValueBuffer(self):
    valueBuf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(valueBuf)]
    if   self.tag == iTag: lib.ElSparseMatrixValueBuffer_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixValueBuffer_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixValueBuffer_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixValueBuffer_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixValueBuffer_z(*args)
    else: DataExcept()
    return valueBuf
  def LockedValueBuffer(self):
    valueBuf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(valueBuf)]
    if   self.tag == iTag: lib.ElSparseMatrixLockedValueBuffer_i(*args)
    elif self.tag == sTag: lib.ElSparseMatrixLockedValueBuffer_s(*args)
    elif self.tag == dTag: lib.ElSparseMatrixLockedValueBuffer_d(*args)
    elif self.tag == cTag: lib.ElSparseMatrixLockedValueBuffer_c(*args)
    elif self.tag == zTag: lib.ElSparseMatrixLockedValueBuffer_z(*args)
    else: DataExcept()
    return valueBuf
