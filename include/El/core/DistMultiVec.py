#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import numpy as np

# DistMultiVec
# ========

lib.ElDistMultiVecCreate_i.argtypes = [POINTER(c_void_p),MPI_Comm]
lib.ElDistMultiVecCreate_i.restype = c_uint
lib.ElDistMultiVecCreate_s.argtypes = [POINTER(c_void_p),MPI_Comm]
lib.ElDistMultiVecCreate_s.restype = c_uint
lib.ElDistMultiVecCreate_d.argtypes = [POINTER(c_void_p),MPI_Comm]
lib.ElDistMultiVecCreate_d.restype = c_uint
lib.ElDistMultiVecCreate_c.argtypes = [POINTER(c_void_p),MPI_Comm]
lib.ElDistMultiVecCreate_c.restype = c_uint
lib.ElDistMultiVecCreate_z.argtypes = [POINTER(c_void_p),MPI_Comm]
lib.ElDistMultiVecCreate_z.restype = c_uint

lib.ElDistMultiVecDestroy_i.argtypes = [c_void_p]
lib.ElDistMultiVecDestroy_i.restype = c_uint
lib.ElDistMultiVecDestroy_s.argtypes = [c_void_p]
lib.ElDistMultiVecDestroy_s.restype = c_uint
lib.ElDistMultiVecDestroy_d.argtypes = [c_void_p]
lib.ElDistMultiVecDestroy_d.restype = c_uint
lib.ElDistMultiVecDestroy_c.argtypes = [c_void_p]
lib.ElDistMultiVecDestroy_c.restype = c_uint
lib.ElDistMultiVecDestroy_z.argtypes = [c_void_p]
lib.ElDistMultiVecDestroy_z.restype = c_uint

lib.ElDistMultiVecEmpty_i.argtypes = [c_void_p]
lib.ElDistMultiVecEmpty_i.restype = c_uint
lib.ElDistMultiVecEmpty_s.argtypes = [c_void_p]
lib.ElDistMultiVecEmpty_s.restype = c_uint
lib.ElDistMultiVecEmpty_d.argtypes = [c_void_p]
lib.ElDistMultiVecEmpty_d.restype = c_uint
lib.ElDistMultiVecEmpty_c.argtypes = [c_void_p]
lib.ElDistMultiVecEmpty_c.restype = c_uint
lib.ElDistMultiVecEmpty_z.argtypes = [c_void_p]
lib.ElDistMultiVecEmpty_z.restype = c_uint

lib.ElDistMultiVecResize_i.argtypes = [c_void_p,iType,iType]
lib.ElDistMultiVecResize_i.restype = c_uint
lib.ElDistMultiVecResize_s.argtypes = [c_void_p,iType,iType]
lib.ElDistMultiVecResize_s.restype = c_uint
lib.ElDistMultiVecResize_d.argtypes = [c_void_p,iType,iType]
lib.ElDistMultiVecResize_d.restype = c_uint
lib.ElDistMultiVecResize_c.argtypes = [c_void_p,iType,iType]
lib.ElDistMultiVecResize_c.restype = c_uint
lib.ElDistMultiVecResize_z.argtypes = [c_void_p,iType,iType]
lib.ElDistMultiVecResize_z.restype = c_uint

lib.ElDistMultiVecSetComm_i.argtypes = [c_void_p,MPI_Comm]
lib.ElDistMultiVecSetComm_i.restype = c_uint
lib.ElDistMultiVecSetComm_s.argtypes = [c_void_p,MPI_Comm]
lib.ElDistMultiVecSetComm_s.restype = c_uint
lib.ElDistMultiVecSetComm_d.argtypes = [c_void_p,MPI_Comm]
lib.ElDistMultiVecSetComm_d.restype = c_uint
lib.ElDistMultiVecSetComm_c.argtypes = [c_void_p,MPI_Comm]
lib.ElDistMultiVecSetComm_c.restype = c_uint
lib.ElDistMultiVecSetComm_z.argtypes = [c_void_p,MPI_Comm]
lib.ElDistMultiVecSetComm_z.restype = c_uint

lib.ElDistMultiVecHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecHeight_i.restype = c_uint
lib.ElDistMultiVecHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecHeight_s.restype = c_uint
lib.ElDistMultiVecHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecHeight_d.restype = c_uint
lib.ElDistMultiVecHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecHeight_c.restype = c_uint
lib.ElDistMultiVecHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecHeight_z.restype = c_uint

lib.ElDistMultiVecWidth_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecWidth_i.restype = c_uint
lib.ElDistMultiVecWidth_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecWidth_s.restype = c_uint
lib.ElDistMultiVecWidth_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecWidth_d.restype = c_uint
lib.ElDistMultiVecWidth_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecWidth_c.restype = c_uint
lib.ElDistMultiVecWidth_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecWidth_z.restype = c_uint

lib.ElDistMultiVecFirstLocalRow_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecFirstLocalRow_i.restype = c_uint
lib.ElDistMultiVecFirstLocalRow_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecFirstLocalRow_s.restype = c_uint
lib.ElDistMultiVecFirstLocalRow_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecFirstLocalRow_d.restype = c_uint
lib.ElDistMultiVecFirstLocalRow_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecFirstLocalRow_c.restype = c_uint
lib.ElDistMultiVecFirstLocalRow_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecFirstLocalRow_z.restype = c_uint

lib.ElDistMultiVecLocalHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecLocalHeight_i.restype = c_uint
lib.ElDistMultiVecLocalHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecLocalHeight_s.restype = c_uint
lib.ElDistMultiVecLocalHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecLocalHeight_d.restype = c_uint
lib.ElDistMultiVecLocalHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecLocalHeight_c.restype = c_uint
lib.ElDistMultiVecLocalHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecLocalHeight_z.restype = c_uint

lib.ElDistMultiVecComm_i.argtypes = [c_void_p,POINTER(MPI_Comm)]
lib.ElDistMultiVecComm_i.restype = c_uint
lib.ElDistMultiVecComm_s.argtypes = [c_void_p,POINTER(MPI_Comm)]
lib.ElDistMultiVecComm_s.restype = c_uint
lib.ElDistMultiVecComm_d.argtypes = [c_void_p,POINTER(MPI_Comm)]
lib.ElDistMultiVecComm_d.restype = c_uint
lib.ElDistMultiVecComm_c.argtypes = [c_void_p,POINTER(MPI_Comm)]
lib.ElDistMultiVecComm_c.restype = c_uint
lib.ElDistMultiVecComm_z.argtypes = [c_void_p,POINTER(MPI_Comm)]
lib.ElDistMultiVecComm_z.restype = c_uint

lib.ElDistMultiVecBlocksize_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecBlocksize_i.restype = c_uint
lib.ElDistMultiVecBlocksize_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecBlocksize_s.restype = c_uint
lib.ElDistMultiVecBlocksize_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecBlocksize_d.restype = c_uint
lib.ElDistMultiVecBlocksize_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecBlocksize_c.restype = c_uint
lib.ElDistMultiVecBlocksize_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMultiVecBlocksize_z.restype = c_uint

lib.ElDistMultiVecGetLocal_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElDistMultiVecGetLocal_i.restype = c_uint
lib.ElDistMultiVecGetLocal_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElDistMultiVecGetLocal_s.restype = c_uint
lib.ElDistMultiVecGetLocal_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElDistMultiVecGetLocal_d.restype = c_uint
lib.ElDistMultiVecGetLocal_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
lib.ElDistMultiVecGetLocal_c.restype = c_uint
lib.ElDistMultiVecGetLocal_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
lib.ElDistMultiVecGetLocal_z.restype = c_uint

lib.ElDistMultiVecSetLocal_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMultiVecSetLocal_i.restype = c_uint
lib.ElDistMultiVecSetLocal_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMultiVecSetLocal_s.restype = c_uint
lib.ElDistMultiVecSetLocal_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMultiVecSetLocal_d.restype = c_uint
lib.ElDistMultiVecSetLocal_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistMultiVecSetLocal_c.restype = c_uint
lib.ElDistMultiVecSetLocal_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistMultiVecSetLocal_z.restype = c_uint

lib.ElDistMultiVecUpdateLocal_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMultiVecUpdateLocal_i.restype = c_uint
lib.ElDistMultiVecUpdateLocal_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMultiVecUpdateLocal_s.restype = c_uint
lib.ElDistMultiVecUpdateLocal_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMultiVecUpdateLocal_d.restype = c_uint
lib.ElDistMultiVecUpdateLocal_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistMultiVecUpdateLocal_c.restype = c_uint
lib.ElDistMultiVecUpdateLocal_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistMultiVecUpdateLocal_z.restype = c_uint

class DistMultiVec(object):
  # Constructors and destructors
  # ============================
  def __init__(self,tag=dTag,comm=MPI_COMM_WORLD(),create=True):
    self.obj = c_void_p()
    self.tag = tag
    CheckTag(tag)
    if create:
      args = [pointer(self.obj),comm]
      if   tag == iTag: lib.ElDistMultiVecCreate_i(*args)
      elif tag == sTag: lib.ElDistMultiVecCreate_s(*args)
      elif tag == dTag: lib.ElDistMultiVecCreate_d(*args)
      elif tag == cTag: lib.ElDistMultiVecCreate_c(*args)
      elif tag == zTag: lib.ElDistMultiVecCreate_z(*args)
      else: DataExcept()
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
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMultiVecEmpty_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecEmpty_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecEmpty_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecEmpty_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecEmpty_z(*args)
    else: DataExcept()
  def Resize(self,height,width):
    args = [self.obj,height,width]
    if   self.tag == iTag: lib.ElDistMultiVecResize_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecResize_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecResize_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecResize_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecResize_z(*args)
    else: DataExcept()
  def SetComm(self,comm):
    args = [self.obj,comm]
    if   self.tag == iTag: lib.ElDistMultiVecSetComm_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecSetComm_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecSetComm_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecSetComm_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecSetComm_z(*args)
    else: DataExcept()
  # Queries
  # =======
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
  def LocalHeight(self):
    localHeight = iType()
    args = [self.obj,pointer(firstLocalRow)]
    if   self.tag == iTag: lib.ElDistMultiVecLocalHeight_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecLocalHeight_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecLocalHeight_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecLocalHeight_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecLocalHeight_z(*args)
    else: DataExcept()
    return localHeight.value
  def Comm(self):
    comm = MPI_Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMultiVecComm_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecComm_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecComm_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecComm_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecComm_z(*args)
    else: DataExcept()
    return comm
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
  # Entrywise manipulation
  # ======================
  def GetLocal(self,iLocal,j):
    value = TagToType(self.tag)()
    args = [self.obj,iLocal,j,pointer(value)]
    if   self.tag == iTag: lib.ElDistMultiVecGetLocal_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecGetLocal_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecGetLocal_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecGetLocal_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecGetLocal_z(*args)
    else: DataExcept()
    return value.value
  def Set(self,iLocal,j,value):
    args = [self.obj,iLocal,j,value]
    if   self.tag == iTag: lib.ElDistMultiVecSetLocal_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecSetLocal_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecSetLocal_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecSetLocal_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecSetLocal_z(*args)
    else: DataExcept()
  def Update(self,iLocal,j,value):
    args = [self.obj,iLocal,j,value]
    if   self.tag == iTag: lib.ElDistMultiVecUpdateLocal_i(*args)
    elif self.tag == sTag: lib.ElDistMultiVecUpdateLocal_s(*args)
    elif self.tag == dTag: lib.ElDistMultiVecUpdateLocal_d(*args)
    elif self.tag == cTag: lib.ElDistMultiVecUpdateLocal_c(*args)
    elif self.tag == zTag: lib.ElDistMultiVecUpdateLocal_z(*args)
    else: DataExcept()
