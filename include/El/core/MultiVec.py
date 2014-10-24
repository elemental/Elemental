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

# MultiVec
# ========

lib.ElMultiVecCreate_i.argtypes = [POINTER(c_void_p)]
lib.ElMultiVecCreate_i.restype = c_uint
lib.ElMultiVecCreate_s.argtypes = [POINTER(c_void_p)]
lib.ElMultiVecCreate_s.restype = c_uint
lib.ElMultiVecCreate_d.argtypes = [POINTER(c_void_p)]
lib.ElMultiVecCreate_d.restype = c_uint
lib.ElMultiVecCreate_c.argtypes = [POINTER(c_void_p)]
lib.ElMultiVecCreate_c.restype = c_uint
lib.ElMultiVecCreate_z.argtypes = [POINTER(c_void_p)]
lib.ElMultiVecCreate_z.restype = c_uint

lib.ElMultiVecDestroy_i.argtypes = [c_void_p]
lib.ElMultiVecDestroy_i.restype = c_uint
lib.ElMultiVecDestroy_s.argtypes = [c_void_p]
lib.ElMultiVecDestroy_s.restype = c_uint
lib.ElMultiVecDestroy_d.argtypes = [c_void_p]
lib.ElMultiVecDestroy_d.restype = c_uint
lib.ElMultiVecDestroy_c.argtypes = [c_void_p]
lib.ElMultiVecDestroy_c.restype = c_uint
lib.ElMultiVecDestroy_z.argtypes = [c_void_p]
lib.ElMultiVecDestroy_z.restype = c_uint

lib.ElMultiVecEmpty_i.argtypes = [c_void_p]
lib.ElMultiVecEmpty_i.restype = c_uint
lib.ElMultiVecEmpty_s.argtypes = [c_void_p]
lib.ElMultiVecEmpty_s.restype = c_uint
lib.ElMultiVecEmpty_d.argtypes = [c_void_p]
lib.ElMultiVecEmpty_d.restype = c_uint
lib.ElMultiVecEmpty_c.argtypes = [c_void_p]
lib.ElMultiVecEmpty_c.restype = c_uint
lib.ElMultiVecEmpty_z.argtypes = [c_void_p]
lib.ElMultiVecEmpty_z.restype = c_uint

lib.ElMultiVecResize_i.argtypes = [c_void_p,iType,iType]
lib.ElMultiVecResize_i.restype = c_uint
lib.ElMultiVecResize_s.argtypes = [c_void_p,iType,iType]
lib.ElMultiVecResize_s.restype = c_uint
lib.ElMultiVecResize_d.argtypes = [c_void_p,iType,iType]
lib.ElMultiVecResize_d.restype = c_uint
lib.ElMultiVecResize_c.argtypes = [c_void_p,iType,iType]
lib.ElMultiVecResize_c.restype = c_uint
lib.ElMultiVecResize_z.argtypes = [c_void_p,iType,iType]
lib.ElMultiVecResize_z.restype = c_uint

lib.ElMultiVecHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecHeight_i.restype = c_uint
lib.ElMultiVecHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecHeight_s.restype = c_uint
lib.ElMultiVecHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecHeight_d.restype = c_uint
lib.ElMultiVecHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecHeight_c.restype = c_uint
lib.ElMultiVecHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecHeight_z.restype = c_uint

lib.ElMultiVecWidth_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecWidth_i.restype = c_uint
lib.ElMultiVecWidth_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecWidth_s.restype = c_uint
lib.ElMultiVecWidth_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecWidth_d.restype = c_uint
lib.ElMultiVecWidth_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecWidth_c.restype = c_uint
lib.ElMultiVecWidth_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElMultiVecWidth_z.restype = c_uint

lib.ElMultiVecGet_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElMultiVecGet_i.restype = c_uint
lib.ElMultiVecGet_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElMultiVecGet_s.restype = c_uint
lib.ElMultiVecGet_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElMultiVecGet_d.restype = c_uint
lib.ElMultiVecGet_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
lib.ElMultiVecGet_c.restype = c_uint
lib.ElMultiVecGet_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
lib.ElMultiVecGet_z.restype = c_uint

lib.ElMultiVecSet_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElMultiVecSet_i.restype = c_uint
lib.ElMultiVecSet_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElMultiVecSet_s.restype = c_uint
lib.ElMultiVecSet_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElMultiVecSet_d.restype = c_uint
lib.ElMultiVecSet_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElMultiVecSet_c.restype = c_uint
lib.ElMultiVecSet_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElMultiVecSet_z.restype = c_uint

lib.ElMultiVecUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElMultiVecUpdate_i.restype = c_uint
lib.ElMultiVecUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElMultiVecUpdate_s.restype = c_uint
lib.ElMultiVecUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElMultiVecUpdate_d.restype = c_uint
lib.ElMultiVecUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElMultiVecUpdate_c.restype = c_uint
lib.ElMultiVecUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElMultiVecUpdate_z.restype = c_uint

class MultiVec(object):
  # Constructors and destructors
  # ============================
  def __init__(self,tag=dTag,create=True):
    self.obj = c_void_p()
    self.tag = tag
    CheckTag(tag)
    if create:
      args = [pointer(self.obj)]
      if   tag == iTag: lib.ElMultiVecCreate_i(*args)
      elif tag == sTag: lib.ElMultiVecCreate_s(*args)
      elif tag == dTag: lib.ElMultiVecCreate_d(*args)
      elif tag == cTag: lib.ElMultiVecCreate_c(*args)
      elif tag == zTag: lib.ElMultiVecCreate_z(*args)
      else: DataExcept()
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElMultiVecDestroy_i(*args)
    elif self.tag == sTag: lib.ElMultiVecDestroy_s(*args)
    elif self.tag == dTag: lib.ElMultiVecDestroy_d(*args)
    elif self.tag == cTag: lib.ElMultiVecDestroy_c(*args)
    elif self.tag == zTag: lib.ElMultiVecDestroy_z(*args)
    else: DataExcept()
  # Assignment and reconfiguration
  # ==============================
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElMultiVecEmpty_i(*args)
    elif self.tag == sTag: lib.ElMultiVecEmpty_s(*args)
    elif self.tag == dTag: lib.ElMultiVecEmpty_d(*args)
    elif self.tag == cTag: lib.ElMultiVecEmpty_c(*args)
    elif self.tag == zTag: lib.ElMultiVecEmpty_z(*args)
    else: DataExcept()
  def Resize(self,height,width):
    args = [self.obj,height,width]
    if   self.tag == iTag: lib.ElMultiVecResize_i(*args)
    elif self.tag == sTag: lib.ElMultiVecResize_s(*args)
    elif self.tag == dTag: lib.ElMultiVecResize_d(*args)
    elif self.tag == cTag: lib.ElMultiVecResize_c(*args)
    elif self.tag == zTag: lib.ElMultiVecResize_z(*args)
    else: DataExcept()
  # Queries
  # =======
  def Height(self):
    height = iType()
    args = [self.obj,pointer(height)]
    if   self.tag == iTag: lib.ElMultiVecHeight_i(*args)
    elif self.tag == sTag: lib.ElMultiVecHeight_s(*args)
    elif self.tag == dTag: lib.ElMultiVecHeight_d(*args)
    elif self.tag == cTag: lib.ElMultiVecHeight_c(*args)
    elif self.tag == zTag: lib.ElMultiVecHeight_z(*args)
    else: DataExcept()
    return height.value
  def Width(self):
    width = iType()
    args = [self.obj,pointer(width)]
    if   self.tag == iTag: lib.ElMultiVecWidth_i(*args)
    elif self.tag == sTag: lib.ElMultiVecWidth_s(*args)
    elif self.tag == dTag: lib.ElMultiVecWidth_d(*args)
    elif self.tag == cTag: lib.ElMultiVecWidth_c(*args)
    elif self.tag == zTag: lib.ElMultiVecWidth_z(*args)
    else: DataExcept()
    return width.value
  # Entrywise manipulation
  # ======================
  def Get(self,i,j):
    value = TagToType(self.tag)()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElMultiVecGet_i(*args)
    elif self.tag == sTag: lib.ElMultiVecGet_s(*args)
    elif self.tag == dTag: lib.ElMultiVecGet_d(*args)
    elif self.tag == cTag: lib.ElMultiVecGet_c(*args)
    elif self.tag == zTag: lib.ElMultiVecGet_z(*args)
    else: DataExcept()
    return value.value
  def Set(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElMultiVecSet_i(*args)
    elif self.tag == sTag: lib.ElMultiVecSet_s(*args)
    elif self.tag == dTag: lib.ElMultiVecSet_d(*args)
    elif self.tag == cTag: lib.ElMultiVecSet_c(*args)
    elif self.tag == zTag: lib.ElMultiVecSet_z(*args)
    else: DataExcept()
  def Update(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElMultiVecUpdate_i(*args)
    elif self.tag == sTag: lib.ElMultiVecUpdate_s(*args)
    elif self.tag == dTag: lib.ElMultiVecUpdate_d(*args)
    elif self.tag == cTag: lib.ElMultiVecUpdate_c(*args)
    elif self.tag == zTag: lib.ElMultiVecUpdate_z(*args)
    else: DataExcept()
