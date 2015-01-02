#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import numpy as np

buffer_from_memory = pythonapi.PyBuffer_FromMemory
buffer_from_memory.restype = ctypes.py_object

buffer_from_memory_RW = pythonapi.PyBuffer_FromReadWriteMemory
buffer_from_memory_RW.restype = ctypes.py_object

# Matrix
# ======

lib.ElMatrixCreate_i.argtypes = [POINTER(c_void_p)]
lib.ElMatrixCreate_i.restype = c_uint
lib.ElMatrixCreate_s.argtypes = [POINTER(c_void_p)]
lib.ElMatrixCreate_s.restype = c_uint
lib.ElMatrixCreate_d.argtypes = [POINTER(c_void_p)]
lib.ElMatrixCreate_d.restype = c_uint
lib.ElMatrixCreate_c.argtypes = [POINTER(c_void_p)]
lib.ElMatrixCreate_c.restype = c_uint
lib.ElMatrixCreate_z.argtypes = [POINTER(c_void_p)]
lib.ElMatrixCreate_z.restype = c_uint

lib.ElMatrixDestroy_i.argtypes = [c_void_p]
lib.ElMatrixDestroy_i.restype = c_uint
lib.ElMatrixDestroy_s.argtypes = [c_void_p]
lib.ElMatrixDestroy_s.restype = c_uint
lib.ElMatrixDestroy_d.argtypes = [c_void_p]
lib.ElMatrixDestroy_d.restype = c_uint
lib.ElMatrixDestroy_c.argtypes = [c_void_p]
lib.ElMatrixDestroy_c.restype = c_uint
lib.ElMatrixDestroy_z.argtypes = [c_void_p]
lib.ElMatrixDestroy_z.restype = c_uint

lib.ElMatrixResize_i.argtypes = [c_void_p,iType,iType]
lib.ElMatrixResize_i.restype = c_uint
lib.ElMatrixResize_s.argtypes = [c_void_p,iType,iType]
lib.ElMatrixResize_s.restype = c_uint
lib.ElMatrixResize_d.argtypes = [c_void_p,iType,iType]
lib.ElMatrixResize_d.restype = c_uint
lib.ElMatrixResize_c.argtypes = [c_void_p,iType,iType]
lib.ElMatrixResize_c.restype = c_uint
lib.ElMatrixResize_z.argtypes = [c_void_p,iType,iType]
lib.ElMatrixResize_z.restype = c_uint

lib.ElMatrixResizeWithLDim_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElMatrixResizeWithLDim_i.restype = c_uint
lib.ElMatrixResizeWithLDim_s.argtypes = [c_void_p,iType,iType,iType]
lib.ElMatrixResizeWithLDim_s.restype = c_uint
lib.ElMatrixResizeWithLDim_d.argtypes = [c_void_p,iType,iType,iType]
lib.ElMatrixResizeWithLDim_d.restype = c_uint
lib.ElMatrixResizeWithLDim_c.argtypes = [c_void_p,iType,iType,iType]
lib.ElMatrixResizeWithLDim_c.restype = c_uint
lib.ElMatrixResizeWithLDim_z.argtypes = [c_void_p,iType,iType,iType]
lib.ElMatrixResizeWithLDim_z.restype = c_uint

lib.ElMatrixEmpty_i.argtypes = [c_void_p]
lib.ElMatrixEmpty_i.restype = c_uint
lib.ElMatrixEmpty_s.argtypes = [c_void_p]
lib.ElMatrixEmpty_s.restype = c_uint
lib.ElMatrixEmpty_d.argtypes = [c_void_p]
lib.ElMatrixEmpty_d.restype = c_uint
lib.ElMatrixEmpty_c.argtypes = [c_void_p]
lib.ElMatrixEmpty_c.restype = c_uint
lib.ElMatrixEmpty_z.argtypes = [c_void_p]
lib.ElMatrixEmpty_z.restype = c_uint

lib.ElMatrixAttach_i.argtypes = [c_void_p,iType,iType,POINTER(iType),iType]
lib.ElMatrixAttach_i.restype = c_uint
lib.ElMatrixAttach_s.argtypes = [c_void_p,iType,iType,POINTER(sType),iType]
lib.ElMatrixAttach_s.restype = c_uint
lib.ElMatrixAttach_d.argtypes = [c_void_p,iType,iType,POINTER(dType),iType]
lib.ElMatrixAttach_d.restype = c_uint
lib.ElMatrixAttach_c.argtypes = [c_void_p,iType,iType,POINTER(cType),iType]
lib.ElMatrixAttach_c.restype = c_uint
lib.ElMatrixAttach_z.argtypes = [c_void_p,iType,iType,POINTER(zType),iType]
lib.ElMatrixAttach_z.restype = c_uint

lib.ElMatrixLockedAttach_i.argtypes = \
  [c_void_p,iType,iType,POINTER(iType),iType]
lib.ElMatrixLockedAttach_i.restype = c_uint
lib.ElMatrixLockedAttach_s.argtypes = \
  [c_void_p,iType,iType,POINTER(sType),iType]
lib.ElMatrixLockedAttach_s.restype = c_uint
lib.ElMatrixLockedAttach_d.argtypes = \
  [c_void_p,iType,iType,POINTER(dType),iType]
lib.ElMatrixLockedAttach_d.restype = c_uint
lib.ElMatrixLockedAttach_c.argtypes = \
  [c_void_p,iType,iType,POINTER(cType),iType]
lib.ElMatrixLockedAttach_c.restype = c_uint
lib.ElMatrixLockedAttach_z.argtypes = \
  [c_void_p,iType,iType,POINTER(zType),iType]
lib.ElMatrixLockedAttach_z.restype = c_uint

lib.ElMatrixControl_i.argtypes = [c_void_p,iType,iType,POINTER(iType),iType]
lib.ElMatrixControl_i.restype = c_uint
lib.ElMatrixControl_s.argtypes = [c_void_p,iType,iType,POINTER(sType),iType]
lib.ElMatrixControl_s.restype = c_uint
lib.ElMatrixControl_d.argtypes = [c_void_p,iType,iType,POINTER(dType),iType]
lib.ElMatrixControl_d.restype = c_uint
lib.ElMatrixControl_c.argtypes = [c_void_p,iType,iType,POINTER(cType),iType]
lib.ElMatrixControl_c.restype = c_uint
lib.ElMatrixControl_z.argtypes = [c_void_p,iType,iType,POINTER(zType),iType]
lib.ElMatrixControl_z.restype = c_uint

lib.ElMatrixHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixHeight_i.restype = c_uint
lib.ElMatrixHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixHeight_s.restype = c_uint
lib.ElMatrixHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixHeight_d.restype = c_uint
lib.ElMatrixHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixHeight_c.restype = c_uint
lib.ElMatrixHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixHeight_z.restype = c_uint

lib.ElMatrixWidth_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixWidth_i.restype = c_uint
lib.ElMatrixWidth_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixWidth_s.restype = c_uint
lib.ElMatrixWidth_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixWidth_d.restype = c_uint
lib.ElMatrixWidth_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixWidth_c.restype = c_uint
lib.ElMatrixWidth_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixWidth_z.restype = c_uint

lib.ElMatrixLDim_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixLDim_i.restype = c_uint
lib.ElMatrixLDim_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixLDim_s.restype = c_uint
lib.ElMatrixLDim_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixLDim_d.restype = c_uint
lib.ElMatrixLDim_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixLDim_c.restype = c_uint
lib.ElMatrixLDim_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixLDim_z.restype = c_uint

lib.ElMatrixMemorySize_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixMemorySize_i.restype = c_uint
lib.ElMatrixMemorySize_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixMemorySize_s.restype = c_uint
lib.ElMatrixMemorySize_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixMemorySize_d.restype = c_uint
lib.ElMatrixMemorySize_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixMemorySize_c.restype = c_uint
lib.ElMatrixMemorySize_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElMatrixMemorySize_z.restype = c_uint

lib.ElMatrixDiagonalLength_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElMatrixDiagonalLength_i.restype = c_uint
lib.ElMatrixDiagonalLength_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElMatrixDiagonalLength_s.restype = c_uint
lib.ElMatrixDiagonalLength_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElMatrixDiagonalLength_d.restype = c_uint
lib.ElMatrixDiagonalLength_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElMatrixDiagonalLength_c.restype = c_uint
lib.ElMatrixDiagonalLength_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElMatrixDiagonalLength_z.restype = c_uint

lib.ElMatrixViewing_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixViewing_i.restype = c_uint
lib.ElMatrixViewing_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixViewing_s.restype = c_uint
lib.ElMatrixViewing_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixViewing_d.restype = c_uint
lib.ElMatrixViewing_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixViewing_c.restype = c_uint
lib.ElMatrixViewing_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixViewing_z.restype = c_uint

lib.ElMatrixFixedSize_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixFixedSize_i.restype = c_uint
lib.ElMatrixFixedSize_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixFixedSize_s.restype = c_uint
lib.ElMatrixFixedSize_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixFixedSize_d.restype = c_uint
lib.ElMatrixFixedSize_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixFixedSize_c.restype = c_uint
lib.ElMatrixFixedSize_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixFixedSize_z.restype = c_uint

lib.ElMatrixLocked_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixLocked_i.restype = c_uint
lib.ElMatrixLocked_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixLocked_s.restype = c_uint
lib.ElMatrixLocked_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixLocked_d.restype = c_uint
lib.ElMatrixLocked_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixLocked_c.restype = c_uint
lib.ElMatrixLocked_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElMatrixLocked_z.restype = c_uint

lib.ElMatrixBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElMatrixBuffer_i.restype = c_uint
lib.ElMatrixBuffer_s.argtypes = [c_void_p,POINTER(POINTER(sType))]
lib.ElMatrixBuffer_s.restype = c_uint
lib.ElMatrixBuffer_d.argtypes = [c_void_p,POINTER(POINTER(dType))]
lib.ElMatrixBuffer_d.restype = c_uint
lib.ElMatrixBuffer_c.argtypes = [c_void_p,POINTER(POINTER(cType))]
lib.ElMatrixBuffer_c.restype = c_uint
lib.ElMatrixBuffer_z.argtypes = [c_void_p,POINTER(POINTER(zType))]
lib.ElMatrixBuffer_z.restype = c_uint

lib.ElMatrixLockedBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElMatrixLockedBuffer_i.restype = c_uint
lib.ElMatrixLockedBuffer_s.argtypes = [c_void_p,POINTER(POINTER(sType))]
lib.ElMatrixLockedBuffer_s.restype = c_uint
lib.ElMatrixLockedBuffer_d.argtypes = [c_void_p,POINTER(POINTER(dType))]
lib.ElMatrixLockedBuffer_d.restype = c_uint
lib.ElMatrixLockedBuffer_c.argtypes = [c_void_p,POINTER(POINTER(cType))]
lib.ElMatrixLockedBuffer_c.restype = c_uint
lib.ElMatrixLockedBuffer_z.argtypes = [c_void_p,POINTER(POINTER(zType))]
lib.ElMatrixLockedBuffer_z.restype = c_uint

lib.ElMatrixGet_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElMatrixGet_i.restype = c_uint
lib.ElMatrixGet_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElMatrixGet_s.restype = c_uint
lib.ElMatrixGet_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElMatrixGet_d.restype = c_uint
lib.ElMatrixGet_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
lib.ElMatrixGet_c.restype = c_uint
lib.ElMatrixGet_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
lib.ElMatrixGet_z.restype = c_uint

lib.ElMatrixGetRealPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElMatrixGetRealPart_c.restype = c_uint
lib.ElMatrixGetRealPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElMatrixGetRealPart_z.restype = c_uint

lib.ElMatrixGetImagPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElMatrixGetImagPart_c.restype = c_uint
lib.ElMatrixGetImagPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElMatrixGetImagPart_z.restype = c_uint

lib.ElMatrixSet_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElMatrixSet_i.restype = c_uint
lib.ElMatrixSet_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElMatrixSet_s.restype = c_uint
lib.ElMatrixSet_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElMatrixSet_d.restype = c_uint
lib.ElMatrixSet_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElMatrixSet_c.restype = c_uint
lib.ElMatrixSet_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElMatrixSet_z.restype = c_uint

lib.ElMatrixSetRealPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElMatrixSetRealPart_c.restype = c_uint
lib.ElMatrixSetRealPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElMatrixSetRealPart_z.restype = c_uint

lib.ElMatrixSetImagPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElMatrixSetImagPart_c.restype = c_uint
lib.ElMatrixSetImagPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElMatrixSetImagPart_z.restype = c_uint

lib.ElMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElMatrixUpdate_i.restype = c_uint
lib.ElMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElMatrixUpdate_s.restype = c_uint
lib.ElMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElMatrixUpdate_d.restype = c_uint
lib.ElMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElMatrixUpdate_c.restype = c_uint
lib.ElMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElMatrixUpdate_z.restype = c_uint

lib.ElMatrixUpdateRealPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElMatrixUpdateRealPart_c.restype = c_uint
lib.ElMatrixUpdateRealPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElMatrixUpdateRealPart_z.restype = c_uint

lib.ElMatrixUpdateImagPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElMatrixUpdateImagPart_c.restype = c_uint
lib.ElMatrixUpdateImagPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElMatrixUpdateImagPart_z.restype = c_uint

lib.ElMatrixMakeReal_c.argtypes = [c_void_p,iType,iType]
lib.ElMatrixMakeReal_c.restype = c_uint
lib.ElMatrixMakeReal_z.argtypes = [c_void_p,iType,iType]
lib.ElMatrixMakeReal_z.restype = c_uint

lib.ElMatrixConjugate_c.argtypes = [c_void_p,iType,iType]
lib.ElMatrixConjugate_c.restype = c_uint
lib.ElMatrixConjugate_z.argtypes = [c_void_p,iType,iType]
lib.ElMatrixConjugate_z.restype = c_uint

lib.ElView_i.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElView_i.restype = c_uint
lib.ElView_s.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElView_s.restype = c_uint
lib.ElView_d.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElView_d.restype = c_uint
lib.ElView_c.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElView_c.restype = c_uint
lib.ElView_z.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElView_z.restype = c_uint

lib.ElLockedView_i.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedView_i.restype = c_uint
lib.ElLockedView_s.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedView_s.restype = c_uint
lib.ElLockedView_d.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedView_d.restype = c_uint
lib.ElLockedView_c.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedView_c.restype = c_uint
lib.ElLockedView_z.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedView_z.restype = c_uint

class Matrix(object):
  def __init__(self,tag=dTag,create=True):
    self.obj = c_void_p()
    CheckTag(tag)
    self.tag = tag
    if create:
      args = [pointer(self.obj)]
      if   tag == iTag: lib.ElMatrixCreate_i(*args)
      elif tag == sTag: lib.ElMatrixCreate_s(*args)
      elif tag == dTag: lib.ElMatrixCreate_d(*args)
      elif tag == cTag: lib.ElMatrixCreate_c(*args)
      elif tag == zTag: lib.ElMatrixCreate_z(*args)
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElMatrixDestroy_i(*args)
    elif self.tag == sTag: lib.ElMatrixDestroy_s(*args)
    elif self.tag == dTag: lib.ElMatrixDestroy_d(*args)
    elif self.tag == cTag: lib.ElMatrixDestroy_c(*args)
    elif self.tag == zTag: lib.ElMatrixDestroy_z(*args)
  def SetType(self,newTag=dTag):
    self.Destroy()
    CheckTag(newTag)
    args = [pointer(self.obj)]
    if   newTag == iTag: lib.ElMatrixCreate_i(*args)
    elif newTag == sTag: lib.ElMatrixCreate_s(*args)
    elif newTag == dTag: lib.ElMatrixCreate_d(*args)
    elif newTag == cTag: lib.ElMatrixCreate_c(*args)
    elif newTag == zTag: lib.ElMatrixCreate_z(*args)
    self.tag = newTag
  def Resize(self,m,n):
    args = [self.obj,m,n]
    if   self.tag == iTag: lib.ElMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElMatrixResize_z(*args)
  def ResizeWithLDim(self,m,n,ldim):
    args = [self.obj,m,n,ldim]
    if   self.tag == iTag: lib.ElMatrixResizeWithLDim_i(*args)
    elif self.tag == sTag: lib.ElMatrixResizeWithLDim_s(*args)
    elif self.tag == dTag: lib.ElMatrixResizeWithLDim_d(*args)
    elif self.tag == cTag: lib.ElMatrixResizeWithLDim_c(*args)
    elif self.tag == zTag: lib.ElMatrixResizeWithLDim_z(*args)
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElMatrixEmpty_z(*args)
  def Attach(self,m,n,buf,ldim):
    args = [self.obj,m,n,buf,ldim]
    if   self.tag == iTag: lib.ElMatrixAttach_i(*args)
    elif self.tag == sTag: lib.ElMatrixAttach_s(*args)
    elif self.tag == dTag: lib.ElMatrixAttach_d(*args)
    elif self.tag == cTag: lib.ElMatrixAttach_c(*args)
    elif self.tag == zTag: lib.ElMatrixAttach_z(*args)
  def LockedAttach(self,m,n,buf,ldim):
    args = [self.obj,m,n,buf,ldim]
    if   self.tag == iTag: lib.ElMatrixLockedAttach_i(*args)
    elif self.tag == sTag: lib.ElMatrixLockedAttach_s(*args)
    elif self.tag == dTag: lib.ElMatrixLockedAttach_d(*args)
    elif self.tag == cTag: lib.ElMatrixLockedAttach_c(*args)
    elif self.tag == zTag: lib.ElMatrixLockedAttach_z(*args)
  def Control(self,m,n,buf,ldim):
    args = [self.obj,m,n,buf,ldim]
    if   self.tag == iTag: lib.ElMatrixControl_i(*args)
    elif self.tag == sTag: lib.ElMatrixControl_s(*args)
    elif self.tag == dTag: lib.ElMatrixControl_d(*args)
    elif self.tag == cTag: lib.ElMatrixControl_c(*args)
    elif self.tag == zTag: lib.ElMatrixControl_z(*args)
  def Height(self):
    m = iType()
    args = [self.obj,pointer(m)]
    if   self.tag == iTag: lib.ElMatrixHeight_i(*args)
    elif self.tag == sTag: lib.ElMatrixHeight_s(*args)
    elif self.tag == dTag: lib.ElMatrixHeight_d(*args)
    elif self.tag == cTag: lib.ElMatrixHeight_c(*args)
    elif self.tag == zTag: lib.ElMatrixHeight_z(*args)
    return m.value
  def Width(self):
    n = iType()
    args = [self.obj,pointer(n)]
    if   self.tag == iTag: lib.ElMatrixWidth_i(*args)
    elif self.tag == sTag: lib.ElMatrixWidth_s(*args)
    elif self.tag == dTag: lib.ElMatrixWidth_d(*args)
    elif self.tag == cTag: lib.ElMatrixWidth_c(*args)
    elif self.tag == zTag: lib.ElMatrixWidth_z(*args)
    return n.value
  def LDim(self):
    ldim = iType()
    args = [self.obj,pointer(ldim)]
    if   self.tag == iTag: lib.ElMatrixLDim_i(*args)
    elif self.tag == sTag: lib.ElMatrixLDim_s(*args)
    elif self.tag == dTag: lib.ElMatrixLDim_d(*args)
    elif self.tag == cTag: lib.ElMatrixLDim_c(*args)
    elif self.tag == zTag: lib.ElMatrixLDim_z(*args)
    return ldim.value
  def MemorySize(self):
    size = iType()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElMatrixMemorySize_i(*args)
    elif self.tag == sTag: lib.ElMatrixMemorySize_s(*args)
    elif self.tag == dTag: lib.ElMatrixMemorySize_d(*args)
    elif self.tag == cTag: lib.ElMatrixMemorySize_c(*args)
    elif self.tag == zTag: lib.ElMatrixMemorySize_z(*args)
    return size.value
  def DiagonalLength(self,offset=0):
    length = iType()
    args = [self.obj,offset,pointer(length)]
    if   self.tag == iTag: lib.ElMatrixDiagonalLength_i(*args)
    if   self.tag == sTag: lib.ElMatrixDiagonalLength_s(*args)
    if   self.tag == dTag: lib.ElMatrixDiagonalLength_d(*args)
    if   self.tag == cTag: lib.ElMatrixDiagonalLength_c(*args)
    if   self.tag == zTag: lib.ElMatrixDiagonalLength_z(*args)
    return length.value
  def Viewing(self):
    viewing = bType()
    args = [self.obj,pointer(viewing)]
    if   self.tag == iTag: lib.ElMatrixViewing_i(*args)
    elif self.tag == sTag: lib.ElMatrixViewing_s(*args)
    elif self.tag == dTag: lib.ElMatrixViewing_d(*args)
    elif self.tag == cTag: lib.ElMatrixViewing_c(*args)
    elif self.tag == zTag: lib.ElMatrixViewing_z(*args)
    return viewing.value
  def FixedSize(self):
    fixed = bType()
    args = [self.obj,pointer(fixed)]
    if   self.tag == iTag: lib.ElMatrixFixedSize_i(*args)
    elif self.tag == sTag: lib.ElMatrixFixedSize_s(*args)
    elif self.tag == dTag: lib.ElMatrixFixedSize_d(*args)
    elif self.tag == cTag: lib.ElMatrixFixedSize_c(*args)
    elif self.tag == zTag: lib.ElMatrixFixedSize_z(*args)
    return fixed.value
  def Locked(self):
    locked = bType()
    args = [self.obj,pointer(locked)]
    if   self.tag == iTag: lib.ElMatrixLocked_i(*args)
    elif self.tag == sTag: lib.ElMatrixLocked_s(*args)
    elif self.tag == dTag: lib.ElMatrixLocked_d(*args)
    elif self.tag == cTag: lib.ElMatrixLocked_c(*args)
    elif self.tag == zTag: lib.ElMatrixLocked_z(*args)
    return locked.value
  def Buffer(self):
    buf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(buf)]
    if   self.tag == iTag: lib.ElMatrixBuffer_i(*args) 
    elif self.tag == sTag: lib.ElMatrixBuffer_s(*args)
    elif self.tag == dTag: lib.ElMatrixBuffer_d(*args)
    elif self.tag == cTag: lib.ElMatrixBuffer_c(*args)
    elif self.tag == zTag: lib.ElMatrixBuffer_z(*args)
    return buf
  def LockedBuffer(self):
    buf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(buf)]
    if   self.tag == iTag: lib.ElMatrixLockedBuffer_i(*args) 
    elif self.tag == sTag: lib.ElMatrixLockedBuffer_s(*args)
    elif self.tag == dTag: lib.ElMatrixLockedBuffer_d(*args)
    elif self.tag == cTag: lib.ElMatrixLockedBuffer_c(*args)
    elif self.tag == zTag: lib.ElMatrixLockedBuffer_z(*args)
    return buf
  def Get(self,i,j):
    value = TagToType(self.tag)()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElMatrixGet_i(*args)
    elif self.tag == sTag: lib.ElMatrixGet_s(*args)
    elif self.tag == dTag: lib.ElMatrixGet_d(*args)
    elif self.tag == cTag: lib.ElMatrixGet_c(*args)
    elif self.tag == zTag: lib.ElMatrixGet_z(*args)
    return value.value
  def GetRealPart(self,i,j):
    if self.tag == cTag:
      value = sType()
      lib.ElMatrixGetRealPart_c(self.obj,i,j,pointer(value))
      return value.value
    elif self.tag == zTag:
      value = dType()
      lib.ElMatrixGetRealPart_z(self.obj,i,j,pointer(value))
      return value.value
    else: return Get(i,j)
  def GetImagPart(self,i,j):
    if   self.tag == iTag: return iType(0).value
    elif self.tag == sTag: return sType(0).value
    elif self.tag == dTag: return dType(0).value
    elif self.tag == cTag:
      value = c_float()
      lib.ElMatrixGetImagPart_c(self.obj,i,j,pointer(value))
      return value.value
    elif self.tag == zTag:
      value = c_double()
      lib.ElMatrixGetImagPart_z(self.obj,i,j,pointer(value))
      return value.value
  def Set(self,i,j,valuePre):
    value = TagToType(self.tag)(valuePre)
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElMatrixSet_i(*args)
    elif self.tag == sTag: lib.ElMatrixSet_s(*args)
    elif self.tag == dTag: lib.ElMatrixSet_d(*args)
    elif self.tag == cTag: lib.ElMatrixSet_c(*args)
    elif self.tag == zTag: lib.ElMatrixSet_z(*args)
  def SetRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetRealPart_z(self.obj,i,j,dType(value))
    else: self.Set(i,j,value)
  def SetImagPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetImagPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetImagPart_z(self.obj,i,j,dType(value))
    else: raise Exception("Datatype does not have an imaginary component")
  def Update(self,i,j,valuePre):
    value = TagToType(self.tag)(valuePre)
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElMatrixUpdate_z(*args)
  def UpdateRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixUpdateRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixUpdateRealPart_z(self.obj,i,j,dType(value))
    else: self.Update(i,j,value)
  def UpdateImagPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixUpdateImagPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixUpdateImagPart_z(self.obj,i,j,dType(value))
    else: raise Exception("Datatype does not have an imaginary component")
  def MakeReal(self,i,j):
    if   self.tag == cTag: lib.ElMatrixMakeReal_c(self.obj,i,j)
    elif self.tag == zTag: lib.ElMatrixMakeReal_z(self.obj,i,j) 
  def Conjugate(self,i,j):
    if   self.tag == cTag: lib.ElMatrixConjugate_c(self.obj,i,j)
    elif self.tag == zTag: lib.ElMatrixConjugate_z(self.obj,i,j)
  def ToNumPy(self):
    m = self.Height()
    n = self.Width()
    ldim = self.LDim()
    locked = self.Locked()
    if   self.tag == iTag:
      # TODO: Switch to 64-bit based upon Elemental's configuration
      entrySize = 4
      bufSize = entrySize*ldim*n
      if locked: buf = buffer_from_memory(self.LockedBuffer(),bufSize)
      else:      buf = buffer_from_memory_RW(self.Buffer(),bufSize)
      return np.ndarray(shape=(m,n),strides=(entrySize,ldim*entrySize),
                        buffer=buf,dtype=np.int32)
    elif self.tag == sTag:
      entrySize = 4
      bufSize = entrySize*ldim*n
      if locked: buf = buffer_from_memory(self.LockedBuffer(),bufSize)
      else:      buf = buffer_from_memory_RW(self.Buffer(),bufSize)
      return np.ndarray(shape=(m,n),strides=(entrySize,ldim*entrySize),
                        buffer=buf,dtype=np.float32)
    elif self.tag == dTag:
      entrySize = 8
      bufSize = entrySize*ldim*n
      if locked: buf = buffer_from_memory(self.LockedBuffer(),bufSize)
      else:      buf = buffer_from_memory_RW(self.Buffer(),bufSize)
      return np.ndarray(shape=(m,n),strides=(entrySize,ldim*entrySize),
                        buffer=buf,dtype=np.float64)
    elif self.tag == cTag: 
      entrySize = 8
      bufSize = entrySize*ldim*n
      if locked: buf = buffer_from_memory(self.LockedBuffer(),bufSize)
      else:      buf = buffer_from_memory_RW(self.Buffer(),bufSize)
      return np.ndarray(shape=(m,n),strides=(entrySize,ldim*entrySize),
                        buffer=buf,dtype=np.complex64)
    elif self.tag == zTag:
      entrySize = 16
      bufSize = entrySize*ldim*n
      if locked: buf = buffer_from_memory(self.LockedBuffer(),bufSize)
      else:      buf = buffer_from_memory_RW(self.Buffer(),bufSize)
      return np.ndarray(shape=(m,n),strides=(entrySize,ldim*entrySize),
                        buffer=buf,dtype=np.complex128)
  def __getitem__(self,indTup):
    iInd, jInd = indTup
    iRan = IndexRange(iInd)
    jRan = IndexRange(jInd)
    ASub = Matrix(self.tag)
    args = [ASub.obj,self.obj,iRan,jRan]
    if self.Locked():
      if   self.tag == iTag: lib.ElLockedView_i(*args)
      elif self.tag == sTag: lib.ElLockedView_s(*args)
      elif self.tag == dTag: lib.ElLockedView_d(*args)
      elif self.tag == cTag: lib.ElLockedView_c(*args)
      elif self.tag == zTag: lib.ElLockedView_z(*args)
      else: raise Exception('Unsupported datatype')
    else:
      if   self.tag == iTag: lib.ElView_i(*args)
      elif self.tag == sTag: lib.ElView_s(*args)
      elif self.tag == dTag: lib.ElView_d(*args)
      elif self.tag == cTag: lib.ElView_c(*args)
      elif self.tag == zTag: lib.ElView_z(*args)
      else: raise Exception('Unsupported datatype')
    return ASub
