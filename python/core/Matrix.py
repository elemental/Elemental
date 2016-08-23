#
#  Copyright (c) 2009-2016, Jack Poulson
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

class Matrix(object):
  # Create an instance
  # ------------------
  lib.ElMatrixCreate_i.argtypes = \
  lib.ElMatrixCreate_s.argtypes = \
  lib.ElMatrixCreate_d.argtypes = \
  lib.ElMatrixCreate_c.argtypes = \
  lib.ElMatrixCreate_z.argtypes = \
    [POINTER(c_void_p)]
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

  # Destroy an instance
  # -------------------
  lib.ElMatrixDestroy_i.argtypes = \
  lib.ElMatrixDestroy_s.argtypes = \
  lib.ElMatrixDestroy_d.argtypes = \
  lib.ElMatrixDestroy_c.argtypes = \
  lib.ElMatrixDestroy_z.argtypes = \
    [c_void_p]
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

  # HERE
  lib.ElMatrixResize_i.argtypes = \
  lib.ElMatrixResize_s.argtypes = \
  lib.ElMatrixResize_d.argtypes = \
  lib.ElMatrixResize_c.argtypes = \
  lib.ElMatrixResize_z.argtypes = \
    [c_void_p,iType,iType]
  def Resize(self,m,n):
    args = [self.obj,m,n]
    if   self.tag == iTag: lib.ElMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElMatrixResize_z(*args)

  lib.ElMatrixResizeWithLDim_i.argtypes = \
  lib.ElMatrixResizeWithLDim_s.argtypes = \
  lib.ElMatrixResizeWithLDim_d.argtypes = \
  lib.ElMatrixResizeWithLDim_c.argtypes = \
  lib.ElMatrixResizeWithLDim_z.argtypes = \
    [c_void_p,iType,iType,iType]
  def ResizeWithLDim(self,m,n,ldim):
    args = [self.obj,m,n,ldim]
    if   self.tag == iTag: lib.ElMatrixResizeWithLDim_i(*args)
    elif self.tag == sTag: lib.ElMatrixResizeWithLDim_s(*args)
    elif self.tag == dTag: lib.ElMatrixResizeWithLDim_d(*args)
    elif self.tag == cTag: lib.ElMatrixResizeWithLDim_c(*args)
    elif self.tag == zTag: lib.ElMatrixResizeWithLDim_z(*args)

  # Set the matrix back to its default (0 x 0) state
  # ------------------------------------------------
  lib.ElMatrixEmpty_i.argtypes = \
  lib.ElMatrixEmpty_s.argtypes = \
  lib.ElMatrixEmpty_d.argtypes = \
  lib.ElMatrixEmpty_c.argtypes = \
  lib.ElMatrixEmpty_z.argtypes = \
    [c_void_p]
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElMatrixEmpty_z(*args)

  # Attach a buffer to the matrix
  # -----------------------------
  lib.ElMatrixAttach_i.argtypes = \
  lib.ElMatrixLockedAttach_i.argtypes = \
    [c_void_p,iType,iType,POINTER(iType),iType]
  lib.ElMatrixAttach_s.argtypes = \
  lib.ElMatrixLockedAttach_s.argtypes = \
    [c_void_p,iType,iType,POINTER(sType),iType]
  lib.ElMatrixAttach_d.argtypes = \
  lib.ElMatrixLockedAttach_d.argtypes = \
    [c_void_p,iType,iType,POINTER(dType),iType]
  lib.ElMatrixAttach_c.argtypes = \
  lib.ElMatrixLockedAttach_c.argtypes = \
    [c_void_p,iType,iType,POINTER(cType),iType]
  lib.ElMatrixAttach_z.argtypes = \
  lib.ElMatrixLockedAttach_z.argtypes = \
    [c_void_p,iType,iType,POINTER(zType),iType]
  def Attach(self,m,n,buf,ldim,locked=False):
    args = [self.obj,m,n,buf,ldim]
    if locked:
      if   self.tag == iTag: lib.ElMatrixLockedAttach_i(*args)
      elif self.tag == sTag: lib.ElMatrixLockedAttach_s(*args)
      elif self.tag == dTag: lib.ElMatrixLockedAttach_d(*args)
      elif self.tag == cTag: lib.ElMatrixLockedAttach_c(*args)
      elif self.tag == zTag: lib.ElMatrixLockedAttach_z(*args)
    else:
      if   self.tag == iTag: lib.ElMatrixAttach_i(*args)
      elif self.tag == sTag: lib.ElMatrixAttach_s(*args)
      elif self.tag == dTag: lib.ElMatrixAttach_d(*args)
      elif self.tag == cTag: lib.ElMatrixAttach_c(*args)
      elif self.tag == zTag: lib.ElMatrixAttach_z(*args)

  lib.ElMatrixControl_i.argtypes = [c_void_p,iType,iType,POINTER(iType),iType]
  lib.ElMatrixControl_s.argtypes = [c_void_p,iType,iType,POINTER(sType),iType]
  lib.ElMatrixControl_d.argtypes = [c_void_p,iType,iType,POINTER(dType),iType]
  lib.ElMatrixControl_c.argtypes = [c_void_p,iType,iType,POINTER(cType),iType]
  lib.ElMatrixControl_z.argtypes = [c_void_p,iType,iType,POINTER(zType),iType]
  def Control(self,m,n,buf,ldim):
    args = [self.obj,m,n,buf,ldim]
    if   self.tag == iTag: lib.ElMatrixControl_i(*args)
    elif self.tag == sTag: lib.ElMatrixControl_s(*args)
    elif self.tag == dTag: lib.ElMatrixControl_d(*args)
    elif self.tag == cTag: lib.ElMatrixControl_c(*args)
    elif self.tag == zTag: lib.ElMatrixControl_z(*args)

  lib.ElMatrixHeight_i.argtypes = \
  lib.ElMatrixHeight_s.argtypes = \
  lib.ElMatrixHeight_d.argtypes = \
  lib.ElMatrixHeight_c.argtypes = \
  lib.ElMatrixHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Height(self):
    m = iType()
    args = [self.obj,pointer(m)]
    if   self.tag == iTag: lib.ElMatrixHeight_i(*args)
    elif self.tag == sTag: lib.ElMatrixHeight_s(*args)
    elif self.tag == dTag: lib.ElMatrixHeight_d(*args)
    elif self.tag == cTag: lib.ElMatrixHeight_c(*args)
    elif self.tag == zTag: lib.ElMatrixHeight_z(*args)
    return m.value

  lib.ElMatrixWidth_i.argtypes = \
  lib.ElMatrixWidth_s.argtypes = \
  lib.ElMatrixWidth_d.argtypes = \
  lib.ElMatrixWidth_c.argtypes = \
  lib.ElMatrixWidth_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def Width(self):
    n = iType()
    args = [self.obj,pointer(n)]
    if   self.tag == iTag: lib.ElMatrixWidth_i(*args)
    elif self.tag == sTag: lib.ElMatrixWidth_s(*args)
    elif self.tag == dTag: lib.ElMatrixWidth_d(*args)
    elif self.tag == cTag: lib.ElMatrixWidth_c(*args)
    elif self.tag == zTag: lib.ElMatrixWidth_z(*args)
    return n.value

  lib.ElMatrixLDim_i.argtypes = \
  lib.ElMatrixLDim_s.argtypes = \
  lib.ElMatrixLDim_d.argtypes = \
  lib.ElMatrixLDim_c.argtypes = \
  lib.ElMatrixLDim_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def LDim(self):
    ldim = iType()
    args = [self.obj,pointer(ldim)]
    if   self.tag == iTag: lib.ElMatrixLDim_i(*args)
    elif self.tag == sTag: lib.ElMatrixLDim_s(*args)
    elif self.tag == dTag: lib.ElMatrixLDim_d(*args)
    elif self.tag == cTag: lib.ElMatrixLDim_c(*args)
    elif self.tag == zTag: lib.ElMatrixLDim_z(*args)
    return ldim.value

  lib.ElMatrixMemorySize_i.argtypes = \
  lib.ElMatrixMemorySize_s.argtypes = \
  lib.ElMatrixMemorySize_d.argtypes = \
  lib.ElMatrixMemorySize_c.argtypes = \
  lib.ElMatrixMemorySize_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def MemorySize(self):
    size = iType()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElMatrixMemorySize_i(*args)
    elif self.tag == sTag: lib.ElMatrixMemorySize_s(*args)
    elif self.tag == dTag: lib.ElMatrixMemorySize_d(*args)
    elif self.tag == cTag: lib.ElMatrixMemorySize_c(*args)
    elif self.tag == zTag: lib.ElMatrixMemorySize_z(*args)
    return size.value

  lib.ElMatrixDiagonalLength_i.argtypes = \
  lib.ElMatrixDiagonalLength_s.argtypes = \
  lib.ElMatrixDiagonalLength_d.argtypes = \
  lib.ElMatrixDiagonalLength_c.argtypes = \
  lib.ElMatrixDiagonalLength_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def DiagonalLength(self,offset=0):
    length = iType()
    args = [self.obj,offset,pointer(length)]
    if   self.tag == iTag: lib.ElMatrixDiagonalLength_i(*args)
    if   self.tag == sTag: lib.ElMatrixDiagonalLength_s(*args)
    if   self.tag == dTag: lib.ElMatrixDiagonalLength_d(*args)
    if   self.tag == cTag: lib.ElMatrixDiagonalLength_c(*args)
    if   self.tag == zTag: lib.ElMatrixDiagonalLength_z(*args)
    return length.value

  lib.ElMatrixViewing_i.argtypes = \
  lib.ElMatrixViewing_s.argtypes = \
  lib.ElMatrixViewing_d.argtypes = \
  lib.ElMatrixViewing_c.argtypes = \
  lib.ElMatrixViewing_z.argtypes = \
    [c_void_p,POINTER(bType)]
  def Viewing(self):
    viewing = bType()
    args = [self.obj,pointer(viewing)]
    if   self.tag == iTag: lib.ElMatrixViewing_i(*args)
    elif self.tag == sTag: lib.ElMatrixViewing_s(*args)
    elif self.tag == dTag: lib.ElMatrixViewing_d(*args)
    elif self.tag == cTag: lib.ElMatrixViewing_c(*args)
    elif self.tag == zTag: lib.ElMatrixViewing_z(*args)
    return viewing.value

  lib.ElMatrixFixedSize_i.argtypes = \
  lib.ElMatrixFixedSize_s.argtypes = \
  lib.ElMatrixFixedSize_d.argtypes = \
  lib.ElMatrixFixedSize_c.argtypes = \
  lib.ElMatrixFixedSize_z.argtypes = \
    [c_void_p,POINTER(bType)]
  def FixedSize(self):
    fixed = bType()
    args = [self.obj,pointer(fixed)]
    if   self.tag == iTag: lib.ElMatrixFixedSize_i(*args)
    elif self.tag == sTag: lib.ElMatrixFixedSize_s(*args)
    elif self.tag == dTag: lib.ElMatrixFixedSize_d(*args)
    elif self.tag == cTag: lib.ElMatrixFixedSize_c(*args)
    elif self.tag == zTag: lib.ElMatrixFixedSize_z(*args)
    return fixed.value

  lib.ElMatrixLocked_i.argtypes = \
  lib.ElMatrixLocked_s.argtypes = \
  lib.ElMatrixLocked_d.argtypes = \
  lib.ElMatrixLocked_c.argtypes = \
  lib.ElMatrixLocked_z.argtypes = \
    [c_void_p,POINTER(bType)]
  def Locked(self):
    locked = bType()
    args = [self.obj,pointer(locked)]
    if   self.tag == iTag: lib.ElMatrixLocked_i(*args)
    elif self.tag == sTag: lib.ElMatrixLocked_s(*args)
    elif self.tag == dTag: lib.ElMatrixLocked_d(*args)
    elif self.tag == cTag: lib.ElMatrixLocked_c(*args)
    elif self.tag == zTag: lib.ElMatrixLocked_z(*args)
    return locked.value

  lib.ElMatrixBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
  lib.ElMatrixBuffer_s.argtypes = [c_void_p,POINTER(POINTER(sType))]
  lib.ElMatrixBuffer_d.argtypes = [c_void_p,POINTER(POINTER(dType))]
  lib.ElMatrixBuffer_c.argtypes = [c_void_p,POINTER(POINTER(cType))]
  lib.ElMatrixBuffer_z.argtypes = [c_void_p,POINTER(POINTER(zType))]
  lib.ElMatrixLockedBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
  lib.ElMatrixLockedBuffer_s.argtypes = [c_void_p,POINTER(POINTER(sType))]
  lib.ElMatrixLockedBuffer_d.argtypes = [c_void_p,POINTER(POINTER(dType))]
  lib.ElMatrixLockedBuffer_c.argtypes = [c_void_p,POINTER(POINTER(cType))]
  lib.ElMatrixLockedBuffer_z.argtypes = [c_void_p,POINTER(POINTER(zType))]
  def Buffer(self,locked=False):
    buf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(buf)]
    if locked:
      if   self.tag == iTag: lib.ElMatrixLockedBuffer_i(*args) 
      elif self.tag == sTag: lib.ElMatrixLockedBuffer_s(*args)
      elif self.tag == dTag: lib.ElMatrixLockedBuffer_d(*args)
      elif self.tag == cTag: lib.ElMatrixLockedBuffer_c(*args)
      elif self.tag == zTag: lib.ElMatrixLockedBuffer_z(*args)
    else:
      if   self.tag == iTag: lib.ElMatrixBuffer_i(*args) 
      elif self.tag == sTag: lib.ElMatrixBuffer_s(*args)
      elif self.tag == dTag: lib.ElMatrixBuffer_d(*args)
      elif self.tag == cTag: lib.ElMatrixBuffer_c(*args)
      elif self.tag == zTag: lib.ElMatrixBuffer_z(*args)
    return buf

  lib.ElMatrixGet_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
  lib.ElMatrixGet_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElMatrixGet_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
  lib.ElMatrixGet_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
  lib.ElMatrixGet_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
  def Get(self,i,j):
    value = TagToType(self.tag)()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElMatrixGet_i(*args)
    elif self.tag == sTag: lib.ElMatrixGet_s(*args)
    elif self.tag == dTag: lib.ElMatrixGet_d(*args)
    elif self.tag == cTag: lib.ElMatrixGet_c(*args)
    elif self.tag == zTag: lib.ElMatrixGet_z(*args)
    return ScalarData(value)

  lib.ElMatrixGetRealPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElMatrixGetRealPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
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

  lib.ElMatrixGetImagPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElMatrixGetImagPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
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

  lib.ElMatrixSet_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElMatrixSet_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElMatrixSet_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElMatrixSet_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElMatrixSet_z.argtypes = [c_void_p,iType,iType,zType]
  def Set(self,i,j,valuePre):
    value = TagToType(self.tag)(valuePre)
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElMatrixSet_i(*args)
    elif self.tag == sTag: lib.ElMatrixSet_s(*args)
    elif self.tag == dTag: lib.ElMatrixSet_d(*args)
    elif self.tag == cTag: lib.ElMatrixSet_c(*args)
    elif self.tag == zTag: lib.ElMatrixSet_z(*args)

  lib.ElMatrixSetRealPart_c.argtypes = [c_void_p,iType,iType,sType]
  lib.ElMatrixSetRealPart_z.argtypes = [c_void_p,iType,iType,dType]
  def SetRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetRealPart_z(self.obj,i,j,dType(value))
    else: self.Set(i,j,value)

  lib.ElMatrixSetImagPart_c.argtypes = [c_void_p,iType,iType,sType]
  lib.ElMatrixSetImagPart_z.argtypes = [c_void_p,iType,iType,dType]
  def SetImagPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetImagPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetImagPart_z(self.obj,i,j,dType(value))
    else: raise Exception("Datatype does not have an imaginary component")

  lib.ElMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
  def Update(self,i,j,valuePre):
    value = TagToType(self.tag)(valuePre)
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElMatrixUpdate_z(*args)

  lib.ElMatrixUpdateRealPart_c.argtypes = [c_void_p,iType,iType,sType]
  lib.ElMatrixUpdateRealPart_z.argtypes = [c_void_p,iType,iType,dType]
  def UpdateRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixUpdateRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixUpdateRealPart_z(self.obj,i,j,dType(value))
    else: self.Update(i,j,value)

  lib.ElMatrixUpdateImagPart_c.argtypes = [c_void_p,iType,iType,sType]
  lib.ElMatrixUpdateImagPart_z.argtypes = [c_void_p,iType,iType,dType]
  def UpdateImagPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixUpdateImagPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixUpdateImagPart_z(self.obj,i,j,dType(value))
    else: raise Exception("Datatype does not have an imaginary component")

  lib.ElMatrixMakeReal_c.argtypes = \
  lib.ElMatrixMakeReal_z.argtypes = \
    [c_void_p,iType,iType]
  def MakeReal(self,i,j):
    if   self.tag == cTag: lib.ElMatrixMakeReal_c(self.obj,i,j)
    elif self.tag == zTag: lib.ElMatrixMakeReal_z(self.obj,i,j) 

  lib.ElMatrixConjugate_c.argtypes = \
  lib.ElMatrixConjugate_z.argtypes = \
    [c_void_p,iType,iType]
  def Conjugate(self,i,j):
    if   self.tag == cTag: lib.ElMatrixConjugate_c(self.obj,i,j)
    elif self.tag == zTag: lib.ElMatrixConjugate_z(self.obj,i,j)

  def ToNumPy(self):
    m = self.Height()
    n = self.Width()
    ldim = self.LDim()
    locked = self.Locked()
    entrySize = TagToSize(self.tag)
    npType = TagToNumpyType(self.tag)
    bufSize = entrySize*ldim*n
    if locked: buf = buffer_from_memory(self.LockedBuffer(),bufSize)
    else:      buf = buffer_from_memory_RW(self.Buffer(),bufSize)
    return np.ndarray(shape=(m,n),strides=(entrySize,ldim*entrySize),
                        buffer=buf,dtype=npType)

  lib.ElView_i.argtypes = \
  lib.ElView_s.argtypes = \
  lib.ElView_d.argtypes = \
  lib.ElView_c.argtypes = \
  lib.ElView_z.argtypes = \
  lib.ElLockedView_i.argtypes = \
  lib.ElLockedView_s.argtypes = \
  lib.ElLockedView_d.argtypes = \
  lib.ElLockedView_c.argtypes = \
  lib.ElLockedView_z.argtypes = \
    [c_void_p,c_void_p,IndexRange,IndexRange]
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
    ASub = Matrix(self.tag)
    args = [ASub.obj,self.obj,iRan,jRan]
    if self.Locked():
      if   self.tag == iTag: lib.ElLockedView_i(*args)
      elif self.tag == sTag: lib.ElLockedView_s(*args)
      elif self.tag == dTag: lib.ElLockedView_d(*args)
      elif self.tag == cTag: lib.ElLockedView_c(*args)
      elif self.tag == zTag: lib.ElLockedView_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElView_i(*args)
      elif self.tag == sTag: lib.ElView_s(*args)
      elif self.tag == dTag: lib.ElView_d(*args)
      elif self.tag == cTag: lib.ElView_c(*args)
      elif self.tag == zTag: lib.ElView_z(*args)
      else: DataExcept()
    return ASub
