#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import ctypes, numpy

from Grid import DefaultGrid

# (Abstract)DistMatrix
# ====================

class DistData(ctypes.Structure):
  _fields_ = [('colDist',ctypes.c_uint),
              ('rowDist',ctypes.c_uint), 
              ('colAlign',iType),
              ('rowAlign',iType),
              ('root',iType),
              ('grid',ctypes.c_void_p)]

class DistMatrix(object):
  def __init__(self,tag=dTag,colDist=MC,rowDist=MR,grid=DefaultGrid()):
    self.obj = ctypes.c_void_p()
    CheckTag(tag)
    colDistVal = ctypes.c_uint(colDist)
    rowDistVal = ctypes.c_uint(rowDist)
    if   tag == iTag: 
      lib.ElDistMatrixCreateSpecific_i \
      (colDistVal,rowDistVal,grid.obj,pointer(self.obj))
    elif tag == sTag: 
      lib.ElDistMatrixCreateSpecific_s \
      (colDistVal,rowDistVal,grid.obj,pointer(self.obj))
    elif tag == dTag: 
      lib.ElDistMatrixCreateSpecific_d \
      (colDistVal,rowDistVal,grid.obj,pointer(self.obj))
    elif tag == cTag: 
      lib.ElDistMatrixCreateSpecific_c \
      (colDistVal,rowDistVal,grid.obj,pointer(self.obj))
    elif tag == zTag: 
      lib.ElDistMatrixCreateSpecific_z \
      (colDistVal,rowDistVal,grid.obj,pointer(self.obj))
    self.tag = tag
  def Destroy(self):
    if   self.tag == iTag: lib.ElDistMatrixDestroy_i(self.obj)
    elif self.tag == sTag: lib.ElDistMatrixDestroy_s(self.obj)
    elif self.tag == dTag: lib.ElDistMatrixDestroy_d(self.obj)
    elif self.tag == cTag: lib.ElDistMatrixDestroy_c(self.obj)
    elif self.tag == zTag: lib.ElDistMatrixDestroy_z(self.obj)
  def Empty(self):
    if   self.tag == iTag: lib.ElDistMatrixEmpty_i(self.obj)
    elif self.tag == sTag: lib.ElDistMatrixEmpty_s(self.obj)
    elif self.tag == dTag: lib.ElDistMatrixEmpty_d(self.obj)
    elif self.tag == cTag: lib.ElDistMatrixEmpty_c(self.obj)
    elif self.tag == zTag: lib.ElDistMatrixEmpty_z(self.obj)
  def EmptyData(self):
    if   self.tag == iTag: lib.ElDistMatrixEmptyData_i(self.obj)
    elif self.tag == sTag: lib.ElDistMatrixEmptyData_s(self.obj)
    elif self.tag == dTag: lib.ElDistMatrixEmptyData_d(self.obj)
    elif self.tag == cTag: lib.ElDistMatrixEmptyData_c(self.obj)
    elif self.tag == zTag: lib.ElDistMatrixEmptyData_z(self.obj)
  def SetGrid(self,grid):
    if   self.tag == iTag: lib.ElDistMatrixSetGrid_i(self.obj,grid.obj)
    elif self.tag == sTag: lib.ElDistMatrixSetGrid_s(self.obj,grid.obj)
    elif self.tag == dTag: lib.ElDistMatrixSetGrid_d(self.obj,grid.obj)
    elif self.tag == cTag: lib.ElDistMatrixSetGrid_c(self.obj,grid.obj)
    elif self.tag == zTag: lib.ElDistMatrixSetGrid_z(self.obj,grid.obj)
  def Resize(self,m,n):
    if   self.tag == iTag: lib.ElDistMatrixResize_i(self.obj,m,n)
    elif self.tag == sTag: lib.ElDistMatrixResize_s(self.obj,m,n)
    elif self.tag == dTag: lib.ElDistMatrixResize_d(self.obj,m,n)
    elif self.tag == cTag: lib.ElDistMatrixResize_c(self.obj,m,n)
    elif self.tag == zTag: lib.ElDistMatrixResize_z(self.obj,m,n)
  def ResizeWithLDim(self,m,n,ldim):
    if   self.tag == iTag: lib.ElDistMatrixResizeWithLDim_i(self.obj,m,n,ldim)
    elif self.tag == sTag: lib.ElDistMatrixResizeWithLDim_s(self.obj,m,n,ldim)
    elif self.tag == dTag: lib.ElDistMatrixResizeWithLDim_d(self.obj,m,n,ldim)
    elif self.tag == cTag: lib.ElDistMatrixResizeWithLDim_c(self.obj,m,n,ldim)
    elif self.tag == zTag: lib.ElDistMatrixResizeWithLDim_z(self.obj,m,n,ldim)
  def MakeConsistent(self,incViewers):
    if   self.tag == iTag: lib.ElDistMatrixMakeConsistent_i(self.obj,incViewers)
    elif self.tag == sTag: lib.ElDistMatrixMakeConsistent_s(self.obj,incViewers)
    elif self.tag == dTag: lib.ElDistMatrixMakeConsistent_d(self.obj,incViewers)
    elif self.tag == cTag: lib.ElDistMatrixMakeConsistent_c(self.obj,incViewers)
    elif self.tag == zTag: lib.ElDistMatrixMakeConsistent_z(self.obj,incViewers)
  def MakeSizeConsistent(self,incViewers):
    if   self.tag == iTag:
      lib.ElDistMatrixMakeSizeConsistent_i(self.obj,incViewers)
    elif self.tag == sTag:
      lib.ElDistMatrixMakeSizeConsistent_s(self.obj,incViewers)
    elif self.tag == dTag:
      lib.ElDistMatrixMakeSizeConsistent_d(self.obj,incViewers)
    elif self.tag == cTag:
      lib.ElDistMatrixMakeSizeConsistent_c(self.obj,incViewers)
    elif self.tag == zTag:
      lib.ElDistMatrixMakeSizeConsistent_z(self.obj,incViewers)
  def Align(self,colAlign,rowAlign,constrain):
    if   self.tag == iTag: 
      lib.ElDistMatrixAlign_i(self.obj,colAlign,rowAlign,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlign_s(self.obj,colAlign,rowAlign,constrain)
    elif self.tag == dTag:
      lib.ElDistMatrixAlign_d(self.obj,colAlign,rowAlign,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlign_c(self.obj,colAlign,rowAlign,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlign_z(self.obj,colAlign,rowAlign,constrain)
  def AlignCols(self,colAlign,constrain):
    if   self.tag == iTag: 
      lib.ElDistMatrixAlignCols_i(self.obj,colAlign,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignCols_s(self.obj,colAlign,constrain)
    elif self.tag == dTag:
      lib.ElDistMatrixAlignCols_d(self.obj,colAlign,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignCols_c(self.obj,colAlign,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignCols_z(self.obj,colAlign,constrain)
  def AlignRows(self,rowAlign,constrain):
    if   self.tag == iTag: 
      lib.ElDistMatrixAlignRows_i(self.obj,rowAlign,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignRows_s(self.obj,rowAlign,constrain)
    elif self.tag == dTag:
      lib.ElDistMatrixAlignRows_d(self.obj,rowAlign,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignRows_c(self.obj,rowAlign,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignRows_z(self.obj,rowAlign,constrain)
  def FreeAlignments(self):
    if   self.tag == iTag: lib.ElDistMatrixFreeAlignments_i(self.obj)
    elif self.tag == sTag: lib.ElDistMatrixFreeAlignments_s(self.obj)
    elif self.tag == dTag: lib.ElDistMatrixFreeAlignments_d(self.obj)
    elif self.tag == cTag: lib.ElDistMatrixFreeAlignments_c(self.obj)
    elif self.tag == zTag: lib.ElDistMatrixFreeAlignments_z(self.obj)
  def SetRoot(self,root):
    if   self.tag == iTag: lib.ElDistMatrixSetRoot_i(self.obj,root)
    elif self.tag == sTag: lib.ElDistMatrixSetRoot_s(self.obj,root)
    elif self.tag == dTag: lib.ElDistMatrixSetRoot_d(self.obj,root)
    elif self.tag == cTag: lib.ElDistMatrixSetRoot_c(self.obj,root)
    elif self.tag == zTag: lib.ElDistMatrixSetRoot_z(self.obj,root)
  def AlignWith(self,distData,constrain):
    if   self.tag == iTag: 
      lib.ElDistMatrixAlignWith_i(self.obj,distData,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignWith_s(self.obj,distData,constrain)
    elif self.tag == dTag:
      lib.ElDistMatrixAlignWith_d(self.obj,distData,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignWith_c(self.obj,distData,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignWith_z(self.obj,distData,constrain)
  def AlignColsWith(self,distData,constrain):
    if   self.tag == iTag: 
      lib.ElDistMatrixAlignColsWith_i(self.obj,distData,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignColsWith_s(self.obj,distData,constrain)
    elif self.tag == dTag:
      lib.ElDistMatrixAlignColsWith_d(self.obj,distData,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignColsWith_c(self.obj,distData,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignColsWith_z(self.obj,distData,constrain)
  def AlignRowsWith(self,distData,constrain):
    if   self.tag == iTag: 
      lib.ElDistMatrixAlignRowsWith_i(self.obj,distData,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignRowsWith_s(self.obj,distData,constrain)
    elif self.tag == dTag:
      lib.ElDistMatrixAlignRowsWith_d(self.obj,distData,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignRowsWith_c(self.obj,distData,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignRowsWith_z(self.obj,distData,constrain)
  def AlignAndResize(self,colAlign,rowAlign,m,n,force,constrain):
    if   self.tag == iTag:
      lib.ElDistMatrixAlignAndResize_i \
      (self.obj,colAlign,rowAlign,m,n,force,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignAndResize_s \
      (self.obj,colAlign,rowAlign,m,n,force,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignAndResize_c \
      (self.obj,colAlign,rowAlign,m,n,force,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignAndResize_z \
      (self.obj,colAlign,rowAlign,m,n,force,constrain)
  def AlignColsAndResize(self,colAlign,m,n,force,constrain):
    if   self.tag == iTag:
      lib.ElDistMatrixAlignColsAndResize_i \
      (self.obj,colAlign,m,n,force,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignColsAndResize_s \
      (self.obj,colAlign,m,n,force,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignColsAndResize_c \
      (self.obj,colAlign,m,n,force,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignColsAndResize_z \
      (self.obj,colAlign,m,n,force,constrain)
  def AlignRowsAndResize(self,rowAlign,m,n,force,constrain):
    if   self.tag == iTag:
      lib.ElDistMatrixAlignRowsAndResize_i \
      (self.obj,rowAlign,m,n,force,constrain)
    elif self.tag == sTag:
      lib.ElDistMatrixAlignRowsAndResize_s \
      (self.obj,rowAlign,m,n,force,constrain)
    elif self.tag == cTag:
      lib.ElDistMatrixAlignRowsAndResize_c \
      (self.obj,rowAlign,m,n,force,constrain)
    elif self.tag == zTag:
      lib.ElDistMatrixAlignRowsAndResize_z \
      (self.obj,rowAlign,m,n,force,constrain)
  def Attach(self,m,n,grid,colAlign,rowAlign,buf,ldim,root):
    EnsureCompatibleBuffer(buf,self.tag)
    if   self.tag == iTag: 
      lib.ElDistMatrixAttach_i \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == sTag:
      lib.ElDistMatrixAttach_s \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == dTag:
      lib.ElDistMatrixAttach_d \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == cTag:
      lib.ElDistMatrixAttach_c \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == zTag:
      lib.ElDistMatrixAttach_z \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
  def LockedAttach(self,m,n,buf,ldim):
    EnsureCompatibleBuffer(buf,self.tag)
    if   self.tag == iTag: 
      lib.ElDistMatrixLockedAttach_i \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == sTag:
      lib.ElDistMatrixLockedAttach_s \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == dTag:
      lib.ElDistMatrixLockedAttach_d \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == cTag:
      lib.ElDistMatrixLockedAttach_c \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
    elif self.tag == zTag:
      lib.ElDistMatrixLockedAttach_z \
      (self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root)
  def Height(self):
    m = iType()
    if   self.tag == iTag: lib.ElDistMatrixHeight_i(self.obj,pointer(m))
    elif self.tag == sTag: lib.ElDistMatrixHeight_s(self.obj,pointer(m))
    elif self.tag == dTag: lib.ElDistMatrixHeight_d(self.obj,pointer(m))
    elif self.tag == cTag: lib.ElDistMatrixHeight_c(self.obj,pointer(m))
    elif self.tag == zTag: lib.ElDistMatrixHeight_z(self.obj,pointer(m))
    return m 
  def Width(self):
    n = iType()
    if   self.tag == iTag: lib.ElDistMatrixWidth_i(self.obj,pointer(n))
    elif self.tag == sTag: lib.ElDistMatrixWidth_s(self.obj,pointer(n))
    elif self.tag == dTag: lib.ElDistMatrixWidth_d(self.obj,pointer(n))
    elif self.tag == cTag: lib.ElDistMatrixWidth_c(self.obj,pointer(n))
    elif self.tag == zTag: lib.ElDistMatrixWidth_z(self.obj,pointer(n))
    return n
  def DiagonalLength(self,offset):
    length = iType()
    if   self.tag == iTag:
      lib.ElDistMatrixDiagonalLength_i(self.obj,pointer(length))
    elif self.tag == sTag:
      lib.ElDistMatrixDiagonalLength_s(self.obj,pointer(length))
    elif self.tag == dTag:
      lib.ElDistMatrixDiagonalLength_d(self.obj,pointer(length))
    elif self.tag == cTag:
      lib.ElDistMatrixDiagonalLength_c(self.obj,pointer(length))
    elif self.tag == zTag:
      lib.ElDistMatrixDiagonalLength_z(self.obj,pointer(length))
    return length
  def Viewing(self):
    viewing = bType()
    if   self.tag == iTag: lib.ElDistMatrixViewing_i(self.obj,pointer(viewing))
    elif self.tag == sTag: lib.ElDistMatrixViewing_s(self.obj,pointer(viewing))
    elif self.tag == dTag: lib.ElDistMatrixViewing_d(self.obj,pointer(viewing))
    elif self.tag == cTag: lib.ElDistMatrixViewing_c(self.obj,pointer(viewing))
    elif self.tag == zTag: lib.ElDistMatrixViewing_z(self.obj,pointer(viewing))
    return viewing
  def Locked(self):
    locked = bType()
    if   self.tag == iTag: lib.ElDistMatrixLocked_i(self.obj,pointer(locked))
    elif self.tag == sTag: lib.ElDistMatrixLocked_s(self.obj,pointer(locked))
    elif self.tag == dTag: lib.ElDistMatrixLocked_d(self.obj,pointer(locked))
    elif self.tag == cTag: lib.ElDistMatrixLocked_c(self.obj,pointer(locked))
    elif self.tag == zTag: lib.ElDistMatrixLocked_z(self.obj,pointer(locked))
    return locked
  def LocalHeight(self):
    mLoc = iType()
    if   self.tag == iTag: lib.ElDistMatrixLocalHeight_i(self.obj,pointer(mLoc))
    elif self.tag == sTag: lib.ElDistMatrixLocalHeight_s(self.obj,pointer(mLoc))
    elif self.tag == dTag: lib.ElDistMatrixLocalHeight_d(self.obj,pointer(mLoc))
    elif self.tag == cTag: lib.ElDistMatrixLocalHeight_c(self.obj,pointer(mLoc))
    elif self.tag == zTag: lib.ElDistMatrixLocalHeight_z(self.obj,pointer(mLoc))
    return mLoc
  def LocalWidth(self):
    nLoc = iType()
    if   self.tag == iTag: lib.ElDistMatrixLocalWidth_i(self.obj,pointer(nLoc))
    elif self.tag == sTag: lib.ElDistMatrixLocalWidth_s(self.obj,pointer(nLoc))
    elif self.tag == dTag: lib.ElDistMatrixLocalWidth_d(self.obj,pointer(nLoc))
    elif self.tag == cTag: lib.ElDistMatrixLocalWidth_c(self.obj,pointer(nLoc))
    elif self.tag == zTag: lib.ElDistMatrixLocalWidth_z(self.obj,pointer(nLoc))
    return nLoc
  def LDim(self):
    ldim = iType()
    if   self.tag == iTag: lib.ElDistMatrixLDim_i(self.obj,pointer(ldim))
    elif self.tag == sTag: lib.ElDistMatrixLDim_s(self.obj,pointer(ldim))
    elif self.tag == dTag: lib.ElDistMatrixLDim_d(self.obj,pointer(ldim))
    elif self.tag == cTag: lib.ElDistMatrixLDim_c(self.obj,pointer(ldim))
    elif self.tag == zTag: lib.ElDistMatrixLDim_z(self.obj,pointer(ldim))
    return ldim
  def Matrix(self):
    # TODO
    raise Exception('Returning the local matrix is not yet supported')
  def LockedMatrix(self):
    # TODO
    raise Exception('Returning the (locked) local matrix is not yet supported')
  def AllocatedMemory(self):
    # TODO: Python analogue of size_t
    raise Exception('Returning the allocated memory is not yet supported')
  def Buffer(self):
    if   self.tag == iTag:
      buf = POINTER(iType)()
      lib.ElDistMatrixBuffer_i(self.obj,pointer(buf))
      return buf
    elif self.tag == sTag:
      buf = POINTER(sType)()
      lib.ElDistMatrixBuffer_s(self.obj,pointer(buf))
      return buf
    elif self.tag == dTag:
      buf = POINTER(dType)()
      lib.ElDistMatrixBuffer_d(self.obj,pointer(buf))
      return buf
    elif self.tag == cTag:
      buf = POINTER(cType)()
      lib.ElDistMatrixBuffer_c(self.obj,pointer(buf))
      return buf
    elif self.tag == zTag:
      buf = POINTER(zType)()
      lib.ElDistMatrixBuffer_z(self.obj,pointer(buf))
      return buf
  def LockedBuffer(self):
    if   self.tag == iTag:
      buf = POINTER(iType)()
      lib.ElDistMatrixLockedBuffer_i(self.obj,pointer(buf))
      return buf
    elif self.tag == sTag:
      buf = POINTER(sType)()
      lib.ElDistMatrixLockedBuffer_s(self.obj,pointer(buf))
      return buf
    elif self.tag == dTag:
      buf = POINTER(dType)()
      lib.ElDistMatrixLockedBuffer_d(self.obj,pointer(buf))
      return buf
    elif self.tag == cTag:
      buf = POINTER(cType)()
      lib.ElDistMatrixLockedBuffer_c(self.obj,pointer(buf))
      return buf
    elif self.tag == zTag:
      buf = POINTER(zType)()
      lib.ElDistMatrixLockedBuffer_z(self.obj,pointer(buf))
      return buf
  def Grid(self):
    # TODO
    raise Exception('Returning the grid is not yet supported')
  def ColConstrained(self):
    colConst = bType() 
    if   self.tag == iTag: 
      lib.ElDistMatrixColConstrained_i(self.obj,pointer(colConst))
    elif self.tag == sTag:
      lib.ElDistMatrixColConstrained_s(self.obj,pointer(colConst))
    elif self.tag == dTag:
      lib.ElDistMatrixColConstrained_d(self.obj,pointer(colConst))
    elif self.tag == cTag:
      lib.ElDistMatrixColConstrained_c(self.obj,pointer(colConst))
    elif self.tag == zTag:
      lib.ElDistMatrixColConstrained_z(self.obj,pointer(colConst))
    return colConst
  def RowConstrained(self):
    rowConst = bType() 
    if   self.tag == iTag: 
      lib.ElDistMatrixRowConstrained_i(self.obj,pointer(rowConst))
    elif self.tag == sTag:
      lib.ElDistMatrixRowConstrained_s(self.obj,pointer(rowConst))
    elif self.tag == dTag:
      lib.ElDistMatrixRowConstrained_d(self.obj,pointer(rowConst))
    elif self.tag == cTag:
      lib.ElDistMatrixRowConstrained_c(self.obj,pointer(rowConst))
    elif self.tag == zTag:
      lib.ElDistMatrixRowConstrained_z(self.obj,pointer(rowConst))
    return rowConst
  # TODO: Finish implementing class
