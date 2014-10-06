#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import numpy

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

lib.ElMatrixGetDiagonal_i.argtypes = [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetDiagonal_i.restype = c_uint
lib.ElMatrixGetDiagonal_s.argtypes = [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetDiagonal_s.restype = c_uint
lib.ElMatrixGetDiagonal_d.argtypes = [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetDiagonal_d.restype = c_uint
lib.ElMatrixGetDiagonal_c.argtypes = [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetDiagonal_c.restype = c_uint
lib.ElMatrixGetDiagonal_z.argtypes = [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetDiagonal_z.restype = c_uint

lib.ElMatrixGetRealPartOfDiagonal_c.argtypes = \
  [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetRealPartOfDiagonal_c.restype = c_uint
lib.ElMatrixGetRealPartOfDiagonal_z.argtypes = \
  [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetRealPartOfDiagonal_z.restype = c_uint

lib.ElMatrixGetImagPartOfDiagonal_i.argtypes = \
  [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetImagPartOfDiagonal_i.restype = c_uint
lib.ElMatrixGetImagPartOfDiagonal_s.argtypes = \
  [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetImagPartOfDiagonal_s.restype = c_uint
lib.ElMatrixGetImagPartOfDiagonal_d.argtypes = \
  [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetImagPartOfDiagonal_d.restype = c_uint
lib.ElMatrixGetImagPartOfDiagonal_c.argtypes = \
  [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetImagPartOfDiagonal_c.restype = c_uint
lib.ElMatrixGetImagPartOfDiagonal_z.argtypes = \
  [c_void_p,iType,POINTER(c_void_p)]
lib.ElMatrixGetImagPartOfDiagonal_z.restype = c_uint

lib.ElMatrixSetDiagonal_i.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetDiagonal_i.restype = c_uint
lib.ElMatrixSetDiagonal_s.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetDiagonal_s.restype = c_uint
lib.ElMatrixSetDiagonal_d.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetDiagonal_d.restype = c_uint
lib.ElMatrixSetDiagonal_c.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetDiagonal_c.restype = c_uint
lib.ElMatrixSetDiagonal_z.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetDiagonal_z.restype = c_uint

lib.ElMatrixSetRealPartOfDiagonal_c.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetRealPartOfDiagonal_c.restype = c_uint
lib.ElMatrixSetRealPartOfDiagonal_z.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetRealPartOfDiagonal_z.restype = c_uint

lib.ElMatrixSetImagPartOfDiagonal_c.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetImagPartOfDiagonal_c.restype = c_uint
lib.ElMatrixSetImagPartOfDiagonal_z.argtypes = [c_void_p,c_void_p,iType]
lib.ElMatrixSetImagPartOfDiagonal_z.restype = c_uint

lib.ElMatrixUpdateDiagonal_i.argtypes = [c_void_p,iType,c_void_p,iType]
lib.ElMatrixUpdateDiagonal_i.restype = c_uint
lib.ElMatrixUpdateDiagonal_s.argtypes = [c_void_p,sType,c_void_p,iType]
lib.ElMatrixUpdateDiagonal_s.restype = c_uint
lib.ElMatrixUpdateDiagonal_d.argtypes = [c_void_p,dType,c_void_p,iType]
lib.ElMatrixUpdateDiagonal_d.restype = c_uint
lib.ElMatrixUpdateDiagonal_c.argtypes = [c_void_p,cType,c_void_p,iType]
lib.ElMatrixUpdateDiagonal_c.restype = c_uint
lib.ElMatrixUpdateDiagonal_z.argtypes = [c_void_p,zType,c_void_p,iType]
lib.ElMatrixUpdateDiagonal_z.restype = c_uint

lib.ElMatrixUpdateRealPartOfDiagonal_c.argtypes = \
  [c_void_p,sType,c_void_p,iType]
lib.ElMatrixUpdateRealPartOfDiagonal_c.restype = c_uint
lib.ElMatrixUpdateRealPartOfDiagonal_z.argtypes = \
  [c_void_p,dType,c_void_p,iType]
lib.ElMatrixUpdateRealPartOfDiagonal_z.restype = c_uint

lib.ElMatrixUpdateImagPartOfDiagonal_c.argtypes = \
  [c_void_p,sType,c_void_p,iType]
lib.ElMatrixUpdateImagPartOfDiagonal_c.restype = c_uint
lib.ElMatrixUpdateImagPartOfDiagonal_z.argtypes = \
  [c_void_p,dType,c_void_p,iType]
lib.ElMatrixUpdateImagPartOfDiagonal_z.restype = c_uint

lib.ElMatrixMakeDiagonalReal_c.argtypes = [c_void_p,iType]
lib.ElMatrixMakeDiagonalReal_c.restype = c_uint
lib.ElMatrixMakeDiagonalReal_z.argtypes = [c_void_p,iType]
lib.ElMatrixMakeDiagonalReal_z.restype = c_uint

lib.ElMatrixConjugateDiagonal_c.argtypes = [c_void_p,iType]
lib.ElMatrixConjugateDiagonal_c.restype = c_uint
lib.ElMatrixConjugateDiagonal_z.argtypes = [c_void_p,iType]
lib.ElMatrixConjugateDiagonal_z.restype = c_uint

# LEFT OFF HERE

class Matrix(object):
  def __init__(self,tag=dTag,create=True):
    self.obj = c_void_p()
    CheckTag(tag)
    self.tag = tag
    if create:
      if   tag == iTag: lib.ElMatrixCreate_i(pointer(self.obj))
      elif tag == sTag: lib.ElMatrixCreate_s(pointer(self.obj))
      elif tag == dTag: lib.ElMatrixCreate_d(pointer(self.obj))
      elif tag == cTag: lib.ElMatrixCreate_c(pointer(self.obj))
      elif tag == zTag: lib.ElMatrixCreate_z(pointer(self.obj))
  def Destroy(self):
    if   self.tag == iTag: lib.ElMatrixDestroy_i(self.obj)
    elif self.tag == sTag: lib.ElMatrixDestroy_s(self.obj)
    elif self.tag == dTag: lib.ElMatrixDestroy_d(self.obj)
    elif self.tag == cTag: lib.ElMatrixDestroy_c(self.obj)
    elif self.tag == zTag: lib.ElMatrixDestroy_z(self.obj)
  def SetType(self,newTag=dTag):
    self.Destroy()
    CheckTag(newTag)
    if   newTag == iTag: lib.ElMatrixCreate_i(pointer(self.obj))
    elif newTag == sTag: lib.ElMatrixCreate_s(pointer(self.obj))
    elif newTag == dTag: lib.ElMatrixCreate_d(pointer(self.obj))
    elif newTag == cTag: lib.ElMatrixCreate_c(pointer(self.obj))
    elif newTag == zTag: lib.ElMatrixCreate_z(pointer(self.obj))
    self.tag = newTag
  def Resize(self,m,n):
    if   self.tag == iTag: lib.ElMatrixResize_i(self.obj,m,n)
    elif self.tag == sTag: lib.ElMatrixResize_s(self.obj,m,n)
    elif self.tag == dTag: lib.ElMatrixResize_d(self.obj,m,n)
    elif self.tag == cTag: lib.ElMatrixResize_c(self.obj,m,n)
    elif self.tag == zTag: lib.ElMatrixResize_z(self.obj,m,n)
  def ResizeWithLDim(self,m,n,ldim):
    if   self.tag == iTag: lib.ElMatrixResizeWithLDim_i(self.obj,m,n,ldim)
    elif self.tag == sTag: lib.ElMatrixResizeWithLDim_s(self.obj,m,n,ldim)
    elif self.tag == dTag: lib.ElMatrixResizeWithLDim_d(self.obj,m,n,ldim)
    elif self.tag == cTag: lib.ElMatrixResizeWithLDim_c(self.obj,m,n,ldim)
    elif self.tag == zTag: lib.ElMatrixResizeWithLDim_z(self.obj,m,n,ldim)
  def Empty(self):
    if   self.tag == iTag: lib.ElMatrixEmpty_i(self.obj)
    elif self.tag == sTag: lib.ElMatrixEmpty_s(self.obj)
    elif self.tag == dTag: lib.ElMatrixEmpty_d(self.obj)
    elif self.tag == cTag: lib.ElMatrixEmpty_c(self.obj)
    elif self.tag == zTag: lib.ElMatrixEmpty_z(self.obj)
  def Attach(self,m,n,buf,ldim):
    if   self.tag == iTag: lib.ElMatrixAttach_i(self.obj,m,n,buf,ldim)
    elif self.tag == sTag: lib.ElMatrixAttach_s(self.obj,m,n,buf,ldim)
    elif self.tag == dTag: lib.ElMatrixAttach_d(self.obj,m,n,buf,ldim)
    elif self.tag == cTag: lib.ElMatrixAttach_c(self.obj,m,n,buf,ldim)
    elif self.tag == zTag: lib.ElMatrixAttach_z(self.obj,m,n,buf,ldim)
  def LockedAttach(self,m,n,buf,ldim):
    if   self.tag == iTag: lib.ElMatrixLockedAttach_i(self.obj,m,n,buf,ldim)
    elif self.tag == sTag: lib.ElMatrixLockedAttach_s(self.obj,m,n,buf,ldim)
    elif self.tag == dTag: lib.ElMatrixLockedAttach_d(self.obj,m,n,buf,ldim)
    elif self.tag == cTag: lib.ElMatrixLockedAttach_c(self.obj,m,n,buf,ldim)
    elif self.tag == zTag: lib.ElMatrixLockedAttach_z(self.obj,m,n,buf,ldim)
  def Control(self,m,n,buf,ldim):
    if   self.tag == iTag: lib.ElMatrixControl_i(self.obj,m,n,buf,ldim)
    elif self.tag == sTag: lib.ElMatrixControl_s(self.obj,m,n,buf,ldim)
    elif self.tag == dTag: lib.ElMatrixControl_d(self.obj,m,n,buf,ldim)
    elif self.tag == cTag: lib.ElMatrixControl_c(self.obj,m,n,buf,ldim)
    elif self.tag == zTag: lib.ElMatrixControl_z(self.obj,m,n,buf,ldim)
  def Height(self):
    m = iType()
    if   self.tag == iTag: lib.ElMatrixHeight_i(self.obj,pointer(m))
    elif self.tag == sTag: lib.ElMatrixHeight_s(self.obj,pointer(m))
    elif self.tag == dTag: lib.ElMatrixHeight_d(self.obj,pointer(m))
    elif self.tag == cTag: lib.ElMatrixHeight_c(self.obj,pointer(m))
    elif self.tag == zTag: lib.ElMatrixHeight_z(self.obj,pointer(m))
    return m
  def Width(self):
    n = iType()
    if   self.tag == iTag: lib.ElMatrixWidth_i(self.obj,pointer(n))
    elif self.tag == sTag: lib.ElMatrixWidth_s(self.obj,pointer(n))
    elif self.tag == dTag: lib.ElMatrixWidth_d(self.obj,pointer(n))
    elif self.tag == cTag: lib.ElMatrixWidth_c(self.obj,pointer(n))
    elif self.tag == zTag: lib.ElMatrixWidth_z(self.obj,pointer(n))
    return n
  def LDim(self):
    ldim = iType()
    if   self.tag == iTag: lib.ElMatrixLDim_i(self.obj,pointer(ldim))
    elif self.tag == sTag: lib.ElMatrixLDim_s(self.obj,pointer(ldim))
    elif self.tag == dTag: lib.ElMatrixLDim_d(self.obj,pointer(ldim))
    elif self.tag == cTag: lib.ElMatrixLDim_c(self.obj,pointer(ldim))
    elif self.tag == zTag: lib.ElMatrixLDim_z(self.obj,pointer(ldim))
    return ldim 
  def MemorySize(self):
    size = iType()
    if   self.tag == iTag: lib.ElMatrixMemorySize_i(self.obj,pointer(size))
    elif self.tag == sTag: lib.ElMatrixMemorySize_s(self.obj,pointer(size))
    elif self.tag == dTag: lib.ElMatrixMemorySize_d(self.obj,pointer(size))
    elif self.tag == cTag: lib.ElMatrixMemorySize_c(self.obj,pointer(size))
    elif self.tag == zTag: lib.ElMatrixMemorySize_z(self.obj,pointer(size))
    return size
  def DiagonalLength(self,offset=iType(0)):
    length = iType()
    if   self.tag == iTag: 
      lib.ElMatrixDiagonalLength_i(self.obj,offset,pointer(length))
    if   self.tag == sTag: 
      lib.ElMatrixDiagonalLength_s(self.obj,offset,pointer(length))
    if   self.tag == dTag: 
      lib.ElMatrixDiagonalLength_d(self.obj,offset,pointer(length))
    if   self.tag == cTag: 
      lib.ElMatrixDiagonalLength_c(self.obj,offset,pointer(length))
    if   self.tag == zTag: 
      lib.ElMatrixDiagonalLength_z(self.obj,offset,pointer(length))
    return length
  def Viewing(self):
    viewing = bType()
    if   self.tag == iTag: lib.ElMatrixViewing_i(self.obj,pointer(viewing))
    elif self.tag == sTag: lib.ElMatrixViewing_s(self.obj,pointer(viewing))
    elif self.tag == dTag: lib.ElMatrixViewing_d(self.obj,pointer(viewing))
    elif self.tag == cTag: lib.ElMatrixViewing_c(self.obj,pointer(viewing))
    elif self.tag == zTag: lib.ElMatrixViewing_z(self.obj,pointer(viewing))
    return viewing 
  def FixedSize(self):
    fixed = bType()
    if   self.tag == iTag: lib.ElMatrixFixedSize_i(self.obj,pointer(fixed))
    elif self.tag == sTag: lib.ElMatrixFixedSize_s(self.obj,pointer(fixed))
    elif self.tag == dTag: lib.ElMatrixFixedSize_d(self.obj,pointer(fixed))
    elif self.tag == cTag: lib.ElMatrixFixedSize_c(self.obj,pointer(fixed))
    elif self.tag == zTag: lib.ElMatrixFixedSize_z(self.obj,pointer(fixed))
    return fixed
  def Locked(self):
    locked = bType()
    if   self.tag == iTag: lib.ElMatrixLocked_i(self.obj,pointer(locked))
    elif self.tag == sTag: lib.ElMatrixLocked_s(self.obj,pointer(locked))
    elif self.tag == dTag: lib.ElMatrixLocked_d(self.obj,pointer(locked))
    elif self.tag == cTag: lib.ElMatrixLocked_c(self.obj,pointer(locked))
    elif self.tag == zTag: lib.ElMatrixLocked_z(self.obj,pointer(locked))
    return locked
  def Buffer(self):
    if   self.tag == iTag:
      buf = POINTER(iType)()
      lib.ElMatrixBuffer_i(self.obj,pointer(buf)) 
      return buf
    elif self.tag == sTag:
      buf = POINTER(sType)()
      lib.ElMatrixBuffer_s(self.obj,pointer(buf))
      return buf
    elif self.tag == dTag: 
      buf = POINTER(dType)()
      lib.ElMatrixBuffer_d(self.obj,pointer(buf))
      return buf
    elif self.tag == cTag:
      buf = POINTER(cType)()
      lib.ElMatrixBuffer_c(self.obj,pointer(buf))
      return buf
    elif self.tag == zTag:
      buf = POINTER(zType)()
      lib.ElMatrixBuffer_z(self.obj,pointer(buf))
      return buf
  def LockedBuffer(self):
    if   self.tag == iTag:
      buf = POINTER(iType)()
      lib.ElMatrixLockedBuffer_i(self.obj,pointer(buf)) 
      return buf
    elif self.tag == sTag:
      buf = POINTER(sType)()
      lib.ElMatrixLockedBuffer_s(self.obj,pointer(buf))
      return buf
    elif self.tag == dTag: 
      buf = POINTER(dType)()
      lib.ElMatrixLockedBuffer_d(self.obj,pointer(buf))
      return buf
    elif self.tag == cTag:
      buf = POINTER(cType)()
      lib.ElMatrixLockedBuffer_c(self.obj,pointer(buf))
      return buf
    elif self.tag == zTag:
      buf = POINTER(zType)()
      lib.ElMatrixLockedBuffer_z(self.obj,pointer(buf))
      return buf
  def Get(self,i,j):
    if   self.tag == iTag:
      value = iType()
      lib.ElMatrixGet_i(self.obj,i,j,pointer(value))
      return value
    elif self.tag == sTag:
      value = sType()
      lib.ElMatrixGet_s(self.obj,i,j,pointer(value))
      return value
    elif self.tag == dTag:
      value = dType()
      lib.ElMatrixGet_d(self.obj,i,j,pointer(value))
    elif self.tag == cTag:
      value = cType()
      lib.ElMatrixGet_c(self.obj,i,j,pointer(value))
      return value
    elif self.tag == zTag:
      value = zType()
      lib.ElMatrixGet_z(self.obj,i,j,pointer(value))
      return value
  def GetRealPart(self,i,j):
    if self.tag == cTag:
      value = sType()
      lib.ElMatrixGetRealPart_c(self.obj,i,j,pointer(value))
      return value
    elif self.tag == zTag:
      value = dType()
      lib.ElMatrixGetRealPart_z(self.obj,i,j,pointer(value))
      return value
    else: return Get(i,j)
  def GetImagPart(self,i,j):
    if   self.tag == iTag: return iType(0)
    elif self.tag == sTag: return sType(0)
    elif self.tag == dTag: return dType(0)
    elif self.tag == cTag:
      value = c_float()
      lib.ElMatrixGetImagPart_c(self.obj,i,j,pointer(value))
      return value
    elif self.tag == zTag:
      value = c_double()
      lib.ElMatrixGetImagPart_z(self.obj,i,j,pointer(value))
      return value
  def Set(self,i,j,value):
    if   self.tag == iTag: lib.ElMatrixSet_i(self.obj,i,j,iType(value))
    elif self.tag == sTag: lib.ElMatrixSet_s(self.obj,i,j,value)
    elif self.tag == dTag: lib.ElMatrixSet_d(self.obj,i,j,value)
    elif self.tag == cTag: lib.ElMatrixSet_c(self.obj,i,j,value)
    elif self.tag == zTag: lib.ElMatrixSet_z(self.obj,i,j,value)
  def SetRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetRealPart_z(self.obj,i,j,dType(value))
    else: Set(i,j,value)
  def SetImagPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixSetImagPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixSetImagPart_z(self.obj,i,j,dType(value))
    else: raise Exception("Datatype does not have an imaginary component")
  def Update(self,i,j,value):
    if   self.tag == iTag: lib.ElMatrixUpdate_i(self.obj,i,j,iType(value))
    elif self.tag == sTag: lib.ElMatrixUpdate_s(self.obj,i,j,value)
    elif self.tag == dTag: lib.ElMatrixUpdate_d(self.obj,i,j,value)
    elif self.tag == cTag: lib.ElMatrixUpdate_c(self.obj,i,j,value)
    elif self.tag == zTag: lib.ElMatrixUpdate_z(self.obj,i,j,value)
  def UpdateRealPart(self,i,j,value):
    if self.tag == cTag: 
      lib.ElMatrixUpdateRealPart_c(self.obj,i,j,sType(value))
    elif self.tag == zTag: 
      lib.ElMatrixUpdateRealPart_z(self.obj,i,j,dType(value))
    else: Update(i,j,value)
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
  def GetDiagonal(self,offset=iType(0)):
    d = Matrix(self.tag,False)
    if   self.tag == iTag:
      lib.ElMatrixGetDiagonal_i(self.obj,offset,pointer(d.obj))
    elif self.tag == sTag:
      lib.ElMatrixGetDiagonal_s(self.obj,offset,pointer(d.obj))
    elif self.tag == dTag:
      lib.ElMatrixGetDiagonal_d(self.obj,offset,pointer(d.obj))
    elif self.tag == cTag:
      lib.ElMatrixGetDiagonal_c(self.obj,offset,pointer(d.obj))
    elif self.tag == zTag:
      lib.ElMatrixGetDiagonal_z(self.obj,offset,pointer(d.obj))
    return d
  def GetRealPartOfDiagonal(self,offset=iType(0)):
    if self.tag == cTag:
      d = Matrix(sTag,False)
      lib.ElMatrixGetRealPartOfDiagonal_c(self.obj,offset,pointer(d.obj))
      return d
    elif self.tag == zTag:
      d = Matrix(dTag,False)
      lib.ElMatrixGetRealPartOfDiagonal_z(self.obj,offset,pointer(d.obj))
      return d
    else: 
      return GetDiagonal(self,offset)
  def GetImagPartOfDiagonal(self,offset=iType(0)):
    if   self.tag == iTag:
      d = Matrix(iTag,False)
      lib.ElMatrixGetImagPartOfDiagonal_i(self.obj,offset,pointer(d.obj))
      return d
    elif self.tag == sTag:
      d = Matrix(sTag,False)
      lib.ElMatrixGetImagPartOfDiagonal_s(self.obj,offset,pointer(d.obj))
      return d
    elif self.tag == dTag:
      d = Matrix(dTag,False)
      lib.ElMatrixGetImagPartOfDiagonal_d(self.obj,offset,pointer(d.obj))
      return d
    elif self.tag == cTag:
      d = Matrix(sTag,False)
      lib.ElMatrixGetImagPartOfDiagonal_c(self.obj,offset,pointer(d.obj))
      return d
    elif self.tag == zTag:
      d = Matrix(dTag,False)
      lib.ElMatrixGetImagPartOfDiagonal_z(self.obj,offset,pointer(d.obj))
      return d
  def SetDiagonal(self,d,offset=iType(0)):
    if type(d) is not Matrix: raise Exception('diagonal must be a Matrix')
    if self.tag != d.tag: raise Exception('Datatypes must match')
    if   self.tag == iTag:
      lib.ElMatrixSetDiagonal_i(self.obj,d.obj,offset)
    elif self.tag == sTag:
      lib.ElMatrixSetDiagonal_s(self.obj,d.obj,offset)
    elif self.tag == dTag:
      lib.ElMatrixSetDiagonal_d(self.obj,d.obj,offset)
    elif self.tag == cTag:
      lib.ElMatrixSetDiagonal_c(self.obj,d.obj,offset)
    elif self.tag == zTag:
      lib.ElMatrixSetDiagonal_z(self.obj,d.obj,offset)
  def SetRealPartOfDiagonal(self,d,offset=iType(0)):
    if type(d) is not Matrix: raise Exception('diagonal must be a Matrix')
    if   self.tag == cTag:
      if d.tag != sTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixSetRealPartOfDiagonal_c(self.obj,d.obj,offset)
    elif self.tag == zTag:
      if d.tag != dTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixSetRealPartOfDiagonal_z(self.obj,d.obj,offset)
    else:
      SetDiagonal(d,offset)
  def SetImagPartOfDiagonal(self,d,offset=iType(0)):
    if type(d) is not Matrix: raise Exception('diagonal must be a Matrix')
    if   self.tag == cTag:
      if d.tag != sTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixSetImagPartOfDiagonal_c(self.obj,d.obj,offset)
    elif self.tag == zTag:
      if d.tag != dTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixSetImagPartOfDiagonal_z(self.obj,d.obj,offset)
    else:
      raise Exception('Cannot set the imaginary part of a real matrix')
  def UpdateDiagonal(self,alpha,d,offset=iType(0)):
    if type(d) is not Matrix: raise Exception('diagonal must be a Matrix')
    if self.tag != d.tag: raise Exception('Datatypes must match')
    if   self.tag == iTag:
      lib.ElMatrixUpdateDiagonal_i(self.obj,alpha,d.obj,offset)
    elif self.tag == sTag:
      lib.ElMatrixUpdateDiagonal_s(self.obj,alpha,d.obj,offset)
    elif self.tag == dTag:
      lib.ElMatrixUpdateDiagonal_d(self.obj,alpha,d.obj,offset)
    elif self.tag == cTag:
      lib.ElMatrixUpdateDiagonal_c(self.obj,alpha,d.obj,offset)
    elif self.tag == zTag:
      lib.ElMatrixUpdateDiagonal_z(self.obj,alpha,d.obj,offset)
  def UpdateRealPartOfDiagonal(self,alpha,d,offset=iType(0)):
    if type(d) is not Matrix: raise Exception('diagonal must be a Matrix')
    if   self.tag == cTag:
      if d.tag != sTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixUpdateRealPartOfDiagonal_c(self.obj,alpha,d.obj,offset)
    elif self.tag == zTag:
      if d.tag != dTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixUpdateRealPartOfDiagonal_z(self.obj,alpha,d.obj,offset)
    else:
      UpdateDiagonal(alpha,d,offset)
  def UpdateImagPartOfDiagonal(self,alpha,d,offset=iType(0)):
    if type(d) is not Matrix: raise Exception('diagonal must be a Matrix')
    if   self.tag == cTag:
      if d.tag != sTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixUpdateImagPartOfDiagonal_c(self.obj,alpha,d.obj,offset)
    elif self.tag == zTag:
      if d.tag != dTag: raise Exception('Datatypes must be compatible')
      lib.ElMatrixUpdateImagPartOfDiagonal_z(self.obj,alpha,d.obj,offset)
    else:
      raise Exception('Cannot update the imaginary part of a real matrix')
  def MakeDiagonalReal(self,offset=iType(0)):
    if   self.tag == cTag: lib.ElMatrixMakeDiagonalReal_c(self.obj,offset)
    elif self.tag == zTag: lib.ElMatrixMakeDiagonalReal_z(self.obj,offset)
  def ConjugateDiagonal(self,offset=iType(0)):
    if   self.tag == cTag: lib.ElMatrixConjugateDiagonal_c(self.obj,offset)
    elif self.tag == zTag: lib.ElMatrixConjugateDiagonal_z(self.obj,offset)
  def GetSubmatrix(self,I,J):
    raise Exception('GetSubmatrix not yet supported by Python interface')
  def GetRealPartOfSubmatrix(self,I,J):
    raise Exception('GetRealPartOfSubmatrix not yet supported in Python')
  def GetImagPartOfSubmatrix(self,I,J):
    raise Exception('GetImagPartOfSubmatrix not yet supported in Python')
  def SetSubmatrix(self,I,J,ASub):
    raise Exception('SetSubmatrix not yet supported by Python interface')
  def SetRealPartOfSubmatrix(self,I,J,ASub):
    raise Exception('SetRealPartOfSubmatrix not yet supported in Python')
  def SetImagPartOfSubmatrix(self,I,J,ASub):
    raise Exception('SetImagPartOfSubmatrix not yet supported in Python')
  def UpdateSubmatrix(self,I,J,ASub):
    raise Exception('UpdateSubmatrix not yet supported by Python interface')
  def UpdateRealPartOfSubmatrix(self,I,J,ASub):
    raise Exception('UpdateRealPartOfSubmatrix not yet supported in Python')
  def UpdateImagPartOfSubmatrix(self,I,J,ASub):
    raise Exception('UpdateImagPartOfSubmatrix not yet supported in Python')
  def MakeSubmatrixReal(self,I,J):
    raise Exception('MakeSubmatrixReal not yet supported by Python interface')
  def ConjugateSubmatrix(self,I,J):
    raise Exception('ConjugateSubmatrix not yet supported by Python interface')
  def ToNumPy(self):
    if   self.tag == iTag:
      # TODO: Switch to 64-bit based upon Elemental's configuration
      buf = buffer_from_memory(Buffer(),4*LDim()*Width())
      return numpy.frombuffer(buf,numpy.int32)
    elif self.tag == sTag:
      buf = buffer_from_memory(Buffer(),4*LDim()*Width())
      return buffer_from_memory(buf,numpy.float32)
    elif self.tag == dTag:
      buf = buffer_from_memory(Buffer(),8*LDim()*Width())
      return numpy.frombuffer(buf,numpy.float64)
    elif self.tag == cTag: 
      buf = buffer_from_memory(Buffer(),8*LDim()*Width())
      return numpy.frombuffer(buf,numpy.complex64)
    elif self.tag == zTag:
      buf = buffer_from_memory(Buffer(),16*LDim()*Width())
      return numpy.frombuffer(buf,numpy.complex128)
