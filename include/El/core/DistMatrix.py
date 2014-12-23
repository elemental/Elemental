#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
from imports     import mpi
import ctypes, numpy

import Matrix as M
import Grid as G

# (Abstract)DistMatrix
# ====================

class DistData(ctypes.Structure):
  _fields_ = [('colDist',c_uint),
              ('rowDist',c_uint), 
              ('colAlign',iType),
              ('rowAlign',iType),
              ('root',iType),
              ('grid',c_void_p)]

lib.ElDistMatrixCreateSpecific_i.argtypes = \
  [c_uint,c_uint,c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixCreateSpecific_i.restype = c_uint
lib.ElDistMatrixCreateSpecific_s.argtypes = \
  [c_uint,c_uint,c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixCreateSpecific_s.restype = c_uint
lib.ElDistMatrixCreateSpecific_d.argtypes = \
  [c_uint,c_uint,c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixCreateSpecific_d.restype = c_uint
lib.ElDistMatrixCreateSpecific_c.argtypes = \
  [c_uint,c_uint,c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixCreateSpecific_c.restype = c_uint
lib.ElDistMatrixCreateSpecific_z.argtypes = \
  [c_uint,c_uint,c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixCreateSpecific_z.restype = c_uint

lib.ElDistMatrixDestroy_i.argtypes = [c_void_p]
lib.ElDistMatrixDestroy_i.restype = c_uint
lib.ElDistMatrixDestroy_s.argtypes = [c_void_p]
lib.ElDistMatrixDestroy_s.restype = c_uint
lib.ElDistMatrixDestroy_d.argtypes = [c_void_p]
lib.ElDistMatrixDestroy_d.restype = c_uint
lib.ElDistMatrixDestroy_c.argtypes = [c_void_p]
lib.ElDistMatrixDestroy_c.restype = c_uint
lib.ElDistMatrixDestroy_z.argtypes = [c_void_p]
lib.ElDistMatrixDestroy_z.restype = c_uint

lib.ElDistMatrixEmpty_i.argtypes = [c_void_p]
lib.ElDistMatrixEmpty_i.restype = c_uint
lib.ElDistMatrixEmpty_s.argtypes = [c_void_p]
lib.ElDistMatrixEmpty_s.restype = c_uint
lib.ElDistMatrixEmpty_d.argtypes = [c_void_p]
lib.ElDistMatrixEmpty_d.restype = c_uint
lib.ElDistMatrixEmpty_c.argtypes = [c_void_p]
lib.ElDistMatrixEmpty_c.restype = c_uint
lib.ElDistMatrixEmpty_z.argtypes = [c_void_p]
lib.ElDistMatrixEmpty_z.restype = c_uint

lib.ElDistMatrixEmptyData_i.argtypes = [c_void_p]
lib.ElDistMatrixEmptyData_i.restype = c_uint
lib.ElDistMatrixEmptyData_s.argtypes = [c_void_p]
lib.ElDistMatrixEmptyData_s.restype = c_uint
lib.ElDistMatrixEmptyData_d.argtypes = [c_void_p]
lib.ElDistMatrixEmptyData_d.restype = c_uint
lib.ElDistMatrixEmptyData_c.argtypes = [c_void_p]
lib.ElDistMatrixEmptyData_c.restype = c_uint
lib.ElDistMatrixEmptyData_z.argtypes = [c_void_p]
lib.ElDistMatrixEmptyData_z.restype = c_uint

lib.ElDistMatrixSetGrid_i.argtypes = [c_void_p,c_void_p]
lib.ElDistMatrixSetGrid_i.restype = c_uint
lib.ElDistMatrixSetGrid_s.argtypes = [c_void_p,c_void_p]
lib.ElDistMatrixSetGrid_s.restype = c_uint
lib.ElDistMatrixSetGrid_d.argtypes = [c_void_p,c_void_p]
lib.ElDistMatrixSetGrid_d.restype = c_uint
lib.ElDistMatrixSetGrid_c.argtypes = [c_void_p,c_void_p]
lib.ElDistMatrixSetGrid_c.restype = c_uint
lib.ElDistMatrixSetGrid_z.argtypes = [c_void_p,c_void_p]
lib.ElDistMatrixSetGrid_z.restype = c_uint

lib.ElDistMatrixResize_i.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixResize_i.restype = c_uint
lib.ElDistMatrixResize_s.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixResize_s.restype = c_uint
lib.ElDistMatrixResize_d.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixResize_d.restype = c_uint
lib.ElDistMatrixResize_c.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixResize_c.restype = c_uint
lib.ElDistMatrixResize_z.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixResize_z.restype = c_uint

lib.ElDistMatrixResizeWithLDim_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMatrixResizeWithLDim_i.restype = c_uint
lib.ElDistMatrixResizeWithLDim_s.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMatrixResizeWithLDim_s.restype = c_uint
lib.ElDistMatrixResizeWithLDim_d.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMatrixResizeWithLDim_d.restype = c_uint
lib.ElDistMatrixResizeWithLDim_c.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMatrixResizeWithLDim_c.restype = c_uint
lib.ElDistMatrixResizeWithLDim_z.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMatrixResizeWithLDim_z.restype = c_uint

lib.ElDistMatrixMakeConsistent_i.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeConsistent_i.restype = c_uint
lib.ElDistMatrixMakeConsistent_s.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeConsistent_s.restype = c_uint
lib.ElDistMatrixMakeConsistent_d.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeConsistent_d.restype = c_uint
lib.ElDistMatrixMakeConsistent_c.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeConsistent_c.restype = c_uint
lib.ElDistMatrixMakeConsistent_z.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeConsistent_z.restype = c_uint

lib.ElDistMatrixMakeSizeConsistent_i.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeSizeConsistent_i.restype = c_uint
lib.ElDistMatrixMakeSizeConsistent_s.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeSizeConsistent_s.restype = c_uint
lib.ElDistMatrixMakeSizeConsistent_d.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeSizeConsistent_d.restype = c_uint
lib.ElDistMatrixMakeSizeConsistent_c.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeSizeConsistent_c.restype = c_uint
lib.ElDistMatrixMakeSizeConsistent_z.argtypes = [c_void_p,bType]
lib.ElDistMatrixMakeSizeConsistent_z.restype = c_uint

lib.ElDistMatrixAlign_i.argtypes = [c_void_p,iType,iType,bType]
lib.ElDistMatrixAlign_i.restype = c_uint
lib.ElDistMatrixAlign_s.argtypes = [c_void_p,iType,iType,bType]
lib.ElDistMatrixAlign_s.restype = c_uint
lib.ElDistMatrixAlign_d.argtypes = [c_void_p,iType,iType,bType]
lib.ElDistMatrixAlign_d.restype = c_uint
lib.ElDistMatrixAlign_c.argtypes = [c_void_p,iType,iType,bType]
lib.ElDistMatrixAlign_c.restype = c_uint
lib.ElDistMatrixAlign_z.argtypes = [c_void_p,iType,iType,bType]
lib.ElDistMatrixAlign_z.restype = c_uint

lib.ElDistMatrixAlignCols_i.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignCols_i.restype = c_uint
lib.ElDistMatrixAlignCols_s.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignCols_s.restype = c_uint
lib.ElDistMatrixAlignCols_d.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignCols_d.restype = c_uint
lib.ElDistMatrixAlignCols_c.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignCols_c.restype = c_uint
lib.ElDistMatrixAlignCols_z.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignCols_z.restype = c_uint

lib.ElDistMatrixAlignRows_i.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignRows_i.restype = c_uint
lib.ElDistMatrixAlignRows_s.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignRows_s.restype = c_uint
lib.ElDistMatrixAlignRows_d.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignRows_d.restype = c_uint
lib.ElDistMatrixAlignRows_c.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignRows_c.restype = c_uint
lib.ElDistMatrixAlignRows_z.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixAlignRows_z.restype = c_uint

lib.ElDistMatrixFreeAlignments_i.argtypes = [c_void_p]
lib.ElDistMatrixFreeAlignments_i.restype = c_uint
lib.ElDistMatrixFreeAlignments_s.argtypes = [c_void_p]
lib.ElDistMatrixFreeAlignments_s.restype = c_uint
lib.ElDistMatrixFreeAlignments_d.argtypes = [c_void_p]
lib.ElDistMatrixFreeAlignments_d.restype = c_uint
lib.ElDistMatrixFreeAlignments_c.argtypes = [c_void_p]
lib.ElDistMatrixFreeAlignments_c.restype = c_uint
lib.ElDistMatrixFreeAlignments_z.argtypes = [c_void_p]
lib.ElDistMatrixFreeAlignments_z.restype = c_uint

lib.ElDistMatrixSetRoot_i.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixSetRoot_i.restype = c_uint
lib.ElDistMatrixSetRoot_s.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixSetRoot_s.restype = c_uint
lib.ElDistMatrixSetRoot_d.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixSetRoot_d.restype = c_uint
lib.ElDistMatrixSetRoot_c.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixSetRoot_c.restype = c_uint
lib.ElDistMatrixSetRoot_z.argtypes = [c_void_p,iType,bType]
lib.ElDistMatrixSetRoot_z.restype = c_uint

lib.ElDistMatrixAlignWith_i.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignWith_i.restype = c_uint
lib.ElDistMatrixAlignWith_s.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignWith_s.restype = c_uint
lib.ElDistMatrixAlignWith_d.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignWith_d.restype = c_uint
lib.ElDistMatrixAlignWith_c.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignWith_c.restype = c_uint
lib.ElDistMatrixAlignWith_z.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignWith_z.restype = c_uint

lib.ElDistMatrixAlignColsWith_i.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignColsWith_i.restype = c_uint
lib.ElDistMatrixAlignColsWith_s.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignColsWith_s.restype = c_uint
lib.ElDistMatrixAlignColsWith_d.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignColsWith_d.restype = c_uint
lib.ElDistMatrixAlignColsWith_c.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignColsWith_c.restype = c_uint
lib.ElDistMatrixAlignColsWith_z.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignColsWith_z.restype = c_uint

lib.ElDistMatrixAlignRowsWith_i.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignRowsWith_i.restype = c_uint
lib.ElDistMatrixAlignRowsWith_s.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignRowsWith_s.restype = c_uint
lib.ElDistMatrixAlignRowsWith_d.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignRowsWith_d.restype = c_uint
lib.ElDistMatrixAlignRowsWith_c.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignRowsWith_c.restype = c_uint
lib.ElDistMatrixAlignRowsWith_z.argtypes = [c_void_p,DistData,bType]
lib.ElDistMatrixAlignRowsWith_z.restype = c_uint

lib.ElDistMatrixAlignAndResize_i.argtypes = \
  [c_void_p,iType,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignAndResize_i.restype = c_uint
lib.ElDistMatrixAlignAndResize_s.argtypes = \
  [c_void_p,iType,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignAndResize_s.restype = c_uint
lib.ElDistMatrixAlignAndResize_d.argtypes = \
  [c_void_p,iType,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignAndResize_d.restype = c_uint
lib.ElDistMatrixAlignAndResize_c.argtypes = \
  [c_void_p,iType,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignAndResize_c.restype = c_uint
lib.ElDistMatrixAlignAndResize_z.argtypes = \
  [c_void_p,iType,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignAndResize_z.restype = c_uint

lib.ElDistMatrixAlignColsAndResize_i.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignColsAndResize_i.restype = c_uint
lib.ElDistMatrixAlignColsAndResize_s.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignColsAndResize_s.restype = c_uint
lib.ElDistMatrixAlignColsAndResize_d.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignColsAndResize_d.restype = c_uint
lib.ElDistMatrixAlignColsAndResize_c.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignColsAndResize_c.restype = c_uint
lib.ElDistMatrixAlignColsAndResize_z.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignColsAndResize_z.restype = c_uint

lib.ElDistMatrixAlignRowsAndResize_i.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignRowsAndResize_i.restype = c_uint
lib.ElDistMatrixAlignRowsAndResize_s.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignRowsAndResize_s.restype = c_uint
lib.ElDistMatrixAlignRowsAndResize_d.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignRowsAndResize_d.restype = c_uint
lib.ElDistMatrixAlignRowsAndResize_c.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignRowsAndResize_c.restype = c_uint
lib.ElDistMatrixAlignRowsAndResize_z.argtypes = \
  [c_void_p,iType,iType,iType,bType,bType]
lib.ElDistMatrixAlignRowsAndResize_z.restype = c_uint

lib.ElDistMatrixAttach_i.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(iType),iType,iType]
lib.ElDistMatrixAttach_i.restype = c_uint
lib.ElDistMatrixAttach_s.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(sType),iType,iType]
lib.ElDistMatrixAttach_s.restype = c_uint
lib.ElDistMatrixAttach_d.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(dType),iType,iType]
lib.ElDistMatrixAttach_d.restype = c_uint
lib.ElDistMatrixAttach_c.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(cType),iType,iType]
lib.ElDistMatrixAttach_c.restype = c_uint
lib.ElDistMatrixAttach_z.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(zType),iType,iType]
lib.ElDistMatrixAttach_z.restype = c_uint

lib.ElDistMatrixLockedAttach_i.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(iType),iType,iType]
lib.ElDistMatrixLockedAttach_i.restype = c_uint
lib.ElDistMatrixLockedAttach_s.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(sType),iType,iType]
lib.ElDistMatrixLockedAttach_s.restype = c_uint
lib.ElDistMatrixLockedAttach_d.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(dType),iType,iType]
lib.ElDistMatrixLockedAttach_d.restype = c_uint
lib.ElDistMatrixLockedAttach_c.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(cType),iType,iType]
lib.ElDistMatrixLockedAttach_c.restype = c_uint
lib.ElDistMatrixLockedAttach_z.argtypes = \
  [c_void_p,iType,iType,c_void_p,iType,iType,POINTER(zType),iType,iType]
lib.ElDistMatrixLockedAttach_z.restype = c_uint

lib.ElDistMatrixHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixHeight_i.restype = c_uint
lib.ElDistMatrixHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixHeight_s.restype = c_uint
lib.ElDistMatrixHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixHeight_d.restype = c_uint
lib.ElDistMatrixHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixHeight_c.restype = c_uint
lib.ElDistMatrixHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixHeight_z.restype = c_uint

lib.ElDistMatrixWidth_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixWidth_i.restype = c_uint
lib.ElDistMatrixWidth_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixWidth_s.restype = c_uint
lib.ElDistMatrixWidth_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixWidth_d.restype = c_uint
lib.ElDistMatrixWidth_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixWidth_c.restype = c_uint
lib.ElDistMatrixWidth_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixWidth_z.restype = c_uint

lib.ElDistMatrixDiagonalLength_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalLength_i.restype = c_uint
lib.ElDistMatrixDiagonalLength_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalLength_s.restype = c_uint
lib.ElDistMatrixDiagonalLength_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalLength_d.restype = c_uint
lib.ElDistMatrixDiagonalLength_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalLength_c.restype = c_uint
lib.ElDistMatrixDiagonalLength_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalLength_z.restype = c_uint

lib.ElDistMatrixViewing_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixViewing_i.restype = c_uint
lib.ElDistMatrixViewing_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixViewing_s.restype = c_uint
lib.ElDistMatrixViewing_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixViewing_d.restype = c_uint
lib.ElDistMatrixViewing_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixViewing_c.restype = c_uint
lib.ElDistMatrixViewing_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixViewing_z.restype = c_uint

lib.ElDistMatrixLocked_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixLocked_i.restype = c_uint
lib.ElDistMatrixLocked_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixLocked_s.restype = c_uint
lib.ElDistMatrixLocked_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixLocked_d.restype = c_uint
lib.ElDistMatrixLocked_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixLocked_c.restype = c_uint
lib.ElDistMatrixLocked_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixLocked_z.restype = c_uint

lib.ElDistMatrixLocalHeight_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalHeight_i.restype = c_uint
lib.ElDistMatrixLocalHeight_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalHeight_s.restype = c_uint
lib.ElDistMatrixLocalHeight_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalHeight_d.restype = c_uint
lib.ElDistMatrixLocalHeight_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalHeight_c.restype = c_uint
lib.ElDistMatrixLocalHeight_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalHeight_z.restype = c_uint

lib.ElDistMatrixLocalWidth_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalWidth_i.restype = c_uint
lib.ElDistMatrixLocalWidth_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalWidth_s.restype = c_uint
lib.ElDistMatrixLocalWidth_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalWidth_d.restype = c_uint
lib.ElDistMatrixLocalWidth_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalWidth_c.restype = c_uint
lib.ElDistMatrixLocalWidth_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLocalWidth_z.restype = c_uint

lib.ElDistMatrixLDim_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLDim_i.restype = c_uint
lib.ElDistMatrixLDim_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLDim_s.restype = c_uint
lib.ElDistMatrixLDim_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLDim_d.restype = c_uint
lib.ElDistMatrixLDim_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLDim_c.restype = c_uint
lib.ElDistMatrixLDim_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixLDim_z.restype = c_uint

lib.ElDistMatrixMatrix_i.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixMatrix_i.restype = c_uint
lib.ElDistMatrixMatrix_s.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixMatrix_s.restype = c_uint
lib.ElDistMatrixMatrix_d.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixMatrix_d.restype = c_uint
lib.ElDistMatrixMatrix_c.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixMatrix_c.restype = c_uint
lib.ElDistMatrixMatrix_z.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixMatrix_z.restype = c_uint

lib.ElDistMatrixLockedMatrix_i.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixLockedMatrix_i.restype = c_uint
lib.ElDistMatrixLockedMatrix_s.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixLockedMatrix_s.restype = c_uint
lib.ElDistMatrixLockedMatrix_d.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixLockedMatrix_d.restype = c_uint
lib.ElDistMatrixLockedMatrix_c.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixLockedMatrix_c.restype = c_uint
lib.ElDistMatrixLockedMatrix_z.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixLockedMatrix_z.restype = c_uint

lib.ElDistMatrixAllocatedMemory_i.argtypes = [c_void_p,POINTER(c_size_t)]
lib.ElDistMatrixAllocatedMemory_i.restype = c_uint
lib.ElDistMatrixAllocatedMemory_s.argtypes = [c_void_p,POINTER(c_size_t)]
lib.ElDistMatrixAllocatedMemory_s.restype = c_uint
lib.ElDistMatrixAllocatedMemory_d.argtypes = [c_void_p,POINTER(c_size_t)]
lib.ElDistMatrixAllocatedMemory_d.restype = c_uint
lib.ElDistMatrixAllocatedMemory_c.argtypes = [c_void_p,POINTER(c_size_t)]
lib.ElDistMatrixAllocatedMemory_c.restype = c_uint
lib.ElDistMatrixAllocatedMemory_z.argtypes = [c_void_p,POINTER(c_size_t)]
lib.ElDistMatrixAllocatedMemory_z.restype = c_uint

lib.ElDistMatrixBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElDistMatrixBuffer_i.restype = c_uint
lib.ElDistMatrixBuffer_s.argtypes = [c_void_p,POINTER(POINTER(sType))]
lib.ElDistMatrixBuffer_s.restype = c_uint
lib.ElDistMatrixBuffer_d.argtypes = [c_void_p,POINTER(POINTER(dType))]
lib.ElDistMatrixBuffer_d.restype = c_uint
lib.ElDistMatrixBuffer_c.argtypes = [c_void_p,POINTER(POINTER(cType))]
lib.ElDistMatrixBuffer_c.restype = c_uint
lib.ElDistMatrixBuffer_z.argtypes = [c_void_p,POINTER(POINTER(zType))]
lib.ElDistMatrixBuffer_z.restype = c_uint

lib.ElDistMatrixLockedBuffer_i.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElDistMatrixLockedBuffer_i.restype = c_uint
lib.ElDistMatrixLockedBuffer_s.argtypes = [c_void_p,POINTER(POINTER(sType))]
lib.ElDistMatrixLockedBuffer_s.restype = c_uint
lib.ElDistMatrixLockedBuffer_d.argtypes = [c_void_p,POINTER(POINTER(dType))]
lib.ElDistMatrixLockedBuffer_d.restype = c_uint
lib.ElDistMatrixLockedBuffer_c.argtypes = [c_void_p,POINTER(POINTER(cType))]
lib.ElDistMatrixLockedBuffer_c.restype = c_uint
lib.ElDistMatrixLockedBuffer_z.argtypes = [c_void_p,POINTER(POINTER(zType))]
lib.ElDistMatrixLockedBuffer_z.restype = c_uint

lib.ElDistMatrixGrid_i.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixGrid_i.restype = c_uint
lib.ElDistMatrixGrid_s.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixGrid_s.restype = c_uint
lib.ElDistMatrixGrid_d.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixGrid_d.restype = c_uint
lib.ElDistMatrixGrid_c.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixGrid_c.restype = c_uint
lib.ElDistMatrixGrid_z.argtypes = [c_void_p,POINTER(c_void_p)]
lib.ElDistMatrixGrid_z.restype = c_uint

lib.ElDistMatrixColConstrained_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixColConstrained_i.restype = c_uint
lib.ElDistMatrixColConstrained_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixColConstrained_s.restype = c_uint
lib.ElDistMatrixColConstrained_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixColConstrained_d.restype = c_uint
lib.ElDistMatrixColConstrained_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixColConstrained_c.restype = c_uint
lib.ElDistMatrixColConstrained_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixColConstrained_z.restype = c_uint

lib.ElDistMatrixRowConstrained_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRowConstrained_i.restype = c_uint
lib.ElDistMatrixRowConstrained_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRowConstrained_s.restype = c_uint
lib.ElDistMatrixRowConstrained_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRowConstrained_d.restype = c_uint
lib.ElDistMatrixRowConstrained_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRowConstrained_c.restype = c_uint
lib.ElDistMatrixRowConstrained_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRowConstrained_z.restype = c_uint

lib.ElDistMatrixRootConstrained_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRootConstrained_i.restype = c_uint
lib.ElDistMatrixRootConstrained_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRootConstrained_s.restype = c_uint
lib.ElDistMatrixRootConstrained_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRootConstrained_d.restype = c_uint
lib.ElDistMatrixRootConstrained_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRootConstrained_c.restype = c_uint
lib.ElDistMatrixRootConstrained_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixRootConstrained_z.restype = c_uint

lib.ElDistMatrixColAlign_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColAlign_i.restype = c_uint
lib.ElDistMatrixColAlign_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColAlign_s.restype = c_uint
lib.ElDistMatrixColAlign_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColAlign_d.restype = c_uint
lib.ElDistMatrixColAlign_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColAlign_c.restype = c_uint
lib.ElDistMatrixColAlign_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColAlign_z.restype = c_uint

lib.ElDistMatrixRowAlign_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowAlign_i.restype = c_uint
lib.ElDistMatrixRowAlign_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowAlign_s.restype = c_uint
lib.ElDistMatrixRowAlign_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowAlign_d.restype = c_uint
lib.ElDistMatrixRowAlign_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowAlign_c.restype = c_uint
lib.ElDistMatrixRowAlign_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowAlign_z.restype = c_uint

lib.ElDistMatrixColShift_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColShift_i.restype = c_uint
lib.ElDistMatrixColShift_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColShift_s.restype = c_uint
lib.ElDistMatrixColShift_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColShift_d.restype = c_uint
lib.ElDistMatrixColShift_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColShift_c.restype = c_uint
lib.ElDistMatrixColShift_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColShift_z.restype = c_uint

lib.ElDistMatrixRowShift_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowShift_i.restype = c_uint
lib.ElDistMatrixRowShift_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowShift_s.restype = c_uint
lib.ElDistMatrixRowShift_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowShift_d.restype = c_uint
lib.ElDistMatrixRowShift_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowShift_c.restype = c_uint
lib.ElDistMatrixRowShift_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowShift_z.restype = c_uint

lib.ElDistMatrixColRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColRank_i.restype = c_uint
lib.ElDistMatrixColRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColRank_s.restype = c_uint
lib.ElDistMatrixColRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColRank_d.restype = c_uint
lib.ElDistMatrixColRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColRank_c.restype = c_uint
lib.ElDistMatrixColRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColRank_z.restype = c_uint

lib.ElDistMatrixRowRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowRank_i.restype = c_uint
lib.ElDistMatrixRowRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowRank_s.restype = c_uint
lib.ElDistMatrixRowRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowRank_d.restype = c_uint
lib.ElDistMatrixRowRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowRank_c.restype = c_uint
lib.ElDistMatrixRowRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowRank_z.restype = c_uint

lib.ElDistMatrixPartialColRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColRank_i.restype = c_uint
lib.ElDistMatrixPartialColRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColRank_s.restype = c_uint
lib.ElDistMatrixPartialColRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColRank_d.restype = c_uint
lib.ElDistMatrixPartialColRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColRank_c.restype = c_uint
lib.ElDistMatrixPartialColRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColRank_z.restype = c_uint

lib.ElDistMatrixPartialRowRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowRank_i.restype = c_uint
lib.ElDistMatrixPartialRowRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowRank_s.restype = c_uint
lib.ElDistMatrixPartialRowRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowRank_d.restype = c_uint
lib.ElDistMatrixPartialRowRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowRank_c.restype = c_uint
lib.ElDistMatrixPartialRowRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowRank_z.restype = c_uint

lib.ElDistMatrixPartialUnionColRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColRank_i.restype = c_uint
lib.ElDistMatrixPartialUnionColRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColRank_s.restype = c_uint
lib.ElDistMatrixPartialUnionColRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColRank_d.restype = c_uint
lib.ElDistMatrixPartialUnionColRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColRank_c.restype = c_uint
lib.ElDistMatrixPartialUnionColRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColRank_z.restype = c_uint

lib.ElDistMatrixPartialUnionRowRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowRank_i.restype = c_uint
lib.ElDistMatrixPartialUnionRowRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowRank_s.restype = c_uint
lib.ElDistMatrixPartialUnionRowRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowRank_d.restype = c_uint
lib.ElDistMatrixPartialUnionRowRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowRank_c.restype = c_uint
lib.ElDistMatrixPartialUnionRowRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowRank_z.restype = c_uint

lib.ElDistMatrixDistRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistRank_i.restype = c_uint
lib.ElDistMatrixDistRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistRank_s.restype = c_uint
lib.ElDistMatrixDistRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistRank_d.restype = c_uint
lib.ElDistMatrixDistRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistRank_c.restype = c_uint
lib.ElDistMatrixDistRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistRank_z.restype = c_uint

lib.ElDistMatrixCrossRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossRank_i.restype = c_uint
lib.ElDistMatrixCrossRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossRank_s.restype = c_uint
lib.ElDistMatrixCrossRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossRank_d.restype = c_uint
lib.ElDistMatrixCrossRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossRank_c.restype = c_uint
lib.ElDistMatrixCrossRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossRank_z.restype = c_uint

lib.ElDistMatrixRedundantRank_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantRank_i.restype = c_uint
lib.ElDistMatrixRedundantRank_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantRank_s.restype = c_uint
lib.ElDistMatrixRedundantRank_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantRank_d.restype = c_uint
lib.ElDistMatrixRedundantRank_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantRank_c.restype = c_uint
lib.ElDistMatrixRedundantRank_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantRank_z.restype = c_uint

lib.ElDistMatrixRoot_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRoot_i.restype = c_uint
lib.ElDistMatrixRoot_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRoot_s.restype = c_uint
lib.ElDistMatrixRoot_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRoot_d.restype = c_uint
lib.ElDistMatrixRoot_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRoot_c.restype = c_uint
lib.ElDistMatrixRoot_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRoot_z.restype = c_uint

lib.ElDistMatrixParticipating_i.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixParticipating_i.restype = c_uint
lib.ElDistMatrixParticipating_s.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixParticipating_s.restype = c_uint
lib.ElDistMatrixParticipating_d.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixParticipating_d.restype = c_uint
lib.ElDistMatrixParticipating_c.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixParticipating_c.restype = c_uint
lib.ElDistMatrixParticipating_z.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistMatrixParticipating_z.restype = c_uint

lib.ElDistMatrixRowOwner_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixRowOwner_i.restype = c_uint
lib.ElDistMatrixRowOwner_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixRowOwner_s.restype = c_uint
lib.ElDistMatrixRowOwner_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixRowOwner_d.restype = c_uint
lib.ElDistMatrixRowOwner_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixRowOwner_c.restype = c_uint
lib.ElDistMatrixRowOwner_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixRowOwner_z.restype = c_uint

lib.ElDistMatrixColOwner_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixColOwner_i.restype = c_uint
lib.ElDistMatrixColOwner_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixColOwner_s.restype = c_uint
lib.ElDistMatrixColOwner_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixColOwner_d.restype = c_uint
lib.ElDistMatrixColOwner_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixColOwner_c.restype = c_uint
lib.ElDistMatrixColOwner_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixColOwner_z.restype = c_uint

lib.ElDistMatrixOwner_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElDistMatrixOwner_i.restype = c_uint
lib.ElDistMatrixOwner_s.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElDistMatrixOwner_s.restype = c_uint
lib.ElDistMatrixOwner_d.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElDistMatrixOwner_d.restype = c_uint
lib.ElDistMatrixOwner_c.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElDistMatrixOwner_c.restype = c_uint
lib.ElDistMatrixOwner_z.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElDistMatrixOwner_z.restype = c_uint

lib.ElDistMatrixLocalRow_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRow_i.restype = c_uint
lib.ElDistMatrixLocalRow_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRow_s.restype = c_uint
lib.ElDistMatrixLocalRow_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRow_d.restype = c_uint
lib.ElDistMatrixLocalRow_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRow_c.restype = c_uint
lib.ElDistMatrixLocalRow_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRow_z.restype = c_uint

lib.ElDistMatrixLocalCol_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalCol_i.restype = c_uint
lib.ElDistMatrixLocalCol_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalCol_s.restype = c_uint
lib.ElDistMatrixLocalCol_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalCol_d.restype = c_uint
lib.ElDistMatrixLocalCol_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalCol_c.restype = c_uint
lib.ElDistMatrixLocalCol_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalCol_z.restype = c_uint

lib.ElDistMatrixLocalRowOffset_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRowOffset_i.restype = c_uint
lib.ElDistMatrixLocalRowOffset_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRowOffset_s.restype = c_uint
lib.ElDistMatrixLocalRowOffset_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRowOffset_d.restype = c_uint
lib.ElDistMatrixLocalRowOffset_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRowOffset_c.restype = c_uint
lib.ElDistMatrixLocalRowOffset_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalRowOffset_z.restype = c_uint

lib.ElDistMatrixLocalColOffset_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalColOffset_i.restype = c_uint
lib.ElDistMatrixLocalColOffset_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalColOffset_s.restype = c_uint
lib.ElDistMatrixLocalColOffset_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalColOffset_d.restype = c_uint
lib.ElDistMatrixLocalColOffset_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalColOffset_c.restype = c_uint
lib.ElDistMatrixLocalColOffset_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixLocalColOffset_z.restype = c_uint

lib.ElDistMatrixGlobalRow_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalRow_i.restype = c_uint
lib.ElDistMatrixGlobalRow_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalRow_s.restype = c_uint
lib.ElDistMatrixGlobalRow_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalRow_d.restype = c_uint
lib.ElDistMatrixGlobalRow_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalRow_c.restype = c_uint
lib.ElDistMatrixGlobalRow_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalRow_z.restype = c_uint

lib.ElDistMatrixGlobalCol_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalCol_i.restype = c_uint
lib.ElDistMatrixGlobalCol_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalCol_s.restype = c_uint
lib.ElDistMatrixGlobalCol_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalCol_d.restype = c_uint
lib.ElDistMatrixGlobalCol_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalCol_c.restype = c_uint
lib.ElDistMatrixGlobalCol_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixGlobalCol_z.restype = c_uint

lib.ElDistMatrixIsLocalRow_i.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalRow_i.restype = c_uint
lib.ElDistMatrixIsLocalRow_s.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalRow_s.restype = c_uint
lib.ElDistMatrixIsLocalRow_d.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalRow_d.restype = c_uint
lib.ElDistMatrixIsLocalRow_c.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalRow_c.restype = c_uint
lib.ElDistMatrixIsLocalRow_z.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalRow_z.restype = c_uint

lib.ElDistMatrixIsLocalCol_i.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalCol_i.restype = c_uint
lib.ElDistMatrixIsLocalCol_s.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalCol_s.restype = c_uint
lib.ElDistMatrixIsLocalCol_d.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalCol_d.restype = c_uint
lib.ElDistMatrixIsLocalCol_c.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalCol_c.restype = c_uint
lib.ElDistMatrixIsLocalCol_z.argtypes = [c_void_p,iType,POINTER(bType)]
lib.ElDistMatrixIsLocalCol_z.restype = c_uint

lib.ElDistMatrixIsLocal_i.argtypes = [c_void_p,iType,iType,POINTER(bType)]
lib.ElDistMatrixIsLocal_i.restype = c_uint
lib.ElDistMatrixIsLocal_s.argtypes = [c_void_p,iType,iType,POINTER(bType)]
lib.ElDistMatrixIsLocal_s.restype = c_uint
lib.ElDistMatrixIsLocal_d.argtypes = [c_void_p,iType,iType,POINTER(bType)]
lib.ElDistMatrixIsLocal_d.restype = c_uint
lib.ElDistMatrixIsLocal_c.argtypes = [c_void_p,iType,iType,POINTER(bType)]
lib.ElDistMatrixIsLocal_c.restype = c_uint
lib.ElDistMatrixIsLocal_z.argtypes = [c_void_p,iType,iType,POINTER(bType)]
lib.ElDistMatrixIsLocal_z.restype = c_uint

lib.ElDistMatrixDistData_i.argtypes = [c_void_p,POINTER(DistData)]
lib.ElDistMatrixDistData_i.restype = c_uint
lib.ElDistMatrixDistData_s.argtypes = [c_void_p,POINTER(DistData)]
lib.ElDistMatrixDistData_s.restype = c_uint
lib.ElDistMatrixDistData_d.argtypes = [c_void_p,POINTER(DistData)]
lib.ElDistMatrixDistData_d.restype = c_uint
lib.ElDistMatrixDistData_c.argtypes = [c_void_p,POINTER(DistData)]
lib.ElDistMatrixDistData_c.restype = c_uint
lib.ElDistMatrixDistData_z.argtypes = [c_void_p,POINTER(DistData)]
lib.ElDistMatrixDistData_z.restype = c_uint

lib.ElDistMatrixDistComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixDistComm_i.restype = c_uint
lib.ElDistMatrixDistComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixDistComm_s.restype = c_uint
lib.ElDistMatrixDistComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixDistComm_d.restype = c_uint
lib.ElDistMatrixDistComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixDistComm_c.restype = c_uint
lib.ElDistMatrixDistComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixDistComm_z.restype = c_uint

lib.ElDistMatrixCrossComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixCrossComm_i.restype = c_uint
lib.ElDistMatrixCrossComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixCrossComm_s.restype = c_uint
lib.ElDistMatrixCrossComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixCrossComm_d.restype = c_uint
lib.ElDistMatrixCrossComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixCrossComm_c.restype = c_uint
lib.ElDistMatrixCrossComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixCrossComm_z.restype = c_uint

lib.ElDistMatrixRedundantComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRedundantComm_i.restype = c_uint
lib.ElDistMatrixRedundantComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRedundantComm_s.restype = c_uint
lib.ElDistMatrixRedundantComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRedundantComm_d.restype = c_uint
lib.ElDistMatrixRedundantComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRedundantComm_c.restype = c_uint
lib.ElDistMatrixRedundantComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRedundantComm_z.restype = c_uint

lib.ElDistMatrixColComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixColComm_i.restype = c_uint
lib.ElDistMatrixColComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixColComm_s.restype = c_uint
lib.ElDistMatrixColComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixColComm_d.restype = c_uint
lib.ElDistMatrixColComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixColComm_c.restype = c_uint
lib.ElDistMatrixColComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixColComm_z.restype = c_uint

lib.ElDistMatrixRowComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRowComm_i.restype = c_uint
lib.ElDistMatrixRowComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRowComm_s.restype = c_uint
lib.ElDistMatrixRowComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRowComm_d.restype = c_uint
lib.ElDistMatrixRowComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRowComm_c.restype = c_uint
lib.ElDistMatrixRowComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixRowComm_z.restype = c_uint

lib.ElDistMatrixPartialColComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialColComm_i.restype = c_uint
lib.ElDistMatrixPartialColComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialColComm_s.restype = c_uint
lib.ElDistMatrixPartialColComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialColComm_d.restype = c_uint
lib.ElDistMatrixPartialColComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialColComm_c.restype = c_uint
lib.ElDistMatrixPartialColComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialColComm_z.restype = c_uint

lib.ElDistMatrixPartialRowComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialRowComm_i.restype = c_uint
lib.ElDistMatrixPartialRowComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialRowComm_s.restype = c_uint
lib.ElDistMatrixPartialRowComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialRowComm_d.restype = c_uint
lib.ElDistMatrixPartialRowComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialRowComm_c.restype = c_uint
lib.ElDistMatrixPartialRowComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialRowComm_z.restype = c_uint

lib.ElDistMatrixPartialUnionColComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionColComm_i.restype = c_uint
lib.ElDistMatrixPartialUnionColComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionColComm_s.restype = c_uint
lib.ElDistMatrixPartialUnionColComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionColComm_d.restype = c_uint
lib.ElDistMatrixPartialUnionColComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionColComm_c.restype = c_uint
lib.ElDistMatrixPartialUnionColComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionColComm_z.restype = c_uint

lib.ElDistMatrixPartialUnionRowComm_i.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionRowComm_i.restype = c_uint
lib.ElDistMatrixPartialUnionRowComm_s.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionRowComm_s.restype = c_uint
lib.ElDistMatrixPartialUnionRowComm_d.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionRowComm_d.restype = c_uint
lib.ElDistMatrixPartialUnionRowComm_c.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionRowComm_c.restype = c_uint
lib.ElDistMatrixPartialUnionRowComm_z.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistMatrixPartialUnionRowComm_z.restype = c_uint

lib.ElDistMatrixColStride_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColStride_i.restype = c_uint
lib.ElDistMatrixColStride_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColStride_s.restype = c_uint
lib.ElDistMatrixColStride_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColStride_d.restype = c_uint
lib.ElDistMatrixColStride_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColStride_c.restype = c_uint
lib.ElDistMatrixColStride_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixColStride_z.restype = c_uint

lib.ElDistMatrixRowStride_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowStride_i.restype = c_uint
lib.ElDistMatrixRowStride_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowStride_s.restype = c_uint
lib.ElDistMatrixRowStride_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowStride_d.restype = c_uint
lib.ElDistMatrixRowStride_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowStride_c.restype = c_uint
lib.ElDistMatrixRowStride_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRowStride_z.restype = c_uint

lib.ElDistMatrixPartialColStride_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColStride_i.restype = c_uint
lib.ElDistMatrixPartialColStride_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColStride_s.restype = c_uint
lib.ElDistMatrixPartialColStride_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColStride_d.restype = c_uint
lib.ElDistMatrixPartialColStride_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColStride_c.restype = c_uint
lib.ElDistMatrixPartialColStride_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialColStride_z.restype = c_uint

lib.ElDistMatrixPartialRowStride_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowStride_i.restype = c_uint
lib.ElDistMatrixPartialRowStride_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowStride_s.restype = c_uint
lib.ElDistMatrixPartialRowStride_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowStride_d.restype = c_uint
lib.ElDistMatrixPartialRowStride_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowStride_c.restype = c_uint
lib.ElDistMatrixPartialRowStride_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialRowStride_z.restype = c_uint

lib.ElDistMatrixPartialUnionColStride_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColStride_i.restype = c_uint
lib.ElDistMatrixPartialUnionColStride_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColStride_s.restype = c_uint
lib.ElDistMatrixPartialUnionColStride_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColStride_d.restype = c_uint
lib.ElDistMatrixPartialUnionColStride_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColStride_c.restype = c_uint
lib.ElDistMatrixPartialUnionColStride_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionColStride_z.restype = c_uint

lib.ElDistMatrixPartialUnionRowStride_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowStride_i.restype = c_uint
lib.ElDistMatrixPartialUnionRowStride_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowStride_s.restype = c_uint
lib.ElDistMatrixPartialUnionRowStride_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowStride_d.restype = c_uint
lib.ElDistMatrixPartialUnionRowStride_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowStride_c.restype = c_uint
lib.ElDistMatrixPartialUnionRowStride_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixPartialUnionRowStride_z.restype = c_uint

lib.ElDistMatrixDistSize_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistSize_i.restype = c_uint
lib.ElDistMatrixDistSize_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistSize_s.restype = c_uint
lib.ElDistMatrixDistSize_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistSize_d.restype = c_uint
lib.ElDistMatrixDistSize_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistSize_c.restype = c_uint
lib.ElDistMatrixDistSize_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixDistSize_z.restype = c_uint

lib.ElDistMatrixCrossSize_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossSize_i.restype = c_uint
lib.ElDistMatrixCrossSize_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossSize_s.restype = c_uint
lib.ElDistMatrixCrossSize_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossSize_d.restype = c_uint
lib.ElDistMatrixCrossSize_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossSize_c.restype = c_uint
lib.ElDistMatrixCrossSize_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixCrossSize_z.restype = c_uint

lib.ElDistMatrixRedundantSize_i.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantSize_i.restype = c_uint
lib.ElDistMatrixRedundantSize_s.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantSize_s.restype = c_uint
lib.ElDistMatrixRedundantSize_d.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantSize_d.restype = c_uint
lib.ElDistMatrixRedundantSize_c.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantSize_c.restype = c_uint
lib.ElDistMatrixRedundantSize_z.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistMatrixRedundantSize_z.restype = c_uint

lib.ElDistMatrixGet_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
lib.ElDistMatrixGet_i.restype = c_uint
lib.ElDistMatrixGet_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElDistMatrixGet_s.restype = c_uint
lib.ElDistMatrixGet_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElDistMatrixGet_d.restype = c_uint
lib.ElDistMatrixGet_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
lib.ElDistMatrixGet_c.restype = c_uint
lib.ElDistMatrixGet_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
lib.ElDistMatrixGet_z.restype = c_uint

lib.ElDistMatrixGetRealPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElDistMatrixGetRealPart_c.restype = c_uint
lib.ElDistMatrixGetRealPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElDistMatrixGetRealPart_z.restype = c_uint

lib.ElDistMatrixGetImagPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
lib.ElDistMatrixGetImagPart_c.restype = c_uint
lib.ElDistMatrixGetImagPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
lib.ElDistMatrixGetImagPart_z.restype = c_uint

lib.ElDistMatrixSet_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMatrixSet_i.restype = c_uint
lib.ElDistMatrixSet_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMatrixSet_s.restype = c_uint
lib.ElDistMatrixSet_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMatrixSet_d.restype = c_uint
lib.ElDistMatrixSet_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistMatrixSet_c.restype = c_uint
lib.ElDistMatrixSet_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistMatrixSet_z.restype = c_uint

lib.ElDistMatrixSetRealPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMatrixSetRealPart_c.restype = c_uint
lib.ElDistMatrixSetRealPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMatrixSetRealPart_z.restype = c_uint

lib.ElDistMatrixSetImagPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMatrixSetImagPart_c.restype = c_uint
lib.ElDistMatrixSetImagPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMatrixSetImagPart_z.restype = c_uint

lib.ElDistMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
lib.ElDistMatrixUpdate_i.restype = c_uint
lib.ElDistMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMatrixUpdate_s.restype = c_uint
lib.ElDistMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMatrixUpdate_d.restype = c_uint
lib.ElDistMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
lib.ElDistMatrixUpdate_c.restype = c_uint
lib.ElDistMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
lib.ElDistMatrixUpdate_z.restype = c_uint

lib.ElDistMatrixUpdateRealPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMatrixUpdateRealPart_c.restype = c_uint
lib.ElDistMatrixUpdateRealPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMatrixUpdateRealPart_z.restype = c_uint

lib.ElDistMatrixUpdateImagPart_c.argtypes = [c_void_p,iType,iType,sType]
lib.ElDistMatrixUpdateImagPart_c.restype = c_uint
lib.ElDistMatrixUpdateImagPart_z.argtypes = [c_void_p,iType,iType,dType]
lib.ElDistMatrixUpdateImagPart_z.restype = c_uint

lib.ElDistMatrixMakeReal_c.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixMakeReal_c.restype = c_uint
lib.ElDistMatrixMakeReal_z.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixMakeReal_z.restype = c_uint

lib.ElDistMatrixConjugate_c.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixConjugate_c.restype = c_uint
lib.ElDistMatrixConjugate_z.argtypes = [c_void_p,iType,iType]
lib.ElDistMatrixConjugate_z.restype = c_uint

lib.ElDistMatrixDiagonalAlignedWith_i.argtypes = \
  [c_void_p,DistData,iType,POINTER(bType)]
lib.ElDistMatrixDiagonalAlignedWith_i.restype = c_uint
lib.ElDistMatrixDiagonalAlignedWith_s.argtypes = \
  [c_void_p,DistData,iType,POINTER(bType)]
lib.ElDistMatrixDiagonalAlignedWith_s.restype = c_uint
lib.ElDistMatrixDiagonalAlignedWith_d.argtypes = \
  [c_void_p,DistData,iType,POINTER(bType)]
lib.ElDistMatrixDiagonalAlignedWith_d.restype = c_uint
lib.ElDistMatrixDiagonalAlignedWith_c.argtypes = \
  [c_void_p,DistData,iType,POINTER(bType)]
lib.ElDistMatrixDiagonalAlignedWith_c.restype = c_uint
lib.ElDistMatrixDiagonalAlignedWith_z.argtypes = \
  [c_void_p,DistData,iType,POINTER(bType)]
lib.ElDistMatrixDiagonalAlignedWith_z.restype = c_uint

lib.ElDistMatrixDiagonalRoot_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalRoot_i.restype = c_uint
lib.ElDistMatrixDiagonalRoot_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalRoot_s.restype = c_uint
lib.ElDistMatrixDiagonalRoot_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalRoot_d.restype = c_uint
lib.ElDistMatrixDiagonalRoot_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalRoot_c.restype = c_uint
lib.ElDistMatrixDiagonalRoot_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalRoot_z.restype = c_uint

lib.ElDistMatrixDiagonalAlign_i.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalAlign_i.restype = c_uint
lib.ElDistMatrixDiagonalAlign_s.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalAlign_s.restype = c_uint
lib.ElDistMatrixDiagonalAlign_d.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalAlign_d.restype = c_uint
lib.ElDistMatrixDiagonalAlign_c.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalAlign_c.restype = c_uint
lib.ElDistMatrixDiagonalAlign_z.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistMatrixDiagonalAlign_z.restype = c_uint

lib.ElDistMatrixGetSubmatrix_i.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetSubmatrix_i.restype = c_uint
lib.ElDistMatrixGetSubmatrix_s.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetSubmatrix_s.restype = c_uint
lib.ElDistMatrixGetSubmatrix_d.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetSubmatrix_d.restype = c_uint
lib.ElDistMatrixGetSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetSubmatrix_c.restype = c_uint
lib.ElDistMatrixGetSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetSubmatrix_z.restype = c_uint

lib.ElDistMatrixGetRealPartOfSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetSubmatrix_c.restype = c_uint
lib.ElDistMatrixGetRealPartOfSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetSubmatrix_z.restype = c_uint

lib.ElDistMatrixGetImagPartOfSubmatrix_i.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfSubmatrix_i.restype = c_uint
lib.ElDistMatrixGetImagPartOfSubmatrix_s.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfSubmatrix_s.restype = c_uint
lib.ElDistMatrixGetImagPartOfSubmatrix_d.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfSubmatrix_d.restype = c_uint
lib.ElDistMatrixGetImagPartOfSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfSubmatrix_c.restype = c_uint
lib.ElDistMatrixGetImagPartOfSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfSubmatrix_z.restype = c_uint

lib.ElDistMatrixSetSubmatrix_i.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetSubmatrix_i.restype = c_uint
lib.ElDistMatrixSetSubmatrix_s.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetSubmatrix_s.restype = c_uint
lib.ElDistMatrixSetSubmatrix_d.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetSubmatrix_d.restype = c_uint
lib.ElDistMatrixSetSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetSubmatrix_c.restype = c_uint
lib.ElDistMatrixSetSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetSubmatrix_z.restype = c_uint

lib.ElDistMatrixSetRealPartOfSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetRealPartOfSubmatrix_c.restype = c_uint
lib.ElDistMatrixSetRealPartOfSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetRealPartOfSubmatrix_z.restype = c_uint

lib.ElDistMatrixSetImagPartOfSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetImagPartOfSubmatrix_c.restype = c_uint
lib.ElDistMatrixSetImagPartOfSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetImagPartOfSubmatrix_z.restype = c_uint

lib.ElDistMatrixUpdateSubmatrix_i.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),iType,c_void_p]
lib.ElDistMatrixUpdateSubmatrix_i.restype = c_uint
lib.ElDistMatrixUpdateSubmatrix_s.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),sType,c_void_p]
lib.ElDistMatrixUpdateSubmatrix_s.restype = c_uint
lib.ElDistMatrixUpdateSubmatrix_d.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),dType,c_void_p]
lib.ElDistMatrixUpdateSubmatrix_d.restype = c_uint
lib.ElDistMatrixUpdateSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),cType,c_void_p]
lib.ElDistMatrixUpdateSubmatrix_c.restype = c_uint
lib.ElDistMatrixUpdateSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),zType,c_void_p]
lib.ElDistMatrixUpdateSubmatrix_z.restype = c_uint

lib.ElDistMatrixUpdateRealPartOfSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),sType,c_void_p]
lib.ElDistMatrixUpdateRealPartOfSubmatrix_c.restype = c_uint
lib.ElDistMatrixUpdateRealPartOfSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),dType,c_void_p]
lib.ElDistMatrixUpdateRealPartOfSubmatrix_z.restype = c_uint

lib.ElDistMatrixUpdateImagPartOfSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),sType,c_void_p]
lib.ElDistMatrixUpdateImagPartOfSubmatrix_c.restype = c_uint
lib.ElDistMatrixUpdateImagPartOfSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),dType,c_void_p]
lib.ElDistMatrixUpdateImagPartOfSubmatrix_z.restype = c_uint

lib.ElDistMatrixMakeSubmatrixReal_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixMakeSubmatrixReal_c.restype = c_uint
lib.ElDistMatrixMakeSubmatrixReal_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixMakeSubmatrixReal_z.restype = c_uint

lib.ElDistMatrixConjugateSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixConjugateSubmatrix_c.restype = c_uint
lib.ElDistMatrixConjugateSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixConjugateSubmatrix_z.restype = c_uint

# TODO: Remove local versions?
lib.ElDistMatrixGetLocalSubmatrix_i.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetLocalSubmatrix_i.restype = c_uint
lib.ElDistMatrixGetLocalSubmatrix_s.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetLocalSubmatrix_s.restype = c_uint
lib.ElDistMatrixGetLocalSubmatrix_d.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetLocalSubmatrix_d.restype = c_uint
lib.ElDistMatrixGetLocalSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixGetLocalSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixGetRealPartOfLocalSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixGetRealPartOfLocalSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixGetImagPartOfLocalSubmatrix_i.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_i.restype = c_uint
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_s.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_s.restype = c_uint
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_d.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_d.restype = c_uint
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType),c_void_p]
lib.ElDistMatrixGetImagPartOfLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixSetLocalSubmatrix_i.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetLocalSubmatrix_i.restype = c_uint
lib.ElDistMatrixSetLocalSubmatrix_s.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetLocalSubmatrix_s.restype = c_uint
lib.ElDistMatrixSetLocalSubmatrix_d.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetLocalSubmatrix_d.restype = c_uint
lib.ElDistMatrixSetLocalSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixSetLocalSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixSetRealPartOfLocalSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetRealPartOfLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixSetRealPartOfLocalSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetRealPartOfLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixSetImagPartOfLocalSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetImagPartOfLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixSetImagPartOfLocalSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),c_void_p]
lib.ElDistMatrixSetImagPartOfLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixUpdateLocalSubmatrix_i.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),iType,c_void_p]
lib.ElDistMatrixUpdateLocalSubmatrix_i.restype = c_uint
lib.ElDistMatrixUpdateLocalSubmatrix_s.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),sType,c_void_p]
lib.ElDistMatrixUpdateLocalSubmatrix_s.restype = c_uint
lib.ElDistMatrixUpdateLocalSubmatrix_d.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),dType,c_void_p]
lib.ElDistMatrixUpdateLocalSubmatrix_d.restype = c_uint
lib.ElDistMatrixUpdateLocalSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),cType,c_void_p]
lib.ElDistMatrixUpdateLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixUpdateLocalSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),zType,c_void_p]
lib.ElDistMatrixUpdateLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixUpdateRealPartOfLocalSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),sType,c_void_p]
lib.ElDistMatrixUpdateRealPartOfLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixUpdateRealPartOfLocalSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),dType,c_void_p]
lib.ElDistMatrixUpdateRealPartOfLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixUpdateImagPartOfLocalSubmatrix_c.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),sType,c_void_p]
lib.ElDistMatrixUpdateImagPartOfLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixUpdateImagPartOfLocalSubmatrix_z.argtypes = \
  [c_void_p,POINTER(iType),POINTER(iType),dType,c_void_p]
lib.ElDistMatrixUpdateImagPartOfLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixMakeLocalSubmatrixReal_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixMakeLocalSubmatrixReal_c.restype = c_uint
lib.ElDistMatrixMakeLocalSubmatrixReal_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixMakeLocalSubmatrixReal_z.restype = c_uint

lib.ElDistMatrixConjugateLocalSubmatrix_c.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixConjugateLocalSubmatrix_c.restype = c_uint
lib.ElDistMatrixConjugateLocalSubmatrix_z.argtypes = \
  [c_void_p,iType,POINTER(iType),iType,POINTER(iType)]
lib.ElDistMatrixConjugateLocalSubmatrix_z.restype = c_uint

lib.ElDistMatrixSumOver_i.argtypes = [c_void_p,mpi.Comm]
lib.ElDistMatrixSumOver_i.restype = c_uint
lib.ElDistMatrixSumOver_s.argtypes = [c_void_p,mpi.Comm]
lib.ElDistMatrixSumOver_s.restype = c_uint
lib.ElDistMatrixSumOver_d.argtypes = [c_void_p,mpi.Comm]
lib.ElDistMatrixSumOver_d.restype = c_uint
lib.ElDistMatrixSumOver_c.argtypes = [c_void_p,mpi.Comm]
lib.ElDistMatrixSumOver_c.restype = c_uint
lib.ElDistMatrixSumOver_z.argtypes = [c_void_p,mpi.Comm]
lib.ElDistMatrixSumOver_z.restype = c_uint

lib.ElViewDist_i.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElViewDist_i.restype = c_uint
lib.ElViewDist_s.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElViewDist_s.restype = c_uint
lib.ElViewDist_d.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElViewDist_d.restype = c_uint
lib.ElViewDist_c.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElViewDist_c.restype = c_uint
lib.ElViewDist_z.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElViewDist_z.restype = c_uint

lib.ElLockedViewDist_i.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedViewDist_i.restype = c_uint
lib.ElLockedViewDist_s.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedViewDist_s.restype = c_uint
lib.ElLockedViewDist_d.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedViewDist_d.restype = c_uint
lib.ElLockedViewDist_c.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedViewDist_c.restype = c_uint
lib.ElLockedViewDist_z.argtypes = [c_void_p,c_void_p,IndexRange,IndexRange]
lib.ElLockedViewDist_z.restype = c_uint

class DistMatrix(object):
  def __init__(self,tag=dTag,colDist=MC,rowDist=MR,grid=G.DefaultGrid()):
    self.obj = c_void_p()
    CheckTag(tag)
    args = [colDist,rowDist,grid.obj,pointer(self.obj)]
    if   tag == iTag: lib.ElDistMatrixCreateSpecific_i(*args)
    elif tag == sTag: lib.ElDistMatrixCreateSpecific_s(*args)
    elif tag == dTag: lib.ElDistMatrixCreateSpecific_d(*args)
    elif tag == cTag: lib.ElDistMatrixCreateSpecific_c(*args)
    elif tag == zTag: lib.ElDistMatrixCreateSpecific_z(*args)
    else: DataExcept()
    self.tag = tag
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMatrixDestroy_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDestroy_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDestroy_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDestroy_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDestroy_z(*args)
    else: DataExcept()
  def Empty(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixEmpty_z(*args)
    else: DataExcept()
  def EmptyData(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMatrixEmptyData_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixEmptyData_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixEmptyData_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixEmptyData_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixEmptyData_z(*args)
    else: DataExcept()
  def SetGrid(self,grid):
    args = [self.obj,grid.obj]
    if   self.tag == iTag: lib.ElDistMatrixSetGrid_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetGrid_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetGrid_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetGrid_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetGrid_z(*args)
    else: DataExcept()
  def Resize(self,m,n):
    args = [self.obj,m,n]
    if   self.tag == iTag: lib.ElDistMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixResize_z(*args)
    else: DataExcept()
  def ResizeWithLDim(self,m,n,ldim):
    args = [self.obj,m,n,ldim]
    if   self.tag == iTag: lib.ElDistMatrixResizeWithLDim_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixResizeWithLDim_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixResizeWithLDim_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixResizeWithLDim_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixResizeWithLDim_z(*args)
    else: DataExcept()
  def MakeConsistent(self,incViewers):
    args = [self.obj,incViewers]
    if   self.tag == iTag: lib.ElDistMatrixMakeConsistent_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixMakeConsistent_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixMakeConsistent_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixMakeConsistent_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeConsistent_z(*args)
    else: DataExcept()
  def MakeSizeConsistent(self,incViewers):
    args = [self.obj,incViewers]
    if   self.tag == iTag: lib.ElDistMatrixMakeSizeConsistent_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixMakeSizeConsistent_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixMakeSizeConsistent_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixMakeSizeConsistent_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeSizeConsistent_z(*args)
    else: DataExcept()
  def Align(self,colAlign,rowAlign,constrain):
    args = [self.obj,colAlign,rowAlign,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlign_z(*args)
    else: DataExcept()
  def AlignCols(self,colAlign,constrain):
    args = [self.obj,colAlign,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignCols_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignCols_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignCols_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignCols_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignCols_z(*args)
    else: DataExcept()
  def AlignRows(self,rowAlign,constrain):
    args = [self.obj,rowAlign,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignRows_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignRows_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignRows_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignRows_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignRows_z(*args)
    else: DataExcept()
  def FreeAlignments(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMatrixFreeAlignments_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixFreeAlignments_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixFreeAlignments_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixFreeAlignments_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixFreeAlignments_z(*args)
    else: DataExcept()
  def SetRoot(self,root,constrain=True):
    args = [self.obj,root,constrain]
    if   self.tag == iTag: lib.ElDistMatrixSetRoot_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetRoot_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetRoot_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetRoot_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetRoot_z(*args)
    else: DataExcept()
  def AlignWith(self,distData,constrain):
    args = [self.obj,distData,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignWith_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignWith_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignWith_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignWith_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignWith_z(*args)
    else: DataExcept()
  def AlignColsWith(self,distData,constrain):
    args = [self.obj,distData,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignColsWith_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignColsWith_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignColsWith_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignColsWith_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignColsWith_z(*args)
    else: DataExcept()
  def AlignRowsWith(self,distData,constrain):
    args = [self.obj,distData,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignRowsWith_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignRowsWith_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignRowsWith_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignRowsWith_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignRowsWith_z(*args)
    else: DataExcept()
  def AlignAndResize(self,colAlign,rowAlign,m,n,force,constrain):
    args = [self.obj,colAlign,rowAlign,m,n,force,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignAndResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignAndResize_s(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignAndResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignAndResize_z(*args)
    else: DataExcept()
  def AlignColsAndResize(self,colAlign,m,n,force,constrain):
    args = [self.obj,colAlign,m,n,force,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignColsAndResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignColsAndResize_s(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignColsAndResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignColsAndResize_z(*args)
    else: DataExcept()
  def AlignRowsAndResize(self,rowAlign,m,n,force,constrain):
    args = [self.obj,rowAlign,m,n,force,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignRowsAndResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignRowsAndResize_s(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignRowsAndResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignRowsAndResize_z(*args)
    else: DataExcept()
  def Attach(self,m,n,grid,colAlign,rowAlign,buf,ldim,root):
    args = [self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root]
    if   self.tag == iTag: lib.ElDistMatrixAttach_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAttach_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAttach_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAttach_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAttach_z(*args)
    else: DataExcept()
  def LockedAttach(self,m,n,buf,ldim):
    args = [self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root]
    if   self.tag == iTag: lib.ElDistMatrixLockedAttach_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLockedAttach_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLockedAttach_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLockedAttach_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLockedAttach_z(*args)
    else: DataExcept()
  def Height(self):
    m = iType()
    args = [self.obj,pointer(m)]
    if   self.tag == iTag: lib.ElDistMatrixHeight_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixHeight_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixHeight_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixHeight_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixHeight_z(*args)
    else: DataExcept()
    return m.value
  def Width(self):
    n = iType()
    args = [self.obj,pointer(n)]
    if   self.tag == iTag: lib.ElDistMatrixWidth_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixWidth_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixWidth_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixWidth_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixWidth_z(*args)
    else: DataExcept()
    return n.value
  def DiagonalLength(self,offset=0):
    length = iType()
    args = [self.obj,offset,pointer(length)]
    if   self.tag == iTag: lib.ElDistMatrixDiagonalLength_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDiagonalLength_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDiagonalLength_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDiagonalLength_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDiagonalLength_z(*args)
    else: DataExcept()
    return length.value
  def Viewing(self):
    viewing = bType()
    args = [self.obj,pointer(viewing)]
    if   self.tag == iTag: lib.ElDistMatrixViewing_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixViewing_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixViewing_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixViewing_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixViewing_z(*args)
    else: DataExcept()
    return viewing.value
  def Locked(self):
    locked = bType()
    args = [self.obj,pointer(locked)]
    if   self.tag == iTag: lib.ElDistMatrixLocked_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLocked_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLocked_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLocked_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLocked_z(*args)
    else: DataExcept()
    return locked.value
  def LocalHeight(self):
    mLoc = iType()
    args = [self.obj,pointer(mLoc)]
    if   self.tag == iTag: lib.ElDistMatrixLocalHeight_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLocalHeight_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLocalHeight_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLocalHeight_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLocalHeight_z(*args)
    else: DataExcept()
    return mLoc.value
  def LocalWidth(self):
    nLoc = iType()
    args = [self.obj,pointer(nLoc)]
    if   self.tag == iTag: lib.ElDistMatrixLocalWidth_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLocalWidth_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLocalWidth_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLocalWidth_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLocalWidth_z(*args)
    else: DataExcept()
    return nLoc.value
  def LDim(self):
    ldim = iType()
    args = [self.obj,pointer(ldim)]
    if   self.tag == iTag: lib.ElDistMatrixLDim_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLDim_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLDim_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLDim_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLDim_z(*args)
    else: DataExcept()
    return ldim.value
  def Matrix(self):
    A = M.Matrix(self.tag,False)
    args = [self.obj,pointer(A.obj)]
    if   self.tag == iTag: lib.ElDistMatrixMatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixMatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixMatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixMatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMatrix_z(*args)
    else: DataExcept()
    return A
  def LockedMatrix(self):
    A = M.Matrix(self.tag,False)
    args = [self.obj,pointer(A.obj)]
    if   self.tag == iTag: lib.ElDistMatrixLockedMatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLockedMatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLockedMatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLockedMatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLockedMatrix_z(*args)
    else: DataExcept()
    return A
  def AllocatedMemory(self):
    allocMem = c_size_t()
    args = [self.obj,pointer(allocMem)]
    if   self.tag == iTag: lib.ElDistMatrixAllocatedMemory_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAllocatedMemory_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAllocatedMemory_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAllocatedMemory_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAllocatedMemory_z(*args)
    else: DataExcept()
  def Buffer(self):
    buf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(buf)]
    if   self.tag == iTag: lib.ElDistMatrixBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixBuffer_z(*args)
    else: DataExcept()
    return buf
  def LockedBuffer(self):
    buf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(buf)]
    if   self.tag == iTag: lib.ElDistMatrixLockedBuffer_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLockedBuffer_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLockedBuffer_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLockedBuffer_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLockedBuffer_z(*args)
    else: DataExcept()
    return buf
  def Grid(self):
    grid = G.Grid()
    args = [self.obj,pointer(grid.obj)]
    if   self.tag == iTag: lib.ElDistMatrixGrid_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGrid_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGrid_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGrid_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGrid_z(*args)
    else: DataExcept()
    return grid 
  def ColConstrained(self):
    colConst = bType() 
    args = [self.obj,pointer(colConst)]
    if   self.tag == iTag: lib.ElDistMatrixColConstrained_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColConstrained_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColConstrained_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColConstrained_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColConstrained_z(*args)
    else: DataExcept()
    return colConst.value
  def RowConstrained(self):
    rowConst = bType() 
    args = [self.obj,pointer(rowConst)]
    if   self.tag == iTag: lib.ElDistMatrixRowConstrained_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowConstrained_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowConstrained_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowConstrained_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowConstrained_z(*args)
    else: DataExcept()
    return rowConst.value
  def RootConstrained(self):
    rootConst = bType() 
    args = [self.obj,pointer(rootConst)]
    if   self.tag == iTag: lib.ElDistMatrixRootConstrained_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRootConstrained_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRootConstrained_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRootConstrained_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRootConstrained_z(*args)
    else: DataExcept()
    return rootConst.value
  def ColAlign(self):
    align = iType()  
    args = [self.obj,pointer(align)]
    if   self.tag == iTag: lib.ElDistMatrixColAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColAlign_z(*args)
    else: DataExcept()
    return align.value
  def RowAlign(self):
    align = iType()  
    args = [self.obj,pointer(align)]
    if   self.tag == iTag: lib.ElDistMatrixRowAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowAlign_z(*args)
    else: DataExcept()
    return align.value
  def ColShift(self):
    shift = iType()
    if   self.tag == iTag: lib.ElDistMatrixColShift_i(self.obj,pointer(shift))
    elif self.tag == sTag: lib.ElDistMatrixColShift_s(self.obj,pointer(shift))
    elif self.tag == dTag: lib.ElDistMatrixColShift_d(self.obj,pointer(shift))
    elif self.tag == cTag: lib.ElDistMatrixColShift_c(self.obj,pointer(shift))
    elif self.tag == zTag: lib.ElDistMatrixColShift_z(self.obj,pointer(shift))
    else: DataExcept()
    return shift.value
  def RowShift(self):
    shift = iType()
    if   self.tag == iTag: lib.ElDistMatrixRowShift_i(self.obj,pointer(shift))
    elif self.tag == sTag: lib.ElDistMatrixRowShift_s(self.obj,pointer(shift))
    elif self.tag == dTag: lib.ElDistMatrixRowShift_d(self.obj,pointer(shift))
    elif self.tag == cTag: lib.ElDistMatrixRowShift_c(self.obj,pointer(shift))
    elif self.tag == zTag: lib.ElDistMatrixRowShift_z(self.obj,pointer(shift))
    else: DataExcept()
    return shift.value
  def ColRank(self):
    rank = iType()
    if   self.tag == iTag: lib.ElDistMatrixColRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: lib.ElDistMatrixColRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: lib.ElDistMatrixColRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: lib.ElDistMatrixColRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: lib.ElDistMatrixColRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def RowRank(self):
    rank = iType()
    if   self.tag == iTag: lib.ElDistMatrixRowRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: lib.ElDistMatrixRowRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: lib.ElDistMatrixRowRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: lib.ElDistMatrixRowRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: lib.ElDistMatrixRowRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def PartialColRank(self):
    rank = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialColRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialColRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialColRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialColRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialColRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def PartialRowRank(self):
    rank = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialRowRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialRowRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialRowRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialRowRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialRowRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def PartialUnionColRank(self):
    rank = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialUnionColRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialUnionColRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialUnionColRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialUnionColRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialUnionColRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def PartialUnionRowRank(self):
    rank = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialUnionRowRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialUnionRowRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialUnionRowRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialUnionRowRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialUnionRowRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def DistRank(self):
    rank = iType()
    if   self.tag == iTag: lib.ElDistMatrixDistRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: lib.ElDistMatrixDistRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: lib.ElDistMatrixDistRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: lib.ElDistMatrixDistRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: lib.ElDistMatrixDistRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def CrossRank(self):
    rank = iType()
    if   self.tag == iTag: lib.ElDistMatrixCrossRank_i(self.obj,pointer(rank))
    elif self.tag == sTag: lib.ElDistMatrixCrossRank_s(self.obj,pointer(rank))
    elif self.tag == dTag: lib.ElDistMatrixCrossRank_d(self.obj,pointer(rank))
    elif self.tag == cTag: lib.ElDistMatrixCrossRank_c(self.obj,pointer(rank))
    elif self.tag == zTag: lib.ElDistMatrixCrossRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def RedundantRank(self):
    rank = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixRedundantRank_i(self.obj,pointer(rank))
    elif self.tag == sTag:
      lib.ElDistMatrixRedundantRank_s(self.obj,pointer(rank))
    elif self.tag == dTag:
      lib.ElDistMatrixRedundantRank_d(self.obj,pointer(rank))
    elif self.tag == cTag:
      lib.ElDistMatrixRedundantRank_c(self.obj,pointer(rank))
    elif self.tag == zTag:
      lib.ElDistMatrixRedundantRank_z(self.obj,pointer(rank))
    else: DataExcept()
    return rank.value
  def Root(self):
    root = iType()
    if   self.tag == iTag: lib.ElDistMatrixRoot_i(self.obj,pointer(root))
    elif self.tag == sTag: lib.ElDistMatrixRoot_s(self.obj,pointer(root))
    elif self.tag == dTag: lib.ElDistMatrixRoot_d(self.obj,pointer(root))
    elif self.tag == cTag: lib.ElDistMatrixRoot_c(self.obj,pointer(root))
    elif self.tag == zTag: lib.ElDistMatrixRoot_z(self.obj,pointer(root))
    else: DataExcept()
    return root.value
  def Participating(self):
    partic = bType()
    if   self.tag == iTag: 
      lib.ElDistMatrixParticipating_i(self.obj,pointer(partic))
    elif self.tag == sTag:
      lib.ElDistMatrixParticipating_s(self.obj,pointer(partic))
    elif self.tag == dTag:
      lib.ElDistMatrixParticipating_d(self.obj,pointer(partic))
    elif self.tag == cTag:
      lib.ElDistMatrixParticipating_c(self.obj,pointer(partic))
    elif self.tag == zTag:
      lib.ElDistMatrixParticipating_z(self.obj,pointer(partic))
    else: DataExcept()
    return partic.value
  def RowOwner(self,i):
    owner = iType()
    if   self.tag == iTag: lib.ElDistMatrixRowOwner_i(self.obj,i,pointer(owner))
    elif self.tag == sTag: lib.ElDistMatrixRowOwner_s(self.obj,i,pointer(owner))
    elif self.tag == dTag: lib.ElDistMatrixRowOwner_d(self.obj,i,pointer(owner))
    elif self.tag == cTag: lib.ElDistMatrixRowOwner_c(self.obj,i,pointer(owner))
    elif self.tag == zTag: lib.ElDistMatrixRowOwner_z(self.obj,i,pointer(owner))
    else: DataExcept()
    return owner.value
  def ColOwner(self,j):
    owner = iType()
    if   self.tag == iTag: lib.ElDistMatrixColOwner_i(self.obj,j,pointer(owner))
    elif self.tag == sTag: lib.ElDistMatrixColOwner_s(self.obj,j,pointer(owner))
    elif self.tag == dTag: lib.ElDistMatrixColOwner_d(self.obj,j,pointer(owner))
    elif self.tag == cTag: lib.ElDistMatrixColOwner_c(self.obj,j,pointer(owner))
    elif self.tag == zTag: lib.ElDistMatrixColOwner_z(self.obj,j,pointer(owner))
    else: DataExcept()
    return owner.value
  def Owner(self,i,j):
    owner = iType()
    if   self.tag == iTag: lib.ElDistMatrixOwner_i(self.obj,i,j,pointer(owner))
    elif self.tag == sTag: lib.ElDistMatrixOwner_s(self.obj,i,j,pointer(owner))
    elif self.tag == dTag: lib.ElDistMatrixOwner_d(self.obj,i,j,pointer(owner))
    elif self.tag == cTag: lib.ElDistMatrixOwner_c(self.obj,i,j,pointer(owner))
    elif self.tag == zTag: lib.ElDistMatrixOwner_z(self.obj,i,j,pointer(owner))
    else: DataExcept()
    return owner.value
  def LocalRow(self,i):
    iLoc = iType()
    if   self.tag == iTag: lib.ElDistMatrixLocalRow_i(self.obj,i,pointer(iLoc))
    elif self.tag == sTag: lib.ElDistMatrixLocalRow_s(self.obj,i,pointer(iLoc))
    elif self.tag == dTag: lib.ElDistMatrixLocalRow_d(self.obj,i,pointer(iLoc))
    elif self.tag == cTag: lib.ElDistMatrixLocalRow_c(self.obj,i,pointer(iLoc))
    elif self.tag == zTag: lib.ElDistMatrixLocalRow_z(self.obj,i,pointer(iLoc))
    else: DataExcept()
    return iLoc.value
  def LocalCol(self,j):
    jLoc = iType()
    if   self.tag == iTag: lib.ElDistMatrixLocalCol_i(self.obj,j,pointer(jLoc))
    elif self.tag == sTag: lib.ElDistMatrixLocalCol_s(self.obj,j,pointer(jLoc))
    elif self.tag == dTag: lib.ElDistMatrixLocalCol_d(self.obj,j,pointer(jLoc))
    elif self.tag == cTag: lib.ElDistMatrixLocalCol_c(self.obj,j,pointer(jLoc))
    elif self.tag == zTag: lib.ElDistMatrixLocalCol_z(self.obj,j,pointer(jLoc))
    else: DataExcept()
    return jLoc.value
  def LocalRowOffset(self,i):
    iLoc = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixLocalRowOffset_i(self.obj,i,pointer(iLoc))
    elif self.tag == sTag:
      lib.ElDistMatrixLocalRowOffset_s(self.obj,i,pointer(iLoc))
    elif self.tag == dTag:
      lib.ElDistMatrixLocalRowOffset_d(self.obj,i,pointer(iLoc))
    elif self.tag == cTag:
      lib.ElDistMatrixLocalRowOffset_c(self.obj,i,pointer(iLoc))
    elif self.tag == zTag:
      lib.ElDistMatrixLocalRowOffset_z(self.obj,i,pointer(iLoc))
    else: DataExcept()
    return iLoc.value
  def LocalColOffset(self,j):
    jLoc = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixLocalColOffset_i(self.obj,j,pointer(jLoc))
    elif self.tag == sTag:
      lib.ElDistMatrixLocalColOffset_s(self.obj,j,pointer(jLoc))
    elif self.tag == dTag:
      lib.ElDistMatrixLocalColOffset_d(self.obj,j,pointer(jLoc))
    elif self.tag == cTag:
      lib.ElDistMatrixLocalColOffset_c(self.obj,j,pointer(jLoc))
    elif self.tag == zTag:
      lib.ElDistMatrixLocalColOffset_z(self.obj,j,pointer(jLoc))
    else: DataExcept()
    return jLoc.value
  def GlobalRow(self,iLoc):
    i = iType()
    if   self.tag == iTag: lib.ElDistMatrixGlobalRow_i(self.obj,iLoc,pointer(i))
    elif self.tag == sTag: lib.ElDistMatrixGlobalRow_s(self.obj,iLoc,pointer(i))
    elif self.tag == dTag: lib.ElDistMatrixGlobalRow_d(self.obj,iLoc,pointer(i))
    elif self.tag == cTag: lib.ElDistMatrixGlobalRow_c(self.obj,iLoc,pointer(i))
    elif self.tag == zTag: lib.ElDistMatrixGlobalRow_z(self.obj,iLoc,pointer(i))
    else: DataExcept()
    return i.value
  def GlobalCol(self,jLoc):
    j = jType()
    if   self.tag == iTag: lib.ElDistMatrixGlobalCol_i(self.obj,jLoc,pointer(j))
    elif self.tag == sTag: lib.ElDistMatrixGlobalCol_s(self.obj,jLoc,pointer(j))
    elif self.tag == dTag: lib.ElDistMatrixGlobalCol_d(self.obj,jLoc,pointer(j))
    elif self.tag == cTag: lib.ElDistMatrixGlobalCol_c(self.obj,jLoc,pointer(j))
    elif self.tag == zTag: lib.ElDistMatrixGlobalCol_z(self.obj,jLoc,pointer(j))
    else: DataExcept()
    return j.value
  def IsLocalRow(self,i):
    isLocal = bType()
    if   self.tag == iTag: 
      lib.ElDistMatrixIsLocalRow_i(self.obj,i,pointer(isLoc))
    elif self.tag == sTag:
      lib.ElDistMatrixIsLocalRow_s(self.obj,i,pointer(isLoc))
    elif self.tag == dTag:
      lib.ElDistMatrixIsLocalRow_d(self.obj,i,pointer(isLoc))
    elif self.tag == cTag:
      lib.ElDistMatrixIsLocalRow_c(self.obj,i,pointer(isLoc))
    elif self.tag == zTag:
      lib.ElDistMatrixIsLocalRow_z(self.obj,i,pointer(isLoc))
    else: DataExcept()
    return isLocal.value
  def IsLocalCol(self,j):
    isLocal = bType()
    if   self.tag == iTag: 
      lib.ElDistMatrixIsLocalCol_i(self.obj,j,pointer(isLoc))
    elif self.tag == sTag:
      lib.ElDistMatrixIsLocalCol_s(self.obj,j,pointer(isLoc))
    elif self.tag == dTag:
      lib.ElDistMatrixIsLocalCol_d(self.obj,j,pointer(isLoc))
    elif self.tag == cTag:
      lib.ElDistMatrixIsLocalCol_c(self.obj,j,pointer(isLoc))
    elif self.tag == zTag:
      lib.ElDistMatrixIsLocalCol_z(self.obj,j,pointer(isLoc))
    else: DataExcept()
    return isLocal.value
  def IsLocal(self,i,j):
    isLocal = bType()
    if   self.tag == iTag: 
      lib.ElDistMatrixIsLocal_i(self.obj,i,j,pointer(isLocal))
    elif self.tag == sTag:
      lib.ElDistMatrixIsLocal_s(self.obj,i,j,pointer(isLocal))
    elif self.tag == dTag:
      lib.ElDistMatrixIsLocal_d(self.obj,i,j,pointer(isLocal))
    elif self.tag == cTag:
      lib.ElDistMatrixIsLocal_c(self.obj,i,j,pointer(isLocal))
    elif self.tag == zTag:
      lib.ElDistMatrixIsLocal_z(self.obj,i,j,pointer(isLocal))
    else: DataExcept()
    return isLocal.value
  def GetDistData(self):
    distData = DistData()
    if   self.tag == iTag: 
      lib.ElDistMatrixDistData_i(self.obj,pointer(distData))
    elif self.tag == sTag:
      lib.ElDistMatrixDistData_s(self.obj,pointer(distData))
    elif self.tag == dTag:
      lib.ElDistMatrixDistData_d(self.obj,pointer(distData))
    elif self.tag == cTag:
      lib.ElDistMatrixDistData_c(self.obj,pointer(distData))
    elif self.tag == zTag:
      lib.ElDistMatrixDistData_z(self.obj,pointer(distData))
    else: DataExcept()
    return distData
  def DistComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: lib.ElDistMatrixDistComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: lib.ElDistMatrixDistComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: lib.ElDistMatrixDistComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: lib.ElDistMatrixDistComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: lib.ElDistMatrixDistComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def CrossComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: lib.ElDistMatrixCrossComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: lib.ElDistMatrixCrossComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: lib.ElDistMatrixCrossComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: lib.ElDistMatrixCrossComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: lib.ElDistMatrixCrossComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def RedundantComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: 
      lib.ElDistMatrixRedundantComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: 
      lib.ElDistMatrixRedundantComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: 
      lib.ElDistMatrixRedundantComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: 
      lib.ElDistMatrixRedundantComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: 
      lib.ElDistMatrixRedundantComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def ColComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: lib.ElDistMatrixColComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: lib.ElDistMatrixColComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: lib.ElDistMatrixColComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: lib.ElDistMatrixColComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: lib.ElDistMatrixColComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def RowComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: lib.ElDistMatrixRowComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: lib.ElDistMatrixRowComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: lib.ElDistMatrixRowComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: lib.ElDistMatrixRowComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: lib.ElDistMatrixRowComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def PartialColComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialColComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialColComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialColComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialColComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialColComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def PartialRowComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialRowComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialRowComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialRowComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialRowComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialRowComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def PartialUnionColComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialUnionColComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialUnionColComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialUnionColComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialUnionColComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialUnionColComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def PartialUnionRowComm(self):
    comm = mpi.Comm()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialUnionRowComm_i(self.obj,pointer(comm))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialUnionRowComm_s(self.obj,pointer(comm))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialUnionRowComm_d(self.obj,pointer(comm))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialUnionRowComm_c(self.obj,pointer(comm))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialUnionRowComm_z(self.obj,pointer(comm))
    else: DataExcept()
    return comm
  def ColStride(self):
    stride = iType()
    if   self.tag == iTag: lib.ElDistMatrixColStride_i(self.obj,pointer(stride))
    elif self.tag == sTag: lib.ElDistMatrixColStride_s(self.obj,pointer(stride))
    elif self.tag == dTag: lib.ElDistMatrixColStride_d(self.obj,pointer(stride))
    elif self.tag == cTag: lib.ElDistMatrixColStride_c(self.obj,pointer(stride))
    elif self.tag == zTag: lib.ElDistMatrixColStride_z(self.obj,pointer(stride))
    else: DataExcept()
    return stride.value
  def RowStride(self):
    stride = iType()
    if   self.tag == iTag: lib.ElDistMatrixRowStride_i(self.obj,pointer(stride))
    elif self.tag == sTag: lib.ElDistMatrixRowStride_s(self.obj,pointer(stride))
    elif self.tag == dTag: lib.ElDistMatrixRowStride_d(self.obj,pointer(stride))
    elif self.tag == cTag: lib.ElDistMatrixRowStride_c(self.obj,pointer(stride))
    elif self.tag == zTag: lib.ElDistMatrixRowStride_z(self.obj,pointer(stride))
    else: DataExcept()
    return stride.value
  def PartialColStride(self):
    stride = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialColStride_i(self.obj,pointer(stride))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialColStride_s(self.obj,pointer(stride))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialColStride_d(self.obj,pointer(stride))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialColStride_c(self.obj,pointer(stride))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialColStride_z(self.obj,pointer(stride))
    else: DataExcept()
    return stride.value
  def PartialRowStride(self):
    stride = iType()
    if   self.tag == iTag: 
      lib.ElDistMatrixPartialRowStride_i(self.obj,pointer(stride))
    elif self.tag == sTag: 
      lib.ElDistMatrixPartialRowStride_s(self.obj,pointer(stride))
    elif self.tag == dTag: 
      lib.ElDistMatrixPartialRowStride_d(self.obj,pointer(stride))
    elif self.tag == cTag: 
      lib.ElDistMatrixPartialRowStride_c(self.obj,pointer(stride))
    elif self.tag == zTag: 
      lib.ElDistMatrixPartialRowStride_z(self.obj,pointer(stride))
    else: DataExcept()
    return stride.value
  # Left off moving upward from HERE
  def PartialUnionColStride(self):
    stride = iType()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionColStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionColStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionColStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionColStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionColStride_z(*args)
    else: DataExcept()
    return stride.value
  def PartialUnionRowStride(self):
    stride = iType()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionRowStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionRowStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionRowStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionRowStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionRowStride_z(*args)
    else: DataExcept()
    return stride.value
  def DistSize(self):
    size = iType()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElDistMatrixDistSize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDistSize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDistSize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDistSize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDistSize_z(*args)
    else: DataExcept()
    return size.value
  def CrossSize(self):
    size = iType()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElDistMatrixCrossSize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixCrossSize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixCrossSize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixCrossSize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixCrossSize_z(*args)
    else: DataExcept()
    return size.value
  def RedundantSize(self):
    size = iType()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElDistMatrixRedundantSize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRedundantSize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRedundantSize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRedundantSize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRedundantSize_z(*args)
    else: DataExcept()
    return size.value
  def Get(self,i,j):
    value = TagToType(self.tag)()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElDistMatrixGet_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGet_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGet_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGet_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGet_z(*args)
    else: DataExcept()
    return value
  def GetRealPart(self,i,j):
    value = TagToType(Base(self.tag))()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElDistMatrixGet_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGet_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGet_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetRealPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetRealPart_z(*args)
    else: DataExcept()
    return value
  def GetImagPart(self,i,j):
    if   self.tag == iTag: return iType(0)
    elif self.tag == sTag: return sType(0)
    elif self.tag == dTag: return dType(0)
    elif self.tag == cTag:
      value = sType()
      lib.ElDistMatrixGetRealPart_c(self.obj,i,j,pointer(value))
      return value
    elif self.tag == zTag:
      value = dType()
      lib.ElDistMatrixGetRealPart_z(self.obj,i,j,pointer(value))
      return value
    else: DataExcept()
  def Set(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElDistMatrixSet_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSet_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSet_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSet_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSet_z(*args)
    else: DataExcept()
  def SetRealPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixSetRealPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetRealPart_z(*args)
    else: self.Set(i,j,value)
  def SetImagPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixSetImagPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetImagPart_z(*args)
    else: raise Exception('Cannot set imaginary part of a real datatype')
  def Update(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElDistMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdate_z(*args)
    else: DataExcept()
  def UpdateRealPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixUpdateRealPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateRealPart_z(*args)
    else: self.Update(i,j,value)
  def UpdateImagPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixUpdateImagPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateImagPart_z(*args)
    else: raise Exception('Cannot update imaginary part of a real datatype')
  def MakeReal(self,i,j):
    args = [self.obj,i,j]
    if   self.tag == cTag: lib.ElDistMatrixMakeReal_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeReal_z(*args)
  def Conjugate(self,i,j):
    args = [self.obj,i,j]
    if   self.tag == cTag: lib.ElDistMatrixConjugate_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixConjugate_z(*args)
  def GetLocal(self,iLoc,jLoc): 
    return self.LockedMatrix().Get(iLoc,jLoc)
  def GetLocalRealPart(self,iLoc,jLoc): 
    return self.LockedMatrix().GetRealPart(iLoc,jLoc)
  def GetLocalImagPart(self,iLoc,jLoc):
    return self.LockedMatrix().GetImagPart(iLoc,jLoc)
  def SetLocal(self,iLoc,jLoc,value):
    self.Matrix().Set(iLoc,jLoc,value)
  def SetLocalRealPart(self,iLoc,jLoc,value):
    self.Matrix().SetRealPart(iLoc,jLoc,value)
  def SetLocalImagPart(self,iLoc,jLoc,value):
    self.Matrix().SetImagPart(iLoc,jLoc,value)
  def UpdateLocal(self,iLoc,jLoc,value):
    self.Matrix().Update(iLoc,jLoc,value)
  def UpdateLocalRealPart(self,iLoc,jLoc,value):
    self.Matrix().UpdateRealPart(iLoc,jLoc,value)
  def UpdateLocalImagPart(self,iLoc,jLoc,value):
    self.Matrix().UpdateImagPart(iLoc,jLoc,value)
  def DiagonalAlignedWith(distData,offset=0):
    aligned = bType()
    args = [self.obj,distData,offset,pointer(aligned)]
    if   self.tag == iTag: lib.ElDistMatrixDiagonalAlignedWith_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDiagonalAlignedWith_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDiagonalAlignedWith_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDiagonalAlignedWith_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDiagonalAlignedWith_z(*args)
    else: DataExcept()
    return aligned.value
  def DiagonalRoot(self,offset=0):
    root = iType()
    args = [self.obj,offset,pointer(root)]
    if   self.tag == iTag: lib.ElDistMatrixDiagonalRoot_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDiagonalRoot_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDiagonalRoot_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDiagonalRoot_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDiagonalRoot_z(*args)
    else: DataExcept()
    return root.value
  def DiagonalAlign(self,offset=0):
    align = iType()
    args = [self.obj,offset,pointer(align)]
    if   self.tag == iTag: lib.ElDistMatrixDiagonalAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDiagonalAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDiagonalAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDiagonalAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDiagonalAlign_z(*args)
    else: DataExcept()
    return align.value
  def GetSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    ASub = DistMatrix(self.tag,STAR,STAR,self.Grid())
    args = [self.obj,numRowInds,rowInd,numColInds,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixGetSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGetSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGetSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetSubmatrix_z(*args)
    else: DataExcept()
    return ASub
  def GetRealPartOfSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    ASub = DistMatrix(Base(self.tag),STAR,STAR,self.Grid())
    args = [self.obj,numRowInds,rowInd,numColInds,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixGetSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGetSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGetSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetRealPartOfSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetRealPartOfSubmatrix_z(*args)
    else: DataExcept()
    return ASub
  def GetImagPartOfSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    ASub = DistMatrix(Base(self.tag),STAR,STAR,self.Grid())
    args = [self.obj,numRowInds,rowInd,numColInds,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixGetImagPartOfSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGetImagPartOfSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGetImagPartOfSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetImagPartOfSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetImagPartOfSubmatrix_z(*args)
    else: DataExcept()
    return ASub
  def SetSubmatrix(self,I,J,ASub):
    if type(ASub) is not DistMatrix: 
      raise Exception('ASub must be a DistMatrix')
    if ASub.tag != self.tag: raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixSetSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetSubmatrix_z(*args)
    else: DataExcept()
  def SetRealPartOfSubmatrix(self,I,J,ASub):
    if type(ASub) is not DistMatrix: 
      raise Exception('ASub must be a DistMatrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixSetSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetRealPartOfSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetRealPartOfSubmatrix_z(*args)
    else: DataExcept()
  def SetImagPartOfSubmatrix(self,I,J,ASub):
    if type(ASub) is not DistMatrix: 
      raise Exception('ASub must be a DistMatrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,ASub.obj]
    if   self.tag == cTag: lib.ElDistMatrixSetImagPartOfSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetImagPartOfSubmatrix_z(*args)
    else: raise Exception('Cannot modify imaginary part of real matrix')
  def UpdateSubmatrix(self,I,J,alphaPre,ASub):
    alpha = TagToType(self.tag)(alphaPre)
    if type(ASub) is not DistMatrix: 
      raise Exception('ASub must be a DistMatrix')
    if ASub.tag != self.tag: raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,alpha,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixUpdateSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixUpdateSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixUpdateSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixUpdateSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateSubmatrix_z(*args)
    else: DataExcept()
  def UpdateRealPartOfSubmatrix(self,I,J,alpha,ASub):
    if type(ASub) is not DistMatrix: 
      raise Exception('ASub must be a DistMatrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,alpha,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixUpdateSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixUpdateSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixUpdateSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixUpdateRealPartOfSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateRealPartOfSubmatrix_z(*args)
    else: DataExcept()
  def UpdateImagPartOfSubmatrix(self,I,J,alpha,ASub):
    if type(ASub) is not DistMatrix: 
      raise Exception('ASub must be a DistMatrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,alpha,ASub.obj]
    if   self.tag == cTag: lib.ElDistMatrixUpdateImagPartOfSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateImagPartOfSubmatrix_z(*args)
    else: raise Exception('Cannot modify imaginary part of real matrix')
  def MakeSubmatrixReal(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,numRowInds,rowInd,numColInds,colInd]
    if   self.tag == cTag: lib.ElDistMatrixMakeSubmatrixReal_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeSubmatrixReal_z(*args)
  def ConjugateSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,numRowInds,rowInd,numColInds,colInd]
    if   self.tag == cTag: lib.ElDistMatrixConjugateSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixConjugateSubmatrix_z(*args)
  #     Get rid of 'Local' member functions?!?
  def GetLocalSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    ASub = Matrix(self.tag)
    args = [self.obj,numRowInds,rowInd,numColInds,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixGetLocalSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGetLocalSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGetLocalSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetLocalSubmatrix_z(*args)
    else: DataExcept()
    return ASub
  def GetRealPartOfLocalSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    ASub = Matrix(Base(self.tag))
    args = [self.obj,numRowInds,rowInd,numColInds,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixGetLocalSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGetLocalSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGetLocalSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetRealPartOfLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetRealPartOfLocalSubmatrix_z(*args)
    else: DataExcept()
    return ASub
  def GetImagPartOfLocalSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    ASub = Matrix(Base(self.tag))
    args = [self.obj,numRowInds,rowInd,numColInds,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixGetImagPartOfLocalSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGetImagPartOfLocalSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGetImagPartOfLocalSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetImagPartOfLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetImagPartOfLocalSubmatrix_z(*args)
    else: DataExcept()
    return ASub
  def SetLocalSubmatrix(self,I,J,ASub):
    if type(ASub) is not Matrix: raise Exception('ASub must be a Matrix')
    if ASub.tag != self.tag: raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixSetLocalSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetLocalSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetLocalSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetLocalSubmatrix_z(*args)
    else: DataExcept()
  def SetRealPartOfLocalSubmatrix(self,I,J,ASub):
    if type(ASub) is not Matrix: raise Exception('ASub must be a Matrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixSetLocalSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetLocalSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetLocalSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetRealPartOfLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetRealPartOfLocalSubmatrix_z(*args)
    else: DataExcept()
  def SetImagPartOfLocalSubmatrix(self,I,J,ASub):
    if type(ASub) is not Matrix: raise Exception('ASub must be a Matrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,ASub.obj]
    if   self.tag == cTag: lib.ElDistMatrixSetImagPartOfLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetImagPartOfLocalSubmatrix_z(*args)
    else: raise Exception('Cannot modify imaginary part of real matrix')
  def UpdateLocalSubmatrix(self,I,J,alphaPre,ASub):
    alpha = TagToType(self.tag)(alphaPre)
    if type(ASub) is not Matrix: raise Exception('ASub must be a Matrix')
    if ASub.tag != self.tag: raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,alpha,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixUpdateLocalSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixUpdateLocalSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixUpdateLocalSubmatrix_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixUpdateLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateLocalSubmatrix_z(*args)
    else: DataExcept()
  def UpdateRealPartOfLocalSubmatrix(self,I,J,alpha,ASub):
    if type(ASub) is not Matrix: raise Exception('ASub must be a Matrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,alpha,ASub.obj]
    if   self.tag == iTag: lib.ElDistMatrixUpdateLocalSubmatrix_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixUpdateLocalSubmatrix_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixUpdateLocalSubmatrix_d(*args)
    elif self.tag == cTag:
      lib.ElDistMatrixUpdateRealPartOfLocalSubmatrix_c(*args)
    elif self.tag == zTag:
      lib.ElDistMatrixUpdateRealPartOfLocalSubmatrix_z(*args)
    else: DataExcept()
  def UpdateImagPartOfLocalSubmatrix(self,I,J,alpha,ASub):
    if type(ASub) is not Matrix: raise Exception('ASub must be a Matrix')
    if ASub.tag != Base(self.tag): raise Exception('Datatypes must conform')
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,rowInd,colInd,alpha,ASub.obj]
    if   self.tag == cTag:
      lib.ElDistMatrixUpdateImagPartOfLocalSubmatrix_c(*args)
    elif self.tag == zTag:
      lib.ElDistMatrixUpdateImagPartOfLocalSubmatrix_z(*args)
    else: raise Exception('Cannot modify imaginary part of real matrix')
  def MakeLocalSubmatrixReal(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,numRowInds,rowInd,numColInds,colInd]
    if   self.tag == cTag: lib.ElDistMatrixMakeLocalSubmatrixReal_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeLocalSubmatrixReal_z(*args)
  def ConjugateLocalSubmatrix(self,I,J):
    numRowInds = len(I)
    numColInds = len(J)
    rowInd = (iType*numRowInds)(*I)
    colInd = (iType*numColInds)(*J)
    args = [self.obj,numRowInds,rowInd,numColInds,colInd]
    if   self.tag == cTag: lib.ElDistMatrixConjugateLocalSubmatrix_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixConjugateLocalSubmatrix_z(*args)
  def SumOver(self,comm):
    args = [self.obj,comm]
    if   self.tag == iTag: lib.ElDistMatrixSumOver_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSumOver_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSumOver_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSumOver_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSumOver_z(*args)
    else: DataExcept()
  def __getitem__(self,indTup):
    iInd, jInd = indTup
    iRan = IndexRange(iInd)
    jRan = IndexRange(jInd)
    distData = self.GetDistData()
    ASub = DistMatrix(self.tag,distData.colDist,distData.rowDist,self.Grid())
    args = [ASub.obj,self.obj,iRan,jRan]
    if self.Locked():
      if   self.tag == iTag: lib.ElLockedViewDist_i(*args)
      elif self.tag == sTag: lib.ElLockedViewDist_s(*args)
      elif self.tag == dTag: lib.ElLockedViewDist_d(*args)
      elif self.tag == cTag: lib.ElLockedViewDist_c(*args)
      elif self.tag == zTag: lib.ElLockedViewDist_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElViewDist_i(*args)
      elif self.tag == sTag: lib.ElViewDist_s(*args)
      elif self.tag == dTag: lib.ElViewDist_d(*args)
      elif self.tag == cTag: lib.ElViewDist_c(*args)
      elif self.tag == zTag: lib.ElViewDist_z(*args)
      else: DataExcept()
    return ASub
