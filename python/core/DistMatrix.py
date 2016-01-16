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
import ctypes, numpy

import Matrix as M
import Grid as G

class DistData(ctypes.Structure):
  _fields_ = [('colDist',c_uint),
              ('rowDist',c_uint), 
              ('colAlign',c_int),
              ('rowAlign',c_int),
              ('root',c_int),
              ('grid',c_void_p)]

class DistMatrix(object):

  # Create an instance with a particular matrix distribution
  # --------------------------------------------------------
  lib.ElDistMatrixCreateSpecific_i.argtypes = \
  lib.ElDistMatrixCreateSpecific_s.argtypes = \
  lib.ElDistMatrixCreateSpecific_d.argtypes = \
  lib.ElDistMatrixCreateSpecific_c.argtypes = \
  lib.ElDistMatrixCreateSpecific_z.argtypes = \
    [c_uint,c_uint,c_void_p,POINTER(c_void_p)]
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

  # Destroy an instance
  # -------------------
  lib.ElDistMatrixDestroy_i.argtypes = \
  lib.ElDistMatrixDestroy_s.argtypes = \
  lib.ElDistMatrixDestroy_d.argtypes = \
  lib.ElDistMatrixDestroy_c.argtypes = \
  lib.ElDistMatrixDestroy_z.argtypes = \
    [c_void_p]
  def Destroy(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMatrixDestroy_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDestroy_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDestroy_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDestroy_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDestroy_z(*args)
    else: DataExcept()

  # Reset the matrix to the default (empty) state
  # ---------------------------------------------
  lib.ElDistMatrixEmpty_i.argtypes = \
  lib.ElDistMatrixEmpty_s.argtypes = \
  lib.ElDistMatrixEmpty_d.argtypes = \
  lib.ElDistMatrixEmpty_c.argtypes = \
  lib.ElDistMatrixEmpty_z.argtypes = \
    [c_void_p,bType]
  def Empty(self,freeMemory=True):
    args = [self.obj,freeMemory]
    if   self.tag == iTag: lib.ElDistMatrixEmpty_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixEmpty_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixEmpty_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixEmpty_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixEmpty_z(*args)
    else: DataExcept()

  # Resize the matrix to 0 x 0, demand ownership, and delete underlying buffer
  # --------------------------------------------------------------------------
  lib.ElDistMatrixEmptyData_i.argtypes = \
  lib.ElDistMatrixEmptyData_s.argtypes = \
  lib.ElDistMatrixEmptyData_d.argtypes = \
  lib.ElDistMatrixEmptyData_c.argtypes = \
  lib.ElDistMatrixEmptyData_z.argtypes = \
    [c_void_p,bType]
  def EmptyData(self,freeMemory=True):
    args = [self.obj,freeMemory]
    if   self.tag == iTag: lib.ElDistMatrixEmptyData_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixEmptyData_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixEmptyData_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixEmptyData_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixEmptyData_z(*args)
    else: DataExcept()

  # Specify the process grid of the matrix (if different, resize to 0x0)
  # --------------------------------------------------------------------
  lib.ElDistMatrixSetGrid_i.argtypes = \
  lib.ElDistMatrixSetGrid_s.argtypes = \
  lib.ElDistMatrixSetGrid_d.argtypes = \
  lib.ElDistMatrixSetGrid_c.argtypes = \
  lib.ElDistMatrixSetGrid_z.argtypes = \
    [c_void_p,c_void_p]
  def SetGrid(self,grid):
    args = [self.obj,grid.obj]
    if   self.tag == iTag: lib.ElDistMatrixSetGrid_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetGrid_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetGrid_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetGrid_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetGrid_z(*args)
    else: DataExcept()

  # Resize the matrix to have the specified height and width
  # --------------------------------------------------------
  lib.ElDistMatrixResize_i.argtypes = \
  lib.ElDistMatrixResize_s.argtypes = \
  lib.ElDistMatrixResize_d.argtypes = \
  lib.ElDistMatrixResize_c.argtypes = \
  lib.ElDistMatrixResize_z.argtypes = \
    [c_void_p,iType,iType]
  def Resize(self,m,n):
    args = [self.obj,m,n]
    if   self.tag == iTag: lib.ElDistMatrixResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixResize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixResize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixResize_z(*args)
    else: DataExcept()

  # Specify the leading dimension in addition to simply resizing
  # ------------------------------------------------------------
  lib.ElDistMatrixResizeWithLDim_i.argtypes = \
  lib.ElDistMatrixResizeWithLDim_s.argtypes = \
  lib.ElDistMatrixResizeWithLDim_d.argtypes = \
  lib.ElDistMatrixResizeWithLDim_c.argtypes = \
  lib.ElDistMatrixResizeWithLDim_z.argtypes = \
    [c_void_p,iType,iType,iType]
  def ResizeWithLDim(self,m,n,ldim):
    args = [self.obj,m,n,ldim]
    if   self.tag == iTag: lib.ElDistMatrixResizeWithLDim_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixResizeWithLDim_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixResizeWithLDim_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixResizeWithLDim_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixResizeWithLDim_z(*args)
    else: DataExcept()

  # Force all processes to share matrix metadata with the root process
  # ------------------------------------------------------------------
  lib.ElDistMatrixMakeConsistent_i.argtypes = \
  lib.ElDistMatrixMakeConsistent_s.argtypes = \
  lib.ElDistMatrixMakeConsistent_d.argtypes = \
  lib.ElDistMatrixMakeConsistent_c.argtypes = \
  lib.ElDistMatrixMakeConsistent_z.argtypes = \
    [c_void_p,bType]
  def MakeConsistent(self,incViewers):
    args = [self.obj,incViewers]
    if   self.tag == iTag: lib.ElDistMatrixMakeConsistent_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixMakeConsistent_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixMakeConsistent_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixMakeConsistent_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeConsistent_z(*args)
    else: DataExcept()

  # Force all processes to set the matrix dimensions to that of the root's
  # ----------------------------------------------------------------------
  lib.ElDistMatrixMakeSizeConsistent_i.argtypes = \
  lib.ElDistMatrixMakeSizeConsistent_s.argtypes = \
  lib.ElDistMatrixMakeSizeConsistent_d.argtypes = \
  lib.ElDistMatrixMakeSizeConsistent_c.argtypes = \
  lib.ElDistMatrixMakeSizeConsistent_z.argtypes = \
    [c_void_p,bType]
  def MakeSizeConsistent(self,incViewers):
    args = [self.obj,incViewers]
    if   self.tag == iTag: lib.ElDistMatrixMakeSizeConsistent_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixMakeSizeConsistent_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixMakeSizeConsistent_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixMakeSizeConsistent_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeSizeConsistent_z(*args)
    else: DataExcept()

  # Set the column and row alignments
  # ---------------------------------
  lib.ElDistMatrixAlign_i.argtypes = \
  lib.ElDistMatrixAlign_s.argtypes = \
  lib.ElDistMatrixAlign_d.argtypes = \
  lib.ElDistMatrixAlign_c.argtypes = \
  lib.ElDistMatrixAlign_z.argtypes = \
    [c_void_p,c_int,c_int,bType]
  def Align(self,colAlign,rowAlign,constrain):
    args = [self.obj,colAlign,rowAlign,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlign_z(*args)
    else: DataExcept()

  # Set the column alignment
  # ------------------------
  lib.ElDistMatrixAlignCols_i.argtypes = \
  lib.ElDistMatrixAlignCols_s.argtypes = \
  lib.ElDistMatrixAlignCols_d.argtypes = \
  lib.ElDistMatrixAlignCols_c.argtypes = \
  lib.ElDistMatrixAlignCols_z.argtypes = \
    [c_void_p,c_int,bType]
  def AlignCols(self,colAlign,constrain):
    args = [self.obj,colAlign,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignCols_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignCols_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignCols_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignCols_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignCols_z(*args)
    else: DataExcept()

  # Set the row alignment
  # ---------------------
  lib.ElDistMatrixAlignRows_i.argtypes = \
  lib.ElDistMatrixAlignRows_s.argtypes = \
  lib.ElDistMatrixAlignRows_d.argtypes = \
  lib.ElDistMatrixAlignRows_c.argtypes = \
  lib.ElDistMatrixAlignRows_z.argtypes = \
    [c_void_p,c_int,bType]
  def AlignRows(self,rowAlign,constrain):
    args = [self.obj,rowAlign,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignRows_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignRows_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignRows_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignRows_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignRows_z(*args)
    else: DataExcept()

  # Unconstrain the row and column alignments for future assignments
  # ----------------------------------------------------------------
  lib.ElDistMatrixFreeAlignments_i.argtypes = \
  lib.ElDistMatrixFreeAlignments_s.argtypes = \
  lib.ElDistMatrixFreeAlignments_d.argtypes = \
  lib.ElDistMatrixFreeAlignments_c.argtypes = \
  lib.ElDistMatrixFreeAlignments_z.argtypes = \
    [c_void_p]
  def FreeAlignments(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMatrixFreeAlignments_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixFreeAlignments_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixFreeAlignments_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixFreeAlignments_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixFreeAlignments_z(*args)
    else: DataExcept()

  # Set the root process within the "cross" communicator
  # ----------------------------------------------------
  lib.ElDistMatrixSetRoot_i.argtypes = \
  lib.ElDistMatrixSetRoot_s.argtypes = \
  lib.ElDistMatrixSetRoot_d.argtypes = \
  lib.ElDistMatrixSetRoot_c.argtypes = \
  lib.ElDistMatrixSetRoot_z.argtypes = \
    [c_void_p,c_int,bType]
  def SetRoot(self,root,constrain=True):
    args = [self.obj,root,constrain]
    if   self.tag == iTag: lib.ElDistMatrixSetRoot_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSetRoot_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSetRoot_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSetRoot_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetRoot_z(*args)
    else: DataExcept()

  # Align this matrix with another distributed matrix as much as possible
  # ---------------------------------------------------------------------
  lib.ElDistMatrixAlignWith_i.argtypes = \
  lib.ElDistMatrixAlignWith_s.argtypes = \
  lib.ElDistMatrixAlignWith_d.argtypes = \
  lib.ElDistMatrixAlignWith_c.argtypes = \
  lib.ElDistMatrixAlignWith_z.argtypes = \
    [c_void_p,DistData,bType]
  def AlignWith(self,distData,constrain):
    args = [self.obj,distData,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignWith_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignWith_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignWith_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignWith_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignWith_z(*args)
    else: DataExcept()

  # Align the column distribution of this matrix with another distributed matrix
  # ----------------------------------------------------------------------------
  lib.ElDistMatrixAlignColsWith_i.argtypes = \
  lib.ElDistMatrixAlignColsWith_s.argtypes = \
  lib.ElDistMatrixAlignColsWith_d.argtypes = \
  lib.ElDistMatrixAlignColsWith_c.argtypes = \
  lib.ElDistMatrixAlignColsWith_z.argtypes = \
    [c_void_p,DistData,bType]
  def AlignColsWith(self,distData,constrain):
    args = [self.obj,distData,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignColsWith_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignColsWith_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignColsWith_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignColsWith_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignColsWith_z(*args)
    else: DataExcept()

  # Align the row distribution of the matrix with another distributed matrix
  # ------------------------------------------------------------------------
  lib.ElDistMatrixAlignRowsWith_i.argtypes = \
  lib.ElDistMatrixAlignRowsWith_s.argtypes = \
  lib.ElDistMatrixAlignRowsWith_d.argtypes = \
  lib.ElDistMatrixAlignRowsWith_c.argtypes = \
  lib.ElDistMatrixAlignRowsWith_z.argtypes = \
    [c_void_p,DistData,bType]
  def AlignRowsWith(self,distData,constrain):
    args = [self.obj,distData,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignRowsWith_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignRowsWith_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAlignRowsWith_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignRowsWith_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignRowsWith_z(*args)
    else: DataExcept()

  # Align this matrix with another and then resize
  # ----------------------------------------------
  lib.ElDistMatrixAlignAndResize_i.argtypes = \
  lib.ElDistMatrixAlignAndResize_s.argtypes = \
  lib.ElDistMatrixAlignAndResize_d.argtypes = \
  lib.ElDistMatrixAlignAndResize_c.argtypes = \
  lib.ElDistMatrixAlignAndResize_z.argtypes = \
    [c_void_p,c_int,c_int,iType,iType,bType,bType]
  def AlignAndResize(self,colAlign,rowAlign,m,n,force,constrain):
    args = [self.obj,colAlign,rowAlign,m,n,force,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignAndResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignAndResize_s(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignAndResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignAndResize_z(*args)
    else: DataExcept()

  # Align the columns with another matrix and then resize
  # -----------------------------------------------------
  lib.ElDistMatrixAlignColsAndResize_i.argtypes = \
  lib.ElDistMatrixAlignColsAndResize_s.argtypes = \
  lib.ElDistMatrixAlignColsAndResize_d.argtypes = \
  lib.ElDistMatrixAlignColsAndResize_c.argtypes = \
  lib.ElDistMatrixAlignColsAndResize_z.argtypes = \
    [c_void_p,c_int,iType,iType,bType,bType]
  def AlignColsAndResize(self,colAlign,m,n,force,constrain):
    args = [self.obj,colAlign,m,n,force,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignColsAndResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignColsAndResize_s(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignColsAndResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignColsAndResize_z(*args)
    else: DataExcept()

  # Align the rows with another matrix and then resize
  # --------------------------------------------------
  lib.ElDistMatrixAlignRowsAndResize_i.argtypes = \
  lib.ElDistMatrixAlignRowsAndResize_s.argtypes = \
  lib.ElDistMatrixAlignRowsAndResize_d.argtypes = \
  lib.ElDistMatrixAlignRowsAndResize_c.argtypes = \
  lib.ElDistMatrixAlignRowsAndResize_z.argtypes = \
    [c_void_p,c_int,iType,iType,bType,bType]
  def AlignRowsAndResize(self,rowAlign,m,n,force,constrain):
    args = [self.obj,rowAlign,m,n,force,constrain]
    if   self.tag == iTag: lib.ElDistMatrixAlignRowsAndResize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAlignRowsAndResize_s(*args)
    elif self.tag == cTag: lib.ElDistMatrixAlignRowsAndResize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAlignRowsAndResize_z(*args)
    else: DataExcept()

  # Attach a buffer to the matrix
  # -----------------------------
  lib.ElDistMatrixAttach_i.argtypes = \
  lib.ElDistMatrixLockedAttach_i.argtypes = \
    [c_void_p,iType,iType,c_void_p,c_int,c_int,POINTER(iType),iType,c_int]
  lib.ElDistMatrixAttach_s.argtypes = \
  lib.ElDistMatrixLockedAttach_s.argtypes = \
    [c_void_p,iType,iType,c_void_p,c_int,c_int,POINTER(sType),iType,c_int]
  lib.ElDistMatrixAttach_d.argtypes = \
  lib.ElDistMatrixLockedAttach_d.argtypes = \
    [c_void_p,iType,iType,c_void_p,c_int,c_int,POINTER(dType),iType,c_int]
  lib.ElDistMatrixAttach_c.argtypes = \
  lib.ElDistMatrixLockedAttach_c.argtypes = \
    [c_void_p,iType,iType,c_void_p,c_int,c_int,POINTER(cType),iType,c_int]
  lib.ElDistMatrixAttach_z.argtypes = \
  lib.ElDistMatrixLockedAttach_z.argtypes = \
    [c_void_p,iType,iType,c_void_p,c_int,c_int,POINTER(zType),iType,c_int]
  def Attach(self,m,n,grid,colAlign,rowAlign,buf,ldim,root,locked=False):
    args = [self.obj,m,n,grid.obj,colAlign,rowAlign,buf,ldim,root]
    if locked:
      if   self.tag == iTag: lib.ElDistMatrixLockedAttach_i(*args)
      elif self.tag == sTag: lib.ElDistMatrixLockedAttach_s(*args)
      elif self.tag == dTag: lib.ElDistMatrixLockedAttach_d(*args)
      elif self.tag == cTag: lib.ElDistMatrixLockedAttach_c(*args)
      elif self.tag == zTag: lib.ElDistMatrixLockedAttach_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistMatrixAttach_i(*args)
      elif self.tag == sTag: lib.ElDistMatrixAttach_s(*args)
      elif self.tag == dTag: lib.ElDistMatrixAttach_d(*args)
      elif self.tag == cTag: lib.ElDistMatrixAttach_c(*args)
      elif self.tag == zTag: lib.ElDistMatrixAttach_z(*args)
      else: DataExcept()

  # Return the height of the matrix
  # -------------------------------
  lib.ElDistMatrixHeight_i.argtypes = \
  lib.ElDistMatrixHeight_s.argtypes = \
  lib.ElDistMatrixHeight_d.argtypes = \
  lib.ElDistMatrixHeight_c.argtypes = \
  lib.ElDistMatrixHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
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

  # Return the width of the matrix
  # ------------------------------
  lib.ElDistMatrixWidth_i.argtypes = \
  lib.ElDistMatrixWidth_s.argtypes = \
  lib.ElDistMatrixWidth_d.argtypes = \
  lib.ElDistMatrixWidth_c.argtypes = \
  lib.ElDistMatrixWidth_z.argtypes = \
    [c_void_p,POINTER(iType)]
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

  # Return the diagonal length of the matrix
  # ----------------------------------------
  lib.ElDistMatrixDiagonalLength_i.argtypes = \
  lib.ElDistMatrixDiagonalLength_s.argtypes = \
  lib.ElDistMatrixDiagonalLength_d.argtypes = \
  lib.ElDistMatrixDiagonalLength_c.argtypes = \
  lib.ElDistMatrixDiagonalLength_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
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

  # Return whether or not the matrix is viewing a buffer (rather than owning it)
  # ----------------------------------------------------------------------------
  lib.ElDistMatrixViewing_i.argtypes = \
  lib.ElDistMatrixViewing_s.argtypes = \
  lib.ElDistMatrixViewing_d.argtypes = \
  lib.ElDistMatrixViewing_c.argtypes = \
  lib.ElDistMatrixViewing_z.argtypes = \
    [c_void_p,POINTER(bType)]
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

  #  Return whether or not the matrix is "locked"/immutable
  # -------------------------------------------------------
  lib.ElDistMatrixLocked_i.argtypes = \
  lib.ElDistMatrixLocked_s.argtypes = \
  lib.ElDistMatrixLocked_d.argtypes = \
  lib.ElDistMatrixLocked_c.argtypes = \
  lib.ElDistMatrixLocked_z.argtypes = \
    [c_void_p,POINTER(bType)]
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

  # Return the local height of the matrix
  # -------------------------------------
  lib.ElDistMatrixLocalHeight_i.argtypes = \
  lib.ElDistMatrixLocalHeight_s.argtypes = \
  lib.ElDistMatrixLocalHeight_d.argtypes = \
  lib.ElDistMatrixLocalHeight_c.argtypes = \
  lib.ElDistMatrixLocalHeight_z.argtypes = \
    [c_void_p,POINTER(iType)]
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

  # Return the local width of the matrix
  # ------------------------------------
  lib.ElDistMatrixLocalWidth_i.argtypes = \
  lib.ElDistMatrixLocalWidth_s.argtypes = \
  lib.ElDistMatrixLocalWidth_d.argtypes = \
  lib.ElDistMatrixLocalWidth_c.argtypes = \
  lib.ElDistMatrixLocalWidth_z.argtypes = \
    [c_void_p,POINTER(iType)]
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

  # Return the leading dimension of the local matrix
  # ------------------------------------------------
  lib.ElDistMatrixLDim_i.argtypes = \
  lib.ElDistMatrixLDim_s.argtypes = \
  lib.ElDistMatrixLDim_d.argtypes = \
  lib.ElDistMatrixLDim_c.argtypes = \
  lib.ElDistMatrixLDim_z.argtypes = \
    [c_void_p,POINTER(iType)]
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

  # Return the local matrix
  # -----------------------
  lib.ElDistMatrixMatrix_i.argtypes = \
  lib.ElDistMatrixMatrix_s.argtypes = \
  lib.ElDistMatrixMatrix_d.argtypes = \
  lib.ElDistMatrixMatrix_c.argtypes = \
  lib.ElDistMatrixMatrix_z.argtypes = \
  lib.ElDistMatrixLockedMatrix_i.argtypes = \
  lib.ElDistMatrixLockedMatrix_s.argtypes = \
  lib.ElDistMatrixLockedMatrix_d.argtypes = \
  lib.ElDistMatrixLockedMatrix_c.argtypes = \
  lib.ElDistMatrixLockedMatrix_z.argtypes = \
    [c_void_p,POINTER(c_void_p)]
  def Matrix(self,locked=False):
    A = M.Matrix(self.tag,False)
    args = [self.obj,pointer(A.obj)]
    if locked:
      if   self.tag == iTag: lib.ElDistMatrixLockedMatrix_i(*args)
      elif self.tag == sTag: lib.ElDistMatrixLockedMatrix_s(*args)
      elif self.tag == dTag: lib.ElDistMatrixLockedMatrix_d(*args)
      elif self.tag == cTag: lib.ElDistMatrixLockedMatrix_c(*args)
      elif self.tag == zTag: lib.ElDistMatrixLockedMatrix_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistMatrixMatrix_i(*args)
      elif self.tag == sTag: lib.ElDistMatrixMatrix_s(*args)
      elif self.tag == dTag: lib.ElDistMatrixMatrix_d(*args)
      elif self.tag == cTag: lib.ElDistMatrixMatrix_c(*args)
      elif self.tag == zTag: lib.ElDistMatrixMatrix_z(*args)
      else: DataExcept()
    return A

  # Return the amount of locally allocated memory
  # ---------------------------------------------
  lib.ElDistMatrixAllocatedMemory_i.argtypes = \
  lib.ElDistMatrixAllocatedMemory_s.argtypes = \
  lib.ElDistMatrixAllocatedMemory_d.argtypes = \
  lib.ElDistMatrixAllocatedMemory_c.argtypes = \
  lib.ElDistMatrixAllocatedMemory_z.argtypes = \
    [c_void_p,POINTER(c_size_t)]
  def AllocatedMemory(self):
    allocMem = c_size_t()
    args = [self.obj,pointer(allocMem)]
    if   self.tag == iTag: lib.ElDistMatrixAllocatedMemory_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixAllocatedMemory_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixAllocatedMemory_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixAllocatedMemory_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixAllocatedMemory_z(*args)
    else: DataExcept()

  # Return the local buffer
  # ----------------------
  lib.ElDistMatrixBuffer_i.argtypes = \
  lib.ElDistMatrixLockedBuffer_i.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  lib.ElDistMatrixBuffer_s.argtypes = \
  lib.ElDistMatrixLockedBuffer_s.argtypes = \
    [c_void_p,POINTER(POINTER(sType))]
  lib.ElDistMatrixBuffer_d.argtypes = \
  lib.ElDistMatrixLockedBuffer_d.argtypes = \
    [c_void_p,POINTER(POINTER(dType))]
  lib.ElDistMatrixBuffer_c.argtypes = \
  lib.ElDistMatrixLockedBuffer_c.argtypes = \
    [c_void_p,POINTER(POINTER(cType))]
  lib.ElDistMatrixBuffer_z.argtypes = \
  lib.ElDistMatrixLockedBuffer_z.argtypes = \
    [c_void_p,POINTER(POINTER(zType))]
  def Buffer(self,locked=False):
    buf = POINTER(TagToType(self.tag))()
    args = [self.obj,pointer(buf)]
    if locked:
      if   self.tag == iTag: lib.ElDistMatrixLockedBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistMatrixLockedBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistMatrixLockedBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistMatrixLockedBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistMatrixLockedBuffer_z(*args)
      else: DataExcept()
    else:
      if   self.tag == iTag: lib.ElDistMatrixBuffer_i(*args)
      elif self.tag == sTag: lib.ElDistMatrixBuffer_s(*args)
      elif self.tag == dTag: lib.ElDistMatrixBuffer_d(*args)
      elif self.tag == cTag: lib.ElDistMatrixBuffer_c(*args)
      elif self.tag == zTag: lib.ElDistMatrixBuffer_z(*args)
      else: DataExcept()
    return buf

  # Return the process grid
  # -----------------------
  lib.ElDistMatrixGrid_i.argtypes = \
  lib.ElDistMatrixGrid_s.argtypes = \
  lib.ElDistMatrixGrid_d.argtypes = \
  lib.ElDistMatrixGrid_c.argtypes = \
  lib.ElDistMatrixGrid_z.argtypes = \
    [c_void_p,POINTER(c_void_p)]
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

  # Return true if the column alignment is pinned at a particular value
  # -------------------------------------------------------------------
  lib.ElDistMatrixColConstrained_i.argtypes = \
  lib.ElDistMatrixColConstrained_s.argtypes = \
  lib.ElDistMatrixColConstrained_d.argtypes = \
  lib.ElDistMatrixColConstrained_c.argtypes = \
  lib.ElDistMatrixColConstrained_z.argtypes = \
    [c_void_p,POINTER(bType)]
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

  # Return true if the row alignment is pinned at a particular value
  # ----------------------------------------------------------------
  lib.ElDistMatrixRowConstrained_i.argtypes = \
  lib.ElDistMatrixRowConstrained_s.argtypes = \
  lib.ElDistMatrixRowConstrained_d.argtypes = \
  lib.ElDistMatrixRowConstrained_c.argtypes = \
  lib.ElDistMatrixRowConstrained_z.argtypes = \
    [c_void_p,POINTER(bType)]
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

  # Return whether or not the root rank is pinned at a particular value
  # -------------------------------------------------------------------
  lib.ElDistMatrixRootConstrained_i.argtypes = \
  lib.ElDistMatrixRootConstrained_s.argtypes = \
  lib.ElDistMatrixRootConstrained_d.argtypes = \
  lib.ElDistMatrixRootConstrained_c.argtypes = \
  lib.ElDistMatrixRootConstrained_z.argtypes = \
    [c_void_p,POINTER(bType)]
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

  # Return the column alignment
  # ---------------------------
  lib.ElDistMatrixColAlign_i.argtypes = \
  lib.ElDistMatrixColAlign_s.argtypes = \
  lib.ElDistMatrixColAlign_d.argtypes = \
  lib.ElDistMatrixColAlign_c.argtypes = \
  lib.ElDistMatrixColAlign_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def ColAlign(self):
    align = c_int()  
    args = [self.obj,pointer(align)]
    if   self.tag == iTag: lib.ElDistMatrixColAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColAlign_z(*args)
    else: DataExcept()
    return align.value

  # Return the row alignment
  # ------------------------
  lib.ElDistMatrixRowAlign_i.argtypes = \
  lib.ElDistMatrixRowAlign_s.argtypes = \
  lib.ElDistMatrixRowAlign_d.argtypes = \
  lib.ElDistMatrixRowAlign_c.argtypes = \
  lib.ElDistMatrixRowAlign_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def RowAlign(self):
    align = c_int()  
    args = [self.obj,pointer(align)]
    if   self.tag == iTag: lib.ElDistMatrixRowAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowAlign_z(*args)
    else: DataExcept()
    return align.value

  # Return the column shift
  # -----------------------
  lib.ElDistMatrixColShift_i.argtypes = \
  lib.ElDistMatrixColShift_s.argtypes = \
  lib.ElDistMatrixColShift_d.argtypes = \
  lib.ElDistMatrixColShift_c.argtypes = \
  lib.ElDistMatrixColShift_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def ColShift(self):
    shift = c_int()
    args = [self.obj,pointer(shift)]
    if   self.tag == iTag: lib.ElDistMatrixColShift_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColShift_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColShift_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColShift_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColShift_z(*args)
    else: DataExcept()
    return shift.value

  # Return the row shift
  # --------------------
  lib.ElDistMatrixRowShift_i.argtypes = \
  lib.ElDistMatrixRowShift_s.argtypes = \
  lib.ElDistMatrixRowShift_d.argtypes = \
  lib.ElDistMatrixRowShift_c.argtypes = \
  lib.ElDistMatrixRowShift_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def RowShift(self):
    shift = c_int()
    args = [self.obj,pointer(shift)]
    if   self.tag == iTag: lib.ElDistMatrixRowShift_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowShift_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowShift_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowShift_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowShift_z(*args)
    else: DataExcept()
    return shift.value

  # Return this process's column rank
  # ---------------------------------
  lib.ElDistMatrixColRank_i.argtypes = \
  lib.ElDistMatrixColRank_s.argtypes = \
  lib.ElDistMatrixColRank_d.argtypes = \
  lib.ElDistMatrixColRank_c.argtypes = \
  lib.ElDistMatrixColRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def ColRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixColRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColRank_z(*args)
    else: DataExcept()
    return rank.value

  # Return this process's row rank
  # ------------------------------
  lib.ElDistMatrixRowRank_i.argtypes = \
  lib.ElDistMatrixRowRank_s.argtypes = \
  lib.ElDistMatrixRowRank_d.argtypes = \
  lib.ElDistMatrixRowRank_c.argtypes = \
  lib.ElDistMatrixRowRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def RowRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixRowRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowRank_z(*args)
    else: DataExcept()
    return rank.value

  # Return this process's partial column rank
  # -----------------------------------------
  lib.ElDistMatrixPartialColRank_i.argtypes = \
  lib.ElDistMatrixPartialColRank_s.argtypes = \
  lib.ElDistMatrixPartialColRank_d.argtypes = \
  lib.ElDistMatrixPartialColRank_c.argtypes = \
  lib.ElDistMatrixPartialColRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialColRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixPartialColRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialColRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialColRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialColRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialColRank_z(*args)
    else: DataExcept()
    return rank.value

  # Return this process's partial row rank
  # --------------------------------------
  lib.ElDistMatrixPartialRowRank_i.argtypes = \
  lib.ElDistMatrixPartialRowRank_s.argtypes = \
  lib.ElDistMatrixPartialRowRank_d.argtypes = \
  lib.ElDistMatrixPartialRowRank_c.argtypes = \
  lib.ElDistMatrixPartialRowRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialRowRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixPartialRowRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialRowRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialRowRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialRowRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialRowRank_z(*args)
    else: DataExcept()
    return rank.value

  # Return this process's partial union column rank
  # -----------------------------------------------
  lib.ElDistMatrixPartialUnionColRank_i.argtypes = \
  lib.ElDistMatrixPartialUnionColRank_s.argtypes = \
  lib.ElDistMatrixPartialUnionColRank_d.argtypes = \
  lib.ElDistMatrixPartialUnionColRank_c.argtypes = \
  lib.ElDistMatrixPartialUnionColRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialUnionColRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionColRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionColRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionColRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionColRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionColRank_z(*args)
    else: DataExcept()
    return rank.value

  # Return this process's partial union row rank
  # --------------------------------------------
  lib.ElDistMatrixPartialUnionRowRank_i.argtypes = \
  lib.ElDistMatrixPartialUnionRowRank_s.argtypes = \
  lib.ElDistMatrixPartialUnionRowRank_d.argtypes = \
  lib.ElDistMatrixPartialUnionRowRank_c.argtypes = \
  lib.ElDistMatrixPartialUnionRowRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialUnionRowRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionRowRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionRowRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionRowRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionRowRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionRowRank_z(*args)
    else: DataExcept()
    return rank.value

  # Return this process's rank in the distribution communicator
  # -----------------------------------------------------------
  lib.ElDistMatrixDistRank_i.argtypes = \
  lib.ElDistMatrixDistRank_s.argtypes = \
  lib.ElDistMatrixDistRank_d.argtypes = \
  lib.ElDistMatrixDistRank_c.argtypes = \
  lib.ElDistMatrixDistRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def DistRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixDistRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDistRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDistRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDistRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDistRank_z(*args)
    else: DataExcept()
    return rank.value

  # HERE
  lib.ElDistMatrixCrossRank_i.argtypes = \
  lib.ElDistMatrixCrossRank_s.argtypes = \
  lib.ElDistMatrixCrossRank_d.argtypes = \
  lib.ElDistMatrixCrossRank_c.argtypes = \
  lib.ElDistMatrixCrossRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def CrossRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixCrossRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixCrossRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixCrossRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixCrossRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixCrossRank_z(*args)
    else: DataExcept()
    return rank.value

  lib.ElDistMatrixRedundantRank_i.argtypes = \
  lib.ElDistMatrixRedundantRank_s.argtypes = \
  lib.ElDistMatrixRedundantRank_d.argtypes = \
  lib.ElDistMatrixRedundantRank_c.argtypes = \
  lib.ElDistMatrixRedundantRank_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def RedundantRank(self):
    rank = c_int()
    args = [self.obj,pointer(rank)]
    if   self.tag == iTag: lib.ElDistMatrixRedundantRank_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRedundantRank_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRedundantRank_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRedundantRank_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRedundantRank_z(*args)
    else: DataExcept()
    return rank.value

  lib.ElDistMatrixRoot_i.argtypes = \
  lib.ElDistMatrixRoot_s.argtypes = \
  lib.ElDistMatrixRoot_d.argtypes = \
  lib.ElDistMatrixRoot_c.argtypes = \
  lib.ElDistMatrixRoot_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def Root(self):
    root = c_int()
    args = [self.obj,pointer(root)]
    if   self.tag == iTag: lib.ElDistMatrixRoot_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRoot_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRoot_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRoot_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRoot_z(*args)
    else: DataExcept()
    return root.value

  lib.ElDistMatrixParticipating_i.argtypes = \
  lib.ElDistMatrixParticipating_s.argtypes = \
  lib.ElDistMatrixParticipating_d.argtypes = \
  lib.ElDistMatrixParticipating_c.argtypes = \
  lib.ElDistMatrixParticipating_z.argtypes = \
    [c_void_p,POINTER(bType)]
  def Participating(self):
    partic = bType()
    args = [self.obj,pointer(partic)]
    if   self.tag == iTag: lib.ElDistMatrixParticipating_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixParticipating_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixParticipating_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixParticipating_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixParticipating_z(*args)
    else: DataExcept()
    return partic.value

  lib.ElDistMatrixRowOwner_i.argtypes = \
  lib.ElDistMatrixRowOwner_s.argtypes = \
  lib.ElDistMatrixRowOwner_d.argtypes = \
  lib.ElDistMatrixRowOwner_c.argtypes = \
  lib.ElDistMatrixRowOwner_z.argtypes = \
    [c_void_p,iType,POINTER(c_int)]
  def RowOwner(self,i):
    owner = c_int()
    args = [self.obj,i,pointer(owner)]
    if   self.tag == iTag: lib.ElDistMatrixRowOwner_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowOwner_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowOwner_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowOwner_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowOwner_z(*args)
    else: DataExcept()
    return owner.value

  lib.ElDistMatrixColOwner_i.argtypes = \
  lib.ElDistMatrixColOwner_s.argtypes = \
  lib.ElDistMatrixColOwner_d.argtypes = \
  lib.ElDistMatrixColOwner_c.argtypes = \
  lib.ElDistMatrixColOwner_z.argtypes = \
    [c_void_p,iType,POINTER(c_int)]
  def ColOwner(self,j):
    owner = c_int()
    args = [self.obj,j,pointer(owner)]
    if   self.tag == iTag: lib.ElDistMatrixColOwner_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColOwner_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColOwner_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColOwner_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColOwner_z(*args)
    else: DataExcept()
    return owner.value

  lib.ElDistMatrixOwner_i.argtypes = \
  lib.ElDistMatrixOwner_s.argtypes = \
  lib.ElDistMatrixOwner_d.argtypes = \
  lib.ElDistMatrixOwner_c.argtypes = \
  lib.ElDistMatrixOwner_z.argtypes = \
    [c_void_p,iType,iType,POINTER(c_int)]
  def Owner(self,i,j):
    owner = c_int()
    args = [self.obj,i,j,pointer(owner)]
    if   self.tag == iTag: lib.ElDistMatrixOwner_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixOwner_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixOwner_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixOwner_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixOwner_z(*args)
    else: DataExcept()
    return owner.value

  lib.ElDistMatrixLocalRow_i.argtypes = \
  lib.ElDistMatrixLocalRow_s.argtypes = \
  lib.ElDistMatrixLocalRow_d.argtypes = \
  lib.ElDistMatrixLocalRow_c.argtypes = \
  lib.ElDistMatrixLocalRow_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def LocalRow(self,i):
    iLoc = iType()
    args = [self.obj,i,pointer(iLoc)]
    if   self.tag == iTag: lib.ElDistMatrixLocalRow_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLocalRow_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLocalRow_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLocalRow_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLocalRow_z(*args)
    else: DataExcept()
    return iLoc.value

  lib.ElDistMatrixLocalCol_i.argtypes = \
  lib.ElDistMatrixLocalCol_s.argtypes = \
  lib.ElDistMatrixLocalCol_d.argtypes = \
  lib.ElDistMatrixLocalCol_c.argtypes = \
  lib.ElDistMatrixLocalCol_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def LocalCol(self,j):
    jLoc = iType()
    args = [self.obj,j,pointer(jLoc)]
    if   self.tag == iTag: lib.ElDistMatrixLocalCol_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLocalCol_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLocalCol_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLocalCol_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLocalCol_z(*args)
    else: DataExcept()
    return jLoc.value

  lib.ElDistMatrixLocalRowOffset_i.argtypes = \
  lib.ElDistMatrixLocalRowOffset_s.argtypes = \
  lib.ElDistMatrixLocalRowOffset_d.argtypes = \
  lib.ElDistMatrixLocalRowOffset_c.argtypes = \
  lib.ElDistMatrixLocalRowOffset_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def LocalRowOffset(self,i):
    iLoc = iType()
    args = [self.obj,i,pointer(iLoc)]
    if   self.tag == iTag: lib.ElDistMatrixLocalRowOffset_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLocalRowOffset_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLocalRowOffset_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLocalRowOffset_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLocalRowOffset_z(*args)
    else: DataExcept()
    return iLoc.value

  lib.ElDistMatrixLocalColOffset_i.argtypes = \
  lib.ElDistMatrixLocalColOffset_s.argtypes = \
  lib.ElDistMatrixLocalColOffset_d.argtypes = \
  lib.ElDistMatrixLocalColOffset_c.argtypes = \
  lib.ElDistMatrixLocalColOffset_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def LocalColOffset(self,j):
    jLoc = iType()
    args = [self.obj,j,pointer(jLoc)]
    if   self.tag == iTag: lib.ElDistMatrixLocalColOffset_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixLocalColOffset_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixLocalColOffset_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixLocalColOffset_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixLocalColOffset_z(*args)
    else: DataExcept()
    return jLoc.value

  lib.ElDistMatrixGlobalRow_i.argtypes = \
  lib.ElDistMatrixGlobalRow_s.argtypes = \
  lib.ElDistMatrixGlobalRow_d.argtypes = \
  lib.ElDistMatrixGlobalRow_c.argtypes = \
  lib.ElDistMatrixGlobalRow_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def GlobalRow(self,iLoc):
    i = iType()
    args = [self.obj,iLoc,pointer(i)]
    if   self.tag == iTag: lib.ElDistMatrixGlobalRow_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGlobalRow_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGlobalRow_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGlobalRow_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGlobalRow_z(*args)
    else: DataExcept()
    return i.value

  lib.ElDistMatrixGlobalCol_i.argtypes = \
  lib.ElDistMatrixGlobalCol_s.argtypes = \
  lib.ElDistMatrixGlobalCol_d.argtypes = \
  lib.ElDistMatrixGlobalCol_c.argtypes = \
  lib.ElDistMatrixGlobalCol_z.argtypes = \
    [c_void_p,iType,POINTER(iType)]
  def GlobalCol(self,jLoc):
    j = iType()
    args = [self.obj,jLoc,pointer(j)]
    if   self.tag == iTag: lib.ElDistMatrixGlobalCol_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGlobalCol_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGlobalCol_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGlobalCol_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGlobalCol_z(*args)
    else: DataExcept()
    return j.value

  lib.ElDistMatrixIsLocalRow_i.argtypes = \
  lib.ElDistMatrixIsLocalRow_s.argtypes = \
  lib.ElDistMatrixIsLocalRow_d.argtypes = \
  lib.ElDistMatrixIsLocalRow_c.argtypes = \
  lib.ElDistMatrixIsLocalRow_z.argtypes = \
    [c_void_p,iType,POINTER(bType)]
  def IsLocalRow(self,i):
    isLocal = bType()
    args = [self.obj,i,pointer(isLoc)]
    if   self.tag == iTag: lib.ElDistMatrixIsLocalRow_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixIsLocalRow_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixIsLocalRow_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixIsLocalRow_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixIsLocalRow_z(*args)
    else: DataExcept()
    return isLocal.value

  lib.ElDistMatrixIsLocalCol_i.argtypes = \
  lib.ElDistMatrixIsLocalCol_s.argtypes = \
  lib.ElDistMatrixIsLocalCol_d.argtypes = \
  lib.ElDistMatrixIsLocalCol_c.argtypes = \
  lib.ElDistMatrixIsLocalCol_z.argtypes = \
    [c_void_p,iType,POINTER(bType)]
  def IsLocalCol(self,j):
    isLocal = bType()
    args = [self.obj,j,pointer(isLoc)]
    if   self.tag == iTag: lib.ElDistMatrixIsLocalCol_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixIsLocalCol_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixIsLocalCol_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixIsLocalCol_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixIsLocalCol_z(*args)
    else: DataExcept()
    return isLocal.value

  lib.ElDistMatrixIsLocal_i.argtypes = \
  lib.ElDistMatrixIsLocal_s.argtypes = \
  lib.ElDistMatrixIsLocal_d.argtypes = \
  lib.ElDistMatrixIsLocal_c.argtypes = \
  lib.ElDistMatrixIsLocal_z.argtypes = \
    [c_void_p,iType,iType,POINTER(bType)]
  def IsLocal(self,i,j):
    isLocal = bType()
    args = [self.obj,i,j,pointer(isLocal)]
    if   self.tag == iTag: lib.ElDistMatrixIsLocal_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixIsLocal_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixIsLocal_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixIsLocal_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixIsLocal_z(*args)
    else: DataExcept()
    return isLocal.value

  lib.ElDistMatrixDistData_i.argtypes = \
  lib.ElDistMatrixDistData_s.argtypes = \
  lib.ElDistMatrixDistData_d.argtypes = \
  lib.ElDistMatrixDistData_c.argtypes = \
  lib.ElDistMatrixDistData_z.argtypes = \
    [c_void_p,POINTER(DistData)]
  def GetDistData(self):
    distData = DistData()
    args = [self.obj,pointer(distData)]
    if   self.tag == iTag: lib.ElDistMatrixDistData_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDistData_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDistData_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDistData_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDistData_z(*args)
    else: DataExcept()
    return distData

  lib.ElDistMatrixDistComm_i.argtypes = \
  lib.ElDistMatrixDistComm_s.argtypes = \
  lib.ElDistMatrixDistComm_d.argtypes = \
  lib.ElDistMatrixDistComm_c.argtypes = \
  lib.ElDistMatrixDistComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def DistComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixDistComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDistComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDistComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDistComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDistComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixCrossComm_i.argtypes = \
  lib.ElDistMatrixCrossComm_s.argtypes = \
  lib.ElDistMatrixCrossComm_d.argtypes = \
  lib.ElDistMatrixCrossComm_c.argtypes = \
  lib.ElDistMatrixCrossComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def CrossComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixCrossComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixCrossComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixCrossComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixCrossComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixCrossComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixRedundantComm_i.argtypes = \
  lib.ElDistMatrixRedundantComm_s.argtypes = \
  lib.ElDistMatrixRedundantComm_d.argtypes = \
  lib.ElDistMatrixRedundantComm_c.argtypes = \
  lib.ElDistMatrixRedundantComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def RedundantComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixRedundantComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRedundantComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRedundantComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRedundantComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRedundantComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixColComm_i.argtypes = \
  lib.ElDistMatrixColComm_s.argtypes = \
  lib.ElDistMatrixColComm_d.argtypes = \
  lib.ElDistMatrixColComm_c.argtypes = \
  lib.ElDistMatrixColComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def ColComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixColComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixRowComm_i.argtypes = \
  lib.ElDistMatrixRowComm_s.argtypes = \
  lib.ElDistMatrixRowComm_d.argtypes = \
  lib.ElDistMatrixRowComm_c.argtypes = \
  lib.ElDistMatrixRowComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def RowComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixRowComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixPartialColComm_i.argtypes = \
  lib.ElDistMatrixPartialColComm_s.argtypes = \
  lib.ElDistMatrixPartialColComm_d.argtypes = \
  lib.ElDistMatrixPartialColComm_c.argtypes = \
  lib.ElDistMatrixPartialColComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def PartialColComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixPartialColComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialColComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialColComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialColComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialColComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixPartialRowComm_i.argtypes = \
  lib.ElDistMatrixPartialRowComm_s.argtypes = \
  lib.ElDistMatrixPartialRowComm_d.argtypes = \
  lib.ElDistMatrixPartialRowComm_c.argtypes = \
  lib.ElDistMatrixPartialRowComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def PartialRowComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixPartialRowComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialRowComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialRowComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialRowComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialRowComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixPartialUnionColComm_i.argtypes = \
  lib.ElDistMatrixPartialUnionColComm_s.argtypes = \
  lib.ElDistMatrixPartialUnionColComm_d.argtypes = \
  lib.ElDistMatrixPartialUnionColComm_c.argtypes = \
  lib.ElDistMatrixPartialUnionColComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def PartialUnionColComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionColComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionColComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionColComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionColComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionColComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixPartialUnionRowComm_i.argtypes = \
  lib.ElDistMatrixPartialUnionRowComm_s.argtypes = \
  lib.ElDistMatrixPartialUnionRowComm_d.argtypes = \
  lib.ElDistMatrixPartialUnionRowComm_c.argtypes = \
  lib.ElDistMatrixPartialUnionRowComm_z.argtypes = \
    [c_void_p,POINTER(mpi.Comm)]
  def PartialUnionRowComm(self):
    comm = mpi.Comm()
    args = [self.obj,pointer(comm)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionRowComm_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionRowComm_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionRowComm_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionRowComm_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionRowComm_z(*args)
    else: DataExcept()
    return comm

  lib.ElDistMatrixColStride_i.argtypes = \
  lib.ElDistMatrixColStride_s.argtypes = \
  lib.ElDistMatrixColStride_d.argtypes = \
  lib.ElDistMatrixColStride_c.argtypes = \
  lib.ElDistMatrixColStride_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def ColStride(self):
    stride = c_int()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixColStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixColStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixColStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixColStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixColStride_z(*args)
    else: DataExcept()
    return stride.value

  lib.ElDistMatrixRowStride_i.argtypes = \
  lib.ElDistMatrixRowStride_s.argtypes = \
  lib.ElDistMatrixRowStride_d.argtypes = \
  lib.ElDistMatrixRowStride_c.argtypes = \
  lib.ElDistMatrixRowStride_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def RowStride(self):
    stride = c_int()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixRowStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRowStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRowStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRowStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRowStride_z(*args)
    else: DataExcept()
    return stride.value

  lib.ElDistMatrixPartialColStride_i.argtypes = \
  lib.ElDistMatrixPartialColStride_s.argtypes = \
  lib.ElDistMatrixPartialColStride_d.argtypes = \
  lib.ElDistMatrixPartialColStride_c.argtypes = \
  lib.ElDistMatrixPartialColStride_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialColStride(self):
    stride = c_int()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixPartialColStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialColStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialColStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialColStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialColStride_z(*args)
    else: DataExcept()
    return stride.value

  lib.ElDistMatrixPartialRowStride_i.argtypes = \
  lib.ElDistMatrixPartialRowStride_s.argtypes = \
  lib.ElDistMatrixPartialRowStride_d.argtypes = \
  lib.ElDistMatrixPartialRowStride_c.argtypes = \
  lib.ElDistMatrixPartialRowStride_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialRowStride(self):
    stride = c_int()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixPartialRowStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialRowStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialRowStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialRowStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialRowStride_z(*args)
    else: DataExcept()
    return stride.value

  lib.ElDistMatrixPartialUnionColStride_i.argtypes = \
  lib.ElDistMatrixPartialUnionColStride_s.argtypes = \
  lib.ElDistMatrixPartialUnionColStride_d.argtypes = \
  lib.ElDistMatrixPartialUnionColStride_c.argtypes = \
  lib.ElDistMatrixPartialUnionColStride_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialUnionColStride(self):
    stride = c_int()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionColStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionColStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionColStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionColStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionColStride_z(*args)
    else: DataExcept()
    return stride.value

  lib.ElDistMatrixPartialUnionRowStride_i.argtypes = \
  lib.ElDistMatrixPartialUnionRowStride_s.argtypes = \
  lib.ElDistMatrixPartialUnionRowStride_d.argtypes = \
  lib.ElDistMatrixPartialUnionRowStride_c.argtypes = \
  lib.ElDistMatrixPartialUnionRowStride_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def PartialUnionRowStride(self):
    stride = c_int()
    args = [self.obj,pointer(stride)]
    if   self.tag == iTag: lib.ElDistMatrixPartialUnionRowStride_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixPartialUnionRowStride_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixPartialUnionRowStride_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixPartialUnionRowStride_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixPartialUnionRowStride_z(*args)
    else: DataExcept()
    return stride.value

  lib.ElDistMatrixDistSize_i.argtypes = \
  lib.ElDistMatrixDistSize_s.argtypes = \
  lib.ElDistMatrixDistSize_d.argtypes = \
  lib.ElDistMatrixDistSize_c.argtypes = \
  lib.ElDistMatrixDistSize_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def DistSize(self):
    size = c_int()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElDistMatrixDistSize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDistSize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDistSize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDistSize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDistSize_z(*args)
    else: DataExcept()
    return size.value

  lib.ElDistMatrixCrossSize_i.argtypes = \
  lib.ElDistMatrixCrossSize_s.argtypes = \
  lib.ElDistMatrixCrossSize_d.argtypes = \
  lib.ElDistMatrixCrossSize_c.argtypes = \
  lib.ElDistMatrixCrossSize_z.argtypes = \
    [c_void_p,POINTER(iType)]
  def CrossSize(self):
    size = c_int()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElDistMatrixCrossSize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixCrossSize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixCrossSize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixCrossSize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixCrossSize_z(*args)
    else: DataExcept()
    return size.value
  
  lib.ElDistMatrixRedundantSize_i.argtypes = \
  lib.ElDistMatrixRedundantSize_s.argtypes = \
  lib.ElDistMatrixRedundantSize_d.argtypes = \
  lib.ElDistMatrixRedundantSize_c.argtypes = \
  lib.ElDistMatrixRedundantSize_z.argtypes = \
    [c_void_p,POINTER(c_int)]
  def RedundantSize(self):
    size = c_int()
    args = [self.obj,pointer(size)]
    if   self.tag == iTag: lib.ElDistMatrixRedundantSize_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixRedundantSize_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixRedundantSize_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixRedundantSize_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixRedundantSize_z(*args)
    else: DataExcept()
    return size.value

  lib.ElDistMatrixGet_i.argtypes = [c_void_p,iType,iType,POINTER(iType)]
  lib.ElDistMatrixGet_s.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElDistMatrixGet_d.argtypes = [c_void_p,iType,iType,POINTER(dType)]
  lib.ElDistMatrixGet_c.argtypes = [c_void_p,iType,iType,POINTER(cType)]
  lib.ElDistMatrixGet_z.argtypes = [c_void_p,iType,iType,POINTER(zType)]
  def Get(self,i,j):
    value = TagToType(self.tag)()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElDistMatrixGet_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGet_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGet_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGet_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGet_z(*args)
    else: DataExcept()
    return ScalarData(value)

  lib.ElDistMatrixGetRealPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElDistMatrixGetRealPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
  def GetRealPart(self,i,j):
    value = TagToType(Base(self.tag))()
    args = [self.obj,i,j,pointer(value)]
    if   self.tag == iTag: lib.ElDistMatrixGet_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixGet_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixGet_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixGetRealPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixGetRealPart_z(*args)
    else: DataExcept()
    return ScalarData(value)

  lib.ElDistMatrixGetImagPart_c.argtypes = [c_void_p,iType,iType,POINTER(sType)]
  lib.ElDistMatrixGetImagPart_z.argtypes = [c_void_p,iType,iType,POINTER(dType)]
  def GetImagPart(self,i,j):
    if   self.tag == iTag: return iType(0).value
    elif self.tag == sTag: return sType(0).value
    elif self.tag == dTag: return dType(0).value
    elif self.tag == cTag:
      value = sType()
      lib.ElDistMatrixGetRealPart_c(self.obj,i,j,pointer(value))
      return value.value
    elif self.tag == zTag:
      value = dType()
      lib.ElDistMatrixGetRealPart_z(self.obj,i,j,pointer(value))
      return value.value
    else: DataExcept()

  lib.ElDistMatrixSet_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistMatrixSet_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMatrixSet_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistMatrixSet_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistMatrixSet_z.argtypes = [c_void_p,iType,iType,zType]
  def Set(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElDistMatrixSet_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixSet_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixSet_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixSet_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSet_z(*args)
    else: DataExcept()

  lib.ElDistMatrixSetRealPart_c.argtypes = \
  lib.ElDistMatrixSetRealPart_z.argtypes = \
    [c_void_p,iType,iType,dType]
  def SetRealPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixSetRealPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetRealPart_z(*args)
    else: self.Set(i,j,value)

  lib.ElDistMatrixSetImagPart_c.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMatrixSetImagPart_z.argtypes = [c_void_p,iType,iType,dType]
  def SetImagPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixSetImagPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixSetImagPart_z(*args)
    else: raise Exception('Cannot set imaginary part of a real datatype')

  lib.ElDistMatrixUpdate_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistMatrixUpdate_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMatrixUpdate_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistMatrixUpdate_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistMatrixUpdate_z.argtypes = [c_void_p,iType,iType,zType]
  def Update(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElDistMatrixUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdate_z(*args)
    else: DataExcept()

  lib.ElDistMatrixUpdateRealPart_c.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMatrixUpdateRealPart_z.argtypes = [c_void_p,iType,iType,dType]
  def UpdateRealPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixUpdateRealPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateRealPart_z(*args)
    else: self.Update(i,j,value)

  lib.ElDistMatrixUpdateImagPart_c.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMatrixUpdateImagPart_z.argtypes = [c_void_p,iType,iType,dType]
  def UpdateImagPart(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == cTag: lib.ElDistMatrixUpdateImagPart_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixUpdateImagPart_z(*args)
    else: raise Exception('Cannot update imaginary part of a real datatype')

  lib.ElDistMatrixMakeReal_c.argtypes = \
  lib.ElDistMatrixMakeReal_z.argtypes = \
    [c_void_p,iType,iType]
  def MakeReal(self,i,j):
    args = [self.obj,i,j]
    if   self.tag == cTag: lib.ElDistMatrixMakeReal_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixMakeReal_z(*args)

  lib.ElDistMatrixConjugate_c.argtypes = \
  lib.ElDistMatrixConjugate_z.argtypes = \
    [c_void_p,iType,iType]
  def Conjugate(self,i,j):
    args = [self.obj,i,j]
    if   self.tag == cTag: lib.ElDistMatrixConjugate_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixConjugate_z(*args)

  lib.ElDistMatrixReserve_i.argtypes = \
  lib.ElDistMatrixReserve_s.argtypes = \
  lib.ElDistMatrixReserve_d.argtypes = \
  lib.ElDistMatrixReserve_c.argtypes = \
  lib.ElDistMatrixReserve_z.argtypes = \
    [c_void_p,iType]
  def Reserve(self,numRemoteEntries):
    args = [self.obj,numRemoteEntries]
    if   self.tag == iTag: lib.ElDistMatrixReserve_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixReserve_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixReserve_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixReserve_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixReserve_z(*args)
    else: DataExcept()

  lib.ElDistMatrixQueueUpdate_i.argtypes = [c_void_p,iType,iType,iType]
  lib.ElDistMatrixQueueUpdate_s.argtypes = [c_void_p,iType,iType,sType]
  lib.ElDistMatrixQueueUpdate_d.argtypes = [c_void_p,iType,iType,dType]
  lib.ElDistMatrixQueueUpdate_c.argtypes = [c_void_p,iType,iType,cType]
  lib.ElDistMatrixQueueUpdate_z.argtypes = [c_void_p,iType,iType,zType]
  def QueueUpdate(self,i,j,value):
    args = [self.obj,i,j,value]
    if   self.tag == iTag: lib.ElDistMatrixQueueUpdate_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixQueueUpdate_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixQueueUpdate_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixQueueUpdate_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixQueueUpdate_z(*args)
    else: DataExcept()

  lib.ElDistMatrixProcessQueues_i.argtypes = \
  lib.ElDistMatrixProcessQueues_s.argtypes = \
  lib.ElDistMatrixProcessQueues_d.argtypes = \
  lib.ElDistMatrixProcessQueues_c.argtypes = \
  lib.ElDistMatrixProcessQueues_z.argtypes = \
    [c_void_p]
  def ProcessQueues(self):
    args = [self.obj]
    if   self.tag == iTag: lib.ElDistMatrixProcessQueues_i(*args) 
    elif self.tag == sTag: lib.ElDistMatrixProcessQueues_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixProcessQueues_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixProcessQueues_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixProcessQueues_z(*args)
    else: DataExcept()

  lib.ElDistMatrixReservePulls_i.argtypes = \
  lib.ElDistMatrixReservePulls_s.argtypes = \
  lib.ElDistMatrixReservePulls_d.argtypes = \
  lib.ElDistMatrixReservePulls_c.argtypes = \
  lib.ElDistMatrixReservePulls_z.argtypes = \
    [c_void_p,iType]
  def ReservePulls(self,numPulls):
    args = [self.obj,numPulls]
    if   self.tag == iTag: lib.ElDistMatrixReservePulls_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixReservePulls_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixReservePulls_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixReservePulls_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixReservePulls_z(*args)
    else: DataExcept()

  lib.ElDistMatrixQueuePull_i.argtypes = \
  lib.ElDistMatrixQueuePull_s.argtypes = \
  lib.ElDistMatrixQueuePull_d.argtypes = \
  lib.ElDistMatrixQueuePull_c.argtypes = \
  lib.ElDistMatrixQueuePull_z.argtypes = \
    [c_void_p,iType,iType]
  def QueuePull(self,i,j):
    args = [self.obj,i,j]
    if   self.tag == iTag: lib.ElDistMatrixQueuePull_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixQueuePull_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixQueuePull_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixQueuePull_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixQueuePull_z(*args)
    else: DataExcept()

  lib.ElDistMatrixProcessPullQueue_i.argtypes = [c_void_p,POINTER(iType)]
  lib.ElDistMatrixProcessPullQueue_s.argtypes = [c_void_p,POINTER(sType)]
  lib.ElDistMatrixProcessPullQueue_d.argtypes = [c_void_p,POINTER(dType)]
  lib.ElDistMatrixProcessPullQueue_c.argtypes = [c_void_p,POINTER(cType)]
  lib.ElDistMatrixProcessPullQueue_z.argtypes = [c_void_p,POINTER(zType)]
  def ProcessPullQueue(self,buf):
    args = [self.obj,buf]
    if   self.tag == iTag: lib.ElDistMatrixProcessPullQueue_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixProcessPullQueue_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixProcessPullQueue_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixProcessPullQueue_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixProcessPullQueue_z(*args)
    else: DataExcept()

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

  lib.ElDistMatrixDiagonalAlignedWith_i.argtypes = \
  lib.ElDistMatrixDiagonalAlignedWith_s.argtypes = \
  lib.ElDistMatrixDiagonalAlignedWith_d.argtypes = \
  lib.ElDistMatrixDiagonalAlignedWith_c.argtypes = \
  lib.ElDistMatrixDiagonalAlignedWith_z.argtypes = \
    [c_void_p,DistData,iType,POINTER(bType)]
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

  lib.ElDistMatrixDiagonalRoot_i.argtypes = \
  lib.ElDistMatrixDiagonalRoot_s.argtypes = \
  lib.ElDistMatrixDiagonalRoot_d.argtypes = \
  lib.ElDistMatrixDiagonalRoot_c.argtypes = \
  lib.ElDistMatrixDiagonalRoot_z.argtypes = \
    [c_void_p,iType,POINTER(c_int)]
  def DiagonalRoot(self,offset=0):
    root = c_int()
    args = [self.obj,offset,pointer(root)]
    if   self.tag == iTag: lib.ElDistMatrixDiagonalRoot_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDiagonalRoot_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDiagonalRoot_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDiagonalRoot_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDiagonalRoot_z(*args)
    else: DataExcept()
    return root.value

  lib.ElDistMatrixDiagonalAlign_i.argtypes = \
  lib.ElDistMatrixDiagonalAlign_s.argtypes = \
  lib.ElDistMatrixDiagonalAlign_d.argtypes = \
  lib.ElDistMatrixDiagonalAlign_c.argtypes = \
  lib.ElDistMatrixDiagonalAlign_z.argtypes = \
    [c_void_p,iType,POINTER(c_int)]
  def DiagonalAlign(self,offset=0):
    align = c_int()
    args = [self.obj,offset,pointer(align)]
    if   self.tag == iTag: lib.ElDistMatrixDiagonalAlign_i(*args)
    elif self.tag == sTag: lib.ElDistMatrixDiagonalAlign_s(*args)
    elif self.tag == dTag: lib.ElDistMatrixDiagonalAlign_d(*args)
    elif self.tag == cTag: lib.ElDistMatrixDiagonalAlign_c(*args)
    elif self.tag == zTag: lib.ElDistMatrixDiagonalAlign_z(*args)
    else: DataExcept()
    return align.value

  lib.ElViewDist_i.argtypes = \
  lib.ElViewDist_s.argtypes = \
  lib.ElViewDist_d.argtypes = \
  lib.ElViewDist_c.argtypes = \
  lib.ElViewDist_z.argtypes = \
  lib.ElLockedViewDist_i.argtypes = \
  lib.ElLockedViewDist_s.argtypes = \
  lib.ElLockedViewDist_d.argtypes = \
  lib.ElLockedViewDist_c.argtypes = \
  lib.ElLockedViewDist_z.argtypes = \
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
