#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
from imports     import mpi

# Grid
# ====

lib.ElDefaultGrid.argtypes = [POINTER(c_void_p)]
lib.ElDefaultGrid.restype = c_uint

lib.ElGridCreate.argtypes = [mpi.Comm,c_uint,POINTER(c_void_p)]
lib.ElGridCreate.restype = c_uint

lib.ElGridCreateSpecific.argtypes = [mpi.Comm,c_int,c_uint,POINTER(c_void_p)]
lib.ElGridCreateSpecific.restype = c_uint

lib.ElGridDestroy.argtypes = [c_void_p]
lib.ElGridDestroy.restype = c_uint

lib.ElGridRow.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridRow.restype = c_uint

lib.ElGridCol.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridCol.restype = c_uint

lib.ElGridRank.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridRank.restype = c_uint

lib.ElGridHeight.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridHeight.restype = c_uint

lib.ElGridWidth.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridWidth.restype = c_uint

lib.ElGridSize.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridSize.restype = c_uint

lib.ElGridOrder.argtypes = [c_void_p,POINTER(c_uint)]
lib.ElGridOrder.restype = c_uint

lib.ElGridColComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridColComm.restype = c_uint

lib.ElGridRowComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridRowComm.restype = c_uint

lib.ElGridComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridComm.restype = c_uint

lib.ElGridMCRank.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridMCRank.restype = c_uint

lib.ElGridMRRank.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridMRRank.restype = c_uint

lib.ElGridVCRank.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridVCRank.restype = c_uint

lib.ElGridVRRank.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridVRRank.restype = c_uint

lib.ElGridMCSize.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridMCSize.restype = c_uint

lib.ElGridMRSize.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridMRSize.restype = c_uint

lib.ElGridVCSize.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridVCSize.restype = c_uint

lib.ElGridVRSize.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridVRSize.restype = c_uint

lib.ElGridMCComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridMCComm.restype = c_uint

lib.ElGridMRComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridMRComm.restype = c_uint

lib.ElGridVCComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridVCComm.restype = c_uint

lib.ElGridVRComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridVRComm.restype = c_uint

lib.ElGridMDComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridMDComm.restype = c_uint

lib.ElGridMDPerpComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridMDPerpComm.restype = c_uint

lib.ElGridGCD.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridGCD.restype = c_uint

lib.ElGridLCM.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridLCM.restype = c_uint

lib.ElGridInGrid.argtypes = [c_void_p,POINTER(bType)]
lib.ElGridInGrid.restype = c_uint

lib.ElGridHaveViewers.argtypes = [c_void_p,POINTER(bType)]
lib.ElGridHaveViewers.restype = c_uint

lib.ElGridOwningRank.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridOwningRank.restype = c_uint

lib.ElGridViewingRank.argtypes = [c_void_p,POINTER(c_int)]
lib.ElGridViewingRank.restype = c_uint

lib.ElGridVCToViewingMap.argtypes = [c_void_p,c_int,POINTER(c_int)]
lib.ElGridVCToViewingMap.restype = c_uint

lib.ElGridOwningGroup.argtypes = [c_void_p,POINTER(mpi.Group)]
lib.ElGridOwningGroup.restype = c_uint

lib.ElGridOwningComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridOwningComm.restype = c_uint

lib.ElGridViewingComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElGridViewingComm.restype = c_uint

lib.ElGridDiagPath.argtypes = [c_void_p,c_int,POINTER(c_int)]
lib.ElGridDiagPath.restype = c_uint

lib.ElGridDiagPathRank.argtypes = [c_void_p,c_int,POINTER(c_int)]
lib.ElGridDiagPathRank.restype = c_uint

lib.ElGridFirstVCRank.argtypes = [c_void_p,c_int,POINTER(c_int)]
lib.ElGridFirstVCRank.restype = c_uint

lib.ElGridFindFactor.argtypes = [c_int,POINTER(c_int)]
lib.ElGridFindFactor.restype = c_uint

class Grid(object):
  def __init__(self,create=True):
    self.obj = c_void_p()
    if create:
      lib.ElDefaultGrid(pointer(self.obj))
  @classmethod
  def FromComm(cls,comm=mpi.COMM_WORLD(),order=COL_MAJOR):
    g = cls(False)
    lib.ElGridCreate(comm,order,pointer(g.obj))
    return g
  @classmethod
  def FromCommSpecific(cls,comm,height,order=COL_MAJOR):
    g = cls(False)
    lib.ElGridCreateSpecific(comm,height,order,pointer(g.obj))
    return g
  def Destroy(self): 
    lib.ElGridDestroy(self.obj)
  def Row(self):
    row = c_int()
    lib.ElGridRow(self.obj,pointer(row))
    return row.value
  def Col(self):
    col = c_int()
    lib.ElGridCol(self.obj,pointer(col))
    return col.value
  def Rank(self):
    rank = c_int()
    lib.ElGridRank(self.obj,pointer(rank))
    return rank.value
  def Height(self):
    height = c_int()
    lib.ElGridHeight(self.obj,pointer(height))
    return height.value
  def Width(self):
    width = c_int()
    lib.ElGridWidth(self.obj,pointer(width))
    return width.value
  def Size(self):
    size = c_int()
    lib.ElGridSize(self.obj,pointer(size))
    return size.value
  def Order(self):
    order = c_uint()
    lib.ElGridOrder(self.obj,pointer(order))
    return order.value
  def ColComm(self):
    colComm = mpi.Comm()
    lib.ElGridColComm(self.obj,pointer(colComm))
    return colComm
  def RowComm(self):
    rowComm = mpi.Comm()
    lib.ElGridColComm(self.obj,pointer(rowComm))
    return rowComm
  def Comm(self):
    comm = mpi.Comm()
    lib.ElGridComm(self.obj,pointer(comm))
    return comm
  def MCRank(self):
    rank = c_int()
    lib.ElGridMCRank(self.obj,pointer(rank))
    return rank.value
  def MRRank(self):
    rank = c_int()
    lib.ElGridMRRank(self.obj,pointer(rank))
    return rank.value
  def VCRank(self):
    rank = c_int()
    lib.ElGridVCRank(self.obj,pointer(rank))
    return rank.value
  def VRRank(self):
    rank = c_int()
    lib.ElGridVRRank(self.obj,pointer(rank))
    return rank.value
  def MCSize(self):
    size = c_int()
    lib.ElGridMCSize(self.obj,pointer(size))
    return size.value
  def MRSize(self):
    size = c_int()
    lib.ElGridMRSize(self.obj,pointer(size))
    return size.value
  def VCSize(self):
    size = c_int()
    lib.ElGridVCSize(self.obj,pointer(size))
    return size.value
  def VRSize(self):
    size = c_int()
    lib.ElGridVRSize(self.obj,pointer(size))
    return size.value
  def MCComm(self):
    comm = mpi.Comm()
    lib.ElGridMCComm(self.obj,pointer(comm))
    return comm
  def MRComm(self):
    comm = mpi.Comm()
    lib.ElGridMRComm(self.obj,pointer(comm))
    return comm
  def VCComm(self):
    comm = mpi.Comm()
    lib.ElGridVCComm(self.obj,pointer(comm))
    return comm
  def VRComm(self):
    comm = mpi.Comm()
    lib.ElGridVRComm(self.obj,pointer(comm))
    return comm
  def MDComm(self):
    comm = mpi.Comm()
    lib.ElGridMDComm(self.obj,pointer(comm))
    return comm
  def MDPerpComm(self):
    comm = mpi.Comm()
    lib.ElGridMDPerpComm(self.obj,pointer(comm))
  def GCD(self):
    gcd = c_int()
    lib.ElGridGCD(self.obj,pointer(gcd))
    return gcd.value
  def LCM(self):
    lcm = c_int()
    lib.ElGridLCM(self.obj,pointer(lcm))
    return lcm.value
  def InGrid(self):
    inGrid = bType()
    lib.ElGridInGrid(self.obj,pointer(inGrid))
    return inGrid.value
  def HaveViewers(self):
    haveViewers = bType() 
    lib.ElGridHaveViewers(self.obj,pointer(haveViewers))
    return haveViewers.value
  def OwningRank(self):
    rank = c_int()
    lib.ElGridOwningRank(self.obj,pointer(rank))
    return rank.value
  def ViewingRank(self):
    rank = c_int()
    lib.ElGridViewingRank(self.obj,pointer(rank))
    return rank.value
  def VCToViewingMap(self,vcRank):
    viewingRank = c_int()
    lib.ElGridVCToViewingMap(self.obj,vcRank,pointer(viewingRank))
    return viewingRank.value
  def OwningGroup(self):
    group = mpi.Group()
    lib.ElGridOwningGroup(self.obj,pointer(group))
    return group
  def OwningComm(self):
    comm = mpi.Comm()
    lib.ElGridOwningComm(self.obj,pointer(comm))
    return comm
  def ViewingComm(self):
    comm = mpi.Comm()
    lib.ElGridViewingComm(self.obj,pointer(comm))
    return comm
  def DiagPath(self,vcRank):
    diagPath = c_int()
    lib.ElGridDiagPath(self.obj,vcRank,pointer(diagPath))
    return diagPath.value
  def DiagPathRank(self,vcRank):
    pathRank = c_int()
    lib.ElGridDiagPathRank(self.obj,vcRank,pointer(pathRank))
    return pathRank.value
  def FirstVCRank(self,diagPath):
    firstVCRank = c_int()
    lib.ElGridFirstVCRank(self.obj,diagPath,pointer(firstVCRank))
    return firstVCRank.value
  # NOTE: The following method is static
  def FindFactor(numProcs):
    factor = c_int()
    lib.ElGridFindFactor(numProcs,pointer(factor))
    return factor.value

def DefaultGrid():
  return Grid()
