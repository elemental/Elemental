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

class Grid(object):

  lib.ElDefaultGrid.argtypes = [POINTER(c_void_p)]
  def __init__(self,create=True):
    self.obj = c_void_p()
    if create:
      lib.ElDefaultGrid(pointer(self.obj))

  @classmethod
  def Default(cls):
    grid = cls(False)
    grid.obj = c_void_p()
    lib.ElDefaultGrid(pointer(grid.obj))
    return grid

  lib.ElTrivialGrid.argtypes = [POINTER(c_void_p)]
  @classmethod
  def Trivial(cls):
    grid = cls(False)
    grid.obj = c_void_p()
    lib.ElTrivialGrid(pointer(grid.obj))
    return grid

  lib.ElGridCreate.argtypes = [mpi.Comm,c_uint,POINTER(c_void_p)]
  @classmethod
  def FromComm(cls,comm=mpi.COMM_WORLD(),order=COL_MAJOR):
    grid = cls(False)
    grid.obj = c_void_p()
    lib.ElGridCreate(comm,order,pointer(grid.obj))
    return grid

  lib.ElGridCreateSpecific.argtypes = [mpi.Comm,c_int,c_uint,POINTER(c_void_p)]
  @classmethod
  def FromCommSpecific(cls,comm,height,order=COL_MAJOR):
    grid = cls(False)
    grid.obj = c_void_p()
    lib.ElGridCreateSpecific(comm,height,order,pointer(grid.obj))
    return grid

  lib.ElGridDestroy.argtypes = [c_void_p]
  def Destroy(self):
    lib.ElGridDestroy(self.obj)

  lib.ElGridRow.argtypes = [c_void_p,POINTER(c_int)]
  def Row(self):
    row = c_int()
    lib.ElGridRow(self.obj,pointer(row))
    return row.value

  lib.ElGridCol.argtypes = [c_void_p,POINTER(c_int)]
  def Col(self):
    col = c_int()
    lib.ElGridCol(self.obj,pointer(col))
    return col.value

  lib.ElGridRank.argtypes = [c_void_p,POINTER(c_int)]
  def Rank(self):
    rank = c_int()
    lib.ElGridRank(self.obj,pointer(rank))
    return rank.value

  lib.ElGridHeight.argtypes = [c_void_p,POINTER(c_int)]
  def Height(self):
    height = c_int()
    lib.ElGridHeight(self.obj,pointer(height))
    return height.value

  lib.ElGridWidth.argtypes = [c_void_p,POINTER(c_int)]
  def Width(self):
    width = c_int()
    lib.ElGridWidth(self.obj,pointer(width))
    return width.value

  lib.ElGridSize.argtypes = [c_void_p,POINTER(c_int)]
  def Size(self):
    size = c_int()
    lib.ElGridSize(self.obj,pointer(size))
    return size.value

  lib.ElGridOrder.argtypes = [c_void_p,POINTER(c_uint)]
  def Order(self):
    order = c_uint()
    lib.ElGridOrder(self.obj,pointer(order))
    return order.value

  lib.ElGridColComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def ColComm(self):
    colComm = mpi.Comm()
    lib.ElGridColComm(self.obj,pointer(colComm))
    return colComm

  lib.ElGridRowComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def RowComm(self):
    rowComm = mpi.Comm()
    lib.ElGridRowComm(self.obj,pointer(rowComm))
    return rowComm

  lib.ElGridComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def Comm(self):
    comm = mpi.Comm()
    lib.ElGridComm(self.obj,pointer(comm))
    return comm

  lib.ElGridMCRank.argtypes = [c_void_p,POINTER(c_int)]
  def MCRank(self):
    rank = c_int()
    lib.ElGridMCRank(self.obj,pointer(rank))
    return rank.value

  lib.ElGridMRRank.argtypes = [c_void_p,POINTER(c_int)]
  def MRRank(self):
    rank = c_int()
    lib.ElGridMRRank(self.obj,pointer(rank))
    return rank.value

  lib.ElGridVCRank.argtypes = [c_void_p,POINTER(c_int)]
  def VCRank(self):
    rank = c_int()
    lib.ElGridVCRank(self.obj,pointer(rank))
    return rank.value

  lib.ElGridVRRank.argtypes = [c_void_p,POINTER(c_int)]
  def VRRank(self):
    rank = c_int()
    lib.ElGridVRRank(self.obj,pointer(rank))
    return rank.value

  lib.ElGridMCSize.argtypes = [c_void_p,POINTER(c_int)]
  def MCSize(self):
    size = c_int()
    lib.ElGridMCSize(self.obj,pointer(size))
    return size.value

  lib.ElGridMRSize.argtypes = [c_void_p,POINTER(c_int)]
  def MRSize(self):
    size = c_int()
    lib.ElGridMRSize(self.obj,pointer(size))
    return size.value

  lib.ElGridVCSize.argtypes = [c_void_p,POINTER(c_int)]
  def VCSize(self):
    size = c_int()
    lib.ElGridVCSize(self.obj,pointer(size))
    return size.value

  lib.ElGridVRSize.argtypes = [c_void_p,POINTER(c_int)]
  def VRSize(self):
    size = c_int()
    lib.ElGridVRSize(self.obj,pointer(size))
    return size.value

  lib.ElGridMCComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def MCComm(self):
    comm = mpi.Comm()
    lib.ElGridMCComm(self.obj,pointer(comm))
    return comm

  lib.ElGridMRComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def MRComm(self):
    comm = mpi.Comm()
    lib.ElGridMRComm(self.obj,pointer(comm))
    return comm

  lib.ElGridVCComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def VCComm(self):
    comm = mpi.Comm()
    lib.ElGridVCComm(self.obj,pointer(comm))
    return comm

  lib.ElGridVRComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def VRComm(self):
    comm = mpi.Comm()
    lib.ElGridVRComm(self.obj,pointer(comm))
    return comm

  lib.ElGridMDComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def MDComm(self):
    comm = mpi.Comm()
    lib.ElGridMDComm(self.obj,pointer(comm))
    return comm

  lib.ElGridMDPerpComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def MDPerpComm(self):
    comm = mpi.Comm()
    lib.ElGridMDPerpComm(self.obj,pointer(comm))
    return comm

  lib.ElGridGCD.argtypes = [c_void_p,POINTER(c_int)]
  def GCD(self):
    gcd = c_int()
    lib.ElGridGCD(self.obj,pointer(gcd))
    return gcd.value

  lib.ElGridLCM.argtypes = [c_void_p,POINTER(c_int)]
  def LCM(self):
    lcm = c_int()
    lib.ElGridLCM(self.obj,pointer(lcm))
    return lcm.value

  lib.ElGridInGrid.argtypes = [c_void_p,POINTER(bType)]
  def InGrid(self):
    inGrid = bType()
    lib.ElGridInGrid(self.obj,pointer(inGrid))
    return inGrid.value

  lib.ElGridHaveViewers.argtypes = [c_void_p,POINTER(bType)]
  def HaveViewers(self):
    haveViewers = bType()
    lib.ElGridHaveViewers(self.obj,pointer(haveViewers))
    return haveViewers.value

  lib.ElGridOwningRank.argtypes = [c_void_p,POINTER(c_int)]
  def OwningRank(self):
    rank = c_int()
    lib.ElGridOwningRank(self.obj,pointer(rank))
    return rank.value

  lib.ElGridViewingRank.argtypes = [c_void_p,POINTER(c_int)]
  def ViewingRank(self):
    rank = c_int()
    lib.ElGridViewingRank(self.obj,pointer(rank))
    return rank.value

  lib.ElGridOwningGroup.argtypes = [c_void_p,POINTER(mpi.Group)]
  def OwningGroup(self):
    group = mpi.Group()
    lib.ElGridOwningGroup(self.obj,pointer(group))
    return group

  lib.ElGridOwningComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def OwningComm(self):
    comm = mpi.Comm()
    lib.ElGridOwningComm(self.obj,pointer(comm))
    return comm

  lib.ElGridViewingComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def ViewingComm(self):
    comm = mpi.Comm()
    lib.ElGridViewingComm(self.obj,pointer(comm))
    return comm

  lib.ElGridDiag.argtypes = [c_void_p,c_int,POINTER(c_int)]
  def Diag(self,vcRank):
    diag = c_int()
    lib.ElGridDiag(self.obj,vcRank,pointer(diag))
    return diag.value

  lib.ElGridDiagRank.argtypes = [c_void_p,c_int,POINTER(c_int)]
  def DiagRank(self,vcRank):
    pathRank = c_int()
    lib.ElGridDiagRank(self.obj,vcRank,pointer(pathRank))
    return pathRank.value

  # TODO(poulson): VCToVR
  # TODO(poulson): VRToVC
  # TODO(poulson): CoordsToVC

  lib.ElGridVCToViewing.argtypes = [c_void_p,c_int,POINTER(c_int)]
  def VCToViewing(self,vcRank):
    viewingRank = c_int()
    lib.ElGridVCToViewing(self.obj,vcRank,pointer(viewingRank))
    return viewingRank.value

  # NOTE: The following method is static
  lib.ElGridDefaultHeight.argtypes = [c_int,POINTER(c_int)]
  def DefaultHeight(numProcs):
    factor = c_int()
    lib.ElGridDefaultHeight(numProcs,pointer(factor))
    return factor.value
