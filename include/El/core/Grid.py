#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *
import ctypes

# Grid
# ====

class Grid(object):
  def __init__(self):
    self.obj = ctypes.c_void_p()
    lib.ElDefaultGrid(pointer(self.obj))
  @classmethod
  def FromComm(self,comm,order):
    if type(comm) is not MPI_Comm: raise Exception('Invalid comm type')
    if type(order) is not ctypes.c_uint: raise Exception('Invalid order type')
    self.obj = ctypes.c_void_p()
    lib.ElGridCreate(comm,order,pointer(self.obj))
  @classmethod
  def FromCommSpecific(self,comm,height,order):
    if type(comm) is not MPI_Comm: raise Exception('Invalid comm type')
    if type(height) is not ctypes.c_int: raise Exception('Invalid height type')
    if type(order) is not ctypes.c_uint: raise Exception('Invalid order type')
    self.obj = ctypes.c_void_p()
    lib.ElGridCreateSpecific(comm,height,order,pointer(self.obj))
  def Destroy(self): lib.ElGridDestroy(self.obj)
  def Row(self):
    row = ctypes.c_int()
    lib.ElGridRow(self.obj,pointer(row))
    return row
  def Col(self):
    col = ctypes.c_int()
    lib.ElGridCol(self.obj,pointer(col))
    return col
  def Rank(self):
    rank = ctypes.c_int()
    lib.ElGridRank(self.obj,pointer(rank))
    return rank
  def Height(self):
    height = ctypes.c_int()
    lib.ElGridHeight(self.obj,pointer(height))
    return height
  def Width(self):
    width = ctypes.c_int()
    lib.ElGridWidth(self.obj,pointer(width))
    return width
  def Size(self):
    size = ctypes.c_int()
    lib.ElGridSize(self.obj,pointer(size))
    return size
  def Order(self):
    order = ctypes.c_uint()
    lib.ElGridOrder(self.obj,pointer(order))
    return order
  def ColComm(self):
    colComm = MPI_Comm()
    lib.ElGridColComm(self.obj,pointer(colComm))
    return colComm
  def RowComm(self):
    rowComm = MPI_Comm()
    lib.ElGridColComm(self.obj,pointer(rowComm))
    return rowComm
  def Comm(self):
    comm = MPI_Comm()
    lib.ElGridComm(self.obj,pointer(comm))
    return comm
  def MCRank(self):
    rank = ctypes.c_int()
    lib.ElGridMCRank(self.obj,pointer(rank))
    return rank
  def MRRank(self):
    rank = ctypes.c_int()
    lib.ElGridMRRank(self.obj,pointer(rank))
    return rank
  def VCRank(self):
    rank = ctypes.c_int()
    lib.ElGridVCRank(self.obj,pointer(rank))
    return rank
  def VRRank(self):
    rank = ctypes.c_int()
    lib.ElGridVRRank(self.obj,pointer(rank))
    return rank
  def MCSize(self):
    size = ctypes.c_int()
    lib.ElGridMCSize(self.obj,pointer(size))
    return size
  def MRSize(self):
    size = ctypes.c_int()
    lib.ElGridMRSize(self.obj,pointer(size))
    return size
  def VCSize(self):
    size = ctypes.c_int()
    lib.ElGridVCSize(self.obj,pointer(size))
    return size
  def VRSize(self):
    size = ctypes.c_int()
    lib.ElGridVRSize(self.obj,pointer(size))
    return size
  def MCComm(self):
    comm = MPI_Comm()
    lib.ElGridMCComm(self.obj,pointer(comm))
    return comm
  def MRComm(self):
    comm = MPI_Comm()
    lib.ElGridMRComm(self.obj,pointer(comm))
    return comm
  def VCComm(self):
    comm = MPI_Comm()
    lib.ElGridVCComm(self.obj,pointer(comm))
    return comm
  def VRComm(self):
    comm = MPI_Comm()
    lib.ElGridVRComm(self.obj,pointer(comm))
    return comm
  def MDComm(self):
    comm = MPI_Comm()
    lib.ElGridMDComm(self.obj,pointer(comm))
    return comm
  def MDPerpComm(self):
    comm = MPI_Comm()
    lib.ElGridMDPerpComm(self.obj,pointer(comm))
  def GCD(self):
    gcd = ctypes.c_int()
    lib.ElGridGCD(self.obj,pointer(gcd))
    return gcd
  def LCM(self):
    lcm = ctypes.c_int()
    lib.ElGridLCM(self.obj,pointer(lcm))
    return lcm
  def InGrid(self):
    inGrid = bType()
    lib.ElGridInGrid(self.obj,pointer(inGrid))
    return inGrid
  def HaveViewers(self):
    haveViewers = bType() 
    lib.ElGridHaveViewers(self.obj,pointer(haveViewers))
    return haveViewers
  def OwningRank(self):
    rank = ctypes.c_int()
    lib.ElGridOwningRank(self.obj,pointer(rank))
    return rank
  def ViewingRank(self):
    rank = ctypes.c_int()
    lib.ElGridViewingRank(self.obj,pointer(rank))
    return rank
  def VCToViewingMap(self,vcRank):
    if type(vcRank) is not ctypes.c_int:
      raise Exception("Expected vcRank to be an integer")
    viewingRank = ctypes.c_int()
    lib.ElGridVCToViewingMap(self.obj,vcRank,pointer(viewingRank))
    return viewingRank
  def OwningGroup(self):
    group = MPI_Group()
    lib.ElGridOwningGroup(self.obj,pointer(group))
    return group
  def OwningComm(self):
    comm = MPI_Comm()
    lib.ElGridOwningComm(self.obj,pointer(comm))
    return comm
  def ViewingComm(self):
    comm = MPI_Comm()
    lib.ElGridViewingComm(self.obj,pointer(comm))
    return comm
  def DiagPath(self,vcRank):
    if type(vcRank) is not ctypes.c_int:
      raise Exception("Expected vcRank to be an integer")
    diagPath = ctypes.c_int()
    lib.ElGridDiagPath(self.obj,vcRank,pointer(diagPath))
    return diagPath
  def DiagPathRank(self,vcRank):
    if type(vcRank) is not ctypes.c_int:
      raise Exception("Expected vcRank to be an integer")
    pathRank = ctypes.c_int()
    lib.ElGridDiagPathRank(self.obj,vcRank,pointer(pathRank))
    return pathRank
  def FirstVCRank(self,diagPath):
    if type(diagPath) is not ctypes.c_int:
      raise Exception("Expected diagPath to be an integer")
    firstVCRank = ctypes.c_int()
    lib.ElGridFirstVCRank(self.obj,diagPath,pointer(firstVCRank))
    return firstVCRank
  # NOTE: The following method is static
  def FindFactor(numProcs):
    if type(numProcs) is not ctypes.c_int:
      raise Exception("Expceted numProcs to be an integer")
    factor = ctypes.c_int()
    lib.ElGridFindFactor(numProcs,pointer(factor))
    return factor

def DefaultGrid():
  return Grid()
