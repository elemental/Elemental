#
# Copyright (c) 2009-2016, Jack Poulson
# All rights reserved.
#
# This file is part of Elemental and is under the BSD 2-Clause License,
# which can be found in the LICENSE file in the root directory, or at
# http://opensource.org/licenses/BSD-2-Clause
#
from ..environment import *

# Query Elemental to determine whether MPI_Comm is an 'int' or a void pointer
# TODO(poulson): Complain if not same size as c_int or c_void_p
commSameSizeAsInteger = bType()
lib.ElMPICommSameSizeAsInteger(pointer(commSameSizeAsInteger))
if commSameSizeAsInteger:
  Comm = c_int
else:
  Comm = c_void_p

# Query Elemental to determine whether MPI_Group is an 'int' or a void pointer
# TODO(poulson): Complain if not same size as c_int or c_void_p
groupSameSizeAsInteger = bType()
lib.ElMPIGroupSameSizeAsInteger(pointer(groupSameSizeAsInteger))
if groupSameSizeAsInteger:
  Group = c_int
else:
  Group = c_void_p

# Return MPI_COMM_WORLD
lib.ElMPICommWorld.argtypes = [POINTER(Comm)]
def COMM_WORLD():
  comm = Comm()
  lib.ElMPICommWorld(pointer(comm))
  return comm

# Return MPI_COMM_SELF
lib.ElMPICommSelf.argtypes = [POINTER(Comm)]
def COMM_SELF():
  comm = Comm()
  lib.ElMPICommSelf(pointer(comm))
  return comm

lib.ElMPICommRank.argtypes = [Comm,POINTER(c_int)]
lib.ElMPIGroupRank.argtypes = [Group,POINTER(c_int)]
def Rank(mpiObj):
  rank = c_int()
  args = [mpiObj,pointer(rank)]
  if type(mpiObj) is Comm:
    lib.ElMPICommRank(*args)
  elif type(mpiObj) is Group:
    lib.ElMPIGroupRank(*args)
  else: TypeExcept()
  return rank.value

lib.ElMPICommSize.argtypes = [Comm,POINTER(c_int)]
lib.ElMPIGroupSize.argtypes = [Group,POINTER(c_int)]
def Size(mpiObj):
  size = c_int()
  args = [mpiObj,pointer(size)]
  if type(mpiObj) is Comm:
    lib.ElMPICommSize(*args)
  elif type(mpiObj) is Group:
    lib.ElMPIGroupSize(*args)
  else: TypeExcept()
  return size.value

lib.ElMPICommFree.argtypes = [POINTER(Comm)]
lib.ElMPIGroupFree.argtypes = [POINTER(Group)]
def Free(mpiObj):
  if type(mpiObj) is POINTER(Comm):
    lib.ElMPICommFree(pointer(mpiObj))
  elif type(mpiObj) is POINTER(Group):
    lib.ElMPIGroupFree(pointer(mpiObj))
  else: DataExcept()

lib.ElMPIWorldRank.argtypes = [POINTER(c_int)]
lib.ElMPIWorldRank.restype = c_uint
def WorldRank():
  rank = c_int()
  lib.ElMPIWorldRank(pointer(rank))
  return rank.value

lib.ElMPIWorldSize.argtypes = [POINTER(c_int)]
lib.ElMPIWorldSize.restype = c_uint
def WorldSize():
  size = c_int()
  lib.ElMPIWorldSize(pointer(size))
  return size.value

lib.ElMPITime.argtypes = [POINTER(c_double)]
lib.ElMPITime.restype = c_uint
def Time():
  walltime = c_double()
  lib.ElMPITime(pointer(walltime))
  return walltime.value
