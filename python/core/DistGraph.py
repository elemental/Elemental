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

import Grid

class DistGraph(object):
  # Constructors and destructors
  # ============================
  lib.ElDistGraphCreate.argtypes = [POINTER(c_void_p),c_void_p]
  def __init__(self,grid=Grid.Grid.Default(),create=True):
    self.obj = c_void_p()
    if create:
      lib.ElDistGraphCreate(pointer(self.obj),grid.obj)

  lib.ElDistGraphDestroy.argtypes = [c_void_p]
  def Destroy(self):
    lib.ElDistGraphDestroy(self.obj)

  # Assignment and reconfiguration
  # ==============================
  lib.ElDistGraphEmpty.argtypes = [c_void_p]
  def Empty(self):
    lib.ElDistGraphEmpty(self.obj)

  lib.ElDistGraphResize.argtypes = [c_void_p,iType,iType]
  def Resize(self,numSources,numTargets=None):
    if numTargets == None:
      numTargets = numSources
    lib.ElDistGraphResize(self.obj,numSources,numTargets)

  lib.ElDistGraphSetGrid.argtypes = [c_void_p,c_void_p]
  def SetGrid(self,grid):
    lib.ElDistGraphSetGrid(self,grid.obj)

  lib.ElDistGraphReserve.argtypes = [c_void_p,iType,iType]
  def Reserve(self,numLocalEdges,numRemoteEdges=0):
    lib.ElDistGraphReserve(self.obj,numLocalEdges,numRemoteEdges)

  lib.ElDistGraphConnect.argtypes = [c_void_p,iType,iType]
  def Connect(self,source,target):
    lib.ElDistGraphConnect(self.obj,source,target)

  lib.ElDistGraphConnectLocal.argtypes = [c_void_p,iType,iType]
  def ConnectLocal(self,localSource,target):
    lib.ElDistGraphConnectLocal(self.obj,localSource,target)

  lib.ElDistGraphDisconnect.argtypes = [c_void_p,iType,iType]
  def Disconnect(self,source,target):
    lib.ElDistGraphDisconnect(self.obj,source,target)

  lib.ElDistGraphDisconnectLocal.argtypes = [c_void_p,iType,iType]
  def DisconnectLocal(self,localSource,target):
    lib.ElDistGraphDisconnectLocal(self.obj,localSource,target)

  lib.ElDistGraphQueueConnection.argtypes = [c_void_p,iType,iType,bType]
  def QueueConnection(self,source,target,passive=False):
    lib.ElDistGraphQueueConnection(self.obj,source,target,passive)

  lib.ElDistGraphQueueLocalConnection.argtypes = [c_void_p,iType,iType]
  def QueueLocalConnection(self,localSource,target):
    lib.ElDistGraphQueueLocalConnection(self.obj,localSource,target)

  lib.ElDistGraphQueueDisconnection.argtypes = [c_void_p,iType,iType,bType]
  def QueueDisconnection(self,source,target,passive=False):
    lib.ElDistGraphQueueDisconnection(self.obj,source,target,passive)

  lib.ElDistGraphQueueLocalDisconnection.argtypes = [c_void_p,iType,iType]
  def QueueLocalDisconnection(self,localSource,target):
    lib.ElDistGraphQueueLocalDisconnection(self.obj,localSource,target)

  lib.ElDistGraphProcessQueues.argtypes = [c_void_p]
  def ProcessQueues(self):
    lib.ElDistGraphProcessQueues(self.obj)

  lib.ElDistGraphProcessLocalQueues.argtypes = [c_void_p]
  def ProcessLocalQueues(self):
    lib.ElDistGraphProcessLocalQueues(self.obj)

  # Queries
  # =======
  lib.ElDistGraphNumSources.argtypes = [c_void_p,POINTER(iType)]
  def NumSources(self):
    numSources = iType()
    lib.ElDistGraphNumSources(self.obj,pointer(numSources))
    return numSources.value

  lib.ElDistGraphNumTargets.argtypes = [c_void_p,POINTER(iType)]
  def NumTargets(self):
    numTargets = iType()
    lib.ElDistGraphNumTargets(self.obj,pointer(numTargets))
    return numTargets.value

  lib.ElDistGraphFirstLocalSource.argtypes = [c_void_p,POINTER(iType)]
  def FirstLocalSource(self):
    firstLocalSource = iType()
    lib.ElDistGraphFirstLocalSource(self.obj,pointer(firstLocalSource))
    return firstLocalSource.value

  lib.ElDistGraphNumLocalSources.argtypes = [c_void_p,POINTER(iType)]
  def NumLocalSources(self):
    numLocalSources = iType()
    lib.ElDistGraphNumLocalSources(self.obj,pointer(numLocalSources))
    return numLocalSources.value

  lib.ElDistGraphNumLocalEdges.argtypes = [c_void_p,POINTER(iType)]
  def NumLocalEdges(self):
    numLocalEdges = iType()
    lib.ElDistGraphNumLocalEdges(self.obj,pointer(numLocalEdges))
    return numLocalEdges.value

  lib.ElDistGraphCapacity.argtypes = [c_void_p,POINTER(iType)]
  def Capacity(self):
    capacity = iType()
    lib.ElDistGraphCapacity(self.obj,pointer(capacity))
    return capacity.value

  lib.ElDistGraphLocallyConsistent.argtypes = [c_void_p,POINTER(bType)]
  def LocallyConsistent(self):
    consistent = bType()
    lib.ElDistGraphLocallyConsistent(self.obj,pointer(consistent))
    return consistent.value

  lib.ElDistGraphGrid.argtypes = [c_void_p,POINTER(c_void_p)]
  def Grid(self):
    grid = Grid.Grid(create=False)
    lib.ElDistGraphGrid(self.obj,pointer(grid.obj))
    return grid

  lib.ElDistGraphBlocksize.argtypes = [c_void_p,POINTER(iType)]
  def Blocksize(self):
    blocksize = iType()
    lib.ElDistGraphBlocksize(self.obj,pointer(blocksize))
    return blocksize.value

  lib.ElDistGraphSource.argtypes = [c_void_p,iType,POINTER(iType)]
  def Source(self,localEdge):
    source = iType()
    lib.ElDistGraphSource(self.obj,localEdge,pointer(source))
    return source.value

  lib.ElDistGraphTarget.argtypes = [c_void_p,iType,POINTER(iType)]
  def Target(self,localEdge):
    target = iType()
    lib.ElDistGraphTarget(self.obj,localEdge,pointer(target))
    return target.value

  lib.ElDistGraphSourceOffset.argtypes = [c_void_p,iType,POINTER(iType)]
  def SourceOffset(self,localSource):
    offset = iType()
    lib.ElDistGraphSourceOffset(self.obj,source,pointer(offset))
    return offset.value

  lib.ElDistGraphOffset.argtypes = [c_void_p,iType,iType,POINTER(iType)]
  def Offset(self,localSource,target):
    offset = iType()
    lib.ElDistGraphOffset(self.obj,source,target,pointer(offset))
    return offset.value

  lib.ElDistGraphNumConnections.argtypes = [c_void_p,iType,POINTER(iType)]
  def NumConnections(self,localSource):
    numConnections = iType()
    lib.ElDistGraphNumConnections(self.obj,localSource,pointer(numConnections))
    return numConnections.value

  lib.ElDistGraphImbalance.argtypes = [c_void_p,POINTER(dType)]
  def Imbalance(self):
    imbalance = dType()
    lib.ElDistGraphImbalance(self.obj,pointer(imbalance))
    return imbalance.value

  lib.ElDistGraphSourceBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
  def SourceBuffer(self):
    sourceBuf = POINTER(iType)()
    lib.ElDistGraphSourceBuffer(self.obj,pointer(sourceBuf))
    return sourceBuf

  lib.ElDistGraphLockedSourceBuffer.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def LockedSourceBuffer(self):
    sourceBuf = POINTER(iType)()
    lib.ElDistGraphLockedSourceBuffer(self.obj,pointer(sourceBuf))
    return sourceBuf

  lib.ElDistGraphTargetBuffer.argtypes = \
  lib.ElDistGraphLockedTargetBuffer.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def TargetBuffer(self,locked=False):
    targetBuf = POINTER(iType)()
    if locked:
      lib.ElDistGraphLockedTargetBuffer(self.obj,pointer(targetBuf))
    else:
      lib.ElDistGraphTargetBuffer(self.obj,pointer(targetBuf))
    return targetBuf
