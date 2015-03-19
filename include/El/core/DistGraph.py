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

class DistGraph(object):
  # Constructors and destructors
  # ============================
  lib.ElDistGraphCreate.argtypes = [POINTER(c_void_p),mpi.Comm]
  def __init__(self,comm=mpi.COMM_WORLD(),create=True):
    self.obj = c_void_p()
    if create:
      lib.ElDistGraphCreate(pointer(self.obj),comm)

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

  lib.ElDistGraphSetComm.argtypes = [c_void_p,mpi.Comm]
  def SetComm(self,comm):
    lib.ElDistGraphSetComm(self,comm)

  lib.ElDistGraphReserve.argtypes = [c_void_p,iType]
  def Reserve(self,numEdges):
    lib.ElDistGraphReserve(self.obj,numEdges)

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

  lib.ElDistGraphQueueConnection.argtypes = [c_void_p,iType,iType]
  def QueueConnection(self,source,target):
    lib.ElDistGraphQueueConnection(self.obj,source,target)

  lib.ElDistGraphQueueLocalConnection.argtypes = [c_void_p,iType,iType]
  def QueueLocalConnection(self,localSource,target):
    lib.ElDistGraphQueueLocalConnection(self.obj,localSource,target)

  lib.ElDistGraphQueueDisconnection.argtypes = [c_void_p,iType,iType]
  def QueueDisconnection(self,source,target):
    lib.ElDistGraphQueueDisconnection(self.obj,source,target)

  lib.ElDistGraphQueueLocalDisconnection.argtypes = [c_void_p,iType,iType]
  def QueueLocalDisconnection(self,localSource,target):
    lib.ElDistGraphQueueLocalDisconnection(self.obj,localSource,target)

  lib.ElDistGraphMakeConsistent.argtypes = [c_void_p]
  def MakeConsistent(self):
    lib.ElDistGraphMakeConsistent(self.obj)

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

  lib.ElDistGraphConsistent.argtypes = [c_void_p,POINTER(bType)]
  def Consistent(self):
    consistent = bType()
    lib.ElDistGraphConsistent(self.obj,pointer(consistent))
    return consistent.value

  lib.ElDistGraphComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
  def Comm(self):
    comm = mpi.Comm()
    lib.ElDistGraphComm(self.obj,pointer(comm))
    return comm

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

  lib.ElDistGraphEdgeOffset.argtypes = [c_void_p,iType,POINTER(iType)]
  def EdgeOffset(self,localSource):
    localEdgeOffset = iType()
    lib.ElDistGraphEdgeOffset(self.obj,source,pointer(localEdgeOffset))
    return localEdgeOffset.value

  lib.ElDistGraphNumConnections.argtypes = [c_void_p,iType,POINTER(iType)]
  def NumConnections(self,localSource):
    numConnections = iType()
    lib.ElDistGraphNumConnections(self.obj,localSource,pointer(numConnections))
    return numConnections.value

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
