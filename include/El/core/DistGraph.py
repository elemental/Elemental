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

# DistGraph
# =========

lib.ElDistGraphCreate.argtypes = [POINTER(c_void_p),mpi.Comm]
lib.ElDistGraphCreate.restype = c_uint

lib.ElDistGraphDestroy.argtypes = [c_void_p]
lib.ElDistGraphDestroy.restype = c_uint

lib.ElDistGraphEmpty.argtypes = [c_void_p]
lib.ElDistGraphEmpty.restype = c_uint

lib.ElDistGraphResize.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphResize.restype = c_uint

lib.ElDistGraphSetComm.argtypes = [c_void_p,mpi.Comm]
lib.ElDistGraphSetComm.restype = c_uint

lib.ElDistGraphReserve.argtypes = [c_void_p,iType]
lib.ElDistGraphReserve.restype = c_uint

lib.ElDistGraphConnect.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphConnect.restype = c_uint

lib.ElDistGraphConnectLocal.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphConnectLocal.restype = c_uint

lib.ElDistGraphDisconnect.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphDisconnect.restype = c_uint

lib.ElDistGraphDisconnectLocal.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphDisconnectLocal.restype = c_uint

lib.ElDistGraphQueueConnection.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphQueueConnection.restype = c_uint

lib.ElDistGraphQueueLocalConnection.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphQueueLocalConnection.restype = c_uint

lib.ElDistGraphQueueDisconnection.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphQueueDisconnection.restype = c_uint

lib.ElDistGraphQueueLocalDisconnection.argtypes = [c_void_p,iType,iType]
lib.ElDistGraphQueueLocalDisconnection.restype = c_uint

lib.ElDistGraphMakeConsistent.argtypes = [c_void_p]
lib.ElDistGraphMakeConsistent.restype = c_uint

lib.ElDistGraphNumSources.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistGraphNumSources.restype = c_uint

lib.ElDistGraphNumTargets.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistGraphNumTargets.restype = c_uint

lib.ElDistGraphFirstLocalSource.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistGraphFirstLocalSource.restype = c_uint

lib.ElDistGraphNumLocalSources.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistGraphNumLocalSources.restype = c_uint

lib.ElDistGraphNumLocalEdges.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistGraphNumLocalEdges.restype = c_uint

lib.ElDistGraphCapacity.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistGraphCapacity.restype = c_uint

lib.ElDistGraphConsistent.argtypes = [c_void_p,POINTER(bType)]
lib.ElDistGraphConsistent.restype = c_uint

lib.ElDistGraphComm.argtypes = [c_void_p,POINTER(mpi.Comm)]
lib.ElDistGraphComm.restype = c_uint

lib.ElDistGraphBlocksize.argtypes = [c_void_p,POINTER(iType)]
lib.ElDistGraphBlocksize.restype = c_uint

lib.ElDistGraphSource.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistGraphSource.restype = c_uint

lib.ElDistGraphTarget.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistGraphTarget.restype = c_uint

lib.ElDistGraphEdgeOffset.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistGraphEdgeOffset.restype = c_uint

lib.ElDistGraphNumConnections.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElDistGraphNumConnections.restype = c_uint

lib.ElDistGraphSourceBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElDistGraphSourceBuffer.restype = c_uint

lib.ElDistGraphLockedSourceBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElDistGraphLockedSourceBuffer.restype = c_uint

lib.ElDistGraphTargetBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElDistGraphTargetBuffer.restype = c_uint

lib.ElDistGraphLockedTargetBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElDistGraphLockedTargetBuffer.restype = c_uint

class DistGraph(object):
  # Constructors and destructors
  # ============================
  def __init__(self,comm=mpi.COMM_WORLD(),create=True):
    self.obj = c_void_p()
    if create:
      lib.ElDistGraphCreate(pointer(self.obj),comm)
  def Destroy(self):
    lib.ElDistGraphDestroy(self.obj)
  # Assignment and reconfiguration
  # ==============================
  def Empty(self):
    lib.ElDistGraphEmpty(self.obj)
  def Resize(self,numSources,numTargets=None):
    if numTargets == None: 
      numTargets = numSources
    lib.ElDistGraphResize(self.obj,numSources,numTargets)
  def SetComm(self,comm):
    lib.ElDistGraphSetComm(self,comm)
  def Reserve(self,numEdges):
    lib.ElDistGraphReserve(self.obj,numEdges)
  def Connect(self,source,target):
    lib.ElDistGraphConnect(self.obj,source,target)
  def ConnectLocal(self,localSource,target):
    lib.ElDistGraphConnectLocal(self.obj,localSource,target)
  def Disconnect(self,source,target):
    lib.ElDistGraphDisconnect(self.obj,source,target)
  def DisconnectLocal(self,localSource,target):
    lib.ElDistGraphDisconnectLocal(self.obj,localSource,target)
  def QueueConnection(self,source,target):
    lib.ElDistGraphQueueConnection(self.obj,source,target)
  def QueueLocalConnection(self,localSource,target):
    lib.ElDistGraphQueueLocalConnection(self.obj,localSource,target)
  def QueueDisconnection(self,source,target):
    lib.ElDistGraphQueueDisconnection(self.obj,source,target)
  def QueueLocalDisconnection(self,localSource,target):
    lib.ElDistGraphQueueLocalDisconnection(self.obj,localSource,target)
  def MakeConsistent(self):
    lib.ElDistGraphMakeConsistent(self.obj)
  # Queries
  # =======
  def NumSources(self):
    numSources = iType()
    lib.ElDistGraphNumSources(self.obj,pointer(numSources))
    return numSources.value
  def NumTargets(self):
    numTargets = iType()
    lib.ElDistGraphNumTargets(self.obj,pointer(numTargets))
    return numTargets.value
  def FirstLocalSource(self):
    firstLocalSource = iType()
    lib.ElDistGraphFirstLocalSource(self.obj,pointer(firstLocalSource))
    return firstLocalSource.value
  def NumLocalSources(self):
    numLocalSources = iType()
    lib.ElDistGraphNumLocalSources(self.obj,pointer(numLocalSources))
    return numLocalSources.value
  def NumLocalEdges(self):
    numLocalEdges = iType()
    lib.ElDistGraphNumLocalEdges(self.obj,pointer(numLocalEdges))
    return numLocalEdges.value
  def Capacity(self):
    capacity = iType()
    lib.ElDistGraphCapacity(self.obj,pointer(capacity))
    return capacity.value
  def Consistent(self):
    consistent = bType()
    lib.ElDistGraphConsistent(self.obj,pointer(consistent))
    return consistent.value
  def Comm(self):
    comm = mpi.Comm()
    lib.ElDistGraphComm(self.obj,pointer(comm))
    return comm
  def Blocksize(self):
    blocksize = iType()
    lib.ElDistGraphBlocksize(self.obj,pointer(blocksize))
    return blocksize.value
  def Source(self,localEdge):
    source = iType()
    lib.ElDistGraphSource(self.obj,localEdge,pointer(source))
    return source.value
  def Target(self,localEdge):
    target = iType()
    lib.ElDistGraphTarget(self.obj,localEdge,pointer(target))
    return target.value
  def EdgeOffset(self,localSource):
    localEdgeOffset = iType()
    lib.ElDistGraphEdgeOffset(self.obj,source,pointer(localEdgeOffset))
    return localEdgeOffset.value
  def NumConnections(self,localSource):
    numConnections = iType()
    lib.ElDistGraphNumConnections(self.obj,localSource,pointer(numConnections))
    return numConnections.value
  def SourceBuffer(self):
    sourceBuf = POINTER(iType)()
    lib.ElDistGraphSourceBuffer(self.obj,pointer(sourceBuf))
    return sourceBuf
  def LockedSourceBuffer(self):
    sourceBuf = POINTER(iType)()
    lib.ElDistGraphLockedSourceBuffer(self.obj,pointer(sourceBuf))
    return sourceBuf
  def TargetBuffer(self):
    targetBuf = POINTER(iType)()
    lib.ElDistGraphTargetBuffer(self.obj,pointer(targetBuf))
    return targetBuf
  def LockedTargetBuffer(self):
    targetBuf = POINTER(iType)()
    lib.ElDistGraphLockedTargetBuffer(self.obj,pointer(targetBuf))
    return targetBuf
