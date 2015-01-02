#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *

# Graph
# ======

lib.ElGraphCreate.argtypes = [POINTER(c_void_p)]
lib.ElGraphCreate.restype = c_uint

lib.ElGraphDestroy.argtypes = [c_void_p]
lib.ElGraphDestroy.restype = c_uint

lib.ElGraphEmpty.argtypes = [c_void_p]
lib.ElGraphEmpty.restype = c_uint

lib.ElGraphResize.argtypes = [c_void_p,iType,iType]
lib.ElGraphResize.restype = c_uint

lib.ElGraphReserve.argtypes = [c_void_p,iType]
lib.ElGraphReserve.restype = c_uint

lib.ElGraphConnect.argtypes = [c_void_p,iType,iType]
lib.ElGraphConnect.restype = c_uint

lib.ElGraphDisconnect.argtypes = [c_void_p,iType,iType]
lib.ElGraphDisconnect.restype = c_uint

lib.ElGraphQueueConnection.argtypes = [c_void_p,iType,iType]
lib.ElGraphQueueConnection.restype = c_uint

lib.ElGraphQueueDisconnection.argtypes = [c_void_p,iType,iType]
lib.ElGraphQueueDisconnection.restype = c_uint

lib.ElGraphMakeConsistent.argtypes = [c_void_p]
lib.ElGraphMakeConsistent.restype = c_uint

lib.ElGraphNumSources.argtypes = [c_void_p,POINTER(iType)]
lib.ElGraphNumSources.restype = c_uint

lib.ElGraphNumTargets.argtypes = [c_void_p,POINTER(iType)]
lib.ElGraphNumTargets.restype = c_uint

lib.ElGraphNumEdges.argtypes = [c_void_p,POINTER(iType)]
lib.ElGraphNumEdges.restype = c_uint

lib.ElGraphCapacity.argtypes = [c_void_p,POINTER(iType)]
lib.ElGraphCapacity.restype = c_uint

lib.ElGraphConsistent.argtypes = [c_void_p,POINTER(bType)]
lib.ElGraphConsistent.restype = c_uint

lib.ElGraphSource.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElGraphSource.restype = c_uint

lib.ElGraphTarget.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElGraphTarget.restype = c_uint

lib.ElGraphEdgeOffset.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElGraphEdgeOffset.restype = c_uint

lib.ElGraphNumConnections.argtypes = [c_void_p,iType,POINTER(iType)]
lib.ElGraphNumConnections.restype = c_uint

lib.ElGraphSourceBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElGraphSourceBuffer.restype = c_uint

lib.ElGraphLockedSourceBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElGraphLockedSourceBuffer.restype = c_uint

lib.ElGraphTargetBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElGraphTargetBuffer.restype = c_uint

lib.ElGraphLockedTargetBuffer.argtypes = [c_void_p,POINTER(POINTER(iType))]
lib.ElGraphLockedTargetBuffer.restype = c_uint

class Graph(object):
  # Constructors and destructors
  # ============================
  def __init__(self,create=True):
    self.obj = c_void_p()
    if create:
      lib.ElGraphCreate(pointer(self.obj))
  def Destroy(self):
    lib.ElGraphDestroy(self.obj)
  # Assignment and reconfiguration
  # ==============================
  def Empty(self):
    lib.ElGraphEmpty(self.obj)
  def Resize(self,numSources,numTargets=None):
    if numTargets == None: 
      numTargets = numSources
    lib.ElGraphResize(self.obj,numSources,numTargets)
  def Reserve(self,numEdges):
    lib.ElGraphReserve(self.obj,numEdges)
  def Connect(self,source,target):
    lib.ElGraphConnect(self.obj,source,target)
  def Disconnect(self,source,target):
    lib.ElGraphDisconnect(self.obj,source,target)
  def QueueConnection(self,source,target):
    lib.ElGraphQueueConnection(self.obj,source,target)
  def QueueDisconnection(self,source,target):
    lib.ElGraphQueueDisconnection(self.obj,source,target)
  def MakeConsistent(self):
    lib.ElGraphMakeConsistent(self.obj)
  # Queries
  # =======
  def NumSources(self):
    numSources = iType()
    lib.ElGraphNumSources(self.obj,pointer(numSources))
    return numSources.value
  def NumTargets(self):
    numTargets = iType()
    lib.ElGraphNumTargets(self.obj,pointer(numTargets))
    return numTargets.value
  def NumEdges(self):
    numEdges = iType()
    lib.ElGraphNumEdges(self.obj,pointer(numEdges))
    return numEdges.value
  def Capacity(self):
    capacity = iType()
    lib.ElGraphCapacity(self.obj,pointer(capacity))
    return capacity.value
  def Consistent(self):
    consistent = bType()
    lib.ElGraphConsistent(self.obj,pointer(consistent))
    return consistent.value
  def Source(self,edge):
    source = iType()
    lib.ElGraphSource(self.obj,edge,pointer(source))
    return source.value
  def Target(self,edge):
    target = iType()
    lib.ElGraphTarget(self.obj,edge,pointer(target))
    return target.value
  def EdgeOffset(self,source):
    edgeOffset = iType()
    lib.ElGraphEdgeOffset(self.obj,source,pointer(edgeOffset))
    return edgeOffset.value
  def NumConnections(self,source):
    numConnections = iType()
    lib.ElGraphNumConnections(self.obj,source,pointer(numConnections))
    return numConnections.value
  def SourceBuffer(self):
    sourceBuf = POINTER(iType)()
    lib.ElGraphSourceBuffer(self.obj,pointer(sourceBuf))
    return sourceBuf
  def LockedSourceBuffer(self):
    sourceBuf = POINTER(iType)()
    lib.ElGraphLockedSourceBuffer(self.obj,pointer(sourceBuf))
    return sourceBuf
  def TargetBuffer(self):
    targetBuf = POINTER(iType)()
    lib.ElGraphTargetBuffer(self.obj,pointer(targetBuf))
    return targetBuf
  def LockedTargetBuffer(self):
    targetBuf = POINTER(iType)()
    lib.ElGraphLockedTargetBuffer(self.obj,pointer(targetBuf))
    return targetBuf
