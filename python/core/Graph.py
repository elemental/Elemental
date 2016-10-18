#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
from environment import *

class Graph(object):
  # Constructors and destructors
  # ============================
  lib.ElGraphCreate.argtypes = [POINTER(c_void_p)]
  def __init__(self,create=True):
    self.obj = c_void_p()
    if create:
      lib.ElGraphCreate(pointer(self.obj))

  lib.ElGraphDestroy.argtypes = [c_void_p]
  def Destroy(self):
    lib.ElGraphDestroy(self.obj)

  # Assignment and reconfiguration
  # ==============================
  lib.ElGraphEmpty.argtypes = [c_void_p]
  def Empty(self):
    lib.ElGraphEmpty(self.obj)

  lib.ElGraphResize.argtypes = [c_void_p,iType,iType]
  def Resize(self,numSources,numTargets=None):
    if numTargets == None: 
      numTargets = numSources
    lib.ElGraphResize(self.obj,numSources,numTargets)

  lib.ElGraphReserve.argtypes = [c_void_p,iType]
  def Reserve(self,numEdges):
    lib.ElGraphReserve(self.obj,numEdges)

  lib.ElGraphConnect.argtypes = [c_void_p,iType,iType]
  def Connect(self,source,target):
    lib.ElGraphConnect(self.obj,source,target)

  lib.ElGraphDisconnect.argtypes = [c_void_p,iType,iType]
  def Disconnect(self,source,target):
    lib.ElGraphDisconnect(self.obj,source,target)

  lib.ElGraphQueueConnection.argtypes = [c_void_p,iType,iType]
  def QueueConnection(self,source,target):
    lib.ElGraphQueueConnection(self.obj,source,target)

  lib.ElGraphQueueDisconnection.argtypes = [c_void_p,iType,iType]
  def QueueDisconnection(self,source,target):
    lib.ElGraphQueueDisconnection(self.obj,source,target)

  lib.ElGraphProcessQueues.argtypes = [c_void_p]
  def ProcessQueues(self):
    lib.ElGraphProcessQueues(self.obj)

  # Queries
  # =======
  lib.ElGraphNumSources.argtypes = [c_void_p,POINTER(iType)]
  def NumSources(self):
    numSources = iType()
    lib.ElGraphNumSources(self.obj,pointer(numSources))
    return numSources.value

  lib.ElGraphNumTargets.argtypes = [c_void_p,POINTER(iType)]
  def NumTargets(self):
    numTargets = iType()
    lib.ElGraphNumTargets(self.obj,pointer(numTargets))
    return numTargets.value

  lib.ElGraphNumEdges.argtypes = [c_void_p,POINTER(iType)]
  def NumEdges(self):
    numEdges = iType()
    lib.ElGraphNumEdges(self.obj,pointer(numEdges))
    return numEdges.value

  lib.ElGraphCapacity.argtypes = [c_void_p,POINTER(iType)]
  def Capacity(self):
    capacity = iType()
    lib.ElGraphCapacity(self.obj,pointer(capacity))
    return capacity.value

  lib.ElGraphConsistent.argtypes = [c_void_p,POINTER(bType)]
  def Consistent(self):
    consistent = bType()
    lib.ElGraphConsistent(self.obj,pointer(consistent))
    return consistent.value

  lib.ElGraphSource.argtypes = [c_void_p,iType,POINTER(iType)]
  def Source(self,edge):
    source = iType()
    lib.ElGraphSource(self.obj,edge,pointer(source))
    return source.value

  lib.ElGraphTarget.argtypes = [c_void_p,iType,POINTER(iType)]
  def Target(self,edge):
    target = iType()
    lib.ElGraphTarget(self.obj,edge,pointer(target))
    return target.value

  lib.ElGraphSourceOffset.argtypes = [c_void_p,iType,POINTER(iType)]
  def SourceOffset(self,source):
    sourceOffset = iType()
    lib.ElGraphSourceOffset(self.obj,source,pointer(sourceOffset))
    return sourceOffset.value

  lib.ElGraphOffset.argtypes = [c_void_p,iType,iType,POINTER(iType)]
  def Offset(self,source,target):
    offset = iType()
    lib.ElGraphOffset(self.obj,source,target,pointer(offset))
    return offset.value

  lib.ElGraphNumConnections.argtypes = [c_void_p,iType,POINTER(iType)]
  def NumConnections(self,source):
    numConnections = iType()
    lib.ElGraphNumConnections(self.obj,source,pointer(numConnections))
    return numConnections.value

  lib.ElGraphSourceBuffer.argtypes = \
  lib.ElGraphLockedSourceBuffer.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def SourceBuffer(self,locked=False):
    sourceBuf = POINTER(iType)()
    if locked:
      lib.ElGraphLockedSourceBuffer(self.obj,pointer(sourceBuf))
    else:
      lib.ElGraphSourceBuffer(self.obj,pointer(sourceBuf))
    return sourceBuf

  lib.ElGraphTargetBuffer.argtypes = \
  lib.ElGraphLockedTargetBuffer.argtypes = \
    [c_void_p,POINTER(POINTER(iType))]
  def TargetBuffer(self):
    targetBuf = POINTER(iType)()
    if locked:
      lib.ElGraphLockedTargetBuffer(self.obj,pointer(targetBuf))
    else:
      lib.ElGraphTargetBuffer(self.obj,pointer(targetBuf))
    return targetBuf

  # TODO: __getitem__
