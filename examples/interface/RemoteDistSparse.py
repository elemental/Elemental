#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time

n0 = 50
n1 = 50
display = True
output = True
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# Stack two 2D finite-difference matrices on top of each other
# and make the last column dense
def RemoteStackedFD2D(N0,N1):
  A = El.DistSparseMatrix()
  height = 2*N0*N1
  width = N0*N1
  A.Resize(height,width)
  customLocalHeight = (height/worldSize) + 1
  A.Reserve(6*customLocalHeight,6*customLocalHeight)
  for s in xrange(worldRank,height,worldSize):
    if s < N0*N1:
      x0 = s % N0
      x1 = s / N0
      A.QueueUpdate( s, s, 1 )
      if x0 > 0:
        A.QueueUpdate( s, s-1, -1 )
      if x0+1 < N0:
        A.QueueUpdate( s, s+1, 2 )
      if x1 > 0:
        A.QueueUpdate( s, s-N0, -3 )
      if x1+1 < N1:
        A.QueueUpdate( s, s+N0, 4 )
    else:
      sRel = s-N0*N1
      x0 = sRel % N0
      x1 = sRel / N0
      A.QueueUpdate( s, sRel, -2 )
      if x0 > 0:
        A.QueueUpdate( s, sRel-1, -1 )
      if x0+1 < N0:
        A.QueueUpdate( s, sRel+1, -2 )
      if x1 > 0:
        A.QueueUpdate( s, sRel-N0, -3 )
      if x1+1 < N1:
        A.QueueUpdate( s, sRel+N0, 3 )

  A.ProcessQueues()
  return A

A = RemoteStackedFD2D(n0,n1)
if display:
  El.Display( A, "A" )
if output:
  El.Print( A, "A" )

El.Finalize()
