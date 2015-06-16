#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time

n0=10
n1=10
display = False
worldSize = El.mpi.WorldSize()
worldRank = El.mpi.WorldRank()

# Place two 2D finite-difference matrices next to each other
# and make the last column dense
def ConcatFD2D(N0,N1):
  A = El.DistSparseMatrix(El.zTag)
  height = N0*N1
  width = 2*N0*N1
  A.Resize(height,width)
  localHeight = A.LocalHeight()
  A.Reserve(11*localHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    x0 = s % N0
    x1 = s / N0
    sRel = s + N0*N1

    A.QueueUpdate( s, s,     11+1j )
    A.QueueUpdate( s, sRel, -20+2j )
    if x0 > 0:
      A.QueueUpdate( s, s-1,    -1+3j  )
      A.QueueUpdate( s, sRel-1, -17+4j )
    if x0+1 < N0:
      A.QueueUpdate( s, s+1,     2+5j  )
      A.QueueUpdate( s, sRel+1, -20+6j )
    if x1 > 0:
      A.QueueUpdate( s, s-N0,    -30+7j )
      A.QueueUpdate( s, sRel-N0, -3+8j  )
    if x1+1 < N1:
      A.QueueUpdate( s, s+N0,    4+9j )
      A.QueueUpdate( s, sRel+N0, 3+10j )

    # The dense last column
    A.QueueUpdate( s, width-1, -10/height );

  A.ProcessLocalQueues()
  return A

A = ConcatFD2D(n0,n1)
b = El.DistMultiVec(El.zTag)
#El.Gaussian( b, n0*n1, 1 )
El.Ones( b, n0*n1, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

ctrl = El.BPCtrl_z()
ctrl.ipmCtrl.mehrotraCtrl.time = True
ctrl.ipmCtrl.mehrotraCtrl.progress = True
ctrl.ipmCtrl.mehrotraCtrl.outerEquil = False
ctrl.ipmCtrl.mehrotraCtrl.innerEquil = True
startBP = time.clock()
x = El.BP( A, b, ctrl )
endBP = time.clock()
if worldRank == 0:
  print "BP time:", endBP-startBP, "seconds"
if display:
  El.Display( x, "x" )

xOneNorm = El.EntrywiseNorm( x, 1 )
e = El.DistMultiVec(El.zTag)
El.Copy( b, e )
El.SparseMultiply( El.NORMAL, -1., A, x, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print "|| x ||_1       =", xOneNorm
  print "|| A x - b ||_2 =", eTwoNorm

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
