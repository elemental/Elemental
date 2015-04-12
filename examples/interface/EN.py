#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El
import time

n0 = 50
n1 = 50
lambda1 = 3
lambda2 = 4
display = True
worldRank = El.mpi.WorldRank()

# Place a 2D finite-difference matrices next to an identity
# and then make the last column dense
#
# NOTE: Increasing the magnitudes of the off-diagonal entries causes the matrix
#       to become significantly more ill-conditioned and stresses the 
#       effective condition number of the SQD solves.
def ConcatFD2D(N0,N1):
  A = El.DistSparseMatrix()
  height = N0*N1
  width = 2*N0*N1
  A.Resize(height,width)
  localHeight = A.LocalHeight()
  A.Reserve(7*localHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    x0 = s % N0
    x1 = s / N0

    # The finite-difference stencil
    A.QueueLocalUpdate( sLoc, s, 15 )
    if x0 > 0:
      A.QueueLocalUpdate( sLoc, s-1,  -1 )
    if x0+1 < N0:
      A.QueueLocalUpdate( sLoc, s+1,   2 )
    if x1 > 0:
      A.QueueLocalUpdate( sLoc, s-N0, -3 )
    if x1+1 < N1:
      A.QueueLocalUpdate( sLoc, s+N0,  4 )

    # The identity
    sRel = s + N0*N1
    A.QueueLocalUpdate( sLoc, sRel,  1  )

    # The dense last column
    A.QueueLocalUpdate( sLoc, width-1, -10/height );

  A.ProcessQueues()
  return A

A = ConcatFD2D(n0,n1)
b = El.DistMultiVec()
El.Gaussian( b, n0*n1, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

ctrl = El.QPAffineCtrl_d()
ctrl.mehrotraCtrl.outerEquil = True
ctrl.mehrotraCtrl.innerEquil = True
ctrl.mehrotraCtrl.scaleTwoNorm = True
ctrl.mehrotraCtrl.basisSize = 15
ctrl.mehrotraCtrl.progress = True
ctrl.mehrotraCtrl.qsdCtrl.progress = True
if worldRank == 0:
  print "lambda1 =", lambda1, "lambda2 =", lambda2
startEN = time.clock()
x = El.EN( A, b, lambda1, lambda2, ctrl )
endEN = time.clock()
if worldRank == 0:
  print "EN time:", endEN-startEN, "seconds"
if display:
  El.Display( x, "x" )

xOneNorm = El.EntrywiseNorm( x, 1 )
xTwoNorm = El.Nrm2( x )
e = El.DistMultiVec()
El.Copy( b, e )
El.SparseMultiply( El.NORMAL, -1., A, x, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print "|| x ||_1       =", xOneNorm
  print "|| x ||_2       =", xTwoNorm
  print "|| A x - b ||_2 =", eTwoNorm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
