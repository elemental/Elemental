#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = 75
n1 = 75
lambda1 = 3
lambda2 = 4
output = False
display = False
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

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
if output:
  El.Print( A, "A" )
  El.Print( b, "b" )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

if worldRank == 0:
  print('lambda1 = {}, lambda2 = {}'.format(lambda1,lambda2))

ctrl = El.QPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True
ctrl.mehrotraCtrl.time = True
ctrl.mehrotraCtrl.solveCtrl.progress = True
ctrl.mehrotraCtrl.solveCtrl.time = True

# Solve *with* resolving the regularization
ctrl.mehrotraCtrl.resolveReg = True
startEN = El.mpi.Time()
x = El.EN( A, b, lambda1, lambda2, ctrl )
endEN = El.mpi.Time()
if worldRank == 0:
  print('EN time (resolve reg.): {} seconds'.format(endEN-startEN))
if display:
  El.Display( x, "x" )

# Solve *without* resolving the regularization
ctrl.mehrotraCtrl.resolveReg = False
startEN = El.mpi.Time()
x = El.EN( A, b, lambda1, lambda2, ctrl )
endEN = El.mpi.Time()
if worldRank == 0:
  print('EN time (no resolve reg.): {} seconds'.format(endEN-startEN))
if display:
  El.Display( x, "x" )

xOneNorm = El.EntrywiseNorm( x, 1 )
xTwoNorm = El.Nrm2( x )
e = El.DistMultiVec()
El.Copy( b, e )
El.Multiply( El.NORMAL, -1., A, x, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print('|| x ||_1       = {}'.format(xOneNorm))
  print('|| x ||_2       = {}'.format(xTwoNorm))
  print('|| A x - b ||_2 = {}'.format(eTwoNorm))

El.Finalize()
