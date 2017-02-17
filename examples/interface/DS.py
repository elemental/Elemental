#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = 25
n1 = 25
numLambdas = 5
startLambda = 0.01
endLambda = 1
display = False
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# Place two 2D finite-difference matrices next to each other
# and make the last column dense
def ConcatFD2D(N0,N1):
  A = El.DistSparseMatrix()
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

    A.QueueLocalUpdate( sLoc, s,     11 )
    A.QueueLocalUpdate( sLoc, sRel, -20 )
    if x0 > 0:
      A.QueueLocalUpdate( sLoc, s-1,    -1  )
      A.QueueLocalUpdate( sLoc, sRel-1, -17 )
    if x0+1 < N0:
      A.QueueLocalUpdate( sLoc, s+1,     2  )
      A.QueueLocalUpdate( sLoc, sRel+1, -20 )
    if x1 > 0:
      A.QueueLocalUpdate( sLoc, s-N0,    -30 )
      A.QueueLocalUpdate( sLoc, sRel-N0, -3  )
    if x1+1 < N1:
      A.QueueLocalUpdate( sLoc, s+N0,    4 )
      A.QueueLocalUpdate( sLoc, sRel+N0, 3 )

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

ctrl = El.LPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True

for j in xrange(0,numLambdas):
  lambd = startLambda + j*(endLambda-startLambda)/(numLambdas-1.)
  if worldRank == 0:
    print('lambda = {}'.format(lambd))

  startDS = El.mpi.Time()
  x = El.DS( A, b, lambd, ctrl )
  endDS = El.mpi.Time()
  if worldRank == 0:
    print('DS time: {} seconds'.format(endDS-startDS))
  if display:
    El.Display( x, "x" )

  xOneNorm = El.EntrywiseNorm( x, 1 )
  r = El.DistMultiVec()
  El.Copy( b, r )
  El.Multiply( El.NORMAL, -1., A, x, 1., r )
  rTwoNorm = El.Nrm2( r )
  t = El.DistMultiVec()
  El.Zeros( t, 2*n0*n1, 1 )
  El.Multiply( El.TRANSPOSE, 1., A, r, 0., t )
  tTwoNorm = El.Nrm2( t )
  tInfNorm = El.MaxNorm( t )
  if display:
    El.Display( r, "r" )
    El.Display( t, "t" )
  if worldRank == 0:
    print('|| x ||_1       = {}'.format(xOneNorm))
    print('|| b - A x ||_2 = {}'.format(rTwoNorm))
    print('|| A^T (b - A x) ||_2 = {}'.format(tTwoNorm))
    print('|| A^T (b - A x) ||_oo = {}'.format(tInfNorm))

El.Finalize()
