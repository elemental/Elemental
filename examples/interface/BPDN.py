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
numLambdas = 3
startLambda = 0
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
  El.Display( A[0:n0*n1,0:n0*n1], "AL" )
  El.Display( A[0:n0*n1,n0*n1:2*n0*n1], "AR" )
  El.Display( b, "b" )

ctrl = El.BPDNCtrl_d()
ctrl.ipmCtrl.mehrotraCtrl.time = True
ctrl.ipmCtrl.mehrotraCtrl.progress = True
ctrl.ipmCtrl.mehrotraCtrl.solveCtrl.progress = True

for j in xrange(0,numLambdas):
  lambd = startLambda + j*(endLambda-startLambda)/(numLambdas-1.)
  if worldRank == 0:
    print('lambda = {}'.format(lambd))

  startBPDN = El.mpi.Time()
  x = El.BPDN( A, b, lambd, ctrl )
  endBPDN = El.mpi.Time()
  if worldRank == 0:
    print('BPDN time: {} seconds'.format(endBPDN-startBPDN))
  if display:
    El.Display( x, "x" )

  xOneNorm = El.EntrywiseNorm( x, 1 )
  e = El.DistMultiVec()
  El.Copy( b, e )
  El.Multiply( El.NORMAL, -1., A, x, 1., e )
  if display:
    El.Display( e, "e" )
  eTwoNorm = El.Nrm2( e )
  if worldRank == 0:
    print('|| x ||_1       = {}'.format(xOneNorm))
    print('|| A x - b ||_2 = {}'.format(eTwoNorm))

El.Finalize()
