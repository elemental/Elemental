#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = 50
n1 = 50
display = False
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# Stack two 2D finite-difference matrices on top of each other
# and make the last column dense
def StackedFD2D(N0,N1):
  A = El.DistSparseMatrix()
  height = 2*N0*N1
  width = N0*N1
  A.Resize(height,width)
  localHeight = A.LocalHeight()
  A.Reserve(6*localHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    if s < N0*N1:
      x0 = s % N0
      x1 = s / N0
      A.QueueLocalUpdate( sLoc, s, 11 )
      if x0 > 0:
        A.QueueLocalUpdate( sLoc, s-1, -1 )
      if x0+1 < N0:
        A.QueueLocalUpdate( sLoc, s+1, 2 )
      if x1 > 0:
        A.QueueLocalUpdate( sLoc, s-N0, -3 )
      if x1+1 < N1:
        A.QueueLocalUpdate( sLoc, s+N0, 4 )
    else:
      sRel = s-N0*N1
      x0 = sRel % N0
      x1 = sRel / N0
      A.QueueLocalUpdate( sLoc, sRel, -2 )
      if x0 > 0:
        A.QueueLocalUpdate( sLoc, sRel-1, -1 )
      if x0+1 < N0:
        A.QueueLocalUpdate( sLoc, sRel+1, -2 )
      if x1 > 0:
        A.QueueLocalUpdate( sLoc, sRel-N0, -3 )
      if x1+1 < N1:
        A.QueueLocalUpdate( sLoc, sRel+N0, 3 )

    # The dense last column
    A.QueueLocalUpdate( sLoc, width-1, -10/height );

  A.ProcessQueues()
  return A

A = StackedFD2D(n0,n1)
b = El.DistMultiVec()
El.Gaussian( b, 2*n0*n1, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

ctrl = El.LPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True
ctrl.mehrotraCtrl.solveCtrl.progress = True
startCP = El.mpi.Time()
x = El.CP( A, b, ctrl )
endCP = El.mpi.Time()
if worldRank == 0:
  print('CP time: {} seconds'.format(endCP-startCP))
if display:
  El.Display( x, "x" )

bTwoNorm = El.Nrm2( b )
bInfNorm = El.MaxNorm( b )
r = El.DistMultiVec()
El.Copy( b, r )
El.Multiply( El.NORMAL, -1., A, x, 1., r )
if display:
  El.Display( r, "r" )
rTwoNorm = El.Nrm2( r )
rInfNorm = El.MaxNorm( r )
if worldRank == 0:
  print('|| b ||_2        = {}'.format(bTwoNorm))
  print('|| b ||_oo       = {}'.format(bInfNorm))
  print('|| A x - b ||_2  = {}'.format(rTwoNorm))
  print('|| A x - b ||_oo = {}'.format(rInfNorm))

startLS = El.mpi.Time()
xLS = El.LeastSquares(A,b)
endLS = El.mpi.Time()
if worldRank == 0:
  print('LS time: {} seconds'.format(endLS-startLS))
if display:
  El.Display( xLS, "x_{LS}" )
rLS = El.DistMultiVec()
El.Copy( b, rLS )
El.Multiply( El.NORMAL, -1., A, xLS, 1., rLS )
if display:
  El.Display( rLS, "A x_{LS} - b" )
rLSTwoNorm = El.Nrm2(rLS)
rLSInfNorm = El.MaxNorm(rLS)
if worldRank == 0:
  print('|| A x_{{LS}} - b ||_2  = {}'.format(rLSTwoNorm))
  print('|| A x_{{LS}} - b ||_oo = {}'.format(rLSInfNorm))

El.Finalize()
