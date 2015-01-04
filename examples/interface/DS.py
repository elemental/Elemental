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

m = 2000
n = 4000
numLambdas = 5
startLambda = 0.01
endLambda = 1
display = True
worldRank = El.mpi.WorldRank()

# Make a sparse matrix with the last column dense
def Rectang(height,width):
  A = El.DistSparseMatrix()
  A.Resize(height,width)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
    A.QueueLocalUpdate( sLoc, s, 11 )
    if s != 0:            A.QueueLocalUpdate( sLoc, s-1,      -1 )
    if s != width-1:      A.QueueLocalUpdate( sLoc, s+1,       2 )
    if s >= height:       A.QueueLocalUpdate( sLoc, s-height, -3 )
    if s <  width-height: A.QueueLocalUpdate( sLoc, s+height,  4 )
    # The dense last column
    A.QueueLocalUpdate( sLoc, width-1, -5/height );

  A.MakeConsistent()
  return A

A = Rectang(m,n)
b = El.DistMultiVec()
El.Gaussian( b, m, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

ctrl = El.LPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True

for j in xrange(0,numLambdas):
  lambd = startLambda + j*(endLambda-startLambda)/(numLambdas-1.)
  if worldRank == 0:
    print "lambda =", lambd

  startDS = time.clock()
  x = El.DS( A, b, lambd, ctrl )
  endDS = time.clock()
  if worldRank == 0:
    print "DS time: ", endDS-startDS

  if display:
    El.Display( x, "x" )

  xOneNorm = El.EntrywiseNorm( x, 1 )
  r = El.DistMultiVec()
  El.Copy( b, r )
  El.SparseMultiply( El.NORMAL, -1., A, x, 1., r )
  rTwoNorm = El.Nrm2( r )
  t = El.DistMultiVec()
  El.Zeros( t, n, 1 )
  El.SparseMultiply( El.TRANSPOSE, 1., A, r, 0., t )
  tTwoNorm = El.Nrm2( t )
  tInfNorm = El.MaxNorm( t )
  if display:
    El.Display( r, "r" )
    El.Display( t, "t" )
  if worldRank == 0:
    print "|| x ||_1       =", xOneNorm
    print "|| b - A x ||_2 =", rTwoNorm
    print "|| A^T (b - A x) ||_2 =", tTwoNorm
    print "|| A^T (b - A x) ||_oo =", tInfNorm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
