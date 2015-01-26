#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time

n = 4000
numLambdas = 5
startLambda = 1.
endLambda = 50.
display = True
worldRank = El.mpi.WorldRank()

def Deriv(height):
  A = El.DistSparseMatrix()
  A.Resize(height-1,height)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(2*localHeight)
  for iLoc in xrange(localHeight):
    i = firstLocalRow + iLoc
    A.QueueLocalUpdate( iLoc, i, 1. )
    A.QueueLocalUpdate( iLoc, i+1, -1. )

  A.MakeConsistent()
  return A

D = Deriv(n)

b = El.DistMultiVec()
El.Gaussian( b, n, 1 )
if display:
  El.Display( b, "b" )

ctrl = El.QPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True

for j in xrange(0,numLambdas):
  lambd = startLambda + j*(endLambda-startLambda)/(numLambdas-1.)
  if worldRank == 0:
    print "lambda =", lambd

  startTV = time.clock()
  x = El.TV( b, lambd, ctrl )
  endTV = time.clock()
  if worldRank == 0:
    print "TV time: ", endTV-startTV

  Dx = El.DistMultiVec() 
  El.Zeros( Dx, n-1, 1 )
  El.SparseMultiply( El.NORMAL, 1., D, x, 0., Dx )
  if display:
    El.Display( x, "x" )
    El.Display( Dx, "Dx" )

  DxOneNorm = El.EntrywiseNorm( Dx, 1 )
  e = El.DistMultiVec()
  El.Copy( b, e )
  El.Axpy( -1., x, e )
  if display:
    El.Display( e, "e" )
  eTwoNorm = El.Nrm2( e )
  if worldRank == 0:
    print "|| D x ||_1   =", DxOneNorm
    print "|| x - b ||_2 =", eTwoNorm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
