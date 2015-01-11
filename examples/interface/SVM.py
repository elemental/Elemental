#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time

# TODO: Add noise so that the soft margin is actually put to use.
#       Perhaps a weighted coin should be flipped based upon the distance to the
#       hyperplane.

m = 4000
n = 2000
numLambdas = 4
startLambda = 1
endLambda = 10
display = True
worldRank = El.mpi.WorldRank()

def Rectang(height,width):
  A = El.DistSparseMatrix()
  A.Resize(height,width)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
    if s < width:
      A.QueueLocalUpdate( sLoc, s,        11 )
    if s >= 1 and s-1 < width:
      A.QueueLocalUpdate( sLoc, s-1,      -1 )
    if s+1 < width:
      A.QueueLocalUpdate( sLoc, s+1,       2 )
    if s >= height and s-height < width:
      A.QueueLocalUpdate( sLoc, s-height, -3 )
    if s+height < width:
      A.QueueLocalUpdate( sLoc, s+height,  4 )
    # The dense last column
    A.QueueLocalUpdate( sLoc, width-1, -5/height );

  A.MakeConsistent()
  return A

# Define a random (affine) hyperplane
wGen = El.DistMultiVec()
El.Gaussian( wGen, n, 1 )
wGenNorm = El.FrobeniusNorm( wGen )
El.Scale( 1./wGenNorm, wGen )
El.Print( wGen, "wGen" )
# TODO: Add support for mpi::Broadcast and randomly generate this
offset = 0.3147

# Define a random set of points
A = Rectang(m,n)

# Label the points based upon their location relative to the hyperplane
d = El.DistMultiVec()
El.Ones( d, m, 1 )
El.SparseMultiply( El.NORMAL, 1., A, wGen, -offset, d )
El.EntrywiseMap( d, lambda alpha : 1. if alpha > 0 else -1. )

El.Print( A, "A" )
El.Print( d, "d" )

if display:
  El.Display( wGen, "wGen" )
  if worldRank == 0:
    print "offset =", offset
  El.Display( A, "A" )
  El.Display( d, "d" )

ctrl = El.QPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True

for j in xrange(0,numLambdas):
  lambd = startLambda + j*(endLambda-startLambda)/(numLambdas-1.)
  if worldRank == 0:
    print "lambda =", lambd

  # TODO: Explicitly return w, beta, and z
  startSVM = time.clock()
  x = El.SVM( A, d, lambd, ctrl )
  endSVM = time.clock()
  if worldRank == 0:
    print "SVM time: ", endSVM-startSVM

  if display:
    El.Display( x, "[w;beta;z]" )

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
