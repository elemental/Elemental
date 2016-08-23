#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

# TODO: Add noise so that the soft margin is actually put to use.
#       Perhaps a weighted coin should be flipped based upon the distance to the
#       hyperplane.

m = 400
n = 200
numLambdas = 4
startLambda = 1
endLambda = 10
display = True
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

def Rectang(height,width):
  A = El.DistMatrix()
  El.Uniform( A, height, width )
  return A

def RectangSparse(height,width):
  A = El.DistMatrix()
  El.Zeros( A, height, width )
  for s in xrange(height):
    if s < width:                        A.Update( s, s,        11 )
    if s >= 1 and s-1 < width:           A.Update( s, s-1,      -1 )
    if s+1 < width:                      A.Update( s, s+1,       2 )
    if s >= height and s-height < width: A.Update( s, s-height, -3 )
    if s+height < width:                 A.Update( s, s+height,  4 )
    # The dense last column
    A.Update( s, width-1, -5/height );    

  return A

# Define a random (affine) hyperplane
wGen = El.DistMatrix()
El.Gaussian( wGen, n, 1 )
wGenNorm = El.FrobeniusNorm( wGen )
El.Scale( 1./wGenNorm, wGen )
El.Print( wGen, "wGen" )
# TODO: Add support for mpi::Broadcast and randomly generate this
offset = 0.3147

# Define a random set of points
A = RectangSparse(m,n)

# Label the points based upon their location relative to the hyperplane
d = El.DistMatrix()
El.Ones( d, m, 1 )
El.Gemv( El.NORMAL, 1., A, wGen, -offset, d )
El.EntrywiseMap( d, lambda alpha : 1. if alpha > 0 else -1. )

El.Print( A, "A" )
El.Print( d, "d" )

if display:
  El.Display( wGen, "wGen" )
  if worldRank == 0:
    print "offset =", offset
  El.Display( A, "A" )
  El.Display( d, "d" )

ctrl = El.SVMCtrl_d()
ctrl.ipmCtrl.mehrotraCtrl.progress = True

for j in xrange(0,numLambdas):
  lambd = startLambda + j*(endLambda-startLambda)/(numLambdas-1.)
  if worldRank == 0:
    print "lambda =", lambd

  startSVM = El.mpi.Time()
  x = El.SVM( A, d, lambd, ctrl )
  endSVM = El.mpi.Time()
  if worldRank == 0:
    print "SVM time: ", endSVM-startSVM

  w = x[0:n,0]
  beta = x.Get(n,0)
  z = x[n+1:n+m+1,0]
  if display:
    El.Display( w, "w" )
    if worldRank==0:
      print "beta=", beta
    El.Display( z, "z" )

  e = El.DistMatrix()
  El.Zeros( e, m, 1 )
  El.Gemv( El.NORMAL, 1., A, w, 1., e )
  El.Shift( e, beta )
  El.DiagonalScale( El.LEFT, El.NORMAL, d, e )
  if display:
    El.Display( e, "diag(d)*(A w + beta)" )
  eTwoNorm = El.Nrm2( e )
  if worldRank == 0:
    print "|| A w + beta - d ||_2 =", eTwoNorm

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
