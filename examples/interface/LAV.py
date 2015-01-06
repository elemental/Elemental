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

m = 4000
n = 2000
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

A = Rectang(m,n)
b = El.DistMultiVec()
El.Gaussian( b, m, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

ctrl = El.LPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True
startLAV = time.clock()
x = El.LAV( A, b, ctrl )
endLAV = time.clock()
if worldRank == 0:
  print "LAV time: ", endLAV-startLAV

if display:
  El.Display( x, "x" )

bTwoNorm = El.Nrm2( b )
bInfNorm = El.MaxNorm( b )
r = El.DistMultiVec()
El.Copy( b, r )
El.SparseMultiply( El.NORMAL, -1., A, x, 1., r )
if display:
  El.Display( r, "r" )
rTwoNorm = El.Nrm2( r )
rOneNorm = El.EntrywiseNorm( r, 1 )
if worldRank == 0:
  print "|| b ||_2       =", bTwoNorm
  print "|| b ||_oo      =", bInfNorm
  print "|| A x - b ||_2 =", rTwoNorm
  print "|| A x - b ||_1 =", rOneNorm

xLS = El.LeastSquares(A,b)
if display:
  El.Display( xLS, "x_{LS}" )
rLS = El.DistMultiVec()
El.Copy( b, rLS )
El.SparseMultiply( El.NORMAL, -1., A, xLS, 1., rLS )
El.Display( rLS, "A x_{LS} - b" )
rLSTwoNorm = El.Nrm2(rLS)
rLSOneNorm = El.EntrywiseNorm(rLS,1)
if worldRank == 0:
  print "|| A x_{LS} - b ||_2 =", rLSTwoNorm
  print "|| A x_{LS} - b ||_1 =", rLSOneNorm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
