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
testMehrotra = True
testIPF = False
manualInit = False
display = False
progress = True
worldRank = El.mpi.WorldRank()

# Make Q a sparse semidefinite matrix
def Semidefinite(height):
  Q = El.DistSparseMatrix()
  Q.Resize(height,height)
  firstLocalRow = Q.FirstLocalRow()
  localHeight = Q.LocalHeight()
  Q.Reserve(localHeight)
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
    Q.QueueLocalUpdate( sLoc, s, 1 );

  Q.MakeConsistent()
  return Q

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

Q = Semidefinite(n)
A = Rectang(m,n)

# Generate a b which implies a primal feasible x
# ==============================================
xGen = El.DistMultiVec()
El.Uniform(xGen,n,1,0.5,0.5)
b = El.DistMultiVec()
El.Zeros( b, m, 1 )
El.SparseMultiply( El.NORMAL, 1., A, xGen, 0., b )

# Generate a c which implies a dual feasible (y,z)
# ================================================
yGen = El.DistMultiVec()
El.Gaussian(yGen,m,1)
c = El.DistMultiVec()
El.Uniform(c,n,1,0.5,0.5)
El.SparseMultiply( El.NORMAL,    -1,  Q, xGen, 1., c )
El.SparseMultiply( El.TRANSPOSE, -1., A, yGen, 1., c )

if display:
  El.Display( Q, "Q" )
  El.Display( A, "A" )
  El.Display( b, "b" )
  El.Display( c, "c" )

# Set up the control structure (and possibly initial guesses)
# ===========================================================
ctrl = El.QPDirectCtrl_d()
xOrig = El.DistMultiVec()
yOrig = El.DistMultiVec()
zOrig = El.DistMultiVec()
if manualInit:
  El.Uniform(xOrig,n,1,0.5,0.4999)
  El.Uniform(yOrig,m,1,0.5,0.4999)
  El.Uniform(zOrig,n,1,0.5,0.4999)
x = El.DistMultiVec()
y = El.DistMultiVec()
z = El.DistMultiVec()

if testMehrotra:
  ctrl.approach = El.QP_MEHROTRA
  ctrl.mehrotraCtrl.primalInitialized = manualInit
  ctrl.mehrotraCtrl.dualInitialized = manualInit
  ctrl.mehrotraCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startMehrotra = time.clock()
  El.QPDirect(Q,A,b,c,x,y,z,ctrl)
  endMehrotra = time.clock()
  if worldRank == 0:
    print "Mehrotra time:", endMehrotra-startMehrotra

  if display:
    El.Display( x, "x Mehrotra" )
    El.Display( y, "y Mehrotra" )
    El.Display( z, "z Mehrotra" )

  d = El.DistMultiVec()
  El.Zeros( d, n, 1 )
  El.SparseMultiply( El.NORMAL, 1., Q, x, 0., d )
  obj = El.Dot(x,d)/2 + El.Dot(c,x)
  if worldRank == 0:
    print "Mehrotra (1/2) x^T Q x + c^T x =", obj

if testIPF:
  ctrl.approach = El.QP_IPF
  ctrl.ipfCtrl.primalInitialized = manualInit
  ctrl.ipfCtrl.dualInitialized = manualInit
  ctrl.ipfCtrl.progress = progress
  ctrl.ipfCtrl.lineSearchCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startIPF = time.clock()
  El.QPDirect(Q,A,b,c,x,y,z,ctrl)
  endIPF = time.clock()
  if worldRank == 0:
    print "IPF time:", endIPF-startIPF

  if display:
    El.Display( x, "x IPF" )
    El.Display( y, "y IPF" )
    El.Display( z, "z IPF" )

  d = El.DistMultiVec()
  El.Zeros( d, n, 1 )
  El.SparseMultiply( El.NORMAL, 1., Q, x, 0., d )
  obj = El.Dot(x,d)/2 + El.Dot(c,x)
  if worldRank == 0:
    print "IPF (1/2) x^T Q x + c^T x =", obj

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
