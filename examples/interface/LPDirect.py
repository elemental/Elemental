#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

m = 2000
n = 4000
testMehrotra = True
testIPF = False
manualInit = False
display = False
progress = True
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# Make a sparse matrix with the last column dense
def Rectang(height,width):
  A = El.DistSparseMatrix()
  A.Resize(height,width)
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    if s < width: 
      A.QueueUpdate( s, s,        11, passive=True )
    if s >= 1 and s-1 < width:
      A.QueueUpdate( s, s-1,      -1, passive=True )
    if s+1 < width:
      A.QueueUpdate( s, s+1,       2, passive=True )
    if s >= height and s-height < width:
      A.QueueUpdate( s, s-height, -3, passive=True )
    if s+height < width: 
      A.QueueUpdate( s, s+height,  4, passive=True )
    # The dense last column
    A.QueueUpdate( s, width-1, -5/height, passive=True );

  A.ProcessQueues()
  return A

A = Rectang(m,n)

# Generate a b which implies a primal feasible x
# ==============================================
xGen = El.DistMultiVec()
El.Uniform(xGen,n,1,0.5,0.5)
b = El.DistMultiVec()
El.Zeros( b, m, 1 )
El.Multiply( El.NORMAL, 1., A, xGen, 0., b )

# Generate a c which implies a dual feasible (y,z)
# ================================================
yGen = El.DistMultiVec()
El.Gaussian(yGen,m,1)
c = El.DistMultiVec()
El.Uniform(c,n,1,0.5,0.5)
El.Multiply( El.TRANSPOSE, -1., A, yGen, 1., c )

if display:
  El.Display( A, "A" )
  El.Display( b, "b" )
  El.Display( c, "c" )

# Set up the control structure (and possibly initial guesses)
# ===========================================================
ctrl = El.LPDirectCtrl_d(isSparse=True)
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
  ctrl.approach = El.LP_MEHROTRA
  ctrl.mehrotraCtrl.primalInit = manualInit
  ctrl.mehrotraCtrl.dualInit = manualInit
  ctrl.mehrotraCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startMehrotra = El.mpi.Time()
  El.LPDirect(A,b,c,x,y,z,ctrl)
  endMehrotra = El.mpi.Time()
  if worldRank == 0:
    print "Mehrotra time:", endMehrotra-startMehrotra

  if display:
    El.Display( x, "x Mehrotra" )
    El.Display( y, "y Mehrotra" )
    El.Display( z, "z Mehrotra" )

  obj = El.Dot(c,x)
  if worldRank == 0:
    print "Mehrotra c^T x =", obj

if testIPF:
  ctrl.approach = El.LP_IPF
  ctrl.ipfCtrl.primalInit = manualInit
  ctrl.ipfCtrl.dualInit = manualInit
  ctrl.ipfCtrl.progress = progress
  ctrl.ipfCtrl.lineSearchCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startIPF = El.mpi.Time()
  El.LPDirect(A,b,c,x,y,z,ctrl)
  endIPF = El.mpi.Time()
  if worldRank == 0:
    print "IPF time:", endIPF-startIPF

  if display:
    El.Display( x, "x IPF" )
    El.Display( y, "y IPF" )
    El.Display( z, "z IPF" )

  obj = El.Dot(c,x)
  if worldRank == 0:
    print "IPF c^T x =", obj

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
