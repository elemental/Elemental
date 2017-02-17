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
k = 3000
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

  A.ProcessQueues()
  return A

A = Rectang(m,n)
G = Rectang(k,n)

# Generate a (b,h) which implies a primal feasible (x,s)
# ======================================================
# b := A xGen
# -----------
xGen = El.DistMultiVec()
El.Gaussian(xGen,n,1)
b = El.DistMultiVec()
El.Zeros( b, m, 1 )
El.Multiply( El.NORMAL, 1., A, xGen, 0., b )
# h := G xGen + sGen
# ------------------
sGen = El.DistMultiVec()
El.Uniform(sGen,k,1,0.5,0.5)
h = El.DistMultiVec()
El.Copy( sGen, h )
El.Multiply( El.NORMAL, 1., G, xGen, 1., h )

# Generate a c which implies a dual feasible (y,z)
# ================================================
yGen = El.DistMultiVec()
El.Gaussian(yGen,m,1)
zGen = El.DistMultiVec()
El.Uniform(zGen,k,1,0.5,0.5)
c = El.DistMultiVec()
El.Zeros(c,n,1)
El.Multiply( El.TRANSPOSE, -1., A, yGen, 1., c )
El.Multiply( El.TRANSPOSE, -1., G, zGen, 1., c )

if display:
  El.Display( A, "A" )
  El.Display( G, "G" )
  El.Display( b, "b" )
  El.Display( c, "c" )
  El.Display( h, "h" )

# Set up the control structure (and possibly initial guesses)
# ===========================================================
ctrl = El.LPAffineCtrl_d()
xOrig = El.DistMultiVec()
yOrig = El.DistMultiVec()
zOrig = El.DistMultiVec()
sOrig = El.DistMultiVec()
if manualInit:
  El.Uniform(xOrig,n,1,0.5,0.4999)
  El.Uniform(yOrig,m,1,0.5,0.4999)
  El.Uniform(zOrig,k,1,0.5,0.4999)
  El.Uniform(sOrig,k,1,0.5,0.4999)
x = El.DistMultiVec()
y = El.DistMultiVec()
z = El.DistMultiVec()
s = El.DistMultiVec()

if testMehrotra:
  ctrl.approach = El.LP_MEHROTRA
  ctrl.mehrotraCtrl.primalInit = manualInit
  ctrl.mehrotraCtrl.dualInit = manualInit
  ctrl.mehrotraCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  El.Copy( sOrig, s )
  startMehrotra = El.mpi.Time()
  El.LPAffine(A,G,b,c,h,x,y,z,s,ctrl)
  endMehrotra = El.mpi.Time()
  if worldRank == 0:
    print('Mehrotra time: {} seconds'.format(endMehrotra-startMehrotra))

  if display:
    El.Display( x, "x Mehrotra" )
    El.Display( y, "y Mehrotra" )
    El.Display( z, "z Mehrotra" )
    El.Display( s, "s Mehrotra" )

  obj = El.Dot(c,x)
  if worldRank == 0:
    print('Mehrotra c^T x = {}'.format(obj))

if testIPF:
  ctrl.approach = El.LP_IPF
  ctrl.ipfCtrl.primalInit = manualInit
  ctrl.ipfCtrl.dualInit = manualInit
  ctrl.ipfCtrl.progress = progress
  ctrl.ipfCtrl.lineSearchCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  El.Copy( sOrig, s )
  startIPF = El.mpi.Time()
  El.LPAffine(A,G,b,c,h,x,y,z,s,ctrl)
  endIPF = El.mpi.Time()
  if worldRank == 0:
    print('IPF time: {} seconds'.format(endIPF-startIPF))

  if display:
    El.Display( x, "x IPF" )
    El.Display( y, "y IPF" )
    El.Display( z, "z IPF" )
    El.Display( s, "s IPF" )

  obj = El.Dot(c,x)
  if worldRank == 0:
    print('IPF c^T x = {}'.format(obj))

El.Finalize()
