#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

m = 1000
n = 2000
testMehrotra = True
testIPF = True
testADMM = False
manualInit = False
display = False
progress = True
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# Make a semidefinite matrix
def Semidefinite(height):
  Q = El.DistMatrix()
  El.Identity( Q, height, height )
  return Q

# Make a dense matrix
def RectangDense(height,width):
  A = El.DistMatrix()
  El.Gaussian( A, height, width )
  return A

Q = Semidefinite(n)
A = RectangDense(m,n)

# Generate a b which implies a primal feasible x
# ==============================================
xGen = El.DistMatrix()
El.Uniform(xGen,n,1,0.5,0.4999)
b = El.DistMatrix()
El.Zeros( b, m, 1 )
El.Gemv( El.NORMAL, 1., A, xGen, 0., b )

# Generate a c which implies a dual feasible (y,z)
# ================================================
yGen = El.DistMatrix()
El.Gaussian(yGen,m,1)
c = El.DistMatrix()
El.Uniform(c,n,1,0.5,0.5)
El.Hemv( El.LOWER,     -1,  Q, xGen, 1., c )
El.Gemv( El.TRANSPOSE, -1., A, yGen, 1., c )

if display:
  El.Display( Q, "Q" )
  El.Display( A, "A" )
  El.Display( b, "b" )
  El.Display( c, "c" )

# Set up the control structure (and possibly initial guesses)
# ===========================================================
ctrl = El.QPDirectCtrl_d()
xOrig = El.DistMatrix()
yOrig = El.DistMatrix()
zOrig = El.DistMatrix()
if manualInit:
  El.Uniform(xOrig,n,1,0.5,0.4999)
  El.Uniform(yOrig,m,1,0.5,0.4999)
  El.Uniform(zOrig,n,1,0.5,0.4999)
x = El.DistMatrix()
y = El.DistMatrix()
z = El.DistMatrix()

if testMehrotra:
  ctrl.approach = El.QP_MEHROTRA
  ctrl.mehrotraCtrl.primalInit = manualInit
  ctrl.mehrotraCtrl.dualInit = manualInit
  ctrl.mehrotraCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startMehrotra = El.mpi.Time()
  El.QPDirect(Q,A,b,c,x,y,z,ctrl)
  endMehrotra = El.mpi.Time()
  if worldRank == 0:
    print('Mehrotra time: {} seconds'.format(endMehrotra-startMehrotra))

  if display:
    El.Display( x, "x Mehrotra" )
    El.Display( y, "y Mehrotra" )
    El.Display( z, "z Mehrotra" )

  d = El.DistMatrix()
  El.Zeros( d, n, 1 )
  El.Hemv( El.LOWER, 1., Q, x, 0., d )
  obj = El.Dot(x,d)/2 + El.Dot(c,x)
  if worldRank == 0:
    print('Mehrotra (1/2) x^T Q x + c^T x = {}'.format(obj))

if testIPF:
  ctrl.approach = El.QP_IPF
  ctrl.ipfCtrl.primalInit = manualInit
  ctrl.ipfCtrl.dualInit = manualInit
  ctrl.ipfCtrl.progress = progress
  ctrl.ipfCtrl.lineSearchCtrl.progress = progress
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startIPF = El.mpi.Time()
  El.QPDirect(Q,A,b,c,x,y,z,ctrl)
  endIPF = El.mpi.Time()
  if worldRank == 0:
    print('IPF time: {} seconds'.format(endIPF-startIPF))

  if display:
    El.Display( x, "x IPF" )
    El.Display( y, "y IPF" )
    El.Display( z, "z IPF" )

  d = El.DistMatrix()
  El.Zeros( d, n, 1 )
  El.Hemv( El.LOWER, 1., Q, x, 0., d )
  obj = El.Dot(x,d)/2 + El.Dot(c,x)
  if worldRank == 0:
    print('IPF c^T x = {}'.format(obj))

El.Finalize()
