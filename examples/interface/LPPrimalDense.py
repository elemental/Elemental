#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El
import time

m = 1000
n = 2000
testMehrotra = True
testIPF = True
testADMM = False
display = False
worldRank = El.mpi.WorldRank()

# Make a sparse matrix with the last column dense
def RectangDense(m,n):
  A = El.DistMatrix()
  El.Zeros( A, m, n )
  for s in xrange(m):
    A.Update( s, s, 11 )
    if s != 0:   A.Update( s, s-1,   1 )
    if s != n-1: A.Update( s, s+1,   2 )
    if s >= m:   A.Update( s, s-m,   3 )
    if s <  n-m: A.Update( s, s+m,   4 )
    # The dense last column
    A.Update( s, n-1, 5./m );

  return A

A = RectangDense(m,n)

# Generate a right-hand side in the positive image
# ================================================
xGen = El.DistMatrix()
El.Uniform(xGen,n,1,0.5,0.4999)
b = El.DistMatrix()
El.Zeros( b, m, 1 )
El.Gemv( El.NORMAL, 1., A, xGen, 0., b )

# Generate a random positive cost function
# ========================================
c = El.DistMatrix()
El.Uniform(c,n,1,0.5,0.4999)

if display:
  El.Display( A,    "A"    )
  El.Display( xGen, "xGen" )
  El.Display( b,    "b"    )
  El.Display( c,    "c"    )

# Generate random initial guesses
# ===============================
xOrig = El.DistMatrix()
yOrig = El.DistMatrix()
zOrig = El.DistMatrix()
El.Uniform(xOrig,n,1,0.5,0.4999)
El.Uniform(yOrig,m,1,0.5,0.4999)
El.Uniform(zOrig,n,1,0.5,0.4999)
x = El.DistMatrix()
y = El.DistMatrix()
z = El.DistMatrix()

if testMehrotra:
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startMehrotra = time.clock()
  El.LPPrimalMehrotra(A,b,c,x,y,z)
  endMehrotra = time.clock()
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
  El.Copy( xOrig, x )
  El.Copy( yOrig, y )
  El.Copy( zOrig, z )
  startIPF = time.clock()
  El.LPPrimalIPF(A,b,c,x,y,z)
  endIPF = time.clock()
  if worldRank == 0:
    print "IPF time:", endIPF-startIPF

  if display:
    El.Display( x, "x IPF" )
    El.Display( y, "y IPF" )
    El.Display( z, "z IPF" )

  obj = El.Dot(c,x)
  if worldRank == 0:
    print "IPF c^T x =", obj

if testADMM:
  startADMM = time.clock()
  x = El.LPPrimalADMM(A,b,c)
  endADMM = time.clock()
  if worldRank == 0:
    print "ADMM time:", endADMM-startADMM

  if display:
    El.Display( x, "x ADMM" )

  obj = El.Dot(c,x)
  if worldRank == 0:
    print "ADMM c^T x =", obj

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
