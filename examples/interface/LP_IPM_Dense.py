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

m = 200
n = 400
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
    A.Update( s, n-1, 5 );

  return A

A = RectangDense(m,n)
El.Display( A, "A" )

# Generate a right-hand side in the positive image
# ================================================
xGen = El.DistMatrix()
El.Uniform(xGen,n,1,0.5,0.4999)
El.Display( xGen, "xGen" )
b = El.DistMatrix()
El.Zeros( b, m, 1 )
El.Gemv( El.NORMAL, 1., A, xGen, 0., b )
El.Display( b, "b" )

# Generate a random positive cost function
# ========================================
c = El.DistMatrix()
El.Uniform(c,n,1,0.5,0.4999)

# Generate random initial guesses
# ===============================
sOrig = El.DistMatrix()
xOrig = El.DistMatrix()
lOrig = El.DistMatrix()
El.Uniform(sOrig,n,1,0.5,0.4999)
El.Uniform(xOrig,n,1,0.5,0.4999)
El.Uniform(lOrig,m,1,0.5,0.4999)
s = El.DistMatrix()
x = El.DistMatrix()
l = El.DistMatrix()

El.Copy( sOrig, s )
El.Copy( xOrig, x )
El.Copy( lOrig, l )
startIPF = time.clock()
El.LinearProgramIPF(A,b,c,s,x,l)
endIPF = time.clock()
if worldRank == 0:
  print "IPF time:", endIPF-startIPF
El.Display( s, "s IPF" )
El.Display( x, "x IPF" )
El.Display( l, "l IPF" )

obj = El.Dot(c,x)
if worldRank == 0:
  print "IPF c^T x =", obj

El.Copy( sOrig, s )
El.Copy( xOrig, x )
El.Copy( lOrig, l )
startMehrotra = time.clock()
El.LinearProgramMehrotra(A,b,c,s,x,l)
endMehrotra = time.clock()
if worldRank == 0:
  print "Mehrotra time:", endMehrotra-startMehrotra
El.Display( s, "s Mehrotra" )
El.Display( x, "x Mehrotra" )
El.Display( l, "l Mehrotra" )

obj = El.Dot(c,x)
if worldRank == 0:
  print "Mehrotra c^T x =", obj

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
