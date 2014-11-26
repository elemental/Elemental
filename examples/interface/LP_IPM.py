#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

m = 200
n = 400

def Rectang(m,n):
  A = El.DistSparseMatrix()
  A.Resize(m,n)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(5*localHeight)
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
    A.QueueLocalUpdate( sLoc, s, 11 )
    if s != 0:   A.QueueLocalUpdate( sLoc, s-1,   1 )
    if s != n-1: A.QueueLocalUpdate( sLoc, s+1,   2 )
    if s >= m:   A.QueueLocalUpdate( sLoc, s-m,   3 )
    if s <  n-m: A.QueueLocalUpdate( sLoc, s+m,   4 )

  A.MakeConsistent()
  return A

A = Rectang(m,n)
El.Display( A, "A" )

# Generate a right-hand side in the positive image
# ================================================
xGen = El.DistMultiVec()
El.Uniform(xGen,n,1,0.5,0.4999)
El.Display( xGen, "xGen" )
b = El.DistMultiVec()
El.Zeros( b, m, 1 )
El.SparseMultiply( El.NORMAL, 1., A, xGen, 0., b )
El.Display( b, "b" )

# Generate a random positive cost function
# ========================================
c = El.DistMultiVec()
El.Uniform(c,n,1,0.5,0.4999)

El.Print( A, "A" )
El.Print( b, "b" )
El.Print( c, "c" )

# Generate random initial guesses
# ===============================
x = El.DistMultiVec()
l = El.DistMultiVec()
s = El.DistMultiVec()
El.Uniform(x,n,1,0.5,0.4999)
El.Uniform(l,m,1,0.5,0.4999)
El.Uniform(s,n,1,0.5,0.4999)

El.LinearProgramIPF(A,b,c,s,x,l)
El.Display( x, "s" )
El.Display( l, "x" )
El.Display( s, "l" )

obj = El.Dot(c,x)
if El.mpi.WorldRank() == 0:
  print "c^T x =", obj

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
