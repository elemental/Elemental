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

m = 2000
n = 4000
testIPF = True
testMehrotra = True
display = False
worldRank = El.mpi.WorldRank()

# Make Q a sparse semidefinite matrix
def Semidefinite(n):
  Q = El.DistSparseMatrix()
  Q.Resize(n,n)
  firstLocalRow = Q.FirstLocalRow()
  localHeight = Q.LocalHeight()
  Q.Reserve(localHeight)
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
    Q.QueueLocalUpdate( sLoc, s, 1 );
  
  Q.MakeConsistent()
  return Q

# Make a sparse matrix with the last column dense
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
    # The dense last column
    A.QueueLocalUpdate( sLoc, n-1, 5./m );

  A.MakeConsistent()
  return A

Q = Semidefinite(n)
A = Rectang(m,n)

# Generate a right-hand side in the positive image
# ================================================
xGen = El.DistMultiVec()
El.Uniform(xGen,n,1,0.5,0.4999)
b = El.DistMultiVec()
El.Zeros( b, m, 1 )
El.SparseMultiply( El.NORMAL, 1., A, xGen, 0., b )

# Generate a random positive cost function
# ========================================
c = El.DistMultiVec()
El.Uniform(c,n,1,0.5,0.4999)

if display:
  El.Display( Q,    "Q"    )
  El.Display( A,    "A"    )
  El.Display( xGen, "xGen" )
  El.Display( b,    "b"    )
  El.Display( c,    "c"    )

# Generate random initial guesses
# ===============================
xOrig = El.DistMultiVec()
lOrig = El.DistMultiVec()
sOrig = El.DistMultiVec()
El.Uniform(xOrig,n,1,0.5,0.4999)
El.Uniform(lOrig,m,1,0.5,0.4999)
El.Uniform(sOrig,n,1,0.5,0.4999)
x = El.DistMultiVec()
l = El.DistMultiVec()
s = El.DistMultiVec()

if testMehrotra:
  El.Copy( sOrig, s )
  El.Copy( xOrig, x )
  El.Copy( lOrig, l )
  startMehrotra = time.clock()
  El.QuadraticProgramMehrotra(Q,A,b,c,s,x,l)
  endMehrotra = time.clock()
  if worldRank == 0:
    print "Mehrotra time:", endMehrotra-startMehrotra

  if display:
    El.Display( x, "s Mehotra" )
    El.Display( l, "x Mehotra" )
    El.Display( s, "l Mehotra" )

  Q_x = El.DistMultiVec()
  El.Zeros( Q_x, n, 1 )
  El.SparseMultiply( El.NORMAL, 1., Q, x, 0., Q_x )
  xTQx = El.Dot(x,Q_x)
  obj = El.Dot(c,x) + xTQx/2
  if worldRank == 0:
    print "Mehrotra primal objective =", obj

if testIPF:
  El.Copy( sOrig, s )
  El.Copy( xOrig, x )
  El.Copy( lOrig, l )
  startIPF = time.clock()
  El.QuadraticProgramIPF(Q,A,b,c,s,x,l)
  endIPF = time.clock()
  if worldRank == 0:
    print "IPF time:", endIPF-startIPF

  if display:
    El.Display( x, "s IPF" )
    El.Display( l, "x IPF" )
    El.Display( s, "l IPF" )

  Q_x = El.DistMultiVec() 
  El.Zeros( Q_x, n, 1 )
  El.SparseMultiply( El.NORMAL, 1., Q, x, 0., Q_x )
  xTQx = El.Dot(x,Q_x)
  obj = El.Dot(c,x) + xTQx/2
  if worldRank == 0:
    print "IPF primal objective =", obj

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
