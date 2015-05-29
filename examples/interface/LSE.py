#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time, random

n0 = n1 = 50
numRowsB = 5
numRHS = 1
display = False
output = False
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

# NOTE: Increasing the magnitudes of the off-diagonal entries by an order of
#       magnitude makes the condition number vastly higher.
def FD2D(N0,N1):
  A = El.DistSparseMatrix()
  height = N0*N1
  width = N0*N1
  A.Resize(height,width)
  localHeight = A.LocalHeight()
  A.Reserve(6*localHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    x0 = s % N0
    x1 = s / N0
    A.QueueLocalUpdate( sLoc, s, 11 )
    if x0 > 0:
      A.QueueLocalUpdate( sLoc, s-1, -1 )
    if x0+1 < N0:
      A.QueueLocalUpdate( sLoc, s+1, 2 )
    if x1 > 0:
      A.QueueLocalUpdate( sLoc, s-N0, -3 )
    if x1+1 < N1:
      A.QueueLocalUpdate( sLoc, s+N0, 4 )

    # The dense last column
    A.QueueLocalUpdate( sLoc, width-1, -10/height );

  A.ProcessQueues()
  return A

def Constraints(numRows,N0,N1):
  B = El.DistSparseMatrix()
  El.Zeros( B, numRows, N0*N1 )
  localHeight = B.LocalHeight()
  B.Reserve( localHeight*N0*N1 )
  for sLoc in xrange(localHeight):
    s = B.GlobalRow(sLoc)
    for j in xrange(N0*N1):
      B.QueueLocalUpdate( sLoc, j, random.uniform(0,1) )

  B.ProcessQueues()
  return B

A = FD2D(n0,n1)
B = Constraints(numRowsB,n0,n1)
if display:
  El.Display( A, "A" )
  El.Display( B, "B" )
if output:
  El.Print( A, "A" )
  El.Print( B, "B" )

C = El.DistMultiVec()
D = El.DistMultiVec()
El.Uniform( C, A.Height(), numRHS )
El.Uniform( D, B.Height(), numRHS )
if display:
  El.Display( C, "C" )
  El.Display( D, "D" )
if output:
  El.Print( C, "C" )
  El.Print( D, "D" )
CNorm = El.FrobeniusNorm( C )
DNorm = El.FrobeniusNorm( D )

baseAlpha = 1e-4
ctrl = El.LeastSquaresCtrl_d()
ctrl.alpha = baseAlpha
ctrl.progress = True
ctrl.equilibrate = True
ctrl.qsdCtrl.relTol = 1e-10
ctrl.qsdCtrl.relTolRefine = 1e-12
ctrl.qsdCtrl.progress = True
startLSE = time.clock()
X = El.LSE(A,B,C,D,ctrl)
endLSE = time.clock()
if worldRank == 0:
  print "LSE time:", endLSE-startLSE, "seconds"
if display:
  El.Display( X, "X" )
if output:
  El.Print( X, "X" )

E = El.DistMultiVec()

El.Copy( C, E )
El.SparseMultiply( El.NORMAL, -1., A, X, 1., E )
residNorm = El.FrobeniusNorm( E )
if display:
  El.Display( E, "C - A X" )
if output:
  El.Print( E, "C - A X" )
if worldRank == 0:
  print "|| C - A X ||_F / || C ||_F =", residNorm/CNorm

El.Copy( D, E )
El.SparseMultiply( El.NORMAL, -1., B, X, 1., E )
equalNorm = El.FrobeniusNorm( E )
if display:
  El.Display( E, "D - B X" )
if output:
  El.Print( E, "D - B X" )
if worldRank == 0:
  print "|| D - B X ||_F / || D ||_F =", equalNorm/DNorm

# Now try solving a weighted least squares problem
# (as lambda -> infinity, the exact solution converges to that of LSE)
def SolveWeighted(A,B,C,D,lambd):
  BScale = El.DistSparseMatrix()
  El.Copy( B, BScale )
  El.Scale( lambd, BScale )

  DScale = El.DistMultiVec()
  El.Copy( D, DScale ) 
  El.Scale( lambd, DScale )

  AEmb = El.VCat(A,BScale)
  CEmb = El.VCat(C,DScale)
  if output:
    El.Print( AEmb, "AEmb" )

  ctrl.alpha = baseAlpha
  if worldRank == 0:
    print "lambda=", lambd, ": ctrl.alpha=", ctrl.alpha
  X=El.LeastSquares(AEmb,CEmb,ctrl)

  El.Copy( C, E )
  El.SparseMultiply( El.NORMAL, -1., A, X, 1., E )
  residNorm = El.FrobeniusNorm( E )
  if display:
    El.Display( E, "C - A X" )
  if output:
    El.Print( E, "C - A X" )
  if worldRank == 0:
    print "lambda=", lambd, ": || C - A X ||_F / || C ||_F =", residNorm/CNorm

  El.Copy( D, E )
  El.SparseMultiply( El.NORMAL, -1., B, X, 1., E )
  equalNorm = El.FrobeniusNorm( E )
  if display:
    El.Display( E, "D - B X" )
  if output:
    El.Print( E, "D - B X" )
  if worldRank == 0:
    print "lambda=", lambd, ": || D - B X ||_F / || D ||_F =", equalNorm/DNorm

SolveWeighted(A,B,C,D,1)
SolveWeighted(A,B,C,D,10)
SolveWeighted(A,B,C,D,100)
SolveWeighted(A,B,C,D,1000)
SolveWeighted(A,B,C,D,10000)
SolveWeighted(A,B,C,D,100000)

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
