#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, random

n0 = n1 = 100
numColsB = 3
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

def Constraints(numCols,N0,N1):
  B = El.DistSparseMatrix()
  El.Zeros( B, N0*N1, numCols )
  localHeight = B.LocalHeight()
  B.Reserve( localHeight*numCols )
  for sLoc in xrange(localHeight):
    s = B.GlobalRow(sLoc)
    for j in xrange(numCols):
      B.QueueLocalUpdate( sLoc, j, random.uniform(0,1) )

  B.ProcessQueues()
  return B

A = FD2D(n0,n1)
B = Constraints(numColsB,n0,n1)
if display:
  El.Display( A, "A" )
  El.Display( B, "B" )
if output:
  El.Print( A, "A" )
  El.Print( B, "B" )

D = El.DistMultiVec()
El.Uniform( D, A.Height(), numRHS )
if display:
  El.Display( D, "D" )
if output:
  El.Print( D, "D" )
DNorm = El.FrobeniusNorm( D )

baseAlpha = 1e-4
ctrl = El.LeastSquaresCtrl_d()
ctrl.alpha = baseAlpha
ctrl.progress = True
ctrl.equilibrate = True
ctrl.solveCtrl.relTol = 1e-10
ctrl.solveCtrl.relTolRefine = 1e-12
ctrl.solveCtrl.progress = True
startGLM = El.mpi.Time()
X,Y = El.GLM(A,B,D,ctrl)
endGLM = El.mpi.Time()
if worldRank == 0:
  print "GLM time:", endGLM-startGLM, "seconds"
if display:
  El.Display( X, "X" )
  El.Display( Y, "Y" )
if output:
  El.Print( X, "X" )
  El.Print( Y, "Y" )

YNorm = El.FrobeniusNorm( Y )
if worldRank == 0:
  print "|| Y ||_F =", YNorm

E = El.DistMultiVec()
El.Copy( D, E )
El.Multiply( El.NORMAL, -1., A, X, 1., E )
El.Multiply( El.NORMAL, -1., B, Y, 1., E )
residNorm = El.FrobeniusNorm( E )
if display:
  El.Display( E, "D - A X - B Y" )
if output:
  El.Print( E, "D - A X - B Y" )
if worldRank == 0:
  print "|| D - A X - B Y ||_F / || D ||_F =", residNorm/DNorm

# Now try solving a weighted least squares problem
# (as lambda -> infinity, the exact solution converges to that of LSE)
def SolveWeighted(A,B,D,lambd):
  AScale = El.DistSparseMatrix()
  El.Copy( A, AScale )
  El.Scale( lambd, AScale )
  AEmb = El.HCat(AScale,B)
  if display:
    El.Display( AEmb, "[lambda*A, B]" )
  if output:
    El.Print( AEmb, "[lambda*A, B]" )

  ctrl.alpha = baseAlpha
  if worldRank == 0:
    print "lambda=", lambd, ": ctrl.alpha=", ctrl.alpha
  XEmb=El.LeastSquares(AEmb,D,ctrl)

  X = XEmb[0:n0*n1,0:numRHS]
  Y = XEmb[n0*n1:n0*n1+numColsB,0:numRHS]
  El.Scale( lambd, X )

  YNorm = El.FrobeniusNorm( Y )
  if worldRank == 0:
    print "lambda=", lambd, ": || Y ||_F =", YNorm

  El.Copy( D, E )
  El.Multiply( El.NORMAL, -1., A, X, 1., E )
  El.Multiply( El.NORMAL, -1., B, Y, 1., E )
  residNorm = El.FrobeniusNorm( E )
  if worldRank == 0:
    print "lambda=", lambd, ": || D - A X - B Y ||_F / || D ||_F =", residNorm/DNorm

SolveWeighted(A,B,D,1)
SolveWeighted(A,B,D,10)
SolveWeighted(A,B,D,100)
SolveWeighted(A,B,D,1000)
SolveWeighted(A,B,D,10000)
SolveWeighted(A,B,D,100000)

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
