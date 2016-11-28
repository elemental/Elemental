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
numRowsB = 5
numRHS = 1
display = False

# NOTE: Increasing the magnitudes of the off-diagonal entries by an order of
#       magnitude makes the condition number vastly higher.
def FD2D(N0,N1):
  A = El.SparseMatrix()
  height = N0*N1
  width = N0*N1
  A.Resize(height,width)
  A.Reserve(6*height)
  for s in xrange(height):
    x0 = s % N0
    x1 = s / N0
    A.QueueUpdate( s, s, 11 )
    if x0 > 0:
      A.QueueUpdate( s, s-1, -1 )
    if x0+1 < N0:
      A.QueueUpdate( s, s+1, 2 )
    if x1 > 0:
      A.QueueUpdate( s, s-N0, -3 )
    if x1+1 < N1:
      A.QueueUpdate( s, s+N0, 4 )

    # The dense last column
    A.QueueUpdate( s, width-1, -10/height );

  A.ProcessQueues()
  return A

def Constraints(numRows,N0,N1):
  B = El.SparseMatrix()
  El.Zeros( B, numRows, N0*N1 )
  B.Reserve( numRows*N0*N1 )
  for s in xrange(numRows):
    for j in xrange(N0*N1):
      B.QueueUpdate( s, j, random.uniform(0,1) )

  B.ProcessQueues()
  return B

A = FD2D(n0,n1)
B = Constraints(numRowsB,n0,n1)
if display:
  El.Display( A, "A" )
  El.Display( B, "B" )

C = El.Matrix()
D = El.Matrix()
El.Uniform( C, A.Height(), numRHS )
El.Uniform( D, B.Height(), numRHS )
if display:
  El.Display( C, "C" )
  El.Display( D, "D" )
CNorm = El.FrobeniusNorm( C )
DNorm = El.FrobeniusNorm( D )

ctrl = El.LeastSquaresCtrl_d()
ctrl.progress = True
ctrl.sqsdCtrl.solveCtrl.relTol = 1e-10
ctrl.sqsdCtrl.solveCtrl.relTolRefine = 1e-12
ctrl.sqsdCtrl.solveCtrl.progress = True

startLSE = El.mpi.Time()
X = El.LSE(A,B,C,D,ctrl)
endLSE = El.mpi.Time()
print('LSE time: {} seconds'.format(endLSE-startLSE))
if display:
  El.Display( X, "X" )

E = El.Matrix()

El.Copy( C, E )
El.Multiply( El.NORMAL, -1., A, X, 1., E )
residNorm = El.FrobeniusNorm( E )
if display:
  El.Display( E, "C - A X" )
print('|| C - A X ||_F / || C ||_F = {}'.format(residNorm/CNorm))

El.Copy( D, E )
El.Multiply( El.NORMAL, -1., B, X, 1., E )
equalNorm = El.FrobeniusNorm( E )
if display:
  El.Display( E, "D - B X" )
print('|| D - B X ||_F / || D ||_F = '.format(equalNorm/DNorm))

# Now try solving a weighted least squares problem
# (as lambda -> infinity, the exact solution converges to that of LSE)
def SolveWeighted(A,B,C,D,lambd):
  BScale = El.SparseMatrix()
  El.Copy( B, BScale )
  El.Scale( lambd, BScale )

  DScale = El.Matrix()
  El.Copy( D, DScale ) 
  El.Scale( lambd, DScale )

  AEmb = El.VCat(A,BScale)
  CEmb = El.VCat(C,DScale)

  X=El.LeastSquares(AEmb,CEmb,ctrl)

  El.Copy( C, E )
  El.Multiply( El.NORMAL, -1., A, X, 1., E )
  residNorm = El.FrobeniusNorm( E )
  if display:
    El.Display( E, "C - A X" )
  print('lambda={}: || C - A X ||_F / || C ||_F = {}'.format(lambd,
    residNorm/CNorm))

  El.Copy( D, E )
  El.Multiply( El.NORMAL, -1., B, X, 1., E )
  equalNorm = El.FrobeniusNorm( E )
  if display:
    El.Display( E, "D - B X" )
  print('lambda={}: || D - B X ||_F / || D ||_F = {}'.format(lambd,
    equalNorm/DNorm))

SolveWeighted(A,B,C,D,1)
SolveWeighted(A,B,C,D,10)
SolveWeighted(A,B,C,D,100)
SolveWeighted(A,B,C,D,1000)
SolveWeighted(A,B,C,D,10000)
SolveWeighted(A,B,C,D,100000)
SolveWeighted(A,B,C,D,1000000)

El.Finalize()
