#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El, time

n0 = n1 = 200
display = False
worldRank = El.mpi.WorldRank()

# A 2D finite-difference matrix with a dense last column
def FD2D(N0,N1):
  A = El.DistSparseMatrix()
  height = N0*N1
  width = N0*N1
  A.Resize(height,width)
  firstLocalRow = A.FirstLocalRow()
  localHeight = A.LocalHeight()
  A.Reserve(6*localHeight)
  for sLoc in xrange(localHeight):
    s = firstLocalRow + sLoc
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
    A.QueueLocalUpdate( sLoc, width-1, -10./height );

  A.MakeConsistent()
  return A

A = FD2D(n0,n1)
x = El.DistMultiVec()
y = El.DistMultiVec()
El.Uniform( x, n0*n1, 1 )
El.Copy( x, y )
if display:
  El.Display( A, "A" )
  El.Display( y, "y" )

yNrm = El.Nrm2(y)
rank = El.mpi.WorldRank()
if rank == 0:
  print "|| y ||_2 =", yNrm

ctrl = El.LeastSquaresCtrl_d()
ctrl.alpha = 0
ctrl.equilibrate = True
ctrl.progress = True
ctrl.qsdCtrl.regPrimal = 1e-6
ctrl.qsdCtrl.regDual = 1e-6
ctrl.qsdCtrl.alg = El.REG_REFINE_FGMRES
ctrl.qsdCtrl.relTol = 1e-9
ctrl.qsdCtrl.relTolRefine = 1e-12
ctrl.qsdCtrl.progress = True

solveStart = time.clock()
El.LinearSolve(A,x,ctrl)
solveStop = time.clock()
if worldRank == 0:
  print "LinearSolve time:", solveStop-solveStart, "seconds"
if display:
  El.Display( x, "x" )
xNrm = El.Nrm2(x)
if rank == 0:
  print "|| x ||_2 =", xNrm

El.SparseMultiply(El.NORMAL,-1.,A,x,1.,y)
if display:
  El.Display( y, "A x - y" )
eNrm = El.Nrm2(y)
if rank == 0:
  print "|| A x - y ||_2 / || y ||_2 =", eNrm/yNrm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
