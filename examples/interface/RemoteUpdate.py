#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

n0 = 50
n1 = 50
display = True
worldSize = El.mpi.WorldSize()
worldRank = El.mpi.WorldRank()

# Stack two 2D finite-difference matrices on top of each other
# and make the last column dense
def StackedFD2D(N0,N1):
  A = El.DistMatrix()
  height = 2*N0*N1
  width = N0*N1
  A.Resize(height,width)
  blocksize = height // worldSize
  myStart = blocksize*worldRank
  if worldRank == worldSize-1:
    myHeight = height - myStart
  else:
    myHeight = blocksize
    
  A.Reserve(6*myHeight)
  for sLoc in xrange(localHeight):
    s = A.GlobalRow(sLoc)
    if s < N0*N1:
      x0 = s % N0
      x1 = s / N0
      A.QueueUpdate( sLoc, s, 11 )
      if x0 > 0:
        A.QueueUpdate( sLoc, s-1, -1 )
      if x0+1 < N0:
        A.QueueUpdate( sLoc, s+1, 2 )
      if x1 > 0:
        A.QueueUpdate( sLoc, s-N0, -30 )
      if x1+1 < N1:
        A.QueueUpdate( sLoc, s+N0, 4 )
    else:
      sRel = s-N0*N1
      x0 = sRel % N0
      x1 = sRel / N0
      A.QueueUpdate( sLoc, sRel, -20 )
      if x0 > 0:
        A.QueueUpdate( sLoc, sRel-1, -17 )
      if x0+1 < N0:
        A.QueueUpdate( sLoc, sRel+1, -20 )
      if x1 > 0:
        A.QueueUpdate( sLoc, sRel-N0, -3 )
      if x1+1 < N1:
        A.QueueUpdate( sLoc, sRel+N0, 3 )

    # The dense last column
    A.QueueUpdate( sLoc, width-1, -10/height );

  A.ProcessQueues()
  return A

A = StackedFD2D(n0,n1)
b = El.DistMatrix()
El.Gaussian( b, 2*n0*n1, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

ctrl = El.LPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True
ctrl.mehrotraCtrl.solveCtrl.relTol = 1e-10
ctrl.mehrotraCtrl.solveCtrl.relTolRefine = 1e-11
ctrl.mehrotraCtrl.solveCtrl.progress = True
startCP = El.mpi.Time()
x = El.CP( A, b, ctrl )
endCP = El.mpi.Time()
if worldRank == 0:
  print "CP time:", endCP-startCP, "seconds"
if display:
  El.Display( x, "x" )

bTwoNorm = El.Nrm2( b )
bInfNorm = El.MaxNorm( b )
r = El.DistMatrix()
El.Copy( b, r )
El.Gemv( El.NORMAL, -1., A, x, 1., r )
if display:
  El.Display( r, "r" )
rTwoNorm = El.Nrm2( r )
rInfNorm = El.MaxNorm( r )
if worldRank == 0:
  print "|| b ||_2        =", bTwoNorm
  print "|| b ||_oo       =", bInfNorm
  print "|| A x - b ||_2  =", rTwoNorm
  print "|| A x - b ||_oo =", rInfNorm

startLS = El.mpi.Time()
xLS = El.LeastSquares(A,b)
endLS = El.mpi.Time()
if worldRank == 0:
  print "LS time:", endLS-startLS, "seconds"
if display:
  El.Display( xLS, "x_{LS}" )
rLS = El.DistMatrix()
El.Copy( b, rLS )
El.Gemv( El.NORMAL, -1., A, xLS, 1., rLS )
if display:
  El.Display( rLS, "A x_{LS} - b" )
rLSTwoNorm = El.Nrm2(rLS)
rLSInfNorm = El.MaxNorm(rLS)
if worldRank == 0:
  print "|| A x_{LS} - b ||_2  =", rLSTwoNorm
  print "|| A x_{LS} - b ||_oo =", rLSInfNorm

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
