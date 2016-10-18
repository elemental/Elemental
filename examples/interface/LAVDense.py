#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

m = 500
n = 250
display = True
worldRank = El.mpi.WorldRank()
worldSize = El.mpi.WorldSize()

def Rectang(height,width):
  A = El.DistMatrix()
  El.Uniform( A, height, width )
  return A

A = Rectang(m,n)
b = El.DistMatrix()
El.Gaussian( b, m, 1 )
if display:
  El.Display( A, "A" )
  El.Display( b, "b" )

ctrl = El.LPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True
startLAV = El.mpi.Time()
x = El.LAV( A, b, ctrl )
endLAV = El.mpi.Time()
if worldRank == 0:
  print "LAV time:", endLAV-startLAV, "seconds"
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
rOneNorm = El.EntrywiseNorm( r, 1 )
if worldRank == 0:
  print "|| b ||_2       =", bTwoNorm
  print "|| b ||_oo      =", bInfNorm
  print "|| A x - b ||_2 =", rTwoNorm
  print "|| A x - b ||_1 =", rOneNorm

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
rLSOneNorm = El.EntrywiseNorm(rLS,1)
if worldRank == 0:
  print "|| A x_{LS} - b ||_2 =", rLSTwoNorm
  print "|| A x_{LS} - b ||_1 =", rLSOneNorm

# Require the user to press a button before the figures are closed
El.Finalize()
if worldSize == 1:
  raw_input('Press Enter to exit')
