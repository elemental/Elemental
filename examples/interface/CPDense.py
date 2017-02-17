#
#  Copyright (c) 2009-2016, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

m = 2000
n = 1000
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
startCP = El.mpi.Time()
x = El.CP( A, b, ctrl )
endCP = El.mpi.Time()
if worldRank == 0:
  print('CP time: {} seconds'.format(endCP-startCP))
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
  print('|| b ||_2        = {}'.format(bTwoNorm))
  print('|| b ||_oo       = {}'.format(bInfNorm))
  print('|| A x - b ||_2  = {}'.format(rTwoNorm))
  print('|| A x - b ||_oo = {}'.format(rInfNorm))

startLS = El.mpi.Time()
xLS = El.LeastSquares(A,b)
endLS = El.mpi.Time()
if worldRank == 0:
  print('LS time: {} seconds'.format(endLS-startLS))
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
  print('|| A x_{{LS}} - b ||_2  = {}'.format(rLSTwoNorm))
  print('|| A x_{{LS}} - b ||_oo = {}'.format(rLSInfNorm))

El.Finalize()
