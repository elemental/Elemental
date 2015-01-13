#
#  Copyright (c) 2009-2015, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El
import time

m = 500
n = 250
display = True
worldRank = El.mpi.WorldRank()

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

startNNLS = time.clock()
x = El.NNLS( A, b )
endNNLS = time.clock()
if worldRank == 0:
  print "NNLS time: ", endNNLS-startNNLS

if display:
  El.Display( x, "x" )

e = El.DistMatrix()
El.Copy( b, e )
El.Gemv( El.NORMAL, -1., A, x, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print "|| A x - b ||_2 =", eTwoNorm

# Note that El.LeastSquares overwrites A by default if A is dense
# (perhaps this should be changed)
ALS = El.DistMatrix()
El.Copy( A, ALS )
startLS = time.clock()
xLS = El.LeastSquares( ALS, b )
endLS = time.clock()
if worldRank == 0:
  print "LS time: ", endLS-startLS
El.Copy( b, e )
El.Gemv( El.NORMAL, -1., A, xLS, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print "|| A x_{LS} - b ||_2 =", eTwoNorm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
