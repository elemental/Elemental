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

startNNLS = El.mpi.Time()
x = El.NNLS( A, b )
endNNLS = El.mpi.Time()
if worldRank == 0:
  print('NNLS time: {} seconds'.format(endNNLS-startNNLS))
if display:
  El.Display( x, "x" )

e = El.DistMatrix()
El.Copy( b, e )
El.Gemv( El.NORMAL, -1., A, x, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print('|| A x - b ||_2 = {}'.format(eTwoNorm))

startLS = El.mpi.Time()
xLS = El.LeastSquares( A, b )
endLS = El.mpi.Time()
if worldRank == 0:
  print('LS time: {} seconds'.format(endLS-startLS))
El.Copy( b, e )
El.Gemv( El.NORMAL, -1., A, xLS, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print('|| A x_{{LS}} - b ||_2 = {}'.format(eTwoNorm))

El.Finalize()
