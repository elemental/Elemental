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

m = 250
n = 500
lambda1 = 3
lambda2 = 4
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

ctrl = El.QPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True

if worldRank == 0:
  print "lambda1 =", lambda1, "lambda2 =", lambda2

startEN = time.clock()
x = El.EN( A, b, lambda1, lambda2, ctrl )
endEN = time.clock()
if worldRank == 0:
  print "EN time: ", endEN-startEN

if display:
  El.Display( x, "x" )

xOneNorm = El.EntrywiseNorm( x, 1 )
xTwoNorm = El.Nrm2( x )
e = El.DistMatrix()
El.Copy( b, e )
El.Gemv( El.NORMAL, -1., A, x, 1., e )
if display:
  El.Display( e, "e" )
eTwoNorm = El.Nrm2( e )
if worldRank == 0:
  print "|| x ||_1       =", xOneNorm
  print "|| x ||_2       =", xTwoNorm
  print "|| A x - b ||_2 =", eTwoNorm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
