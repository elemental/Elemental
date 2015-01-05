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

m = 2000
n = 1000
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

ctrl = El.LPAffineCtrl_d()
ctrl.mehrotraCtrl.progress = True
startCP = time.clock()
x = El.CP( A, b, ctrl )
endCP = time.clock()
if worldRank == 0:
  print "CP time: ", endCP-startCP

if display:
  El.Display( x, "x" )

bTwoNorm = El.Nrm2( b )
bInfNorm = El.MaxNorm( b )
r = El.DistMatrix()
El.Copy( b, r )
El.Gemv( El.NORMAL, -1., A, x, 1., r )
if display:
  El.Display( r, "r" )
rInfNorm = El.MaxNorm( r )
if worldRank == 0:
  print "|| b ||_2        =", bTwoNorm
  print "|| b ||_oo       =", bInfNorm
  print "|| A x - b ||_oo =", rInfNorm

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
