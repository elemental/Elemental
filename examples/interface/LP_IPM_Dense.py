#
#  Copyright (c) 2009-2014, Jack Poulson
#  All rights reserved.
#
#  This file is part of Elemental and is under the BSD 2-Clause License, 
#  which can be found in the LICENSE file in the root directory, or at 
#  http://opensource.org/licenses/BSD-2-Clause
#
import El

m = 200
n = 400

def RectangDense(m,n):
  A = El.DistMatrix()
  El.Zeros( A, m, n )
  for s in xrange(m):
    A.Update( s, s, 11 )
    if s != 0:   A.Update( s, s-1,   1 )
    if s != n-1: A.Update( s, s+1,   2 )
    if s >= m:   A.Update( s, s-m,   3 )
    if s <  n-m: A.Update( s, s+m,   4 )

  return A

A = RectangDense(m,n)
El.Display( A, "A" )

# Generate a right-hand side in the positive image
# ================================================
xGen = El.DistMatrix()
El.Uniform(xGen,n,1,0.5,0.4999)
El.Display( xGen, "xGen" )
b = El.DistMatrix()
El.Zeros( b, m, 1 )
El.Gemv( El.NORMAL, 1., A, xGen, 0., b )
El.Display( b, "b" )

# Generate a random positive cost function
# ========================================
c = El.DistMatrix()
El.Uniform(c,n,1,0.5,0.4999)

# Generate random initial guesses
# ===============================
x = El.DistMatrix()
l = El.DistMatrix()
s = El.DistMatrix()
El.Uniform(x,n,1,0.5,0.4999)
El.Uniform(l,m,1,0.5,0.4999)
El.Uniform(s,n,1,0.5,0.4999)

muTol = 1e-10
rbTol = 1e-10
rcTol = 1e-10
maxIts = 1000
sigma = 0.9
gamma = 1e-3
beta = 1.5
psi = 100
progress = True
El.LinearProgramIPF(A,b,c,s,x,l,muTol,rbTol,rcTol,maxIts,
                    sigma,gamma,beta,psi,progress)
El.Display( s, "s" )
El.Display( x, "x" )
El.Display( l, "l" )

obj = El.Dot(c,x)
if El.mpi.WorldRank() == 0:
  print "c^T x =", obj

# Require the user to press a button before the figures are closed
commSize = El.mpi.Size( El.mpi.COMM_WORLD() )
El.Finalize()
if commSize == 1:
  raw_input('Press Enter to exit')
